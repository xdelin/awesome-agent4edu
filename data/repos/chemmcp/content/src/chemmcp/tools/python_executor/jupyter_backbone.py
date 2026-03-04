import os
import time
import uuid
import json
import logging
import tempfile
import requests
import atexit
import signal
import re
import docker
import queue
from websocket import WebSocketApp
import threading
import subprocess

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Dockerfile for building executor image
DOCKERFILE = """\
FROM continuumio/miniconda3:latest
RUN apt-get update && apt-get install -y \
    cmake \
    pkg-config \
    build-essential \
    && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir jupyter jupyter_kernel_gateway rdkit matplotlib
"""


def has_gpu() -> bool:
    try:
        # this will return non-zero and raise if no NVIDIA driver / no GPUs
        subprocess.run(
            ["nvidia-smi"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

def strip_ansi(text: str) -> str:
    """Remove ANSI escape sequences for clean output."""
    ansi_re = re.compile(r'\x1B\[[0-?]*[ -/]*[@-~]')
    return ansi_re.sub('', text)

class KernelSession:
    """
    Simple Jupyter kernel session over WebSocketApp.
    """
    def __init__(self, base_url: str, convid: str, kernel_name: str = "python3"):
        self.base_url = base_url.rstrip('/')
        self.ws_url = self.base_url.replace('http://','ws://').replace('https://','wss://')
        self.convid = convid
        self.kernel_name = kernel_name
        self.msg_queue = queue.Queue()

        # HTTP connectivity check
        try:
            resp = requests.get(f"{self.base_url}/api", timeout=5)
            resp.raise_for_status()
        except requests.RequestException as e:
            raise ConnectionError(f"Cannot reach Jupyter Gateway at {self.base_url}: {e}")

        # Start kernel
        self.kernel_id = self._start_kernel()
        # Initialize WebSocketApp
        self._init_ws_app()
        # Allow connection
        time.sleep(1)
        # Disable ANSI coloring
        self.execute("%colors nocolor")

    def _start_kernel(self) -> str:
        resp = requests.post(f"{self.base_url}/api/kernels", json={"name": self.kernel_name})
        resp.raise_for_status()
        kid = resp.json()["id"]
        logger.info(f"Started kernel {kid}")
        return kid

    def _init_ws_app(self):
        """
        Create WebSocketApp to push incoming messages into a queue.
        """
        # Tear down existing WS if any
        try:
            self.ws_app.close()
            self._ws_thread.join(timeout=1)
        except Exception:
            pass
        # Clear queue
        try:
            self.msg_queue.queue.clear()
        except Exception:
            pass

        endpoint = f"{self.ws_url}/api/kernels/{self.kernel_id}/channels"
        def on_message(ws, message):
            self.msg_queue.put(message)
        def on_error(ws, error):
            logger.error(f"WebSocket error: {error}")
        def on_close(ws, code, reason):
            logger.info(f"WebSocket closed: {code} {reason}")
        def on_open(ws):
            logger.info("WebSocket connection opened")

        self.ws_app = WebSocketApp(
            endpoint,
            on_message=on_message,
            on_error=on_error,
            on_close=on_close,
            on_open=on_open,
        )
        self._ws_thread = threading.Thread(
            target=self.ws_app.run_forever,
            daemon=True
        )
        self._ws_thread.start()

    def execute(self, code: str, timeout: int = 60) -> list:
        """
        Execute code and return list of messages (stream, result, data, error).
        """
        msg_id = uuid.uuid4().hex
        header = {"msg_id": msg_id, "username": "", "session": self.convid,
                  "msg_type": "execute_request", "version": "5.0"}
        msg = {"header": header, "parent_header": {}, "metadata": {},
               "content": {"code": code, "silent": False,
                           "store_history": False, "user_expressions": {},
                           "allow_stdin": False},
               "channel": "shell"}
        # Send request
        self.ws_app.send(json.dumps(msg))

        messages = []
        start = time.time()
        while True:
            elapsed = time.time() - start
            if elapsed > timeout:
                raise TimeoutError(f"Execution timed out after {timeout}s")
            try:
                raw = self.msg_queue.get(timeout=timeout - elapsed)
            except Exception:
                continue
            try:
                m = json.loads(raw)
            except ValueError:
                continue
            if m.get('parent_header', {}).get('msg_id') != msg_id:
                continue
            mtype = m.get('msg_type')
            content = m.get('content', {})
            if mtype in ('stream','execute_result','display_data','error'):
                messages.append({'msg_type': mtype, 'content': content})
            if mtype in ('execute_reply','error'):
                break
        return messages

    def shutdown(self):
        try:
            self.ws_app.close()
            self._ws_thread.join(timeout=1)
        except Exception:
            pass
        try:
            requests.delete(f"{self.base_url}/api/kernels/{self.kernel_id}")
            logger.info(f"Shutdown kernel {self.kernel_id}")
        except Exception as e:
            logger.debug(f"Error shutting down kernel: {e}")

class JupyterBackbone:
    """
    Manages Jupyter Kernel Gateway in Docker.
    """
    DEFAULT_IMAGE = "chemmcp-python-executor:latest"
    JUPYTER_PORT = 8888

    def __init__(self, kernel_url: str = None,
                 image: str = None, mem_limit: str = "8g", cpus: float = 2.0):
        self.conv_id = str(uuid.uuid4())
        self.image = image or self.DEFAULT_IMAGE
        self.mem_limit = mem_limit
        self.cpus = cpus
        self.container = None
        self.session = None

        # Cleanup on exit/signals
        atexit.register(self.close)
        for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
            try:
                signal.signal(sig, lambda s, f: self.close())
            except Exception:
                pass

        if kernel_url:
            self.kernel_url = kernel_url.rstrip('/')
        else:
            client = docker.from_env()
            # Ensure or build image
            try:
                client.images.get(self.image)
                logger.info(f"Docker image {self.image} found")
            except docker.errors.ImageNotFound:
                logger.info(f"Building Docker image {self.image}")
                tmp = tempfile.mkdtemp()
                with open(os.path.join(tmp, 'Dockerfile'), 'w') as df:
                    df.write(DOCKERFILE)
                client.images.build(path=tmp, tag=self.image)
                logger.info("Image built locally")
            # Launch container
            ports = {f"{self.JUPYTER_PORT}/tcp": None}
            device_req = docker.types.DeviceRequest(count=-1, capabilities=[['gpu']])
            cmd = (
                f"jupyter kernelgateway --ip=0.0.0.0 --port={self.JUPYTER_PORT}"
            )
            run_kwargs = dict(
                image=self.image,
                command=cmd,
                detach=True,
                remove=True,
                ports=ports,
                mem_limit=self.mem_limit,
                nano_cpus=int(self.cpus * 1e9),
            )
            if has_gpu():
                device_req = docker.types.DeviceRequest(count=-1, capabilities=[["gpu"]])
                run_kwargs["device_requests"] = [device_req]
            self.container = client.containers.run(**run_kwargs)
            self.container.reload()
            binding = self.container.attrs['NetworkSettings']['Ports'][f"{self.JUPYTER_PORT}/tcp"][0]
            port = int(binding['HostPort'])
            self.kernel_url = f"http://localhost:{port}"
            logger.info(f"Gateway container started at {self.kernel_url}")

        self.session = KernelSession(self.kernel_url, self.conv_id)

    def execute_jupyter(self, code: str, timeout: int = 600) -> list:
        logger.info(f"Executing (conv_id={self.conv_id}): {code}")
        msgs = self.session.execute(code, timeout)
        return msgs

    def close(self):
        if self.session:
            self.session.shutdown()
            self.session = None
        if self.container:
            try:
                self.container.stop()
            except Exception as e:
                logger.debug(f"Error stopping container: {e}")
            self.container = None

    def __del__(self):
        self.close()
