import asyncio


class BaseAppleScriptOperations:
    """Base class with common AppleScript execution functionality."""

    @staticmethod
    async def execute_applescript(script: str) -> str:
        """Execute AppleScript and return result."""
        process = await asyncio.create_subprocess_exec(
            "osascript",
            "-e",
            script,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
        )
        stdout, stderr = await process.communicate()

        if process.returncode != 0:
            raise RuntimeError(f"AppleScript error: {stderr.decode()}")

        return stdout.decode().strip()
