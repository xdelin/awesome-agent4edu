import os
import json
import stat
import datetime
import mimetypes
import subprocess
import uuid
from pathlib import Path
from typing import List, Optional, Dict, Any, Union
from fastapi import FastAPI, HTTPException, Query, Path as PathParam, Body, Response
from fastapi.responses import JSONResponse, FileResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field
from fastapi_mcp import FastApiMCP

app = FastAPI(title="Docker Filesystem API", description="API for exploring the Docker container filesystem")

# Mount static files directory for direct file access
app.mount("/static", StaticFiles(directory="/"), name="static")

# Define allowed base directories for security
ALLOWED_BASE_DIRS = [
    "/manim",
    "/app",
    "/media",
    "/usr/local",
    "/tmp"
]

@app.get("/")
def read_root():
    """API root endpoint"""
    return {
        "message": "Docker Filesystem API is running",
        "endpoints": [
            {
                "path": "/list-files",
                "description": "List files and directories in the container filesystem"
            }
        ]
    }

@app.get("/list-files")
def list_files(
    directory: str = Query("/", description="Directory path to list (must be under an allowed base directory)"),
    recursive: bool = Query(False, description="Whether to list files recursively"),
    show_hidden: bool = Query(False, description="Whether to show hidden files (starting with .)"),
    max_depth: int = Query(1, description="Maximum depth for recursive listing"),
    pattern: Optional[str] = Query(None, description="Filter files by pattern (glob syntax)")
):
    """
    List files and directories in the Docker container filesystem.
    
    Security restrictions:
    - Only allowed base directories can be accessed
    - Path traversal attempts are blocked
    - Certain system directories are inaccessible
    """
    
    # Normalize the path
    directory = os.path.normpath(directory)
    
    # Security check: Prevent path traversal attacks
    if ".." in directory.split(os.sep):
        raise HTTPException(
            status_code=403, 
            detail="Path traversal attempts are not allowed"
        )
    
    # Security check: Ensure the directory is under an allowed base directory
    if not any(directory == base or directory.startswith(f"{base}/") for base in ALLOWED_BASE_DIRS):
        raise HTTPException(
            status_code=403,
            detail=f"Access is only allowed to these base directories: {', '.join(ALLOWED_BASE_DIRS)}"
        )
    
    # Check if directory exists
    if not os.path.exists(directory):
        raise HTTPException(
            status_code=404,
            detail=f"Directory not found: {directory}"
        )
    
    # Check if path is a directory
    if not os.path.isdir(directory):
        raise HTTPException(
            status_code=400,
            detail=f"Path is not a directory: {directory}"
        )
    
    # List files and directories
    results = []
    
    if recursive:
        for root, dirs, files in os.walk(directory):
            # Check depth
            relative_path = os.path.relpath(root, directory)
            depth = len(relative_path.split(os.sep)) if relative_path != "." else 0
            
            if depth > max_depth:
                continue
                
            # Filter hidden files/dirs if needed
            if not show_hidden:
                dirs[:] = [d for d in dirs if not d.startswith('.')]
                files = [f for f in files if not f.startswith('.')]
                
            # Pattern filtering
            if pattern:
                import fnmatch
                files = [f for f in files if fnmatch.fnmatch(f, pattern)]
                
            # Get file information
            for filename in files:
                file_path = os.path.join(root, filename)
                results.append(get_file_info(file_path, directory))
                
            # Add directories with trailing slash
            for dirname in dirs:
                dir_path = os.path.join(root, dirname)
                results.append(get_file_info(dir_path, directory))
    else:
        # Non-recursive listing
        try:
            entries = os.listdir(directory)
            
            # Filter hidden files if needed
            if not show_hidden:
                entries = [entry for entry in entries if not entry.startswith('.')]
                
            # Pattern filtering
            if pattern:
                import fnmatch
                entries = [entry for entry in entries if fnmatch.fnmatch(entry, pattern)]
                
            # Get file information
            for entry in entries:
                entry_path = os.path.join(directory, entry)
                results.append(get_file_info(entry_path, directory))
                
        except PermissionError:
            raise HTTPException(
                status_code=403,
                detail=f"Permission denied: {directory}"
            )
    
    # Sort results: directories first, then files alphabetically
    results.sort(key=lambda x: (not x["is_dir"], x["name"]))
    
    return {
        "directory": directory,
        "parent_directory": os.path.dirname(directory) if directory != "/" else None,
        "count": len(results),
        "results": results
    }

def get_file_info(file_path: str, base_dir: str) -> Dict[str, Any]:
    """Get detailed information about a file or directory"""
    try:
        stat_info = os.stat(file_path)
        is_dir = os.path.isdir(file_path)
        rel_path = os.path.relpath(file_path, base_dir) if base_dir != file_path else ""
        
        # Get file size (0 for directories)
        size = stat_info.st_size if not is_dir else 0
        
        # Format timestamps
        mtime = datetime.datetime.fromtimestamp(stat_info.st_mtime).isoformat()
        ctime = datetime.datetime.fromtimestamp(stat_info.st_ctime).isoformat()
        
        # Get file permissions
        mode = stat_info.st_mode
        perms = ""
        for who in "USR", "GRP", "OTH":
            for what in "R", "W", "X":
                perms += what if mode & getattr(stat, f"S_I{what}{who}") else "-"
        
        return {
            "name": os.path.basename(file_path),
            "path": file_path,
            "relative_path": rel_path,
            "is_dir": is_dir,
            "size": size,
            "size_human": format_size(size),
            "modified_time": mtime,
            "created_time": ctime,
            "permissions": perms,
            "static_url": f"/static{file_path}" if not is_dir else None
        }
    except (FileNotFoundError, PermissionError):
        # Return minimal info if we can't access the file
        return {
            "name": os.path.basename(file_path),
            "path": file_path,
            "relative_path": os.path.relpath(file_path, base_dir) if base_dir != file_path else "",
            "is_dir": os.path.isdir(file_path) if os.path.exists(file_path) else None,
            "error": "Permission denied or file not found"
        }

def format_size(size_bytes: int) -> str:
    """Format file size in human-readable format"""
    if size_bytes == 0:
        return "0B"
    
    size_names = ("B", "KB", "MB", "GB", "TB")
    i = 0
    while size_bytes >= 1024 and i < len(size_names) - 1:
        size_bytes /= 1024
        i += 1
    
    return f"{size_bytes:.2f}{size_names[i]}"

# Create a model for the write file request
class WriteFileRequest(BaseModel):
    content: str = Field(..., description="Content to write to the file")
    
@app.post("/write-file")
def write_file(
    filepath: str = Query(..., description="Path where the file should be written, including filename"),
    overwrite: bool = Query(False, description="Whether to overwrite existing files"),
    create_dirs: bool = Query(True, description="Whether to create parent directories if they don't exist"),
    file_data: WriteFileRequest = Body(..., description="File content to write")
):
    """
    Write content to a file in the Docker container filesystem.
    
    Security restrictions:
    - Only allowed base directories can be written to
    - Path traversal attempts are blocked
    - Existing files are not overwritten unless explicitly specified
    """
    
    # Get the content from the request
    content = file_data.content
    
    # Normalize the path
    filepath = os.path.normpath(filepath)
    
    # Security check: Prevent path traversal attacks
    if ".." in filepath.split(os.sep):
        raise HTTPException(
            status_code=403, 
            detail="Path traversal attempts are not allowed"
        )
    
    # Security check: Ensure the file is under an allowed base directory
    if not any(filepath == base or filepath.startswith(f"{base}/") for base in ALLOWED_BASE_DIRS):
        raise HTTPException(
            status_code=403,
            detail=f"Writing is only allowed to these base directories: {', '.join(ALLOWED_BASE_DIRS)}"
        )
    
    # Get the directory part of the path
    directory = os.path.dirname(filepath)
    
    # Check if the directory exists, create it if it doesn't and create_dirs is True
    if not os.path.exists(directory):
        if create_dirs:
            try:
                os.makedirs(directory, exist_ok=True)
            except PermissionError:
                raise HTTPException(
                    status_code=403,
                    detail=f"Permission denied: Cannot create directory {directory}"
                )
            except Exception as e:
                raise HTTPException(
                    status_code=500,
                    detail=f"Failed to create directory {directory}: {str(e)}"
                )
        else:
            raise HTTPException(
                status_code=404,
                detail=f"Directory not found: {directory}"
            )
    
    # Check if the file already exists
    if os.path.exists(filepath) and not overwrite:
        raise HTTPException(
            status_code=409,
            detail=f"File already exists: {filepath}. Use overwrite=true to replace it."
        )
    
    # Try to write the file
    try:
        with open(filepath, 'w') as file:
            file.write(content)
            
        # Get file info for the response
        file_info = get_file_info(filepath, os.path.dirname(filepath))
        
        return {
            "status": "success",
            "message": "File written successfully",
            "file": file_info,
            "static_url": f"/static{filepath}"
        }
    except PermissionError:
        raise HTTPException(
            status_code=403,
            detail=f"Permission denied: Cannot write to {filepath}"
        )
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to write file: {str(e)}"
        )
@app.get("/download-file")
def download_file(
    filepath: str = Query(..., description="Path to the file to download"),
    as_attachment: bool = Query(True, description="Whether to download as attachment or display in browser"),
    filename: Optional[str] = Query(None, description="Optional filename to use for the download (defaults to original filename)")
):
    """
    Download a file from the Docker container filesystem with proper content disposition headers.
    
    Security restrictions:
    - Only allowed base directories can be accessed
    - Path traversal attempts are blocked
    """
    
    # Normalize the path
    filepath = os.path.normpath(filepath)
    
    # Security check: Prevent path traversal attacks
    if ".." in filepath.split(os.sep):
        raise HTTPException(
            status_code=403, 
            detail="Path traversal attempts are not allowed"
        )
    
    # Security check: Ensure the file is under an allowed base directory
    if not any(filepath == base or filepath.startswith(f"{base}/") for base in ALLOWED_BASE_DIRS):
        raise HTTPException(
            status_code=403,
            detail=f"Access is only allowed to these base directories: {', '.join(ALLOWED_BASE_DIRS)}"
        )
    
    # Check if the file exists
    if not os.path.exists(filepath):
        raise HTTPException(
            status_code=404,
            detail=f"File not found: {filepath}"
        )
    
    # Check if the file is a regular file
    if not os.path.isfile(filepath):
        raise HTTPException(
            status_code=400,
            detail=f"Path is not a file: {filepath}"
        )
    
    # Get filename if not provided
    download_filename = filename or os.path.basename(filepath)
    
    # Return file response with appropriate content disposition
    return FileResponse(
        path=filepath,
        filename=download_filename,
        media_type=mimetypes.guess_type(filepath)[0],
        content_disposition_type="attachment" if as_attachment else "inline"
    ) 

@app.post("/run-manim")
def run_manim(
    filepath: str = Query(..., description="Path to the Python file with Manim scenes (must be under an allowed base directory)"),
    scene_name: str = Query(..., description="Name of the scene class to render"),
    quality: str = Query("medium_quality", 
                        description="Quality setting that affects rendering resolution and speed: "
                                   "low_quality (-ql): 480p, 15fps - Fast rendering, good for previews; "
                                   "medium_quality (-qm): 720p, 30fps - Good balance of quality and speed; "
                                   "high_quality (-qh): 1080p, 60fps - High quality but slower rendering; "
                                   "production_quality (-qk): 1440p, 60fps - Highest quality, slowest rendering"),
    preview: bool = Query(False, description="Whether to automatically open the output file after rendering (adds -p flag)"),
    format: Optional[str] = Query(None, description="Output format (e.g., 'png', 'gif', 'mp4', 'webm') - adds --format flag"),
    transparent: bool = Query(False, description="Renders the background as transparent, if possible (adds -t flag)"),
    save_last_frame: bool = Query(False, description="Only renders the last frame of the scene (adds --save_last_frame flag)"),
    from_animation_number: Optional[int] = Query(None, description="Start rendering from a specific animation number (adds -n flag)"),
    upto_animation_number: Optional[int] = Query(None, description="End rendering at a specific animation number (adds -n flag)"),
    resolution: Optional[str] = Query(None, description="Resolution in the format WIDTHxHEIGHT (e.g., '1920x1080') - adds -r flag"),
    frame_rate: Optional[int] = Query(None, description="Frame rate in frames per second - adds -f flag"),
    color: Optional[str] = Query(None, description="Background color of the scene (e.g., '#ffffff', 'WHITE') - adds -c flag"),
    additional_args: Optional[List[str]] = Query(None, description="Additional command-line arguments to pass to Manim")
):
    """
    Run Manim on a Python file in the Docker container to generate an animation.
    
    Manim is a mathematical animation engine that allows you to create explanatory math videos.
    This endpoint runs Manim with the specified options and returns links to download the generated videos/images.
    
    Security restrictions:
    - Only allowed base directories can be accessed
    - Path traversal attempts are blocked
    """
    
    # Normalize the path
    filepath = os.path.normpath(filepath)
    
    # Security check: Prevent path traversal attacks
    if ".." in filepath.split(os.sep):
        raise HTTPException(
            status_code=403, 
            detail="Path traversal attempts are not allowed"
        )
    
    # Security check: Ensure the file is under an allowed base directory
    if not any(filepath == base or filepath.startswith(f"{base}/") for base in ALLOWED_BASE_DIRS):
        raise HTTPException(
            status_code=403,
            detail=f"Access is only allowed to these base directories: {', '.join(ALLOWED_BASE_DIRS)}"
        )
    
    # Check if the file exists
    if not os.path.exists(filepath):
        raise HTTPException(
            status_code=404,
            detail=f"File not found: {filepath}"
        )
    
    # Check if the path is a file
    if not os.path.isfile(filepath):
        raise HTTPException(
            status_code=400,
            detail=f"Path is not a file: {filepath}"
        )
    
    # Generate a unique job ID
    job_id = str(uuid.uuid4())
    
    # Prepare output directory
    output_dir = f"/manim/output/{job_id}"
    os.makedirs(output_dir, exist_ok=True)
    
    # Build the command
    cmd = ["python3", "-m", "manim"]
    
    # Add quality flag
    if quality == "low_quality":
        cmd.append("-ql")
    elif quality == "medium_quality":
        cmd.append("-qm")
    elif quality == "high_quality":
        cmd.append("-qh")
    elif quality == "production_quality":
        cmd.append("-qk")
    
    # Add preview flag if requested
    if preview:
        cmd.append("-p")
    
    # Add transparent flag if requested
    if transparent:
        cmd.append("-t")
    
    # Add save_last_frame flag if requested
    if save_last_frame:
        cmd.append("--save_last_frame")
    
    # Add format flag if provided
    if format:
        cmd.extend(["--format", format])
    
    # Add resolution flag if provided
    if resolution:
        cmd.extend(["-r", resolution])
    
    # Add frame rate flag if provided
    if frame_rate:
        cmd.extend(["-f", str(frame_rate)])
    
    # Add color flag if provided
    if color:
        cmd.extend(["-c", color])
    
    # Add animation number range if provided
    if from_animation_number is not None or upto_animation_number is not None:
        anim_range = []
        if from_animation_number is not None:
            anim_range.append(str(from_animation_number))
        else:
            anim_range.append("")
            
        if upto_animation_number is not None:
            anim_range.append(str(upto_animation_number))
            
        cmd.extend(["-n", ",".join(anim_range)])
    
    # Add output directory
    cmd.extend(["--output_file", output_dir])
    
    # Add the file path and scene name
    cmd.append(filepath)
    cmd.append(scene_name)
    
    # Add any additional arguments
    if additional_args:
        cmd.extend(additional_args)
    
    try:
        # Run the command
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        # Check both the output directory and the parent output directory for files
        # Manim sometimes places files directly in the output directory with the job ID as part of the filename
        files = []
        file_info = []
        
        # Check in the job-specific output directory
        if os.path.exists(output_dir):
            dir_files = [f for f in os.listdir(output_dir) if os.path.isfile(os.path.join(output_dir, f))]
            for file in dir_files:
                file_path = os.path.join(output_dir, file)
                info = get_file_info(file_path, output_dir)
                info["download_url"] = f"/download-file?filepath={file_path}"
                info["static_url"] = f"/static{file_path}"
                file_info.append(info)
                files.append(file)
        
        # Check for files in the parent output directory that match the job_id pattern
        parent_output_dir = "/manim/output"
        if os.path.exists(parent_output_dir):
            parent_files = [f for f in os.listdir(parent_output_dir) 
                          if os.path.isfile(os.path.join(parent_output_dir, f)) and job_id in f]
            for file in parent_files:
                file_path = os.path.join(parent_output_dir, file)
                info = get_file_info(file_path, parent_output_dir)
                info["download_url"] = f"/download-file?filepath={file_path}"
                info["static_url"] = f"/static{file_path}"
                file_info.append(info)
                files.append(file)
        
        return {
            "job_id": job_id,
            "status": "success",
            "command": " ".join(cmd),
            "output": result.stdout,
            "files": file_info,
            "output_directory": output_dir
        }
    
    except subprocess.CalledProcessError as e:
        return {
            "job_id": job_id,
            "status": "error",
            "command": " ".join(cmd),
            "error": e.stderr,
            "returncode": e.returncode
        }

# Create the MCP server and mount it to the FastAPI app
mcp = FastApiMCP(app)
mcp.mount() 