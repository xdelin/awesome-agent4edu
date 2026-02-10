import requests
import zipfile
import os
import tempfile
from tqdm import tqdm

def download_and_extract_zenodo_zip(zenodo_url: str, output_folder: str):
    """
    Download a zip file from Zenodo to the system temp folder, extract it to the specified folder, and delete the zip file.
    Shows a download progress bar.
    
    Args:
        zenodo_url (str): Direct URL to the zip file hosted on Zenodo.
        output_folder (str): Path to the folder where files will be extracted.
    """
    os.makedirs(output_folder, exist_ok=True)
    
    # Create a temp file path
    with tempfile.NamedTemporaryFile(delete=False, suffix=".zip") as tmp_file:
        temp_zip_path = tmp_file.name

    # Download with progress bar
    with requests.get(zenodo_url, stream=True) as r:
        r.raise_for_status()
        total_size = int(r.headers.get('content-length', 0))
        with open(temp_zip_path, 'wb') as f, tqdm(
            total=total_size, unit='B', unit_scale=True, desc="Downloading"
        ) as pbar:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))
    
    # Extract
    with zipfile.ZipFile(temp_zip_path, 'r') as zip_ref:
        zip_ref.extractall(output_folder)
    
    # Clean up
    os.remove(temp_zip_path)
