import subprocess
import gzip
from pathlib import Path
import shutil

def run_stride(pdb_file):
    # Define the relative path to the stride executable
    stride_executable = Path("stride/stride")
    
    # Ensure the stride executable exists
    if not stride_executable.exists():
        print(f"Error: Stride executable not found at {stride_executable}")
        return None

    # Create output directory if it doesn't exist
    output_folder = Path("stride_output")
    output_folder.mkdir(exist_ok=True)

    # Decompress .gz file to a temporary .pdb file
    decompressed_pdb_file = pdb_file.with_suffix('')
    with gzip.open(pdb_file, 'rb') as gz_file:
        with open(decompressed_pdb_file, 'wb') as decompressed_file:
            shutil.copyfileobj(gz_file, decompressed_file)
    
    # Define the output file path
    output_file = output_folder / f"{decompressed_pdb_file.stem}_stride_output.txt"

    # Run Stride on the decompressed PDB file and save the output
    with open(output_file, 'w') as out_file:
        try:
            subprocess.run(
                [str(stride_executable), str(decompressed_pdb_file)],  # Ensure paths are strings
                stdout=out_file,                          # Redirect output to file
                stderr=subprocess.PIPE,                   # Capture errors
                check=True                               # Raise error on failure
            )
            print(f"Stride output saved to: {output_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error running Stride on {decompressed_pdb_file}: {e.stderr.decode()}")

    # Clean up the decompressed temporary file
    decompressed_pdb_file.unlink()

    return output_file


if __name__ == "__main__":
    # Define the directory containing your PDB files
    pdb_dir = Path("UP000002311_559292_YEAST_v4")
    
    # Check if the directory exists
    if not pdb_dir.exists():
        print(f"Error: Directory {pdb_dir} does not exist.")
        exit(1)
    
    # Process .pdb.gz files only
    pdb_files = list(pdb_dir.glob("*.pdb.gz"))
    
    if not pdb_files:
        print(f"No PDB files found in directory: {pdb_dir}")
        exit(1)
    
    for pdb_file in pdb_files:
        print(f"Processing PDB file: {pdb_file.name}...")
        run_stride(pdb_file)
        print(f"Finished processing {pdb_file.name}\n")
    
    print("All files processed!")
