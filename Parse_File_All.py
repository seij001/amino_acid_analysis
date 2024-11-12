import gzip
from pathlib import Path

def decompress_pdb_gz(pdb_gz_file):
    """Decompress a .pdb.gz file and return its content as a string."""
    with gzip.open(pdb_gz_file, 'rt') as gz_file:
        return gz_file.read()

def parse_stride_file(stride_file):
    """Parse the Stride output file to extract beta sheet regions."""
    beta_sheets = []
    with open(stride_file, 'r') as file:
        for line in file:
            if line.startswith("ASG") and "Strand" in line:
                parts = line.split()
                if len(parts) >= 9:
                    start_res = (parts[1], parts[3], parts[2])  # (Residue type, number, chain)
                    beta_sheets.append((start_res))
    return beta_sheets

def parse_pdb_content(pdb_content):
    """Parse the PDB content to extract residues and their confidence scores."""
    residues = []
    for line in pdb_content.splitlines():
        if line.startswith("ATOM"):
            try:
                residue_type = line[17:20].strip()
                chain_id = line[21].strip()
                residue_number = int(line[22:26].strip())
                confidence_score = float(line[60:66].strip())

                if confidence_score >= 85:
                    residues.append((residue_type, chain_id, residue_number, confidence_score))
            except ValueError:
                continue
    return residues

def extract_beta_sequences(pdb_residues, beta_sheets):
    """Extract sequences within beta sheet regions that meet the confidence score criteria."""
    filtered_sequences = set()
    
    for start_res in beta_sheets:
        start_type, start_num, start_chain = start_res
        start_num = int(start_num)

        for (res_type, chain_id, res_num, conf_score) in pdb_residues:
            if chain_id == start_chain and start_num <= res_num:
                filtered_sequences.add((res_type, chain_id, res_num, conf_score))
    
    return list(filtered_sequences)

def save_sequences_to_file(sequences, output_file):
    """Save the filtered sequences to a text file."""
    output_file.parent.mkdir(exist_ok=True, parents=True)
    with open(output_file, 'w') as file:
        file.write("Residue | Chain | Residue Number | Confidence Score\n")
        file.write("-------------------------------------------------\n")
        for res_type, chain_id, res_num, conf_score in sequences:
            file.write(f"{res_type} {chain_id} {res_num} | {conf_score:.3f}\n")

def process_all_pdb_files(directory_path):
    """Process all PDB files in the directory and save individual results for each file."""
    directory_path = Path(directory_path)
    if not directory_path.exists() or not directory_path.is_dir():
        print(f"Error: Directory '{directory_path}' not found or is not a directory.")
        return
    
    for pdb_gz_file in directory_path.glob("*.pdb.gz"):
        print(f"Processing {pdb_gz_file.name}...")
        
        pdb_content = decompress_pdb_gz(pdb_gz_file)
        pdb_stem = pdb_gz_file.stem.replace(".pdb", "")

        stride_output_folder = Path("stride_output")
        stride_file = stride_output_folder / f"{pdb_stem}_stride_output.txt"
        
        if not stride_file.exists():
            print(f"Warning: Stride output file '{stride_file}' not found. Skipping this file.")
            continue

        beta_sheets = parse_stride_file(stride_file)
        pdb_residues = parse_pdb_content(pdb_content)
        
        filtered_sequences = extract_beta_sequences(pdb_residues, beta_sheets)

        # Save each PDB file's results in a separate output file
        output_file = Path("extracted") / f"{pdb_stem}_filtered.txt"
        save_sequences_to_file(filtered_sequences, output_file)
        print(f"Results saved to: {output_file}")

# Usage Example
if __name__ == "__main__":
    input_directory = "D:/abroad/Duke/CBB520/assginment3/amino_acid_analysis/UP000002311_559292_YEAST_v4"
    process_all_pdb_files(input_directory)
