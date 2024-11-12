import gzip
from pathlib import Path
import re

def decompress_pdb_gz(pdb_gz_file):
    """Decompress a .pdb.gz file and return its content as a string."""
    with gzip.open(pdb_gz_file, 'rt') as gz_file:
        return gz_file.read()

def parse_stride_file(stride_file):
    """Parse the Stride output file to extract beta sheet regions."""
    beta_sheets = []
    with open(stride_file, 'r') as file:
        for line in file:
            # Extract only lines that indicate beta sheets ("LOC" and "Strand")
            if line.startswith("ASG") and "Strand" in line:
                parts = line.split()
                print(parts)
                if len(parts) >= 9:
                    start_res = (parts[1], parts[3], parts[2])  # (Residue type, number, chain)
                    #end_res = (parts[6], parts[7], parts[8])    # (Residue type, number, chain)
                    beta_sheets.append((start_res))
    print(beta_sheets)
    return beta_sheets


def parse_pdb_content(pdb_content):
    """Parse the PDB content to extract residues and their confidence scores."""
    residues = []
    for line in pdb_content.splitlines():
        if line.startswith("ATOM"):
            try:
                residue_type = line[17:20].strip()  # Residue type (e.g., GLU)
                chain_id = line[21].strip()        # Chain ID (e.g., A)
                residue_number = int(line[22:26].strip())  # Residue number (e.g., 1072)
                confidence_score = float(line[60:66].strip())  # Confidence score (e.g., 85.0)

                # Only include residues with a confidence score >= 85
                if confidence_score >= 85:
                    residues.append((residue_type, chain_id, residue_number, confidence_score))
            except ValueError:
                continue  # Skip lines where parsing fails
    return residues

def extract_beta_sequences(pdb_residues, beta_sheets):
    """Extract sequences within beta sheet regions that meet the confidence score criteria."""
    filtered_sequences = []
    
    for start_res in beta_sheets:
        start_type, start_num, start_chain = start_res
        #end_type, end_num, end_chain = end_res
        start_num = int(start_num)
        #end_num = int(end_num)

        # Filter residues that fall within the beta sheet region and have a confidence score >= 85
        sequence = [
            (res_type, chain_id, res_num, conf_score)
            for (res_type, chain_id, res_num, conf_score) in pdb_residues
            if chain_id == start_chain and start_num <= res_num 
        ]

        if sequence:
            filtered_sequences.append(sequence)
    
    return filtered_sequences

def save_sequences_to_file(sequences, output_file):
    """Save the filtered sequences to a text file."""
    output_file.parent.mkdir(exist_ok=True, parents=True)
    with open(output_file, 'w') as file:
        file.write("Residue | Chain | Residue Number | Confidence Score\n")
        file.write("-------------------------------------------------\n")
        for sequence in sequences:
            for res_type, chain_id, res_num, conf_score in sequence:
                file.write(f"{res_type} {chain_id} {res_num} | {conf_score:.3f}\n")
            file.write("\n")

def process_pdb_and_stride(pdb_gz_path):
    """Main function to process the PDB and corresponding Stride file."""
    pdb_gz_path = Path(pdb_gz_path)
    
    if not pdb_gz_path.exists() or not pdb_gz_path.suffix == ".gz":
        print(f"Error: File '{pdb_gz_path}' not found or is not a .pdb.gz file.")
        return

    # Load the PDB content from the .pdb.gz file
    pdb_content = decompress_pdb_gz(pdb_gz_path)

    # Remove the ".pdb" part from the stem to match your Stride file naming convention
    pdb_stem = pdb_gz_path.stem.replace(".pdb", "")

    # Locate the corresponding Stride output file in the "stride_output" folder
    stride_output_folder = Path("stride_output")
    stride_file = stride_output_folder / f"{pdb_stem}_stride_output.txt"
    
    if not stride_file.exists():
        print(f"Error: Stride output file '{stride_file}' not found.")
        return

    # Parse the Stride output and PDB content
    beta_sheets = parse_stride_file(stride_file)
    pdb_residues = parse_pdb_content(pdb_content)
    
    # Extract sequences in beta sheets with confidence score >= 85
    filtered_sequences = extract_beta_sequences(pdb_residues, beta_sheets)

    # Define the output file path in the "extracted" folder
    output_file = Path("extracted") / f"{pdb_stem}_filtered.txt"
    save_sequences_to_file(filtered_sequences, output_file)

    print(f"Processing completed. Results saved to: {output_file}")

# Usage Example
if __name__ == "__main__":
    pdb_gz_path = "AF-Q99394-F1-model_v4.pdb.gz"
    process_pdb_and_stride(pdb_gz_path)