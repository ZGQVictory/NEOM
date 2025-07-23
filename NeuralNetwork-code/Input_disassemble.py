import json
import os
import sys

def process_peptide_file(filename):
    """
    Process input.peptide file to extract 9-mers into a list.
    
    :param filename: Path to the peptide file.
    :return: A list of 9-mers.
    """
    peptides = []
    with open(filename, "r") as file:
        for line in file:
            line = line.strip()
            if len(line) == 9:  # Validate peptide length
                peptides.append(line)
            else:
                raise ValueError(f"Invalid peptide length in line: {line}")
    return peptides

def process_fasta_file(filename):
    """
    Process input.fasta file to extract sequences and split them into 9-mers.
    
    :param filename: Path to the fasta file.
    :return: A list of 9-mers.
    """
    peptides = []
    number = 0
    with open(filename, "r") as file:
        sequence = ""
        for line in file:
            if line.startswith(">"):  # Skip header lines
                number += 1
                continue
            sequence += line.strip()  # Concatenate sequence lines
        assert number == 1 # Assert that only 1 long peptide in this file

        # Split the sequence into 9-mers
        for i in range(0, len(sequence) - 8):  # Sliding window of 9
            peptides.append(sequence[i:i+9])
    return peptides

def update_initialization_json(start_peptide, template_file, output_file):
    """
    Update the 'start peptide' default value in Initialization_template.json.
    
    :param start_peptide: The starting peptide sequence to set.
    :param template_file: Path to the Initialization_template.json file.
    :param output_file: Path to save the updated JSON file.
    :return: None
    """
    with open(template_file, "r") as infile:
        data = json.load(infile)
    
    # Update the 'start peptide' field
    for param in data["parameters"]:
        if param["name"] == "start peptide":
            param["default"] = start_peptide
            break
    
    # Save the updated JSON to a new file
    with open(output_file, "w") as outfile:
        json.dump(data, outfile, indent=4)
    print(f"Updated JSON saved to {output_file}")

def main(filename, template_file="Initialization_template.json", output_dir=None):
    """
    Main function to handle input.peptide or input.fasta files and update JSON.
    
    :param filename: The input file (either input.peptide or input.fasta).
    :param template_file: Path to the Initialization_template.json file.
    :param output_file: Path to save the updated JSON file.
    :return: None
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"The file {filename} does not exist.")
    
    if filename.endswith(".peptide"):
        peptides = process_peptide_file(filename)
    elif filename.endswith(".fasta"):
        peptides = process_fasta_file(filename)
    else:
        raise ValueError("Unsupported file type. Please provide a .peptide or .fasta file.")
    
    if not peptides:
        raise ValueError("No valid peptides were extracted from the file.")
    
    # Use the first peptide in the list as the 'start peptide'
    start_peptides = list(set(peptides))

    for start_peptide in start_peptides:
        print(f"Selected start peptide: {start_peptide}")

        # Update the JSON file
        output_file = f"Initialization_{start_peptide}.json"
        if output_dir != None:
            output_file = os.path.join(output_dir, output_file)

        update_initialization_json(start_peptide, template_file, output_file)

# Example usage
if __name__ == "__main__":
    # Replace with your file paths
    input_file = sys.argv[1]  # "input.peptide" or "input.fasta"
    output_dir = sys.argv[2] 
    main(input_file, output_dir=output_dir)
