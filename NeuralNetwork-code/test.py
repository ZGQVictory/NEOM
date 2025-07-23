import pandas as pd
import os

def process_accepted_peptides(logfile_path, output_excel="accepted_peptides_results.xlsx"):
    """
    Process the Output_accepted_peptides.log file to extract the starting peptide and accepted peptides,
    and save them into an Excel file.

    :param logfile_path: Path to the Output_accepted_peptides.log file.
    :param output_excel: Path to save the extracted information as an Excel file.
    :return: None
    """
    try:
        with open(logfile_path, "r") as file:
            lines = file.readlines()
        
        if len(lines) < 3:
            raise ValueError("The log file does not contain the expected structure.")

        # Extract the starting peptide label and sequence
        starting_label = lines[0].strip()  # "starting peptide":
        if not starting_label.startswith('"starting peptide":'):
            raise ValueError("The first line does not indicate the starting peptide.")
        
        starting_peptide = lines[1].strip()  # The actual starting peptide sequence

        # Ensure the third line is the accepted peptides label
        accepted_label = lines[2].strip()  # "accepted peptides":
        if not accepted_label.startswith('"accepted peptides":'):
            raise ValueError("The third line does not indicate accepted peptides.")

        # Extract the accepted peptides
        accepted_peptides = [line.strip() for line in lines[3:] if line.strip()]  # Remaining lines

        # Create a DataFrame
        data = {
            "Type": ["Starting Peptide"] + ["Accepted Peptide"] * len(accepted_peptides),
            "Peptide Sequence": [starting_peptide] + accepted_peptides,
        }
        df = pd.DataFrame(data)

        # Save to Excel
        df.to_excel(output_excel, index=False)
        print(f"Processed peptides saved to {output_excel}")

    except Exception as e:
        print(f"An error occurred while processing {logfile_path}: {e}")

output_results_path = "D:\\AllData\\Science\\pmhc-I\\Network-version\\NeuralNetwork-code\\usr\\initialization_VMNILLQQQ_output_20241121_215633"

accepted_peptides_logfile = os.path.join(output_results_path, "Output_accepted_peptides.log")
if os.path.exists(accepted_peptides_logfile):
    process_accepted_peptides(accepted_peptides_logfile, output_excel=f"{output_results_path}/accepted_peptides_results.xlsx")
else:
    print(f"{accepted_peptides_logfile} not found. Skipping peptide processing.")