import numpy as np

#############################################
# This script is used to calculate RMSD values of each conformer related to the
# first one, which is found to be the lowest-energy one from CREST run. This kind
# of analysis would be helpful to study the flexibility of the interested intermediate/TS,
# as well as identifying the conformers with significant/slight strcutural changes.
#############################################

def read_xyz(file_path):
    """
    Reads an XYZ file and extracts atomic coordinates for each conformer.
    :param file_path: Path to the XYZ file.
    :return: List of conformers, where each conformer is a list of coordinates (Nx3 numpy array).
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()

    conformers = []
    i = 0
    while i < len(lines):
        num_atoms = int(lines[i].strip())
        i += 2  # Skip the atom count and comment line

        coords = []
        for _ in range(num_atoms):
            parts = lines[i].strip().split()
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            coords.append([x, y, z])
            i += 1

        conformers.append(np.array(coords))

    return conformers

def adjust_relative_to_first_atom(coords):
    """
    Adjusts atomic coordinates relative to the first atom in the molecule.
    :param coords: Numpy array of shape (N, 3) representing atomic coordinates.
    :return: Adjusted numpy array of coordinates.
    """
    return coords - coords[0]  # Subtract the coordinates of the first atom

def calculate_rmsd(coords1, coords2):
    """
    Calculates the RMSD between two sets of atomic coordinates.
    :param coords1: Numpy array of shape (N, 3) for the first conformer.
    :param coords2: Numpy array of shape (N, 3) for the second conformer.
    :return: RMSD value.
    """
    diff = coords1 - coords2
    return np.sqrt(np.sum(diff ** 2) / len(coords1))

def write_rmsd_to_file(rmsd_values, output_file):
    """
    Writes the RMSD values to a file and classifies conformers into categories.
    :param rmsd_values: List of RMSD values.
    :param output_file: Path to the output text file.
    """
    
    # manually designed ranges, subject to changes based on needs
    categories = {
        'smaller_than_0.5': [],
        '0.5_to_1.0': [],
        '1.0_to_2.0': [],
        'above_2.0': []
    }
    
    # if the above ranges change, the corresponding values below should also be modified
    for i, rmsd in enumerate(rmsd_values):
        if rmsd < 0.5:
            categories['smaller_than_0.5'].append(i + 1)
        elif 0.5 <= rmsd < 1.0:
            categories['0.5_to_1.0'].append(i + 1)
        elif 1.0 <= rmsd < 2.0:
            categories['1.0_to_2.0'].append(i + 1)
        else:
            categories['above_2.0'].append(i + 1)

    with open(output_file, 'w') as f:
        f.write("RMSD Values and Classification:\n")
        for i, rmsd in enumerate(rmsd_values):
            f.write(f"Conformer {i+1}: RMSD = {rmsd:.4f} \u00c5\n")
        
        f.write("\nClassification:\n")
        for category, conformers in categories.items():
            f.write(f"{category.replace('_', ' ').title()}: {', '.join(map(str, conformers))}\n")

def main(file_path, output_file):
    # Step 1: Read all conformers from the XYZ file
    conformers = read_xyz(file_path)
    
    # Step 2: Adjust all conformers relative to their first atom
    adjusted_conformers = [adjust_relative_to_first_atom(c) for c in conformers]

    # Step 3: Calculate RMSD relative to the first conformer
    reference_conformer = adjusted_conformers[0]
    rmsd_values = []
    for i, conformer in enumerate(adjusted_conformers):
        rmsd = calculate_rmsd(reference_conformer, conformer)
        rmsd_values.append(rmsd)

    # Step 4: Write RMSD values and classification to file
    write_rmsd_to_file(rmsd_values, output_file)
    print(f"RMSD values and classifications written to {output_file}")

if __name__ == "__main__":
    file_path = "crest_conformers.xyz"  # Replace with your actual file path
    output_file = "rmsd_output.txt"  # Replace with desired output file path
    main(file_path, output_file)

