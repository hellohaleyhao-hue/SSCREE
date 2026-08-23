from pymol import cmd, stored
import math
import time
import re

class StrainCalculator:
    amino_acids = (
        'ALA','ARG','ASP','ASN','CYS','GLU','GLN','GLY','HIS','ILE',
        'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'
    )
    def __init__(self, pathway: str):
        self.pathway = pathway
        self.include_hydrogen_atoms = False
        
        # --- Calculation Parameters ---
        # Computational cutoff for steric clashes (in Angstroms)
        self.search_distance = 5.0

        # Standard distance tolerance is 0.4 Angstroms due to flexible backbones
        self.buffer = 0.4

        print("NOTE: Starting this calculate_strain_value script reinitializes everything.\n")

    def help(self):
        print()
        print("calculate_strain <mutation (e.g. Val30Met)>, <chain (e.g. chain C, default: chain A)>\n"
            "     -> calculate strain value for that amino acid.")
        print("set_h <on/off>\n"
            "     -> include or exclude hydrogen atoms in the steric clash calculations.")
        print("help\n"
            "     -> display this help message.")

    def reset(self):
        cmd.reinitialize()
        cmd.load(self.pathway)

        if self.include_hydrogen_atoms:
            cmd.h_add()

    def bool_conversion(self, value: str) -> bool:
        if value.lower() in ("yes", "true", "on", "1"):
            return True
        elif value.lower() in ("no", "false", "off", "0"):
            return False
        else:
            raise ValueError("Invalid input. Please use 'on' or 'off'.")

    def calculate_strain_value(self, mutation: str, chain: str = "chain A") -> list:
        start_time = time.perf_counter()

        self.reset()
        strain_value_data = []
        h_filter = "" if self.include_hydrogen_atoms else " and not elem H"

        # --- Input Validation ---
        match = re.fullmatch(r"([A-Za-z]{3})(\d+)([A-Za-z]{3})", mutation)
        if not match:
            print(f"[!] Error: Invalid format for mutation '{mutation}'. Please use the format <Wild-Type Residue><Residue Position><Mutant Residue> "
                "(e.g., Val30Met) with residue names being the three letter code, case insensitive")
            return []

        position = match.group(2)
        wild_type_residue = match.group(1).upper()
        mutant_residue = match.group(3).upper()
        residue = f"{chain} and resi {position}"

        if wild_type_residue not in self.amino_acids or mutant_residue not in self.amino_acids:
            print(f"[!] Error: Invalid amino acid in mutation '{mutation}'. Please use three-letter codes for standard amino acids.")
            return []

        stored.residue_name = []
        cmd.iterate(f"{residue} and name CA", "stored.residue_name.append(resn)")
        residue_name = stored.residue_name[0] if stored.residue_name else None

        if len(stored.residue_name) > 1:
            print(f"[!] Error: Something went wrong. Program detected {len(stored.residue_name)} CA atoms at position {position} in {chain}.")
            return []
        if not residue_name:
            print(f"[!] Error: Resi {position} not found in {chain}.")
            return []
        if residue_name != wild_type_residue:
            print(f"[!] Error: Wild-type residue mismatch at position {position} in {chain}. Expected {wild_type_residue}, but found {residue_name}.")
            return []

        # --- Rotamer Validation ---
        cmd.wizard("mutagenesis")
        cmd.get_wizard().do_select(residue)
        cmd.get_wizard().set_mode(mutant_residue)
        num_rotamers = cmd.count_states("mutation")
        cmd.set_wizard()

        if num_rotamers < 1:
            print(f"[!] Error: No rotamers found for mutation '{mutation}'.")
            return []

        for rotamer_number in range(1, num_rotamers + 1):

            cmd.wizard("mutagenesis")
            cmd.get_wizard().do_select(residue)
            cmd.get_wizard().set_mode(mutant_residue)
            cmd.frame(rotamer_number)
            Statistical_frequency = cmd.get_title("mutation", rotamer_number)
            cmd.get_wizard().apply()
            cmd.set_wizard()

            Total_protein_strain = 0
            Total_nucleic_strain = 0
            Total_organic_strain = 0   # (Ligands/Cofactors)
            Total_inorganic_strain = 0 # (Ions/Metals)
            Total_solvent_strain = 0   # (Water)
            Total_other_strain = 0     # (Other compounds)

            # Exclude HA (PyMOL sidechain selection includes backbone H-alpha)
            model_resi = cmd.get_model(f"{residue} and sidechain and not name HA {h_filter}")

            for atom_resi in model_resi.atom:
                x_i, y_i, z_i = atom_resi.coord

                cur_atom = f"model {atom_resi.model} and index {atom_resi.index}"
                model_surr = cmd.get_model(
                    f"({cur_atom} around {self.search_distance}) "
                    f"and not (byres ({cur_atom})) {h_filter}"
                )
                for atom_sur in model_surr.atom:
                    x_j, y_j, z_j = atom_sur.coord

                    Distance = math.sqrt((x_i - x_j) ** 2 + (y_i - y_j) ** 2 + (z_i - z_j) ** 2)
                    Overlap = max(0, atom_resi.vdw + atom_sur.vdw - Distance - self.buffer)
                    
                    # Quadratic strain model for potential repulsive energy estimation
                    Strain = Overlap ** 2
                    
                    target_atom = f"model {atom_sur.model} and index {atom_sur.index}"

                    if cmd.count_atoms(f"{target_atom} and polymer.protein"):
                        Total_protein_strain += Strain
                    elif cmd.count_atoms(f"{target_atom} and polymer.nucleic"):
                        Total_nucleic_strain += Strain
                    elif cmd.count_atoms(f"{target_atom} and organic"):
                        Total_organic_strain += Strain
                    elif cmd.count_atoms(f"{target_atom} and inorganic"):
                        Total_inorganic_strain += Strain
                    elif cmd.count_atoms(f"{target_atom} and solvent"):
                        Total_solvent_strain += Strain
                    else:
                        Total_other_strain += Strain

            Total_strain = (
                Total_protein_strain
                + Total_nucleic_strain
                + Total_organic_strain
                + Total_inorganic_strain
                + Total_solvent_strain
                + Total_other_strain
            )
            strain_value_data.append([
                rotamer_number,
                Statistical_frequency,
                Total_strain,
                Total_protein_strain, 
                Total_nucleic_strain, 
                Total_organic_strain, 
                Total_inorganic_strain, 
                Total_solvent_strain,
                Total_other_strain
            ])
            self.reset()

        cmd.sync()
        cmd.refresh()
        end_time = time.perf_counter()

        print(f"\nStrain results for {mutation} (run_time: {end_time - start_time:.2f}s)")
        print("=" * 111)
        print(
            f"{'Rotamer':<8} | {'Frequency':<10} | {'Total':<10} | "
            f"{'Protein':<10} | {'Nucleic':<10} | {'Organic':<10} | "
            f"{'Inorganic':<10} | {'Solvent':<10} | {'Other':<10}"
        )
        print("-" * 111)
        for row in strain_value_data:
            print(
                f"{row[0]:<8} | {row[1]:<10} | {row[2]:<10.4f} | "
                f"{row[3]:<10.4f} | {row[4]:<10.4f} | {row[5]:<10.4f} | "
                f"{row[6]:<10.4f} | {row[7]:<10.4f} | {row[8]:<10.4f}"
            )
        print("=" * 111)

        return strain_value_data

    def set_include_hydrogen_atoms(self, include: str):
        self.Include_hydrogen_atoms = self.bool_conversion(include)

        if self.Include_hydrogen_atoms:
            print("Hydrogen atoms will be included in the steric clash calculations.")
        else:
            print("Hydrogen atoms will be excluded from the steric clash calculations.")

# --- Register commands in PyMOL ---
model = StrainCalculator("/.../example_protein.pse")

cmd.extend("help", model.help)
cmd.extend("calculate_strain", model.calculate_strain_value)
cmd.extend("set_h", model.set_include_hydrogen_atoms)
