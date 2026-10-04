"""
SSCREE: Steric Strain Calculation and Repulsive Energy Estimation
-----------------------------------------------------------------

A PyMOL-integrated Python package for evaluating steric overlap 
and calculating repulsive strain across amino acid rotamers.
"""

import math
import re
import time
from dataclasses import dataclass
from pymol import cmd, stored


@dataclass
class StrainParameters:
    # --- User-configurable parameters for steric strain calculations ---

    # Default setting for including hydrogen atoms in steric clash calculations
    include_hydrogen_atoms: bool = False

    # Computational cutoff for steric clashes (in Angstroms)
    search_distance: float = 5.0 # 5 Angstroms is a reasonable distance to capture potential steric clashes without excessive computation time

    # Standard distance tolerance is 0.4 Angstroms due to flexible backbones
    buffer: float = 0.4 # 0.4 Angstroms is standard for steric clash calculations to account for minor deviations in atomic positions and flexibility in the protein structure


class StrainCalculator:
    amino_acids = (
        'ALA','ARG','ASP','ASN','CYS','GLU','GLN','GLY','HIS','ILE',
        'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'
    )

    def __init__(self):
        self.params = StrainParameters()
        self.pathway = None

        print("-" * 20)
        print("NOTE: Starting this calculate_strain_value script reinitializes everything.\n")
        print("Type 'help' for instructions on how to use this script.")

    def config(self, pathway: str):
        self.pathway = pathway
        print(f"Protein structure file pathway set to: {self.pathway}")

    def help(self):
        """Displays comprehensive usage instructions for SSCREE commands in PyMOL."""
        help_text = """
===================================================================================================
                          SSCREE: Steric Strain & Repulsive Energy Estimator                       
===================================================================================================

COMMAND OVERVIEW
----------------
  1. config <pathway>
     Sets the target PDB structure file path for calculations. Must be run before calculating strain.
     
     Examples:
       PyMOL> config data/TTR_alphafold_homotetramer.pdb
       PyMOL> config /path/to/structure.pdb

  2. calculate_strain <mutation>, [chain]
     Calculates steric overlap and repulsive strain for every rotamer of the specified mutation.
     
     Parameters:
       - mutation : Standard 3-letter code variant string in <WT><Pos><MUT> format (e.g., Val30Met).
       - chain    : Target chain selection (default: chain A).
     
     Examples:
       PyMOL> calculate_strain Val30Met
       PyMOL> calculate_strain Val30Met, chain C
       PyMOL> calculate_strain Leu55Pro, chain A

  3. set_h <on/off>
     Enables or disables hydrogen atoms in steric strain distance evaluations (default: off).
     
     Examples:
       PyMOL> set_h on
       PyMOL> set_h off

  4. help
     Displays this documentation message.

===================================================================================================
"""
        print(help_text)

    def reset(self):
        if not self.pathway:
            raise ValueError("No PDB file pathway configured. Call the config function first.")

        cmd.reinitialize()
        cmd.load(self.pathway)

        if self.params.include_hydrogen_atoms:
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
        h_filter = "" if self.params.include_hydrogen_atoms else " and not elem H"

        # --- Input Validation ---
        match = re.fullmatch(r"([A-Za-z]{3})(\d+)([A-Za-z]{3})", mutation)
        if not match:
            raise ValueError(
                f"Invalid format for mutation '{mutation}'. "
                "Expected format: <WT><Position><Mutant> (e.g., Val30Met)."
            )
        
        position = match.group(2)
        wild_type_residue = match.group(1).upper()
        mutant_residue = match.group(3).upper()
        residue = f"{chain} and resi {position}"

        if wild_type_residue not in self.amino_acids or mutant_residue not in self.amino_acids:
            raise ValueError(
                f"Invalid amino acid code in '{mutation}'. "
                "Please use standard three-letter amino acid codes."
            )
        
        stored.residue_name = []
        cmd.iterate(f"{residue} and name CA", "stored.residue_name.append(resn)")
        residue_name = stored.residue_name[0] if stored.residue_name else None

        if len(stored.residue_name) > 1:
            raise RuntimeError(
                f"Ambiguous selection: detected {len(stored.residue_name)} C-alpha atoms "
                f"at position {position} in {chain}."
            )
        
        if not residue_name:
            raise KeyError(f"Residue position {position} not found in {chain}.")
        
        if residue_name != wild_type_residue:
            raise ValueError(
                f"Sequence mismatch at position {position} in {chain}: "
                f"expected {wild_type_residue}, but structure contains {residue_name}."
            )
        
        # --- Rotamer Validation ---
        cmd.wizard("mutagenesis")
        cmd.get_wizard().do_select(residue)
        cmd.get_wizard().set_mode(mutant_residue)
        num_rotamers = cmd.count_states("mutation")
        cmd.set_wizard()

        for rotamer_number in range(1, num_rotamers + 1):

            cmd.wizard("mutagenesis")
            cmd.get_wizard().do_select(residue)
            cmd.get_wizard().set_mode(mutant_residue)
            cmd.frame(rotamer_number)
            raw_frequency = cmd.get_title("mutation", rotamer_number)
            Statistical_frequency = raw_frequency.strip() if raw_frequency else "N/A"
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
                    f"({cur_atom} around {self.params.search_distance}) "
                    f"and not (byres ({cur_atom})) {h_filter}"
                )
                for atom_sur in model_surr.atom:
                    x_j, y_j, z_j = atom_sur.coord

                    Distance = math.sqrt((x_i - x_j) ** 2 + (y_i - y_j) ** 2 + (z_i - z_j) ** 2)
                    Overlap = max(0, atom_resi.vdw + atom_sur.vdw - Distance - self.params.buffer)
                    
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
        print(f"\nStrain results for {mutation} on {chain} (run_time: {end_time - start_time:.2f}s)")
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
        self.params.include_hydrogen_atoms = self.bool_conversion(include)

        if self.params.include_hydrogen_atoms:
            print("Hydrogen set to on (included in calculation)")
        else:
            print("Hydrogen set to off (excluded from calculation)")

# --- Register commands in PyMOL ---
model = StrainCalculator()

cmd.extend("help", model.help)
cmd.extend("calculate_strain", model.calculate_strain_value)
cmd.extend("set_h", model.set_include_hydrogen_atoms)
cmd.extend("config", model.config)
