import lib.sscree as sc

sc.model.config("data/TTR_alphafold_homotetramer.pdb")
sc.model.calculate_strain_value("Val30Met", chain="chain A")
