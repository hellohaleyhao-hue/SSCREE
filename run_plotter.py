import lib.strain_plotter as sp

# =================== EXAMPLE USAGE ====================

sp.run_pipeline(
    variant_file_path="data/variant_list.txt",
    structure_path="data/TTR_alphafold_homotetramer.pdb",
    output_path="/Users/apple/Desktop/beeswarm_plot.png",
    strain_category="protein",
    default_chain="chain A"
)

# ======================================================

# ===================== TEMPLATE =======================
# sp.run_pipeline(
#     variant_file_path="data/variant_list.txt",
#     structure_path="data/TTR_alphafold_homotetramer.pdb",
#     output_path=f"outputs/beeswarm_plot.png",
#     strain_category="protein",
#     default_chain="chain A"
# )