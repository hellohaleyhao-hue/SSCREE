# SSCREE

**Steric Strain Calculation and Repulsive Energy Estimation (SSCREE)**

SSCREE is a PyMOL-integrated Python tool for evaluating local steric overlap caused by an amino-acid substitution. It enumerates available mutant side-chain rotamers through PyMOL's Mutagenesis Wizard, calculates a geometric steric-strain score for each rotamer, separates the score by the type of surrounding atom, and returns the results as structured Python data and a formatted PyMOL-console table.

The repository also contains a batch-analysis and visualization module that applies SSCREE to a list of variants and generates comparative categorical beeswarm plots.

> **Scientific scope:** SSCREE produces a geometric steric-overlap score. The score is useful as one structural evidence feature when comparing variants, including variants of uncertain significance (VUS), but it is not itself a clinical pathogenicity classifier or a calibrated physical energy calculation.

---

## Table of Contents

- [SSCREE](#sscree)
  - [Table of Contents](#table-of-contents)
  - [Project Overview](#project-overview)
    - [Motivation](#motivation)
    - [What the project contains](#what-the-project-contains)
  - [How SSCREE Works](#how-sscree-works)
- [Repository Structure](#repository-structure)
    - [File descriptions](#file-descriptions)
- [Requirements](#requirements)
  - [Core SSCREE](#core-sscree)
    - [PyMOL](#pymol)
  - [Plotting module](#plotting-module)
- [Installation and Setup](#installation-and-setup)
  - [1. Clone the repository](#1-clone-the-repository)
  - [2. Put the scripts in the expected locations](#2-put-the-scripts-in-the-expected-locations)
  - [3. Make the plotting dependencies available](#3-make-the-plotting-dependencies-available)
  - [4. Launch PyMOL](#4-launch-pymol)
- [Protein Structure Requirements](#protein-structure-requirements)
    - [Chain selection](#chain-selection)
    - [Residue numbering](#residue-numbering)
- [Variant Format](#variant-format)
    - [Chain specification in batch files](#chain-specification-in-batch-files)
- [PyMOL Command Reference](#pymol-command-reference)
  - [`config`](#config)
  - [`calculate_strain`](#calculate_strain)
  - [`set_h`](#set_h)
  - [`help`](#help)
- [Python API](#python-api)
    - [Configure a structure](#configure-a-structure)
    - [Calculate one variant](#calculate-one-variant)
- [Basic Usage](#basic-usage)
  - [Option 1: Use the PyMOL command line](#option-1-use-the-pymol-command-line)
  - [Option 2: Use `run_sscree.py`](#option-2-use-run_sscreepy)
- [Batch Analysis and Plotting](#batch-analysis-and-plotting)
- [Variant Input File](#variant-input-file)
- [Running the Plotting Pipeline](#running-the-plotting-pipeline)
- [Selecting the Strain Category](#selecting-the-strain-category)
- [Plot Interpretation and Configuration](#plot-interpretation-and-configuration)
  - [Point size](#point-size)
  - [Frequency handling](#frequency-handling)
  - [Y-axis](#y-axis)
- [Output Data Format](#output-data-format)
    - [Column definitions](#column-definitions)
- [Steric-Strain Calculation Method](#steric-strain-calculation-method)
  - [Which atoms are evaluated?](#which-atoms-are-evaluated)
- [Parameters](#parameters)
    - [Hydrogen atoms](#hydrogen-atoms)
    - [Search distance](#search-distance)
    - [Buffer](#buffer)
- [TTR Example Dataset](#ttr-example-dataset)
- [Troubleshooting](#troubleshooting)
  - [`ModuleNotFoundError: No module named 'lib'`](#modulenotfounderror-no-module-named-lib)
  - [`No PDB file pathway configured`](#no-pdb-file-pathway-configured)
  - [`Invalid format for mutation`](#invalid-format-for-mutation)
  - [`Invalid amino acid code`](#invalid-amino-acid-code)
  - [`Residue position ... not found`](#residue-position--not-found)
  - [`Sequence mismatch`](#sequence-mismatch)
  - [`Ambiguous selection`](#ambiguous-selection)
  - [Rotamer frequency is `N/A`](#rotamer-frequency-is-na)
  - [The plot was not generated](#the-plot-was-not-generated)
  - [Plotting imports fail](#plotting-imports-fail)
- [Reproducibility](#reproducibility)
- [Limitations and Interpretation](#limitations-and-interpretation)
  - [What a higher score means](#what-a-higher-score-means)
  - [What the score does not directly provide](#what-the-score-does-not-directly-provide)
  - [Dependence on structure and settings](#dependence-on-structure-and-settings)
- [Development Notes](#development-notes)
  - [Core class](#core-class)
  - [Calculation reset behavior](#calculation-reset-behavior)
  - [Rotamer workflow](#rotamer-workflow)
- [Batch Plotter Architecture](#batch-plotter-architecture)
  - [`parse_variant_line()`](#parse_variant_line)
  - [`get_strain_data()`](#get_strain_data)
  - [`plot_strain_data()`](#plot_strain_data)
- [Example End-to-End Workflow](#example-end-to-end-workflow)
    - [Single-variant test](#single-variant-test)
    - [Batch pipeline](#batch-pipeline)
- [Extending SSCREE](#extending-sscree)
- [License and Citation](#license-and-citation)
  - [Summary](#summary)


---

## Project Overview

### Motivation

Missense variants can alter protein structure and function through several mechanisms, including changes in packing, steric compatibility, bonding, conformational preferences, solvent exposure, and interactions with other residues or molecules.

A **steric clash** occurs when atoms are geometrically closer than permitted by the sum of their Van der Waals radii, after applying a user-defined tolerance. SSCREE quantifies this local geometric overlap around a mutated residue.

For a given variant, SSCREE evaluates every rotamer supplied by PyMOL's Mutagenesis Wizard rather than selecting one conformation manually. This preserves the rotamer-frequency information reported by the wizard and allows the structural strain associated with alternative side-chain conformations to be compared.

### What the project contains

The repository has two main computational layers:

1. **SSCREE core (`lib/sscree.py`)**
   - Loads a user-specified structure.
   - Validates a variant.
   - Uses PyMOL's Mutagenesis Wizard to enumerate mutant rotamers.
   - Calculates steric overlap for each rotamer.
   - Separates strain into molecular categories.
   - Prints and returns the results.

2. **Strain plotting module (`lib/strain_plotter.py`)**
   - Reads a variant list.
   - Supports a chain specified on each variant line.
   - Runs SSCREE across the variants.
   - Extracts one selected strain category.
   - Generates two-panel categorical beeswarm plots.
   - Scales point size according to rotamer statistical frequency.
   - Marks missing (`N/A`) frequencies separately.
   - Saves a high-resolution PNG.

---

## How SSCREE Works

The core workflow is:

```text
Input protein structure
        │
        ▼
Choose variant + chain
        │
        ▼
Validate mutation format
        │
        ├── Validate amino-acid codes
        ├── Locate residue
        └── Verify WT residue matches structure
        │
        ▼
Open PyMOL Mutagenesis Wizard
        │
        ▼
Enumerate available mutant rotamers
        │
        ▼
For each rotamer:
        │
        ├── Apply rotamer
        ├── Read PyMOL rotamer frequency
        ├── Identify mutated side-chain atoms
        ├── Find surrounding atoms within cutoff
        ├── Calculate Van der Waals overlap
        ├── Square the overlap
        └── Classify surrounding atoms
        │
        ▼
Return strain values
        │
        ▼
Optional batch plotting
```

The program therefore measures **local geometric strain around the substituted side chain**, rather than performing molecular dynamics, quantum chemistry, free-energy calculations, or force-field minimization.

---

# Repository Structure

The scripts are designed to be organized approximately as follows:

```text
SSCREE/
│
├── README.md
│
├── lib/
│   ├── sscree.py
│   └── strain_plotter.py
│
├── data/
│   ├── TTR_alphafold_homotetramer.pdb
│   └── variant_list.txt
│
├── outputs/
│   └── beeswarm_plot.png
│
├── run_sscree.py
└── run_plotter.py
```

### File descriptions

| File | Purpose |
|---|---|
| `lib/sscree.py` | Core SSCREE calculation engine and PyMOL command registration |
| `lib/strain_plotter.py` | Batch analysis, data extraction, and beeswarm visualization |
| `run_sscree.py` | Example Python runner for one variant |
| `run_plotter.py` | Example runner for the batch plotting pipeline |
| `data/variant_list.txt` | Variant + chain input file |
| `data/TTR_alphafold_homotetramer.pdb` | Example TTR protein structure |
| `outputs/` | Directory for generated figures |

The runner scripts use imports of the form:

```python
import lib.sscree as sc
```

and

```python
import lib.strain_plotter as sp
```

Therefore, run commands from the **repository root** so that the `lib` package can be resolved correctly.

---

# Requirements

## Core SSCREE

### PyMOL

SSCREE was developed and tested with:

- **PyMOL 3.1.8**
- Schrödinger, LLC

The original project documentation specifies compatibility with **PyMOL 2.x and 3.x**, with PyMOL 3.1.8 being the tested version.

The core module uses:

- `pymol.cmd`
- `pymol.stored`
- PyMOL's Mutagenesis Wizard
- PyMOL atom properties such as Van der Waals radius

No third-party Python packages are required by the **core `sscree.py` module** beyond the Python environment supplied by PyMOL.

The Python standard-library modules used by the core are:

```text
dataclasses
math
re
time
```

## Plotting module

`lib/strain_plotter.py` additionally imports:

```text
matplotlib
pandas
seaborn
```

Therefore, these packages must be available in the Python environment in which PyMOL is running the plotting code.

The plotting module also imports:

```python
import lib.sscree as sc
```

so it still requires PyMOL even though it is a Python module.

---

# Installation and Setup

## 1. Clone the repository

```bash
git clone <YOUR-REPOSITORY-URL>
cd SSCREE
```

## 2. Put the scripts in the expected locations

Make sure the repository contains:

```text
lib/sscree.py
lib/strain_plotter.py
run_sscree.py
run_plotter.py
```

and that the required input files are under:

```text
data/
```

Create the output directory before running the example plotting script:

```text
outputs/
```

The current `run_plotter.py` writes directly to:

```text
outputs/beeswarm_plot.png
```

and the plotting module does not create the output directory automatically.

## 3. Make the plotting dependencies available

The core script does not require Matplotlib, pandas, or seaborn. The plotting pipeline does.

Install the plotting dependencies into the Python environment used by PyMOL:

```bash
python -m pip install matplotlib pandas seaborn
```

The exact installation command may depend on how PyMOL was installed. The important requirement is that these packages are importable from the same Python environment that executes the PyMOL-integrated code.

## 4. Launch PyMOL

SSCREE is designed to run inside PyMOL because the calculation engine directly calls the PyMOL API and Mutagenesis Wizard.

---

# Protein Structure Requirements

A structure file must be supplied through the SSCREE `config` command or Python API.

The project documentation and plotting module support protein structures supplied as:

- `.pdb`
- `.mmcif`
- PyMOL session files such as `.pse`

The structure should contain all molecular components that are intentionally meant to participate in the steric analysis.

Examples include:

- protein chains
- ligands
- cofactors
- ions
- metals
- solvent

At the same time, **extraneous objects should be removed**. SSCREE evaluates surrounding atoms from the loaded PyMOL state, so an unrelated object can contribute to the calculated score if it is present in the calculation environment.

### Chain selection

The default target is:

```text
chain A
```

A different chain can be specified for an individual mutation.

For example:

```text
calculate_strain Val30Met, chain C
```

The chain string is used directly as part of a PyMOL selection.

### Residue numbering

SSCREE does not perform sequence-number conversion or automatic mapping between numbering schemes.

The position in:

```text
Val30Met
```

is interpreted as the PyMOL residue identifier:

```text
resi 30
```

within the selected chain.

Therefore, the residue numbering in the input mutation list must match the numbering in the structure file.

This is especially important when working with proteins that have:

- signal peptides
- propeptides
- truncated constructs
- crystallographic numbering
- mature-protein numbering
- engineered residues
- missing residues

---

# Variant Format

SSCREE expects a standard three-letter amino-acid mutation:

```text
<WT><Position><Mutant>
```

Examples:

```text
Val30Met
Thr5Ile
Ser23Arg
Leu55Pro
```

The implementation accepts alphabetic amino-acid codes case-insensitively, but the recommended format is the conventional three-letter capitalization shown above.

Both the wild-type and mutant residue must be one of the standard 20 amino acids:

```text
ALA ARG ASP ASN CYS GLU GLN GLY HIS ILE
LEU LYS MET PHE PRO SER THR TRP TYR VAL
```

### Chain specification in batch files

The batch parser supports all of the following forms:

```text
Thr5Ile
Thr5Ile A
Thr5Ile, A
Thr5Ile, chain A
```

When no chain is supplied, the plotting pipeline uses its configured default chain, which is:

```text
chain A
```

The supplied TTR variant list uses:

```text
Thr5Ile     A
Gly6Ser     A
Cys10Arg    A
...
```

---

# PyMOL Command Reference

When `lib/sscree.py` is loaded, it registers four PyMOL commands.

## `config`

Sets the protein structure file used for subsequent calculations.

```pymol
config <pathway>
```

Example:

```pymol
config data/TTR_alphafold_homotetramer.pdb
```

The structure path must be configured before calling `calculate_strain`.

---

## `calculate_strain`

Calculates steric strain for every available rotamer of a mutation.

```pymol
calculate_strain <mutation>, [chain]
```

Examples:

```pymol
calculate_strain Val30Met
```

```pymol
calculate_strain Val30Met, chain C
```

```pymol
calculate_strain Leu55Pro, chain A
```

The default chain is:

```text
chain A
```

---

## `set_h`

Turns hydrogen-atom inclusion on or off.

```pymol
set_h on
```

or:

```pymol
set_h off
```

Default:

```text
off
```

When hydrogen inclusion is enabled, SSCREE calls:

```python
cmd.h_add()
```

after resetting and loading the structure.

When hydrogen inclusion is disabled, hydrogen atoms are excluded through the atom-selection filter used in the calculation.

---

## `help`

Displays the built-in SSCREE command documentation directly in the PyMOL console.

```pymol
help
```

---

# Python API

The core module creates a global calculator instance:

```python
model = StrainCalculator()
```

The PyMOL commands are registered using:

```python
cmd.extend("help", model.help)
cmd.extend("calculate_strain", model.calculate_strain_value)
cmd.extend("set_h", model.set_include_hydrogen_atoms)
cmd.extend("config", model.config)
```

Therefore, the Python equivalents are:

```python
sc.model.config(...)
```

and:

```python
sc.model.calculate_strain_value(...)
```

### Configure a structure

```python
import lib.sscree as sc

sc.model.config("data/TTR_alphafold_homotetramer.pdb")
```

### Calculate one variant

```python
results = sc.model.calculate_strain_value(
    "Val30Met",
    chain="chain A"
)
```

The function returns the complete strain dataset for all available mutant rotamers.

---

# Basic Usage

## Option 1: Use the PyMOL command line

Load the core script:

```pymol
run lib/sscree.py
```

Configure the structure:

```pymol
config data/TTR_alphafold_homotetramer.pdb
```

Calculate a variant:

```pymol
calculate_strain Val30Met, chain A
```

Change hydrogen handling when needed:

```pymol
set_h on
```

or:

```pymol
set_h off
```

---

## Option 2: Use `run_sscree.py`

The supplied runner contains:

```python
import lib.sscree as sc

sc.model.config("data/TTR_alphafold_homotetramer.pdb")
sc.model.calculate_strain_value("Val30Met", chain="chain A")
```

This provides a minimal Python-level example for reproducing one calculation.

The key point is that `calculate_strain_value` is called through:

```python
sc.model
```

rather than as a top-level module function.

---

# Batch Analysis and Plotting

The plotting layer is implemented in:

```text
lib/strain_plotter.py
```

It provides a higher-level pipeline:

```python
sp.run_pipeline(...)
```

The pipeline:

1. Configures SSCREE with the selected structure.
2. Reads the variant list.
3. Splits the list into the first 15 and second 15 entries.
4. Calculates SSCREE values for every variant.
5. Extracts one strain category.
6. Builds two stacked beeswarm panels sharing the y-axis.
7. Scales point size according to rotamer frequency.
8. Saves the figure at 300 DPI.

---

# Variant Input File

The default example uses:

```text
data/variant_list.txt
```

Each non-empty line contains a variant and a chain.

The supplied 30-entry TTR dataset is:

```text
Thr5Ile     A
Gly6Ser     A
Cys10Arg    A
Cys10Tyr    A
Pro11Leu    A
Leu12Val    A
Asp18Asn    A
Asp18Gly    A
Asp18Tyr    A
Ala19Asp    A
Ala19Gly    A
Ser23Arg    A
Pro24Leu    A
Val30Met    A
His31Asn    A
Arg34Ser    A
Lys35Glu    A
Ala37Asp    A
Lys48Thr    A
Ser50Asn    A
Glu61Ala    A
Glu66Gln    A
Tyr78His    A
Gly83Asp    A
Ala97Gly    A
Pro102Arg   A
Tyr105His   A
Ala120Ser   A
Asn124Ser   A
Lys126Arg   A
```

The plotting code specifically uses:

```python
variant_list[:15]
```

for panel 1 and:

```python
variant_list[15:30]
```

for panel 2.

Thus, the current plotting design is built around a 30-variant dataset divided into two groups of 15.

---

# Running the Plotting Pipeline

The supplied `run_plotter.py` contains:

```python
import lib.strain_plotter as sp

sp.run_pipeline(
    variant_file_path="data/variant_list.txt",
    structure_path="data/TTR_alphafold_homotetramer.pdb",
    output_path="outputs/beeswarm_plot.png",
    strain_category="protein",
    default_chain="chain A"
)
```

This can be executed through the PyMOL Python environment.

The default example plots:

```text
protein
```

strain.

---

# Selecting the Strain Category

`run_pipeline()` accepts the following values for `strain_category`:

| Key | Description | Result column |
|---|---|---:|
| `total` | Total strain | 2 |
| `protein` | Protein strain | 3 |
| `nucleic` | Nucleic-acid strain | 4 |
| `organic` | Organic compound strain | 5 |
| `inorganic` | Inorganic compound strain | 6 |
| `solvent` | Solvent strain | 7 |
| `other` | Other-object strain | 8 |

Example:

```python
sp.run_pipeline(
    variant_file_path="data/variant_list.txt",
    structure_path="data/TTR_alphafold_homotetramer.pdb",
    output_path="outputs/protein_strain.png",
    strain_category="protein",
    default_chain="chain A"
)
```

For total strain:

```python
sp.run_pipeline(
    variant_file_path="data/variant_list.txt",
    structure_path="data/TTR_alphafold_homotetramer.pdb",
    output_path="outputs/total_strain.png",
    strain_category="total",
    default_chain="chain A"
)
```

---

# Plot Interpretation and Configuration

The plotter creates two vertically stacked panels with a shared y-axis:

```text
Panel 1: variants 1–15
Panel 2: variants 16–30
```

Each point corresponds to one mutant rotamer.

## Point size

Point size is determined from the rotamer statistical frequency reported by PyMOL.

The implementation maps:

```text
frequency = 0%  → size 20
frequency = 100% → size 500
```

using linear interpolation.

A rotamer whose frequency is reported as `N/A` receives a fixed marker size of:

```text
80
```

The legend shows representative frequencies of:

```text
5%
20%
40%
60%
```

plus:

```text
N/A
```

## Frequency handling

A frequency string such as:

```text
16.2%
```

is converted to a numerical fraction:

```text
0.162
```

for marker-size scaling.

If a frequency is `N/A`, it is treated as missing for numerical scaling and displayed using the dedicated `N/A` marker convention.

## Y-axis

The strain axis uses a symmetric logarithmic (`symlog`) scale so that values near zero can be displayed while retaining logarithmic behavior over larger values.

The current implementation uses:

```text
linthresh = 0.001
linscale = 0.3
```

and explicitly displays major ticks at:

```text
0
0.001
0.01
0.1
1
10
100
```

The current plotting range is:

```text
-0.0008 to 150
```

Because SSCREE strain is based on squared distances/overlaps, the plotted strain unit is represented as:

```text
Å²
```

---

# Output Data Format

`calculate_strain_value()` returns one row for each mutant rotamer.

Each row has the form:

```python
[
    rotamer_number,
    statistical_frequency,
    total_strain,
    protein_strain,
    nucleic_strain,
    organic_strain,
    inorganic_strain,
    solvent_strain,
    other_strain
]
```

For example, conceptually:

```python
[
    1,
    "16.2%",
    total,
    protein,
    nucleic,
    organic,
    inorganic,
    solvent,
    other
]
```

The actual numerical values are calculated from the structure and are not hard-coded.

### Column definitions

| Index | Name | Meaning |
|---:|---|---|
| `0` | Rotamer | PyMOL rotamer number |
| `1` | Frequency | Statistical frequency reported by the Mutagenesis Wizard |
| `2` | Total | Sum across all surrounding-object categories |
| `3` | Protein | Contribution from protein atoms |
| `4` | Nucleic | Contribution from nucleic-acid atoms |
| `5` | Organic | Contribution from organic compounds such as ligands/cofactors |
| `6` | Inorganic | Contribution from inorganic compounds such as ions/metals |
| `7` | Solvent | Contribution from solvent/water |
| `8` | Other | Contribution from everything not captured by the preceding categories |

The PyMOL console also prints these values in a formatted table.

---

# Steric-Strain Calculation Method

SSCREE uses a geometric overlap model based on the Van der Waals radii assigned to atoms in PyMOL.

For each relevant pair of atoms:

- `i` = atom in the mutant side chain
- `j` = surrounding atom

the interatomic distance is:

```text
dᵢⱼ = sqrt((xᵢ-xⱼ)² + (yᵢ-yⱼ)² + (zᵢ-zⱼ)²)
```

The raw geometric overlap is:

```text
Overlapᵢⱼ = max(0, vdwᵢ + vdwⱼ - dᵢⱼ - buffer)
```

where:

- `vdwᵢ` is the Van der Waals radius stored for atom `i`
- `vdwⱼ` is the Van der Waals radius stored for atom `j`
- `dᵢⱼ` is the measured distance
- `buffer` is the configured tolerance

Only positive overlap contributes to strain.

The contribution of a pair is then:

```text
Strainᵢⱼ = Overlapᵢⱼ²
```

The total score is the sum of all pair contributions:

```text
Total Strain = Σ (Overlapᵢⱼ²)
```

The category-specific scores are summed independently and then combined:

```text
Total
= Protein
+ Nucleic
+ Organic
+ Inorganic
+ Solvent
+ Other
```

All distances are measured in **Angstroms (Å)**, and the squared-overlap score therefore has units of **Å²**.

## Which atoms are evaluated?

For each rotamer, SSCREE selects atoms from the mutated residue's side chain.

The implementation explicitly excludes:

```text
HA
```

because PyMOL's side-chain selection can include the backbone H-alpha atom.

Hydrogen atoms are additionally excluded by default.

For every mutant side-chain atom, SSCREE searches for surrounding atoms within the configured search distance and excludes atoms belonging to the same residue.

---

# Parameters

The calculation parameters are stored in:

```python
StrainParameters
```

The current defaults are:

| Parameter | Default | Meaning |
|---|---:|---|
| `include_hydrogen_atoms` | `False` | Whether hydrogen atoms participate |
| `search_distance` | `5.0 Å` | Radius used to search for surrounding atoms |
| `buffer` | `0.4 Å` | Tolerance subtracted from the VDW-overlap calculation |

### Hydrogen atoms

Default:

```python
include_hydrogen_atoms = False
```

Hydrogens can be enabled through:

```pymol
set_h on
```

and disabled with:

```pymol
set_h off
```

### Search distance

Default:

```python
search_distance = 5.0
```

This defines the maximum distance for identifying surrounding atoms that could contribute to local strain.

### Buffer

Default:

```python
buffer = 0.4
```

This provides a tolerance for small positional deviations and structural flexibility.

Changing this parameter changes the numerical strain scale and therefore should be kept consistent when comparing variants.

---

# TTR Example Dataset

The repository's example structure is:

```text
data/TTR_alphafold_homotetramer.pdb
```

The supplied model contains:

- **4 chains:** A, B, C, and D
- **127 residues per chain**
- residue numbering from **1 through 127** in each chain

The supplied variant list contains **30 variants**, all targeting:

```text
chain A
```

The variants are distributed throughout the TTR sequence and include examples such as:

```text
Thr5Ile
Cys10Arg
Asp18Gly
Ser23Arg
Val30Met
Arg34Ser
Glu66Gln
Tyr78His
Gly83Asp
Pro102Arg
Tyr105His
Ala120Ser
Lys126Arg
```

The included model and variant file are therefore already configured to work with the example runner scripts.

---

# Troubleshooting

## `ModuleNotFoundError: No module named 'lib'`

Make sure you are running from the repository root.

The expected layout is:

```text
SSCREE/
├── lib/
├── data/
├── outputs/
├── run_sscree.py
└── run_plotter.py
```

The runner imports:

```python
import lib.sscree as sc
```

so the repository root must be on Python's import path.

---

## `No PDB file pathway configured`

The calculation must be configured first.

Run:

```pymol
config data/TTR_alphafold_homotetramer.pdb
```

before:

```pymol
calculate_strain Val30Met
```

The Python equivalent is:

```python
sc.model.config("data/TTR_alphafold_homotetramer.pdb")
```

---

## `Invalid format for mutation`

SSCREE expects:

```text
<WT><Position><Mutant>
```

For example:

```text
Val30Met
```

Incorrect formats include strings that do not contain a three-letter wild-type code, a numeric position, and a three-letter mutant code.

---

## `Invalid amino acid code`

Both residues must use standard three-letter amino-acid codes.

Examples:

```text
Val
Met
Thr
Asp
Lys
```

---

## `Residue position ... not found`

The specified residue position does not exist in the selected chain.

Check:

- the chain identifier
- the residue number
- whether the structure uses a different numbering scheme
- whether the structure contains missing residues

---

## `Sequence mismatch`

SSCREE verifies that the wild-type residue specified in the mutation string agrees with the residue actually present in the structure.

For:

```text
Val30Met
```

the structure must contain:

```text
VAL
```

at residue 30 of the selected chain.

This check helps prevent accidentally applying a mutation to the wrong position or numbering system.

---

## `Ambiguous selection`

SSCREE raises an error if more than one C-alpha atom is detected at the requested position in the selected chain.

This can occur when the loaded structure contains duplicate or multiple structural states/models that make the residue selection ambiguous.

Inspect the structure in PyMOL and make sure the intended state/model and residue selection are being used.

---

## Rotamer frequency is `N/A`

The script retrieves the frequency through:

```python
cmd.get_title("mutation", rotamer_number)
```

If PyMOL does not provide a usable title for a rotamer, SSCREE stores:

```text
N/A
```

The strain calculation can still proceed; the plotting module gives `N/A` points a fixed marker size.

---

## The plot was not generated

The plotting pipeline intentionally skips figure generation if **every calculated value in the selected strain category is zero across all variants and rotamers**.

The console reports:

```text
NOTICE: ALL CALCULATED ... VALUES ARE ZERO ACROSS ALL VARIANTS AND ROTAMERS.
Skipping plot generation.
```

Also confirm that:

```text
outputs/
```

exists before saving the image.

---

## Plotting imports fail

If `matplotlib`, `pandas`, or `seaborn` cannot be imported, install them in the same Python environment used by PyMOL.

Remember that the plotting module still depends on PyMOL through:

```python
import lib.sscree as sc
```

so a normal standalone Python interpreter is not sufficient unless it has access to the PyMOL Python API.

---

# Reproducibility

For meaningful comparisons between variants, keep the following constant unless the change is intentionally part of the experiment:

- PyMOL version
- protein structure
- residue numbering
- target chain
- hydrogen setting
- search distance
- buffer
- selected strain category
- variant parsing convention

The current SSCREE defaults are:

```text
PyMOL:                 3.1.8 (tested)
Hydrogens:             off
Search distance:       5.0 Å
Buffer:                0.4 Å
Default chain:         chain A
```

When reporting results, it is recommended to record the structure file and calculation settings alongside the resulting strain values.

---

# Limitations and Interpretation

SSCREE is intentionally a **geometric steric-overlap estimator**, not a full physical simulation.

## What a higher score means

Under identical calculation settings, a higher SSCREE score means that the evaluated mutant rotamer contains more or larger geometric overlaps according to the implemented Van der Waals-radius model.

This can provide evidence of increased local steric incompatibility.

## What the score does not directly provide

SSCREE does not directly calculate:

- experimental stability change (`ΔΔG`)
- free energy
- molecular-dynamics trajectories
- quantum-mechanical energy
- folding probability
- aggregation kinetics
- clinical pathogenicity probability

A high strain score should therefore be interpreted as **one structural feature**, not as a standalone diagnosis or pathogenicity determination.

Similarly, a low score does not prove that a variant is functionally neutral.

## Dependence on structure and settings

Because the calculation is based on coordinates and PyMOL atom properties, the numerical result depends on:

- the selected structure
- structural resolution/model quality
- residue numbering
- the selected chain
- the presence or absence of ligands, solvent, ions, and other objects
- hydrogen inclusion
- search distance
- buffer value
- the rotamers generated by PyMOL

For this reason, absolute values from different structural models or parameter settings should not automatically be treated as directly interchangeable.

---

# Development Notes

## Core class

The central class is:

```python
StrainCalculator
```

It stores:

```python
self.params
self.pathway
```

and exposes methods for:

- configuration
- help
- resetting/reloading the structure
- Boolean conversion
- strain calculation
- hydrogen-atom configuration

## Calculation reset behavior

Each `calculate_strain_value()` call begins by resetting PyMOL and reloading the configured structure.

Each rotamer is also followed by a reset before the next rotamer calculation.

This is intentional: each rotamer calculation starts from a freshly loaded structure rather than accumulating mutations or state changes from previous rotamers.

Because of this design, repeated calculations can involve substantial PyMOL initialization overhead.

## Rotamer workflow

For each rotamer:

1. SSCREE launches the Mutagenesis Wizard.
2. The target residue is selected.
3. The mutant amino acid is specified.
4. The requested PyMOL state/rotamer is selected.
5. The statistical frequency is retrieved.
6. The mutation is applied.
7. Mutant side-chain atoms are selected.
8. Surrounding atoms are identified.
9. Pairwise overlap is calculated.
10. Strain is accumulated by molecular category.
11. The results are stored.
12. PyMOL is reset before the next rotamer.

---

# Batch Plotter Architecture

The plotting module contains three major functions.

## `parse_variant_line()`

Converts input text into:

```python
(variant, chain_selection)
```

and normalizes a chain such as:

```text
A
```

into:

```text
chain A
```

## `get_strain_data()`

Loops through the supplied variants and calls:

```python
sc.model.calculate_strain_value(...)
```

It then extracts:

```python
[row[0], row[1], row[col_idx]]
```

for each rotamer.

The result is a dictionary of:

```text
variant → rotamer/frequency/selected-strain data
```

## `plot_strain_data()`

Converts the results to a pandas DataFrame and creates:

- two categorical panels
- swarm-distributed points
- frequency-dependent point sizes
- an `N/A` marker convention
- a shared symlog y-axis
- a frequency legend
- a 300-DPI PNG output

---

# Example End-to-End Workflow

From the repository root:

```text
1. Start PyMOL
2. Load SSCREE
3. Configure the structure
4. Test one mutation
5. Confirm the output
6. Run the batch pipeline
7. Inspect the generated figure
```

### Single-variant test

```pymol
run lib/sscree.py
config data/TTR_alphafold_homotetramer.pdb
set_h off
calculate_strain Val30Met, chain A
```

### Batch pipeline

Use the supplied runner:

```python
import lib.strain_plotter as sp

sp.run_pipeline(
    variant_file_path="data/variant_list.txt",
    structure_path="data/TTR_alphafold_homotetramer.pdb",
    output_path="outputs/beeswarm_plot.png",
    strain_category="protein",
    default_chain="chain A"
)
```

---

# Extending SSCREE

The current implementation is intentionally modular.

Potential areas for future extension include:

- adding additional scoring models
- exporting calculations directly to CSV/JSON
- adding weighted or frequency-averaged variant scores
- supporting alternative structural models
- adding additional interaction categories
- making group sizes configurable instead of fixed at 15 + 15
- exposing plotting parameters as function arguments
- adding automated validation of the input dataset
- adding automated test cases for residue/chain validation

Any extension that changes the numerical strain definition should be documented explicitly so that results remain reproducible.

---

# License and Citation

No formal license or citation metadata is specified in the current source files.

Before publishing or distributing the repository, add the license you intend to use and, if this code is used in a paper, provide an appropriate software citation.

A useful repository citation should identify:

- the project name: **SSCREE**
- the repository URL
- the version or commit
- the author(s)
- the date accessed, when appropriate

If SSCREE is used as part of a research workflow, describe it as a **steric-overlap/strain calculation tool** and report the calculation parameters used.

---

## Summary

SSCREE is a PyMOL-based computational workflow for systematically quantifying local steric overlap caused by missense substitutions.

Its core calculation:

```text
mutant rotamer
      ↓
side-chain / neighborhood geometry
      ↓
Van der Waals overlap
      ↓
squared overlap
      ↓
category-specific strain
      ↓
total strain
```

The companion plotting module extends this into a batch workflow for comparing multiple variants and visualizing their rotamer-specific strain distributions.

The current repository is configured around a 30-variant TTR example dataset and a four-chain, 127-residue-per-chain TTR structural model, while the underlying code is designed to accept other protein structures, variants, and chains as well.
