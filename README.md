# Descriptors-Formulation
## Overview
A general-purpose Python script that resolves chemical names to SMILES structures and generates configurable molecular descriptor datasets for QSAR/QSPR workflows. Descriptor families are toggleable via JSON configuration, and results are exported to a structured Excel workbook. Built on RDKit with optional Mordred and thermo integration.

## Features
`
- Chemical name to SMILES resolution
- RDKit descriptor generation
- Optional Mordred descriptor calculation
- Optional experimental property enrichment (PubChem PUG View and thermo)
- Configurable descriptor families via JSON
- Parallel processing
- Excel export with summary statistics and descriptors split by descriptor library

## Requirements

The list of dependencies can be found in [requirements.txt](`requirements.txt`).

Mandatory dependencies:

- pandas
- requests
- rdkit

Optional dependencies:

- tqdm
- mordred
- thermo

The requirements for this script can be installed using pip with the command `pip install -r requirements.txt`.

## Basic usage

The script can be ran with the following command in the terminal:

```bash
python "Descriptors Formulation Script.py" --input Input.csv --output features.xlsx
```

Optional flags include:

```
--use-mordred              Enable Mordred descriptor calculation
--clear-cache              Reset SMILES resolution cache
--properties-config        Provide external property configuration via the .JSON file
--no-prefer-input-smiles   Ignore provided SMILES column
--output                   Specify the location of the output file
```

## Input, output, and configuration

The input file can either be an Excel spreadsheet (.xlsx) or a Comma Separated Values file (.csv). Inside the input file, there must be a column indicating the names of the chemicals. Optionally, there is the ability to specify the SMILES format and/or the CAS number if known. See the example `Input.csv` file provided for a clear example of this.

The output file is an Excel spreadsheet containing multiple pages, including a summary sheet, the full feature table, and separate tables containing RDKit, mordred, and external properties respectively. In all cases, they also contain PubChem metadata for reference. 

Configuration of descriptor families and optional data enrichment can be configured by editing the `properties_config.json` file. There, you are able to specify which specific descriptors from RDKit, mordred, PubChem, or thermo to include or exclude at your discretion. It also allows the enabling and disabling of specific descriptor families as well, but be aware this may not work as intended all the time. See the `properties_config.json` file for an example configuration to begin with. To run the script with these modified descriptors, add `--properties-config properties_config.json` to the terminal command.

## Citation

See the [citation.cff](`citation.cff`) file for more information on how to cite this repository and to export it to different formats (APA, BibTeX, etc)

## License

See the [LICENSE](LICENSE) file for details.
