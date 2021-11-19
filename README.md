![gust logo](misc/gust.svg)
An easy breezy snp-based whole genome phylogenetic pipeline

## Installation
### 0. Pre-requisites
1. Any kind of conda environment (like `anaconda`, `miniconda`, `mamba` (recommended!), `micromamba`)
2. A `git` installation
3. Linux system (macOS probably works too, untested)

### 1. Clone this repositroy
Everything you'll need is right here in this repository, which will also be your project directory. You can clone
this repository for every new project (it's very lightweight)
```bash
git clone https://github.com/pdimens/gust.git
```
Feel free to rename this folder whatever is relevant for your project.

### 2. Initiate the `gust` conda environment 
**note**: replace `conda` with the appropriate environment framework, (`mamba`, `micromamba`, etc.)
Use the provided `gust.yaml` conda configuration file to create a new conda environment
```bash
cd gust   # enter the gust directory
conda env create -f misc/gust.yaml
```
This will create a new environment called `gust`, which can be activated by
```bash
conda activate gust
```

## Usage
1. Put all of the genomes you want included in analysis in the `genomes` folder
2. Specify the name of the genome you want to use as the reference in `config.yml`
    - the genome should be in the `genomes/` folder
    - use **just** the name, not the full path
    - e.g. `"bostauros.fasta"` ✅  vs `"genomes/bostauros.fasta"` ❌
3. Specify any other additional parameters in `config.yml` if you want them (bwa kmer size, fragment size, etc.)
