![gust logo](.other/gust.svg)
An easy breezy snp-based whole genome phylogenetic pipeline 🌪️

## Installation
### Pre-requisites
1. Any kind of conda environment (like `anaconda`, `miniconda`, `mamba` (recommended!), `micromamba`)
2. A `git` installation
3. Linux system (macOS probably works too, untested)

### 1. Clone this repositroy
Everything you'll need is right here in this repository, which will also be your project directory. You'll need to clone
this repository for every new project (it's very lightweight). Feel free to rename this folder `gust` to whatever is relevant for your project.
```bash
git clone https://github.com/pdimens/gust.git
```
### 2. Initiate the `gust` conda environment 
🌪️ **if you already have a `gust` conda environment, skip this step** 🌪️

**note**: replace `conda` with the appropriate environment framework, (`mamba`, `micromamba`, etc.)

Use the provided `gust.yaml` conda configuration file to create a new conda environment
```bash
cd gust   # enter the gust directory
conda env create -f .other/gust_env.yaml
```
This will create a new environment called `gust` containing all the software dependencies. Activate the environment with:
```bash
conda activate gust
```

## Usage
### Preparation
It's minimal, I swear.
1. Put all of the genomes you want included in analysis in the `genomes` folder
    - **make sure the files end with `.fasta`**
    - I haven't figured out how to make snakemake more flexible with this (PR's welcome!)
2. Specify the name of the genome you want to use as the reference in `config.yml`
    - the genome should be in the `genomes/` folder
    - use **just** the name, not the full path
    - e.g. `"bostauros.fasta"` ✅  vs `"genomes/bostauros.fasta"` ❌
3. Specify any other additional parameters in `config.yml` if you want them (bwa kmer size, fragment size, etc.)

### Running
If you call `gust` without arguments you will see help text, otherwise:
```bash
./gust numThreads configFile
```
where `numThreads` is the number of threads to use (>1, no default) and
`configFile` is the name of the configuration file (optional, defaults to `config.yml`).
#### examples
Running `gust` with 25 threads using the configuration file `config.yml`
```bash
./gust 25
```
Running `gust` with 10 threads using the configuration file `frag250.yml`
```bash
./gust 10 frag250.yml
```
