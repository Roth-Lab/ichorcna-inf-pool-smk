# Evaluate ichorCNA on simulated cfDNA samples
A [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow that generates in-silico mixtures and runs [ichorCNA](https://github.com/Roth-Lab/myichorcna).

# Getting started
This pipeline requires that [conda](https://github.com/conda-forge/miniforge) and [Snakemake](https://snakemake.readthedocs.io/en/stable/) be installed; the [Bioconda](https://bioconda.github.io/#usage) package channel must also be configured.

**Dependencies**
* [conda](https://github.com/conda-forge/miniforge), version >24.7.1
* [Snakemake](https://snakemake.readthedocs.io/en/stable/), version >=9.14.8

**Environment**
1. Ensure that you have a working `conda` installation, you can do this by installing [Miniforge](https://github.com/conda-forge/miniforge#install).
2. Configure the [Bioconda channel](https://bioconda.github.io/#usage) and set strict channel priority:
   ```
   conda config --add channels bioconda
   conda config --add channels conda-forge
   conda config --set channel_priority strict
   ```
3. Install [Snakemake](https://snakemake.readthedocs.io/en/stable/):
   ```
   conda create -c conda-forge -c bioconda --name snakemake snakemake'>=9.12'
   ```

**Workflow**

1. Create a working directory for the workflow:
   ```
   mkdir -p path/to/project-workdir
   cd path/to/project-workdir
   ```
2. Clone the workflow repository through git:
      ```
      git clone --depth 1 https://github.com/Roth-Lab/ichorcna-inf-pool-tf-smk.git
      ```
# Usage

**Configuration**

For a full description of all available pipeline options, please refer to the pipeline [schema](schemas/config.schema.yaml). Modify the configuration file, [config.yaml](configs/example.yaml) to suit your dataset.

**Run Workflow**
1. Navigate to the project directory and activate the snakemake environment:
   ```
   cd path/to/project-workdir
   conda activate snakemake
   ```
2. Run a dry-run of the pipeline to confirm the ruleset and outputs are as you expect:
   ```
   snakemake --cores <number-of-CPU-cores-to-use> --configfile <path/to/config-file> -n 
   ```
3. Run the pipeline:
   ```
   snakemake --cores <number-of-CPU-cores-to-use> --configfile <path/to/config-file>
   ```

---------
# Output
The main outputs of the pipeline are posterior distributions on tumour fraction and clone prevalences.
More on the contents of these output files can be found in the [cfClone repository](https://github.com/Roth-Lab/cfClone).

**Example workflow output folder structure:**
```
.
├── <out-dir>
│   ├── config.yaml
│   └── input
│       └── clone_cn
│           └── clone_cn.tsv.gz
└── <pipeline_dir>
    └── log
        └── output
            └── input
                └── clone_cn
                    └── clone_cn.tsv.log

```