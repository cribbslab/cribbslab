# Botnar Computational Group

This is a collection of ruffus pipelines and R scripts used for the cribbs group
at the botnar research centre Oxford.

## Installation

### Conda environment

I recommend using miniconda to install a conda environment followed by mamba:

    conda install mamba
    mamba env create -f conda/environments/cribbslab.yml

### Manual installation

The repository can also be installed manually, but dependancies will need to be installed seperately::

    python setup.py develop
    cribbslab --help

## Usage

Run ``cribbslab --help`` to list available workflows.

To run a pipeline, create a project directory, generate a configuration file,
and run a target:

```bash
mkdir my_run && cd my_run
cribbslab <workflow> config    # creates pipeline.yml
# edit pipeline.yml (paths, kit, threads, ...)
cribbslab <workflow> make full --local -j 8
```

### Single-cell long-read RNA-seq (`sclong`)

The 10x ONT single-cell pipeline has a dedicated guide covering input
preparation (FASTQ, reference, GTF), configuration, run targets, and optional
fusion/SNV steps:

**[cribbslab/pipeline_sclong/README.md](cribbslab/pipeline_sclong/README.md)**
