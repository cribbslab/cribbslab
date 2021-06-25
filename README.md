# Botnar Computational Group

This is a collection of pipelines and R scripts used for the cribbs group
at the botnar research centre Oxford.

## Installation

### Conda environment

I recommend using miniconda to install a conda environment followed by mamba:

    conda install mamba
    mamba env create -f conda/environments/cgat-core.yml

### Manual installation

The repository can also be installed manually, but dependancies will need to be installed seperately::

    python setup.py develop
    cribbslab --help

## Usage

    Run the ``cribbslab --help`` command view the help documentation for how to run the available workflows.

    To run a pipeline first generate a configuration file::

        cribbslab <workflow> config

    Then run the pipeline::

        cribbslab <workflow> make full -v5
