# Installation

You don't have to build a local environment for each pipeline! All pipelines are executable in a [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/) container with [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) so that you can run the pipelines on your local machine, HPC, or cloud environment.

## Docker

To install Docker, please follow the instructions [here](https://docs.docker.com/install/).

## Singularity

Singualrity is a container platform that is more suitable for HPC environments. To install Singularity, please follow the instructions [here](https://docs.sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps). Note that you need to have root access to install Singularity.

## Snakemake

For executing any snakefiles ending with `.smk`, you need to install [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) and [Singularity](https://sylabs.io/singularity/). Please follow the instructions [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) to install Snakemake.

## Clone the repository

You can pull the GitHub repository containing all the pipelines by running the following command:

``` bash
# Clone the repository
git clone https://github.com/NaotoKubota/SnakeNgs
```
