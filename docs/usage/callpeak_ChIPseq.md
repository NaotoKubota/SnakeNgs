# callpeak_ChIPseq.smk

Snakemake workflow for peak calling of ChIP-seq data.

!!! note

    Please make sure that you have [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html) and [Snakemake](https://snakemake.readthedocs.io/en/stable/) installed on your system and cloned the [SnakeNgs](https://github.com/NaotoKubota/SnakeNgs) repository.

## Workflow

<figure markdown="span">
	![callpeak_ChIPseq.smk rulegraph](https://github.com/NaotoKubota/SnakeNgs/blob/develop/img/callpeak_ChIPseq_rulegraph.svg?raw=true){ width="1000" align="center" }
</figure>

<span style="font-size: 0.8em; color: rgba(0, 0, 0, 0.4);">The rulegraph was created by [snakevision](https://github.com/OpenOmics/snakevision).</span>

1. Peak calling using [MACS](https://github.com/macs3-project/MACS) with the parameter specified in the `config.yaml`.
2. Convert bedgraph files to bigwig files using [bedgraphtobigwig](https://genome.ucsc.edu/goldenPath/help/bigWig.html).
3. Make summary statistics using [MultiQC](https://multiqc.info/).

## Usage

``` bash
snakemake -s /path/to/SnakeNgs/snakefile/callpeak_ChIPseq.smk \
--configfile /path/to/config.yaml \
--cores <int> \
--use-singularity \
--rerun-incomplete
```

`config.yaml` should contain the following information:

``` yaml
general:
  workdir: /path/to/output
  experiment_table: /path/to/experiment_table.tsv

macs2:
  broad: False # Boolean
  broad_cutoff: # Float when broad is True
  genomesize: hs # ["hs", "mm", "ce", "dm"]
  fileformat: BAM # ["BED", "BAM", "SAM", "BEDPE", "BAMPE"]
  qvalue: 0.05

bedgraphtobigwig:
  assembly: hg38
```

- `experiment_table.tsv` should contain the following information:

``` text
sample	target	control
sample1	/path/to/sample1_target.bam	/path/to/control_A.bam
sample2	/path/to/sample2_target.bam	/path/to/control_B.bam
sample3	/path/to/sample3_target.bam	/path/to/control_A.bam
```

`target` and `control` are paths to the BAM files for the target and control (input) samples, respectively.

- `macs2` contains the parameters for MACS2 peak calling. Please refer to the [MACS2 documentation](https://macs3-project.github.io/MACS/).

- `bedgraphtobigwig` contains the assembly information for the conversion of bedgraph files to bigwig files.

## Docker image used in the workflow

- [quay.io/biocontainers/macs2:2.2.9.1--py39hf95cd2a_0](https://quay.io/repository/biocontainers/macs2)
- [quay.io/biocontainers/ucsc-bedgraphtobigwig:445--h954228d_0](https://quay.io/repository/biocontainers/ucsc-bedgraphtobigwig)
- [multiqc/multiqc:v1.25](https://hub.docker.com/r/multiqc/multiqc)
