# MAJIQ.smk

Snakemake workflow for differential RNA splicing analysis using [MAJIQ](https://majiq.biociphers.org/).

!!! note

    Please make sure that you have [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html) and [Snakemake](https://snakemake.readthedocs.io/en/stable/) installed on your system and cloned the [SnakeNgs](https://github.com/NaotoKubota/SnakeNgs) repository.

## Workflow

<figure markdown="span">
	![MAJIQ.smk rulegraph](https://github.com/NaotoKubota/SnakeNgs/blob/develop/img/MAJIQ_rulegraph.svg?raw=true){ width="1000" align="center" }
</figure>

<span style="font-size: 0.8em; color: rgba(0, 0, 0, 0.4);">The rulegraph was created by [snakevision](https://github.com/OpenOmics/snakevision).</span>

1. Parse the experiment table to create a setting file for MAJIQ.
2. Build splicegraph using [MAJIQ](https://majiq.biociphers.org/) `build`.
3. Differential splicing analysis using [MAJIQ](https://majiq.biociphers.org/) `heterogen` and `deltapsi`.
4. Make a tab-separated summary of the results using [MAJIQ](https://majiq.biociphers.org/) `voila tsv`.
5. Modulize the results using [MAJIQ](https://majiq.biociphers.org/) `voila modulize`.

## Usage

``` bash
snakemake -s /path/to/SnakeNgs/snakefile/MAJIQ.smk \
--configfile /path/to/config.yaml \
--cores <int> \
--use-singularity \
--rerun-incomplete
```

`config.yaml` should contain the following information:

``` yaml
workdir: /path/to/output
container: /path/to/singularity_image
experiment_table: /path/to/experiment_table.tsv
gff: /path/to/reference_transcriptome.gff3
genome: hg19 # hg19, hg38, mm10, mm39, etc.
strandness: None # [None, forward, reverse]
```

- `/path/to/output` is the directory where the output files will be saved.
- `/path/to/singularity_image` is the path to the Singularity image file for MAJIQ.

!!! note

	As of now, there is no public Docker/Singularity image available for MAJIQ. You need to build the image by yourself. Please refer to the [MAJIQ installation guide](https://majiq.biociphers.org/app_download/) for more information.

- `/path/to/experiment_table.tsv` is a tab-separated file, which is same as the one used in [Shiba](https://github.com/NaotoKubota/Shiba).

``` text
sample	bam	group
sample1	/path/to/sample1.bam	Ref
sample2	/path/to/sample2.bam	Ref
sample3	/path/to/sample3.bam	Ref
sample4	/path/to/sample4.bam	Alt
sample5	/path/to/sample5.bam	Alt
sample6	/path/to/sample6.bam	Alt
```

The `group` column must be specified as `Ref` or `Alt`. This workflow will perform the differential splicing analysis between the two groups.

- `/path/to/reference_transcriptome.gff3` is the reference transcriptome in GFF3 format (e.g. `Homo_sapiens.GRCh38.106.gff3` for human transcriptome).
- `genome` is the genome version used in the analysis (e.g. `hg19`, `hg38`, `mm10`, `mm39`, etc.).
- `strandness` is the strandness of the RNA-seq data. (`None`, `forward`, `reverse`)

## Docker image used in the workflow

No public Docker image is available for MAJIQ.
