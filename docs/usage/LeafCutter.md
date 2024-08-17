# LeafCutter.smk

Snakemake workflow for differential RNA splicing analysis using [LeafCutter](https://davidaknowles.github.io/leafcutter/).

!!! note

    Please make sure that you have [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html) and [Snakemake](https://snakemake.readthedocs.io/en/stable/) installed on your system and cloned the [SnakeNgs](https://github.com/NaotoKubota/SnakeNgs) repository.

## Workflow

<figure markdown="span">
	![LeafCutter.smk rulegraph](https://github.com/NaotoKubota/SnakeNgs/blob/develop/img/LeafCutter_rulegraph.svg?raw=true){ width="1000" align="center" }
</figure>

<span style="font-size: 0.8em; color: rgba(0, 0, 0, 0.4);">The rulegraph was created by [snakevision](https://github.com/OpenOmics/snakevision).</span>

1. Extract junction reads using [regtools](https://regtools.readthedocs.io/en/latest/) `junctions extract`.
2. Cluster introns using [LeafCutter](https://davidaknowles.github.io/leafcutter/) `leafcutter_cluster_regtools.py`.
3. Extract exon information using [LeafCutter](https://davidaknowles.github.io/leafcutter/) `gtf_to_exons.R`.
4. Differential splicing analysis using [LeafCutter](https://davidaknowles.github.io/leafcutter/) `leafcutter_ds.R`.
5. Plot splice junctions using [LeafCutter](https://davidaknowles.github.io/leafcutter/) `ds_plots.R`.
6. Make annotation codes using [LeafCutter](https://davidaknowles.github.io/leafcutter/) `gtf2leafcutter.pl`.
7. Prepare results for visualization in LeafViz using [LeafCutter](https://davidaknowles.github.io/leafcutter/) `prepare_results.R`.
8. Classify clusters using [LeafCutter](https://davidaknowles.github.io/leafcutter/) `classify_clusters.R`.

## Usage

``` bash
snakemake -s /path/to/SnakeNgs/snakefile/LeafCutter.smk \
--configfile <path to config.yaml> \
--cores <int> \
--use-singularity \
--rerun-incomplete
```

`config.yaml` should contain the following information:

``` yaml
workdir: /path/to/output
experiment_table: /path/to/experiment_table.tsv
gtf: /path/to/reference_transcriptome.gtf

# Junction read filtering
minimum_anchor_length: 6 # Minimum anchor length
minimum_intron_length: 70 # Minimum intron length
maximum_intron_length: 500000 # Maximum intron length
strand: XS # Strand information in BAM file

# Intron clustering
minimum_reads: 10

# Differential splicing analysis
min_coverage: 10 # Minimum coverage
min_samples_per_intron: 1 # Minimum samples per intron
min_samples_per_group: 1 # Minimum samples per group
FDR: 0.05 # False discovery rate
```

- `/path/to/output` is the directory where the output files will be saved.
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

- `/path/to/reference_transcriptome.gtf` is the reference transcriptome in GTF format (e.g. `Homo_sapiens.GRCh38.106.gtf` for human transcriptome).

- `minimum_anchor_length` is the minimum anchor length for the junction reads.
- `minimum_intron_length` is the minimum intron length for the junction reads.
- `maximum_intron_length` is the maximum intron length for the junction reads.
- `strand` is the strand information in the BAM file, where `XS` is used for unstranded data. Please refer to the [regtools documentation](https://regtools.readthedocs.io/en/latest/commands/junctions-extract/) for more information.
- `minimum_reads` is the minimum number of reads required to cluster the introns.
- `min_coverage` is the minimum coverage required for the intron to be considered.
- `min_samples_per_intron` is the minimum number of samples required for the intron to be considered.
- `min_samples_per_group` is the minimum number of samples required for the group to be considered.
- `FDR` is the false discovery rate for the differential splicing analysis.

## Docker image used in the workflow

- [griffithlab/regtools:release-1.0.0](https://hub.docker.com/r/griffithlab/regtools)
- [naotokubota/leafcutter:0.2.9](https://hub.docker.com/r/naotokubota/leafcutter)
