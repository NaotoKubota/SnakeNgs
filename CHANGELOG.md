# Change log

All notable changes to this SnakeNgs project will be documented in this file.

## [v0.4.0] - 2026-??-??

### Added

- `preprocessing_RNAseq_long.smk`: Long-read RNA-seq preprocessing pipeline.
- `preprocessing_CLIPseq.smk`: Unified CLIP-seq preprocessing pipeline (HITS-CLIP, iCLIP-seq, etc.) with `fastp_args`, `outFilterMultimapNmax`, and `normalize_bigwig` config parameters.
- `ONT_plasmid_assembly.smk`: ONT sequencing data analysis for plasmid assembly.
- `snakefile/common/containers.smk`: Centralized Docker container definitions for all pipelines.
- `snakefile/common/functions.smk`: Shared utility functions (e.g., `check_chromosome_length`).
- GitHub Actions CI workflow for dry-run testing all Snakefiles on push/PR to `develop`.

### Changed

- `preprocessing_RNAseq.smk`:
  - Merged `preprocessing_RNAseq.smk` and `preprocessing_RNAseq_single.smk` into a single file with `layout` config parameter (`"paired"` or `"single"`).
  - Extracted common rules (`sort`, `bigwig`, `makeRefFlat`, `makeRibosomalInterval`, `CollectRnaSeqMetrics`) into `snakefile/common/`.
  - Check chromosome lengths before indexing the bam file and make .csi index if a chromosome is over 512 Mbp.
- `preprocessing_ChIPseq.smk`:
  - Merged `preprocessing_ChIPseq.smk` and `preprocessing_ChIPseq_single.smk` into a single file with `layout` config parameter (`"paired"` or `"single"`).
  - Extracted common rules (`sort`, `markdup`, `index`, `plotFingerprint`, `bigwig`) into `snakefile/common/`.
- All Snakefiles: Centralized container definitions into `snakefile/common/containers.smk` using `CONTAINERS["key"]` references instead of hardcoded strings.
- All Snakefiles: Moved per-rule `wildcard_constraints` to workflow-level.

### Removed

- `preprocessing_RNAseq_single.smk`: Merged into `preprocessing_RNAseq.smk`.
- `preprocessing_ChIPseq_single.smk`: Merged into `preprocessing_ChIPseq.smk`.
- `preprocessing_HITSCLIP.smk`: Merged into `preprocessing_CLIPseq.smk`.
- `preprocessing_iCLIPseq.smk`: Merged into `preprocessing_CLIPseq.smk`.

## [v0.3.2] - 2025-03-30

### Changed

- Use `multiqc/multiqc:v1.28` image for MultiQC.

## [v0.3.1] - 2025-03-30

### Deprecated

- `ngsFetch` is deprecated and will be removed in the next version. Please use [ngsfetch](https://github.com/NaotoKubota/ngsfetch) instead.

## [v0.3.0] - 2025-03-07

### Added

- `footprinting_timeseries_ATACseq.smk`: ATAC-seq footprinting time series analysis by TOBIAS.
- `STARsolo.smk`: gene count quantification from single-cell/nucleus RNA-seq data by STARsolo.
- `cellranger_count.smk`: gene count quantification from single-cell/nucleus RNA-seq data by Cell Ranger.

### Changed

- `callpeak_ATACseq.smk`: Remove the `bedgraphtobigwig` rule.
- `preprocessing_ChIPseq.smk`: Add `CollectInsertSizeMetrics` rule to collect metrics.
- `kb-nac.smk`: Add rules for summarizing filtered BUS files and MultiQC.

## [v0.2.0] - 2024-10-08

### Added

### Changed

- `MAJIQ.smk`:
  - Use local singularity image for MAJIQ as there is no public image available.
  - Add a rule to run `majiq deltapsi` module.
- `Whippet.smk`: Use single thread to allow parallel execution of the `quant` rule.
- `LeafCutter.smk`: Remove rule `plottingSpliceJunctions`.
- `preprocessing_RNAseq.smk`: Add rule `CollectRnaSeqMetrics` and `CollectInsertSizeMetrics` to collect metrics.
- Use `multiqc/multiqc:v1.25` image for MultiQC.

### Fixed

## [v0.1.4] - 2024-08-17

### Added

### Changed

- Update README.md.

### Fixed

- Fix documentation.

## [v0.1.3] - 2024-08-17

### Added

### Changed

### Fixed

- Fix documentation.

## [v0.1.2] - 2024-08-17

### Added

### Changed

### Fixed

- Fix the bug in the `Dockerfile`.

## [v0.1.1] - 2024-08-17

### Added

- Initial release

### Changed

### Fixed
