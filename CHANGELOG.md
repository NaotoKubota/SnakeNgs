# Change log

All notable changes to this SnakeNgs project will be documented in this file.

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
