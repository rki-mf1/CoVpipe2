# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.6] 2022-07-27

### Fixed

- report can handle Kraken2 output with no human, unclassified or tax_id reads

## [0.2.5] 2022-07-19

### Changed

- changed `bcftools consensus` parameter from `-I` to `-H I`
  - does not change behaviour
- genotype adjustment for all heterzygote variants
- updated default pangolin version to 4.0.6
- renamed parameter `pangolin_skip_scorpio` to `pangolin_scorpio` (default behavior of skipping scorpio does not change)
- Kraken2 read output is not compressed - changed the file naming accordingly

### Removed

- removed unsued script
- removed unused conda yamls

## [0.2.4] 2022-06-14

### Added

- added sc2rf 
  - strict params to find potential recombinants
- added skip-scorpio parameter for pangolin
  - activeted per defaul
- added this changelog

### Changed

- separate result tabled get published
- better primer/reference id check
  - fist fasta id has to be at least once in the primer file (some primer kits have a non-sc2 amplicon)

## [0.2.3] 2022-06-01

### Changed

- kraken2 database can be set by user
- hard filtered variants get published in the variant calling result directory

### Removed

- removed `?` removal since deletion adjustment is removed
- update processes for pangolin and nextclade

### Fixed

- fixed container execution for pangolin and nextclade
- fixed disabling/enabling of `var_sap` parameter

## [0.2.2] 2022-06-01

### Changed

- updated workflow figure
- more stable automated conda update for pangolin and nextclade directly from anaconda.org
- updated help/readme

### Removed

- removed deletion adjustment since bcftools works as expected

## [0.2.1] 2022-05-09

### Added

- LCS update functionality

## [0.2.0] 2022-05-05

### Added

- LCS fro mixed sample detection

### Changed

- updated versions of samtools, liftoff, kraken, freebayes, bcftools, fastp, bedtools

### Fixed

- `bedtools subtract` instead of `intersect` for low coverage mask generation

## [0.1.2] 2022-04-13

### Changed

- updated help and readme
- disabel `var-sap` in default
- removed workflow version from consensus fasta header

## [0.1.1] 2022-04-07

- fixed krona fails with mamba/conda execution
- fixes pangolin update
- update nextclade

## [0.1.1-alpha] 2022-04-07

### Added

- DESH support

## [0.1.0-alpha] 2022-04-06

- initial commit