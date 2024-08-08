# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## [0.1.0] - 2024-08-07

### Added
- Started this CHANGELOG.
- Defined `AveChunk.__getitem__` and `AveChunk.__len__` functions for custom indexing and length from `AveChunk` objects.
- Type annotations and docstrings to the public functions of `AveChunk`.

### Changed
- Renamed the repository to `lammps-post` and the package to `lmp_post`. Resetting the version number back to 0.1.0.
- Replaced the setup.py script with a pyproject.toml file for pip installation.
- Updated the README, adding additional information and sections.
- Changed module and class names to mirror the corresponding LAMMPS `fix ave/chunk` command: `chunkavg.py` -> `avechunk.py` and `ChunkAvg` -> `AveChunk`
- Applied Black formatting to `avechunk.py` and `thermoout.py` modules.
- Updated the notebooks in `jupyter-notebooks` with new names of the package and `AveChunk` class, as well as polished the text a bit.


## [Unreleased] - yyyy-mm-dd

N/A

### Added

### Changed

### Fixed