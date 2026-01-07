# Changelog
All notable changes to this project will be documented in this file.

The format is based on Keep a Changelog, and this project adheres to Semantic Versioning.

## [0.1.0] - 2026-01-07
### Added
- Packaged a CLI tool (`clump-history`) with subcommands:
  - `run`: forward-model Δ47 along a thermal history and generate a 2×3 scenario figure (PDF/SVG).
  - `ufit`: apply a constrained U-shaped peak adjustment to a thermal history and export an adjusted CSV.
- Added `--outdir` to write outputs into a specified directory (auto-created if missing).
- Refactored notebook-style workflow into a maintainable Python package structure (`src/` layout) to support `python -m clump_history ...` usage.
- Saved figures as both `.pdf` and `.svg` for publication-quality output.

### Notes
- Default input paths assume a workspace layout containing `./datasets/`.
  If you run the CLI from another directory, please pass `--thermal` / `--test` with explicit paths.
