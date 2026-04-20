# Local Model Fallback Directory

Put supplementary model files here when the repository-level `models/` directory is incomplete.

Reproduce workflows resolve model files in this order:
1. `models/<model_filename>`
2. Any nested subdirectory under `models/` (recursive filename match)
3. `reproduce/models/<model_filename>`
4. Any nested subdirectory under `reproduce/models/` (recursive filename match)

This allows local reproduction scripts to run even if some model files are missing from the top-level models directory.
