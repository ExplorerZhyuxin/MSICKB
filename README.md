# MSICKB analysis

This repository contains analysis code and reproducible workflows for the MSICKB manuscript.

## Download
- **Database snapshot (v1.0.0)**: https://github.com/ExplorerZhyuxin/MSICKB/releases/tag/v1.0.0
- **Latest releases**: https://github.com/ExplorerZhyuxin/MSICKB/releases
- **MSICKB website**: <YOUR_WEBSITE_URL>

## Contents
- `scripts/`: analysis scripts used to generate network statistics, hub analyses, sensitivity analyses, and TCGA/cBioPortal validation.
- `results/figures/`: generated figures used in the manuscript and supplement.
- `results/tables/`: generated tables used in the manuscript and supplement.

## Requirements
- R (>=4.3) and/or Python (>=3.10)
- See `environment.yml` (or `requirements.txt` / `renv.lock`) for dependencies.

## Reproducibility
1. Clone this repository
2. Install dependencies
3. Run scripts in order:
   - `scripts/01_build_network.*`
   - `scripts/02_hub_analysis.*`
   - `scripts/03_assay_sensitivity.*`
   - `scripts/05_tcga_cbioportal_validation.*`

## Data availability
Due to licensing and size constraints, raw literature-curated data are not fully redistributed in this repository.
A curated snapshot (tables used in the manuscript) is available at: [LINK_TO_DOWNLOAD_PAGE_OR_RELEASE]
All analyses are scripted to run from the provided processed tables.

## License
MIT License.
