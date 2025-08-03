# Decoding the interconnected splicing patterns of hepatitis B virus and host using large language and deep learning models

This repository contains the full pipeline for [our manuscript](https://doi.org/10.1101/2025.07.28.667110), covering data preprocessing,  statistical analysis, and embedding extraction, dimensional reduction, and clustering, and splice site predictions. All steps are documented within the provided Jupyter notebooks.


## Setup and Installation
Here are the steps to re-create conda environments used in this study. For simplicity, these steps will install Miniforge3 in the home directory.
1. Download and install Miniforge3
```
wget https://github.com/conda-forge/miniforge/releases/download/25.3.0-3/Miniforge3-25.3.0-3-Linux-x86_64.sh
bash Miniforge3-25.3.0-3-Linux-x86_64.sh -b -p $HOME/miniforge3 # can be installed elsewhere
eval "$(/$HOME/miniforge3/bin/conda shell.bash hook)"
```
2. Download this GitHub repository
```
git clone https://github.com/lcscs12345/HBV_splicing_paper_2025.git
```
3. Create conda environments
```
# conda environment with various utilities installed
conda env create --file environment_files/environment_utils.yml -p $HOME/miniforge3/envs/utils
# Fix pyfasta to make it compatible with python 3.13.5 and numpy 2.3.1
sh environment_files/fix_pyfasta.sh

# conda environment with OpenSpliceAI and dependencies installed
conda env create --file environment_files/environment_openspliceai.yml -p $HOME/miniforge3/envs/openspliceai

# conda environment for SpliceBERT dependencies
conda env create --file environment_files/environment_llm.yml -p $HOME/miniforge3/envs/llm
```
4. Download project files and SpliceBERT and OpenSpliceAI models from [Zenodo](https://doi.org/10.5281/zenodo.16730945) and unzip them within this repository. The complete directory structure should look like:
```
HBV_splicing_paper_2025/
└── data/
└── environment_files/
└── jupyter_notebooks/
└── ref/ 
└── results/
└── scripts/
└── src/
```


## Data Description
The `data/processed_files/` directory contains intermediate data files generated
during the various stages of the analysis pipeline. These files are typically derived
from raw data through complex transformations and computations performed by Bash
commands and Python scripts executed within the Jupyter notebooks.

Key output files include:
- **Cleaned Data Files (`.csv`, `.pkl.gz`):** Data files in results/data.
- **Figure Files  (`.pdf`, '.png'):** Plots in results/figures.

For a comprehensive list and descriptions of all automatically detected project files, please refer to [`./PROJECT_FILES.md`](./PROJECT_FILES.md).

## Notebooks
### [01-splicing_efficiency.ipynb: Quantification of spliced HBV RNAs across tissues and cell lines](https://github.com/lcscs12345/HBV_splicing_paper_2025/jupyter_notebooks/01-splicing_efficiency.ipynb)
Compares the proportions of spliced HBV RNAs and evaluates splice site-level splicing efficiency across liver biopsy tissues and cultured cell lines. It highlights that the splice site-level measurement may serve as a stronger biomarker than the overall proportions of HBV splicing.
- **Input Files:**
  - `../data/huh7/SraRunTable.csv`: Comma-separated values file with tabular data.
  - `../data/map.txt`: Tab-delimited text file with processed results.
  - `../data/processed_files/cosi.hbv.txt`: Tab-delimited text file with processed results.
  - `../data/processed_files/cosi_long.pkl.gz`: Long-format coSI scores for HBV splice donor and acceptor sites.
  - `../data/processed_files/mgen-7-492-s002.xlsx`: Supplementary file from our previous study
  - `../data/processed_files/track.lol.txt`: Tab-delimited text file with processed results.
  - `../data/tcons.txt`: Tab-delimited text file with processed results.
- **Output Files:**
  - This notebook generates 4 'Generated plot or figure output.' files, 2 'Tab-delimited text file with processed results.' files, 1 'Final proportions of spliced HBV RNAs used in results.' files, 6 other files. For a comprehensive list, see the [Project Files README](./PROJECT_FILES.md).

### [02-splicebert_embeddings.ipynb: Nucleotide embedding analysis for HBV and host splicing patterns](https://github.com/lcscs12345/HBV_splicing_paper_2025/jupyter_notebooks/02-splicebert_embeddings.ipynb)
Extracts nucleotide embeddings, applies dimensionality reduction, and performs clustering to uncover splice site sequence patterns in both HBV and host genomes.
- **Input Files:**
  - This notebook uses 3 'Coordinate mapping file in MAFFT mapout format.' files, 3 'Splice site statistics across clusters.' files, 2 'Splice donor/acceptor site coordinates.' files, 7 other files. For a comprehensive list, see the [Project Files README](./PROJECT_FILES.md).
- **Output Files:**
  - This notebook generates 10 'Splice donor/acceptor site coordinates.' files, 8 'Genomic features in BED format.' files, 7 'Generated plot or figure output.' files, 19 other files. For a comprehensive list, see the [Project Files README](./PROJECT_FILES.md).

### [03a-splicebert_prediction-donor.ipynb: Splicing propensity and donor site classification using SpliceBERT](https://github.com/lcscs12345/HBV_splicing_paper_2025/jupyter_notebooks/03a-splicebert_prediction-donor.ipynb)
Analyses splicing propensity and classifies HBV splice donor sites using SpliceBERT, leveraging transformer-based sequence representations to identify predictive features. Includes performance metrics evaluation and splice donor site conservation analysis to assess model accuracy and evolutionary constraints, respectively.
- **Input Files:**
  - `../data/processed_files/cosi_long.pkl.gz`: Long-format coSI scores for HBV splice donor and acceptor sites.
- **Output Files:**
  - `../ref/hbvdb/pgrna/pgrna_flank200.txt`: Tab-delimited text file with processed results.
  - `../results/figures/fig4/umap.logit_donors.leiden.png`: Generated plot or figure output.
  - `../results/figures/fig4/umap.logit_donors.png`: Generated plot or figure output.
  - `../results/figures/fig5/conservation_donors.png`: Generated plot or figure output.
  - `../results/figures/fig5/logit_donors.png`: Generated plot or figure output.

### [03b-splicebert_prediction-acceptor.ipynb: Splicing propensity and acceptor site classification using SpliceBERT](https://github.com/lcscs12345/HBV_splicing_paper_2025/jupyter_notebooks/03b-splicebert_prediction-acceptor.ipynb)
Analyses splicing propensity and classifies HBV splice acceptor sites using SpliceBERT, leveraging transformer-based sequence representations to identify predictive features. Includes performance metrics evaluation and splice acceptor site conservation analysis to assess model accuracy and evolutionary constraints, respectively.
- **Input Files:**
  - `../data/processed_files/cosi_long.pkl.gz`: Long-format coSI scores for HBV splice donor and acceptor sites.
- **Output Files:**
  - `../results/figures/fig4/umap.logit_acceptors.leiden.png`: Generated plot or figure output.
  - `../results/figures/fig4/umap.logit_acceptors.png`: Generated plot or figure output.
  - `../results/figures/fig5/conservation_acceptors.png`: Generated plot or figure output.
  - `../results/figures/fig5/logit_acceptors.png`: Generated plot or figure output.

### [04-spliceai_prediction.ipynb: Splice site classification with OpenSpliceAI](https://github.com/lcscs12345/HBV_splicing_paper_2025/jupyter_notebooks/04-spliceai_prediction.ipynb)
Integrates a deep learning-based framework into the study by performing splice site classification with OpenSpliceAI. Includes performance metrics evaluation to assess model accuracy.
- **Input Files:**
  - `../data/processed_files/cosi_long.pkl.gz`: Long-format coSI scores for HBV splice donor and acceptor sites.
  - `../ref/hbvdb/pgrna/pgrna_flank200.txt`: Tab-delimited text file with processed results.
- **Output Files:**
  - `../results/figures/fig5/openspliceai_acceptors.png`: Generated plot or figure output.
  - `../results/figures/fig5/openspliceai_donors.png`: Generated plot or figure output.

### [05-sequence_logo.ipynb: Plot splice site motif frequencies](https://github.com/lcscs12345/HBV_splicing_paper_2025/jupyter_notebooks/05-sequence_logo.ipynb)
Helps generate the frequency of splice site motifs across genomic sequences.
- **Input Files:**
  - `../data/processed_files/cosi_long.pkl.gz`: Long-format coSI scores for HBV splice donor and acceptor sites.
- **Output Files:**
  - This notebook generates 10 'FASTA file containing nucleotide sequences.' files, 6 'Genomic features in BED format.' files, 6 'Mapped coSI scores to exonic splice sites.' files, 2 other files. For a comprehensive list, see the [Project Files README](./PROJECT_FILES.md).

## Scripts
- [scripts/common.py](https://github.com/lcscs12345/HBV_splicing_paper_2025scripts/common.py): Utility functions for reading fasta sequences, motif extraction, and calculate performance metrics for a classifier.
- [scripts/dicts.py](https://github.com/lcscs12345/HBV_splicing_paper_2025scripts/dicts.py): Dictionaries for HBV and human splice sites.
- [scripts/generate_readme.py](https://github.com/lcscs12345/HBV_splicing_paper_2025scripts/generate_readme.py): Functions to generate README files automatically.
- [scripts/openspliceai_helpers.py](https://github.com/lcscs12345/HBV_splicing_paper_2025scripts/openspliceai_helpers.py): Wrapper for OpenSpliceAI to generate and format splicing predictions from input sequences using PyTorch and pandas.
- [scripts/splicebert_helpers.py](https://github.com/lcscs12345/HBV_splicing_paper_2025scripts/splicebert_helpers.py): Wrapper for SpliceBERT splice site classification using HuggingFace Transformers with a sliding window approach, alongside utilities for identifying non-splice sites and calculating normalised mutual information (NMI).
- [scripts/track.r](https://github.com/lcscs12345/HBV_splicing_paper_2025scripts/track.r): R functions to generate track plots for HBV splice variants and splicing efficiency at each splice site in the viral genome.

