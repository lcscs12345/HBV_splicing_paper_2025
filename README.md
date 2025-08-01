##### Lim CS, Brown CM (2025). Decoding the interconnected splicing patterns of hepatitis B virus and host using large language and deep learning models. BioRxiv. doi: https://doi.org/10.1101/2025.07.28.667110

This repository provides Jupyter Notebooks and scripts used in this study.

Here are the steps to re-create conda environments used in this study. For simplicity, these steps will install Miniforge3 in the home directory.
1. Download and install Miniforge3
```
wget https://github.com/conda-forge/miniforge/releases/download/25.3.0-3/Miniforge3-25.3.0-3-Linux-x86_64.sh
bash Miniforge3-25.3.0-3-Linux-x86_64.sh -b -p $HOME/miniforge3 # can be installed elsewhere
eval "$(/$HOME/miniforge3/bin/conda shell.bash hook)"
```
2. Create conda environments using yml files
```
# conda environment for various utilities
conda env create --file environment_files/environment_utils.yml -p /$HOME/miniforge3/envs/utils
# Fix pyfasta to make it compatible with python 3.13.5 and numpy 2.3.1
sh environment_files/fix_pyfasta.sh

# conda environment for OpenSpliceAI dependencies
# OpenSpliceAI will be installed within the environment
conda env create --file environment_files/environment_openspliceai.yml -p /$HOME/miniforge3/envs/openspliceai

# conda environment for SpliceBERT dependencies
# SpliceBERT will need to be downloaded separately from https://github.com/chenkenbio/SpliceBERT/ or Zenodo
conda env create --file environment_files/environment_llm.yml -p /$HOME/miniforge3/envs/llm
```


