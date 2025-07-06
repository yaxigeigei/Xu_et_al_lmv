# np_babble

# Setup Python Environment

1. Create a conda environment with python 3.11 and install main packages:
    ```bash
    conda create -n babble python=3.11 jupyterlab pandas dill matplotlib
    ```
2. Install pynwb from conda-forge:
    ```bash 
    conda install -c conda-forge pynwb
    ```
3. Install local packages in editable mode:
    ```bash
    cd MSessionExplorer
    pip3 install -e .

    cd sca
    pip3 install -e .
    ```
4. Install PyTorch with CUDA 12.6:
    ```bash
    pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu126
    ```
5. Install graph visualization packages:
    Download graphviz from https://graphviz.org/download/
    and install it in the system.
    ```bash
    pip3 install torchview torchviz graphviz
    ```
