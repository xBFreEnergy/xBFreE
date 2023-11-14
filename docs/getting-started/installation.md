---
template: main.html
title: Installation
---
# **xBFreE** Installation

Currently, **xBFreE** can be installed using three ways:

- [Manual installation](#manual-installation)
- Docker Package (Coming soon!)
- Singularity Package (Coming soon!)

## Manual installation
You can carry out the installation of **xBFreE** in three ways: 

[`conda environment`](#gs_gs_conda_environment){#gs_gs_conda_environment} (Recommended for PC installation)
:   The conda environment provides a clean and efficient way of installing **xBFreE**. It also allows having 
    different versions of xBFreE in isolated environments, thus reducing the possibility of incompatibility with 
    other packages. Installation time is also less since it does not require the compilation of AmberTools or GROMACS.
    To install xBFreE in a conda environment, please, follows these steps:
    
    ??? "1-Miniconda Installation (skip this step if conda is already installed in your PC)"
    
        Download and install [Miniconda]
    
        <div class="termy">
    
        ```bash
        $ curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
        ---> 100%
    
        $ chmod +x Miniconda3-latest-Linux-x86_64.sh
    
        $ ./Miniconda3-latest-Linux-x86_64.sh
        ---> 100%
    
        Successful miniconda intallation
        ```
    
        </div>
    
        ??? note "Copy described intructions"     
    
            ``` bash 
            curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh   # (1)!
            chmod +x Miniconda3-latest-Linux-x86_64.sh                                      # (2)!
            ./Miniconda3-latest-Linux-x86_64.sh                                             # (3)!
           
            ```
        
            1. Download Miniconda installer
            2. Change permissions for the installer
            3. Execute and install miniconda

    !!! -n "2-Create conda environment and install libraries"
    
        !!! Info "Important"
            Make sure to have conda installed in your computer. 
    
        === "`*.yml file`"    
            Installing **xBFreE** using a requierement yml file. 
    
            :material-file-download-outline:{:.heart } Download **[env.yml](../env.yml)** file
    
            <div class="termy">
    
            ```console
            // Create a new environment and use the *.yml file to install dependencies
            $ conda env create -n xbfree-env --file env.yml
    
            // To use xBFreE, just activate the environment
            $ conda activate xbfree-env
            ```
                
            </div>
    
            !!! important
                The latest release of **xBFreE** available in PYPI will be installed
    
            ??? note "Copy described intructions"     
    
                ``` bash 
                conda env create -n xbfree-env --file env.yml                                    # (1)!
                conda activate xbfree-env                                                        # (2)!
               
                ```
            
                1. Create the `xbfree-env` environment and use the *.yml file to install dependencies
                2. Activate `xbfree-env` environment
    
    
        === "`pip`"    
            Installing dependencies
    
            <div class="termy">
    
            ```console
            // Update conda
            $ conda update conda
            
            // Create a new environment and activate it
            $ conda create -n xbfree-env python=3.9 -y -q 
            $ conda activate xbfree-env
            
            // Install mpi4py and AmberTools
            $ conda install -c conda-forge mpi4py=3.1.3 ambertools=21.12 compilers=1.2.0 -y -q
            
            // Install updated version of ParmEd
            $ python -m pip install git+https://github.com/Valdes-Tresanco-MS/ParmEd.git@v3.4
            
            // Install PyQt5 required to use the GUI analyzer tool (xBFreE-Analyzer). Not needed for HPC
            $ python -m pip install pyqt5
    
            // (Optional) Install GROMACS
            $ conda install -c conda-forge gromacs==2022.4 -y -q
            
            // Install the latest version available in PYPI and xBFreE-Analyzer
            $ python -m pip install xbfree[xbfree-analyzer]
            ```                
            </div>
    
            ??? note "Copy described intructions"     
    
                ``` bash 
                conda update conda
                conda create -n xbfree-env python=3.9 -y -q                                      # (1)
                conda activate xbfree-env                                                        # (2)
                conda install -c conda-forge mpi4py=3.1.3 ambertools=21.12 compilers=1.2.0 -y -q      # (3)
                python -m pip install git+https://github.com/Valdes-Tresanco-MS/ParmEd.git@v3.4 # (4)
                python -m pip install pyqt5                                                     # (5)
                # Optional
                conda install -c conda-forge gromacs==2022.4 -y -q                                 # (6)
               
                ```
            
                1. Create `xbfree-env` environment
                2. Activate `xbfree-env` environment
                3. Install dependencies
                4. Install ParmEd
                5. Install PyQt5 if you will use xBFreE-Analyzer
                6. (Optional) Install GROMACS if GROMACS is not installed in your machine
    
        [Miniconda]: https://docs.conda.io/en/latest/miniconda.html

[`virtual environment`](#gs_gs_virtual_environment){#gs_gs_virtual_environment} (Recommended for HPC installation)
:   The virtual environment provides a clean and efficient way of installing **xBFreE**. It also allows having 
different versions of **xBFreE** in isolated environment, thus reducing the possibility of incompatibility with 
other packages. Differently to `conda environment` this environment only contain python packages and not 
pre-compiled libraries.

:   In HPC, admins recommend using a virtual environment (_virtualenv_) instead of _conda_ to take advantage of optimized 
compiled programs. In this case, you will require AmberTools, Gromacs, and other dependencies previously compiled 
(generally as modules). The required libraries for each dependency will depend on the HPC the user is working on. A 
good example and how the _virtualenv_ works in an HPC can be 
found [here](https://docs.alliancecan.ca/wiki/GROMACS#gmx_MMPBSA).
    
[`AmberTools compilation`](#gs_gs_AmberTools_compilation){#gs_gs_AmberTools_compilation} (Recommended when you need to modify sander or compile AmberTools in a specific way)
:   This way, we assume that you have AmberTools compiled on your machine and that you want to do an installation 
without worrying about enabling or disabling conda environments. It also involves user compilation of GROMACS, which 
takes considerable installation time. This way also requires installed packages to be compatible and installation 
errors are more frequent.

    !!! AmberTools compilation
        [Follow the oficial AmberTools installation according to your OS](https://ambermd.org/Installation.php)
        !!! note
            We asume that AmberTools and their shell environment are correctly configured
    
        === "Rolling/stable release"
            **INSTALLATION**
            <div class="termy">
            ```console
            // Install updated ParmEd
            $ amber.python -m pip install git+https://github.com/Valdes-Tresanco-MS/ParmEd.git@v3.4
            // Install xBFreE
            $ amber.python -m pip install xBFreE                                               
            ```
            </div>
    
            **UPDATE**
            <div class="termy">
            ```console
            // Update xBFreE
            $ amber.python -m pip install xBFreE -U
            ```
            </div>
            
        === "development version" 
            **INSTALLATION**
            <div class="termy">
            ```bash
            // Install updated ParmEd
            $ amber.python -m pip install git+https://github.com/Valdes-Tresanco-MS/ParmEd.git@v3.4
            // Install xBFreE
            $ amber.python -m pip install git+https://github.com/Valdes-Tresanco-MS/xBFreE     
            ```
            </div>
    
            **UPDATE**
            <div class="termy">
            ```bash
            amber.python -m pip install git+https://github.com/Valdes-Tresanco-MS/xBFreE -U
            ```
            </div>
    
            !!! warning
                Install/update xBFreE from the master branch of GitHub repository. This version is only recommended 
                to test a new version or to try temporary solutions to reported bugs.
    
        !!! danger
            If you get an error related to installing `mpi4py`, you may want to install this package manually from 
            `conda-forge` as follows:
    
            ```
            amber.conda install -c conda-forge mpi4py=3.1.3
            ```
            
            If you get an error related to `pip`, you may want to install this package manually as follows:
            
            ```
            amber.conda install pip
            ```

## Docker Package
Coming soon!

## Singularity Package
Coming soon!

## Testing the operation of xBFreE
After preparing everything to run `xBFreE`, it only remains to check its correct operation. To know how to do it, 
consult the documentation of [`xBFreE_test`](../examples/gmx_MMPBSA_test.md#running-gmx_mmpbsa_test)
