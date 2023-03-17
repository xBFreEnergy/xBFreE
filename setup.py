# ##############################################################################
#                           GPLv3 LICENSE INFO                                 #
#                                                                              #
#  Copyright (C) 2020  Mario S. Valdes-Tresanco and Mario E. Valdes-Tresanco   #
#  Copyright (C) 2014  Jason Swails, Bill Miller III, and Dwight McGee         #
#                                                                              #
#   Project: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA                  #
#                                                                              #
#   This program is free software; you can redistribute it and/or modify it    #
#  under the terms of the GNU General Public License version 3 as published    #
#  by the Free Software Foundation.                                            #
#                                                                              #
#  This program is distributed in the hope that it will be useful, but         #
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  #
#  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License    #
#  for more details.                                                           #
# ##############################################################################

from setuptools import setup, find_packages
import versioneer
import sys

with open("README.md", "r") as f:
    LONG_DESCRIPTION = f.read()

if sys.version_info[:2] < (3, 8):
    raise RuntimeError("seaborn requires python >= 3.8.")

setup(
    name='BFEnergy',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    package_data={"BFEnergy": ["data/*", 'data/gmxMMPBSA/*', 'data/xvv_files/*', 'analyzer/style/*', 'GMXMMPBSA.sh']},
    license='GPLv3',
    author='Mario S. Valdes-Tresanco and Mario E. Valdes-Tresanco ',
    author_email='mariosergiovaldes145@gmail.com',
    maintainer='Mario S. Valdes-Tresanco',
    maintainer_email='mariosergiovaldes145@gmail.com',
    url='https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA',
    description="gmx_MMPBSA is a new tool based on AMBER's MMPBSA.py aiming to perform end-state free energy  "
                "calculations with GROMACS files",
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    keywords=['BFEnergy', 'MMPBSA', 'MMGBSA', 'LIE', 'GROMACS', 'AMBER', 'NAMD', 'AmberTools'],
    install_requires=['pandas>=1.2.2', 'seaborn>=0.11.2', 'mpi4py>=3.1.3', 'scipy>=1.6.1', 'matplotlib>=3.5.1',
                      'h5py'],
    entry_points={
        "console_scripts": [
            "xbfree=GMXMMPBSA.app:gmxmmpbsa",
            "xBFreE=GMXMMPBSA.app:gmxmmpbsa",
            "namd_MMPBSA=GMXMMPBSA.app:gmxmmpbsa",
            "gmx_LIE=GMXMMPBSA.app:gmxmmpbsa",
            "amber_LIE=GMXMMPBSA.app:gmxmmpbsa",
            "namd_LIE=GMXMMPBSA.app:gmxmmpbsa",
            "BFEnergyAna=GMXMMPBSA.app:gmxmmpbsa_ana",
            "BFEnergyTest=GMXMMPBSA.app:gmxmmpbsa_test"]}
)
