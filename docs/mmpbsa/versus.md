---
template: main.html
title: xBFreE vs other programs
---

# Comparison of **xBFreE** vs other programs
This comparison is based on the documentation of the different programs


### Calculation features
| Feature                      |                [xBFreE][0]                |              [gmx_MMPBSA][1]              |              [MMPBSA.py][4]               |                 [g_mmpbsa][2]                  |       [CaFe][3]       |
|:-----------------------------|:-----------------------------------------:|:-----------------------------------------:|:-----------------------------------------:|:----------------------------------------------:|:---------------------:|
| **Binding free energies**    |                 PB and GB                 |                 PB and GB                 |                 PB and GB                 |                       PB                       |          PB           |
| * PB models                  |           Linear and Non-Linear           |           Linear and Non-Linear           |         Linear and Non-Linear[^1]         |             Linear and Non-Linear              | Linear and Non-Linear |
| * GB models                  |          1, 2, 5, 7, 8 and NSR6           |          1, 2, 5, 7, 8 and NSR6           |             1, 2, 5, 7 and 8              |                                                |                       |
| **Stability**                | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |                                                |                       |
| **Alanine scanning**         | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} [^2] |                       |
| **Entropy corrections** [^3] |           NMODE, QH, IE, and C2           |           NMODE, QH, IE, and C2           |               NMODE and QH                |                                                |                       |
| **Decomposition schemes**    |         Per-Residues and Per-Wise         |         Per-Residues and Per-Wise         |         Per-Residues and Per-Wise         |                  Per-Residues                  |                       |
| **QM/MMGBSA**                | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |                                                |                       |
| **MM/3D-RISM**               | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |                                                |                       |
| **Support Membrane Protein** | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |                                                |                       |
| **Approximations**           |                 ST and MT                 |                 ST and MT                 |                 ST and MT                 |                       ST                       |          ST           |


### Technical features
| Feature                   |                [xBFreE][0]                |              [gmx_MMPBSA][1]              |              [MMPBSA.py][4]               |        [g_mmpbsa][2]         |      [CaFe][3]      |
|:--------------------------|:-----------------------------------------:|:-----------------------------------------:|:-----------------------------------------:|:----------------------------:|:-------------------:|
| **MD Program**            |     GROMACS, AMBER, NAMD, CHARMM [^7]     |                  GROMACS                  |                   AMBER                   |           GROMACS            |     NAMD, AMBER     |
| * GROMACS Version         |              5.x and 20xx.x               |            4.x, 5.x and 20xx.x            |                    ---                    |   4.x, 5.x and 2016+ [^6]    |                     |
| **Dependencies**          |        AmberTools >=20 and Delphi         |              AmberTools >=20              |                AmberTools                 | APBS (1.2.x, 1.3.x or 1.4.x) | VMD, APBS or Delphi |
| **Parallel computation**  | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |       Depends on APBS        |                     |
| **Steps for:**            |                                           |                                           |                                           |                              |                     |
| * Calculation and Summary |                    One                    |                    One                    |                    One                    |           Multiple           |         One         |
| * Analysis                |                    One                    |                    One                    |                 Multiple                  |           Multiple           |      Multiple       |


## Analysis features
Please, check the **xBFreE-Analyzer** documentation for more information.

| Feature                         |                [xBFreE][0]                |              [gmx_MMPBSA][1]              |              [MMPBSA.py][4]               |                [g_mmpbsa][2]                 | [CaFe][3] |
|:--------------------------------|:-----------------------------------------:|:-----------------------------------------:|:-----------------------------------------:|:--------------------------------------------:|:---------:|
| **API**                         | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |                                              |           |
| **Analyzer Tool**               | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |                   [^4]                    |                                              |           |
| * Multiple systems at same time | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |                                           |                                              |           |
| * Correlation between systems   | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |                                           |  :material-check-bold:{.scale_icon_medium}   |           |
| * Per-residue energies to PDB   | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |                                           |  :material-check-bold:{.scale_icon_medium}   |           |
| * Interactive visualization     | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |                                           |                                              |           |
| ** _3D Molecular Visualization_ |                   PyMOL                   |                   PyMOL                   |                                           |                                              |           |
| ** _Interactive Charts_         | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |                                           |                 static image                 |           |
| * Plotting tool                 |              xBFreE-Analyzer              |              gmx_MMPBSA_ana               |       API and graphics library [^5]       |                internal tools                |           |
| * Energetic Terms charts        |                    All                    |                    All                    |                                           | ΔG~polar~, ΔG~nonpolar~, ΔE~MM~ and ΔG~bind~ |           |
| * Export data to CSV file       | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |                                              |           |
| ** _Energy Summary_             | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |                                              |           |
| ** _Individual Energetic Terms_ | :material-check-bold:{.scale_icon_medium} | :material-check-bold:{.scale_icon_medium} |                                           |                                              |           |


  [^1]: Requires the user to modify manually the *.mdin input files 
  [^2]: Without documentation
  [^3]: NMODE = Normal modes approximation, QH = Quasic-Harmony approximation, IE = Interaction Entropy
approximation, and C2 = C2 Entropy
  [^4]: We plan to extend xBFreE-Analyzer compatibility to MMPBSA.py's results
  [^5]: Currently there is a repository ([AmberUtils][5]) for analysing the results
  [^6]: GROMACS 20xx.x is not officially supported. There is a Pull Request that offers a minimum compatibility 
with versions higher than 2016.x one, but still with limitations
  [^7]: xBFreE currently supports a variety of topology/trajectory formats including but not limited to: *.top, *.prmtop,
  *.psf, and *.xtc, *.mdcrd, *.dcd. These formats can be used by other MD engines, and hence they are supported in xBFreE 
  as well.


  [0]: https://github.com/xBFreEnergy/xBFreE
  [1]: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA
  [2]: https://github.com/RashmiKumari/g_mmpbsa
  [3]: https://github.com/huiliucode/cafe_plugin
  [4]: https://ambermd.org/doc12/Amber21.pdf#chapter.36
  [5]: https://github.com/williamdlees/AmberUtils

  