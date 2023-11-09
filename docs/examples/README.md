---
template: main.html
title: Examples
---

Here you can find a representation of almost all the types of calculations and analyses that you can perform 
with xBFreE. Although each example focuses on specific cases, you can use xBFreE on systems that combine 
a number of different components (_i.e._ metalloprotein-ligand complex, Protein-DNA-ligand, etc.). In addition, 
several types of calculations (_e.g._ GB, Alanine scanning and Per-residue decomposition; PB, Interaction Entropy, 
and Per-wise decomposition) can be also performed in the same run for a specific system.

## Systems

This is a representation of the systems that can be processed and analyzed with xBFreE. Our program employs a robust 
method to process input structures. Even if your system is not represented here, still have a go, you won't 
be disappointed! 😀

* [Protein-protein](Protein_protein/README.md)[^1][^2][^3]
* [Protein-ligand](Protein_ligand/ST/README.md)[^1][^2]  

[//]: # (* [Protein-DNA]&#40;Protein_DNA/README.md&#41;[^1][^2][^3])

[//]: # (* [Protein-glycan]&#40;Protein_glycan/README.md&#41;[^1][^2][^3])

[//]: # (* [MMPBSA with membrane proteins]&#40;Protein_membrane/README.md&#41;[^1][^2]  )

[//]: # (* [Metalloprotein-ligand]&#40;Metalloprotein_ligand/README.md&#41;[^1][^2])

[//]: # (* [Multicomponent system &#40;Protein-DNA-RNA-Ions-Ligand&#41;]&#40;Comp_receptor/README.md&#41;[^1][^2][^3])

[//]: # (* COVID-19 related proteins)

[//]: # (    * [Info]&#40;COVID-19_related_proteins/README.md&#41;)

[//]: # (    * [Main protease]&#40;COVID-19_related_proteins/Main_protease_7l5d/README.md&#41;)

[//]: # (    * [Papain-like protease]&#40;COVID-19_related_proteins/Papain-like_protease_7koj/README.md&#41;)

[//]: # (    * [S1-ACE2 complex]&#40;COVID-19_related_proteins/S1-ACE2_complex_7dmu/README.md&#41;)

[//]: # (    * [S1 RBD with antibody]&#40;COVID-19_related_proteins/S1_RBD_with_antibody_6zlr/README.md&#41;)

[//]: # ()
[//]: # (## CHARMMff support)

[//]: # ()
[//]: # (This section focuses more on how to work with systems prepared with CHARMM force fields. We only show few examples )

[//]: # (for better clarity.)

[//]: # ()
[//]: # (* [Protein-ligand]&#40;Protein_ligand_CHARMMff/README.md&#41;[^1][^2])

[//]: # (* [Protein-ligand complex embedded in membrane]&#40;Protein_membrane_CHARMMff/README.md&#41;[^1])

[//]: # (* [Protein-ligand with LPH atoms]&#40;Protein_ligand_LPH_atoms_CHARMMff/README.md&#41;)

[//]: # ()
[//]: # (## OPLSff support)

[//]: # ()
[//]: # (This section focuses more on how to work with systems prepared with OPLS force fields. We only show few examples )

[//]: # (for better clarity.)

[//]: # ()
[//]: # (* [Protein-protein]&#40;OPLS/protein_protein/README.md&#41;)

## Analysis

This section focuses on the analysis that can be performed with xBFreE. Although each example focuses on specific 
cases, you can use xBFreE to perform several types of calculations (_e.g._ GB, Alanine scanning and Per-residue 
decomposition; PB, Interaction Entropy, and Per-residue decomposition) in the same run for a specific system.

* [Single Trajectory Protocol](Protein_ligand/ST/README.md)[^1][^2][^3]
* [Multiple Trajectory Protocol](Protein_ligand/MT/README.md)[^1]
* Binding free energy calculations
    * [Binding free energy calculation with GB](Protein_ligand/ST/README.md)
    * [Binding free energy calculation with GBNSR6](GBNSR6/README.md)
    * [Binding free energy calculation with linear PB (LPBE)](Linear_PB_solver/README.md)
    * [Binding free energy calculation with NonLinear PB (non-LPBE)](NonLinear_PB_solver/README.md)  
    * [Binding free energy calculation with 3D-RISM model](3D-RISM/README.md)[^1]
* [Alanine scanning](Alanine_scanning/README.md)[^1][^2][^3]
* [Decomposition analysis](Decomposition_analysis/README.md)[^1][^2][^3]
* Entropy
    * [Interaction Entropy calculations](Entropy_calculations/Interaction_Entropy/README.md)[^1][^2][^3]
    * [NMODE Entropy calculations](Entropy_calculations/nmode/README.md)[^1]
    * [C2 Entropy calculations](Entropy_calculations/C2_Entropy/README.md)
* [Stability calculations](Stability/README.md)[^1][^2][^3]
* [QM/MMGBSA calculations](QM_MMGBSA/README.md)
* [Correlation](Correlation/README.md)

[//]: # (## Support for psf_dcd files)

[//]: # ()
[//]: # (This section focuses on how to work with psf-dcd files. These files are used for several MD simulation )

[//]: # (programs such as NAMD, OpenMM or GENESIS. We plan to add more examples in the near future.)

[//]: # ()
[//]: # (* [Protein-protein binding free energy calculations]&#40;psf_dcd/protein_protein/README.md&#41;)

[//]: # (* [Protein-ligand binding free energy calculations]&#40;psf_dcd/protein_ligand/README.md&#41;)

[//]: # (* [Binding free energy calculations in multicomponent systems]&#40;psf_dcd/multicomponent_system/README.md&#41;)

 [^1]: It is part of the `All` set defined with `-t 0` in `gmx_MMPBSA_test`
 [^2]: It is part of the `Minimal` set defined with `-t 1` in `gmx_MMPBSA_test`
 [^3]: It is part of the `Fast` set defined with `-t 2` in `gmx_MMPBSA_test`
