## **`&rism` namelist variables**

!!! note "Keep in mind"
    * A default 3drism input file can be created as follows:

        === "GROMACS"
            ```
            xbfree gmx_MMPBSA --create_input rism
            ```
        === "AMBER"
            ```
            xbfree amber_MMPBSA --create_input rism
            ```
        === "NAMD"
            ```
            xbfree namd_MMPBSA --create_input rism
            ```
        === "CHARMM"
            ```
            xbfree charmm_MMPBSA --create_input rism
            ```
    
    * `3D-RISM` calculations are performed with the `rism3d.snglpnt` program built with AmberTools, written by Tyler 
    Luchko. It is the most expensive, yet most statistical mechanically rigorous solvation model. See 
        * [Introduction to RISM](https://ambermd.org/doc12/Amber21.pdf#section.7.1) for a thorough description RISM 
        theory.
        * [General workflow for using 3D-RISM](https://ambermd.org/doc12/Amber21.pdf#section.7.3)
        * Practical considerations on:
            * [Computational Requirements and Parallel Scaling of RISM](https://ambermd.org/doc12/Amber21.pdf#subsection.7.2.1)
            * [Numerical Accuracy of RISM](https://ambermd.org/doc12/Amber21.pdf#subsection.7.2.3)
            * [Convergence issues](https://ambermd.org/doc12/Amber21.pdf#subsection.7.3.1)
    * A sample 3drism input file is shown [here](input_file.md#mm3d-rism)
    * A tutorial on binding free energy calculation with 3D-RISM is available [here](examples/3D-RISM/README.md)
    * We have included more variables in 3D-RISM calculations than the ones available in the MMPBSA.py original code. 
    That way, users can be more in control and tackle various issues (_e.g._, convergence issues).
    * One advantage of `3D-RISM` is that an arbitrary solvent can be chosen; you just need to change the `xvvfile` 
    specified on the command line (see `-xvvfile` flag in [gmx_MMPBSA command line](gmx_MMPBSA_command-line.md). The 
    default solvent is `$AMBERHOME/AmberTools/test/rism1d/tip3p-kh/tip3p.xvv.save`. In case this file 
    doesn't exist, a copy `path_to_GMXMMPBSA/data/xvv_files/tip3p.xvv` is used. You can find examples of precomputed 
    `.xvv` files for SPC/E and TIP3P water in `$AMBERHOME/AmberTools/test/rism1d` or 
    `path_to_GMXMMPBSA/data/xvv_files` folders.

  [7]: https://ambermd.org/doc12/Amber21.pdf#chapter.7
  [8]: https://ambermd.org/doc12/Amber21.pdf#subsection.36.3.2

[`xvv`](#mmpbsa_ifv_xvv){#mmpbsa_ifv_xvv} (Default = "tip3p")
:   Define the selected solvent for 3D-RISM. These solvent xvv files are contained in xBFreE, but you can define a 
new one simply adding the file path.

    * tip3p
    * spc
    * spc-nacl-3
    * spc_mmpbsa_py


### **Closure approximations**

[`closure`](#mmpbsa_ifv_closure){#mmpbsa_ifv_closure} (Default = "kh")
:   Comma separate list of closure approximations. If more than one closure is provided, the 3D-RISM solver will use 
the closures in order to obtain a solution for the last closure in the list when no previous solutions are available.
The solution for the last closure in the list is used for all output. The use of several closures combined with 
different tolerances can be useful to overcome convergence issues (see [§7.3.1](https://ambermd.org/doc12/Amber21.
pdf#subsection.7.3.1))

    * "kh": Kovalenko-Hirata
    * "hnc": Hyper-netted chain equation
    * "psen": Partial Series Expansion of order-n where “n” is a positive integer (_e.g._, "pse3")

    !!! example "Examples"
        === "One closure"
                 closure = pse3
        === "Several closures"
                 closure = kh, pse3

### **Solvation free energy corrections**

[`gfcorrection`](#mmpbsa_ifv_gfcorrection){#mmpbsa_ifv_gfcorrection} (Default = 0)
:    Compute the Gaussian fluctuation excess chemical potential functional. 
See [§7.1.2](https://ambermd.org/doc12/Amber21.pdf#subsection.7.1.2)

    * 0: Off
    * 1: On

[`pcpluscorrection`](#mmpbsa_ifv_pcpluscorrection){#mmpbsa_ifv_pcpluscorrection} (Default = 0)
:    Compute the PC+/3D-RISM excess chemical potential functional.
See [§7.2.4](https://ambermd.org/doc12/Amber21.pdf#subsection.7.2.4)

    * 0: Off
    * 1: On

### **Long-range asymptotics**

!!! info
    Long-range asymptotics are used to analytically account for solvent distribution beyond the solvent box. 
    Long-range asymptotics are always used when calculating a solution but can be omitted for
    the subsequent thermodynamic calculations, though it is not recommended.

[`noasympcorr`](#mmpbsa_ifv_noasympcorr){#mmpbsa_ifv_noasympcorr} (Default = 1) 
:   Use long-range asymptotic corrections for thermodynamic calculations.

    * 0: Do not use long-range corrections
    * 1: Use the long-range corrections


[`treeDCF`](#mmpbsa_ifv_treeDCF){#mmpbsa_ifv_treeDCF} (Default = 1)
:   Use direct sum, or the treecode approximation to calculate the direct correlation function long-range asymptotic 
correction.

    * 0: Use direct sum
    * 1: Use treecode approximation

[`treeTCF`](#mmpbsa_ifv_treeTCF){#mmpbsa_ifv_treeTCF} (Default = 1)
:   Use direct sum, or the treecode approximation to calculate the total correlation function long-range asymptotic 
correction.

    * 0: Use direct sum
    * 1: Use treecode approximation

[`treeCoulomb`](#mmpbsa_ifv_treeCoulomb){#mmpbsa_ifv_treeCoulomb} (Default = 1)
:   Use direct sum, or the treecode approximation to calculate the Coulomb potential energy.

    * 0: Use direct sum
    * 1: Use treecode approximation


[`treeDCFMAC`](#mmpbsa_ifv_treeDCFMAC){#mmpbsa_ifv_treeDCFMAC} (Default = 0.1)
:   Treecode multipole acceptance criterion for the direct correlation function long-range asymptotic correction.


[`treeTCFMAC`](#mmpbsa_ifv_treeTCFMAC){#mmpbsa_ifv_treeTCFMAC} (Default = 0.1)
:   Treecode multipole acceptance criterion for the total correlation function long-range asymptotic correction.


[`treeCoulombMAC`](#mmpbsa_ifv_treeCoulombMAC){#mmpbsa_ifv_treeCoulombMAC} (Default = 0.1)
:   Treecode multipole acceptance criterion for the Coulomb potential energy.


[`treeDCFOrder`](#mmpbsa_ifv_treeDCFOrder){#mmpbsa_ifv_treeDCFOrder} (Default = 2)
:   Treecode Taylor series order for the direct correlation function long-range asymptotic correction.


[`treeTCFOrder`](#mmpbsa_ifv_treeTCFOrder){#mmpbsa_ifv_treeTCFOrder} (Default = 2)
:   Treecode Taylor series order for the total correlation function long-range asymptotic correction. Note that the 
Taylor expansion used does not converge exactly to the TCF long-range asymptotic correction, so a very high order 
will not necessarily increase accuracy.


[`treeCoulombOrder`](#mmpbsa_ifv_treeCoulombOrder){#mmpbsa_ifv_treeCoulombOrder} (Default = 2)
:   Treecode Taylor series order for the Coulomb potential energy.


[`treeDCFN0`](#mmpbsa_ifv_treeDCFN0){#mmpbsa_ifv_treeDCFN0} (Default = 500)
:   Maximum number of grid points contained within the treecode leaf clusters for the direct correlation function 
long-range asymptotic correction. This sets the depth of the hierarchical octtree.


[`treeTCFN0`](#mmpbsa_ifv_treeTCFN0){#mmpbsa_ifv_treeTCFN0} (Default = 500)
:   Maximum number of grid points contained within the treecode leaf clusters for the total correlation function 
long-range asymptotic correction. This sets the depth of the hierarchical octtree.


[`treeCoulombN0`](#mmpbsa_ifv_treeCoulombN0){#mmpbsa_ifv_treeCoulombN0} (Default = 500)
:   Maximum number of grid points contained within the treecode leaf clusters for the Coulomb potential energy. This 
sets the depth of the hierarchical octtree.


### **Solvation box**

!!! info
    The non-periodic solvation box super-cell can be defined as variable or fixed in size. When a
    variable box size is used, the box size will be adjusted to maintain a minimum buffer distance between the atoms
    of the solute and the box boundary. This has the advantage of maintaining the smallest possible box size while
    adapting to change of solute shape and orientation. Alternatively, the box size can be specified at run-time. This
    box size will be used for the duration of the sander calculation. Solvent box dimensions have a strong effect on 
    the numerical precision of 3D-RISM. See [§7.2.3](https://ambermd.org/doc12/Amber21.pdf#subsection.7.2.3) for 
    recommendation on selecting an appropriate box size and resolution.

##### **Variable box size**

!!! info "Keep in mind"
    It is recommended to avoid specifying a large, prime number of processes (≥ 7) when using a variable solvation 
    box size.

[`buffer`](#mmpbsa_ifv_buffer){#mmpbsa_ifv_buffer} (Default = 14)
:   Minimum distance (in Å) between solute and edge of solvation box. Specify this with `grdspc` below. Mutually 
exclusive with `ng` and `solvbox`. See [§7.2.3](https://ambermd.org/doc12/Amber21.pdf#subsection.7.2.3) for details on 
how this affects numerical accuracy and how this interacts with `ljTolerance`, and `tolerance`

    * when < 0: Use fixed box size (see `ng` and `solvbox` below)
    * when >= 0: Use `buffer` distance

[`grdspc`](#mmpbsa_ifv_grdspc){#mmpbsa_ifv_grdspc}(Default = 0.5,0.5,0.5)
:   Grid spacing (in Å) of the solvation box. Specify this with `buffer` above. Mutually exclusive with `ng` and 
`solvbox`.

##### **Fixed box size**

[`ng`](#mmpbsa_ifv_ng){#mmpbsa_ifv_ng} (Default = -1,-1,-1)
:   Comma separated number of grid points to use in the x, y, and z directions. Used only if buffer < 0. Mutually 
exclusive with `buffer` and `grdspc` above, and paired with `solvbox` below.

    !!! warning 
        No default, this must be set if buffer < 0. As a general requirement, the number of grids points in each 
        dimension must be divisible by two, and the number of grid points in the z-axis must be divisible by the 
        number of processes.

        As an example: define like `ng=1000,1000,1000`, where all numbers are divisible by two 
        and you can use 1, 2, 4, 5, 8, 10... pocessors, all divisors of 1000 (value in the z-axis).

        Take into account that at a certain level, running RISM in 
        parallel may actually hurt performance, since previous solutions are used 
        as an initial guess for the next frame, hastening convergence. Running in parallel loses this advantage. Also, 
        due to the overhead involved in which each thread is required to load every topology file when calculating 
        energies, parallel scaling will begin to fall off as the number of threads reaches the number of frames. 

[`solvbox`](#mmpbsa_ifv_solvbox){#mmpbsa_ifv_solvbox} (Default = -1,-1,-1)
:    Sets the size in Å of the fixed size solvation box. Used only if `buffer` < 0. Mutually exclusive with `buffer` 
and `grdspc` above, and paired with `ng` above. 

    !!! warning 
        No default, this must be set if buffer < 0. Define like `solvbox=20,20,20`

[`solvcut`](#mmpbsa_ifv_solvcut){#mmpbsa_ifv_solvcut}  (Default = 14)
:   Cutoff used for solute-solvent interactions. The default value is that of buffer. Therefore, if you set `buffer` < 
0 and specify `ng` and `solvbox` instead, you must set `solvcut` to a nonzero value; otherwise the program will quit in 
error.

### **Solution convergence**

[`tolerance`](#mmpbsa_ifv_tolerance){#mmpbsa_ifv_tolerance} (Default = 0.00001)
:   A comma-separated list of maximum residual values for solution convergence. This has a strong effect on the 
cost of 3D-RISM calculations (smaller value for tolerance -> more computation). When used in combination with a list 
of closures it is possible to define different tolerances for each of the closures. This can be useful for difficult 
to converge calculations (see [§7.4.1](https://ambermd.org/doc12/Amber21.pdf#page=120&zoom=100,96,798)). For the sake of 
efficiency, it is best to use as high a tolerance as possible for all but the last closure. 
See [§7.2.3](https://ambermd.org/doc12/Amber21.pdf#subsection.7.2.3) for details on how this affects numerical 
accuracy and how this interacts with `ljTolerance`, `buffer`, and `solvbox`. Three formats of list are possible:

    * one tolerance: All closures but the last use a tolerance of 1. The last tolerance in the list is used
    by the last closure. In practice this is the most efficient.
    * two tolerances: All closures but the last use the first tolerance in the list. The last tolerance in the
    list is used by the last closure.
    * n tolerances: Tolerances from the list are assigned to the closure list in order.

    !!! example "Examples"
        === "One closure/One tolerance"
                closure = pse3, tolerance=0.00001
            
            A tolerance of `0.00001` will be used for clousure `pse3`
        === "Several closures/One tolerance"
                 closure = kh, pse3, tolerance=0.00001

            A tolerance of `1` will be used for clousure `kh`, while `0.00001` will be used for clousure `pse3`. 
            Equivalent to `closure = kh, pse3, tolerance=1,0.00001`
        === "Several closures/Two tolerances"
                 closure = kh, pse2, pse3, tolerance=0.01, 0.00001

            A tolerance of `0.01` will be used for clousures `kh` and `pse2`, while `0.00001` will be used for clousure 
            `pse3`. Equivalent to `closure = kh, pse2, pse3, tolerance=0.01,0.01,0.00001`
        === "Several closures/Several tolerances"
                 closure = kh,pse2,pse3, tolerance=0.1,0.01,0.00001

            A tolerance of `0.1` will be used for clousure `kh`, `0.01` will be used for clousure `pse2`, while `0.00001` 
            will be used for clousure `pse3`.


[`ljTolerance`](#mmpbsa_ifv_ljTolerance){#mmpbsa_ifv_ljTolerance} (Default = -1)
:   Lennard-Jones accuracy (Optional.) Determines the Lennard-Jones cutoff distance based on the desired accuracy of 
the calculation. See [§7.2.3](https://ambermd.org/doc12/Amber21.pdf#subsection.7.2.3) for details on how this affects 
numerical accuracy and how this interacts with `tolerance`, `buffer`, and `solvbox`.


[`asympKSpaceTolerance`](#mmpbsa_ifv_asympKSpaceTolerance){#mmpbsa_ifv_asympKSpaceTolerance} (Default = -1)
:   Tolerance reciprocal space long range asymptotics accuracy (Optional.) Determines the reciprocal space long 
range asymptotic cutoff distance based on the desired accuracy of the calculation. 
See [§7.2.3](https://ambermd.org/doc12/Amber21.pdf#subsection.7.2.3) for details on how this affects numerical 
accuracy. Possible values are:

    * `when < 0`: asympKSpaceTolerance = tolerance/10
    * `when = 0`: no cutoff
    * `when > 0`: given value determines the maximum error in the reciprocal-space long range asymptotics calculations


[`mdiis_del`](#mmpbsa_ifv_mdiis_del){#mmpbsa_ifv_mdiis_del} (Default = 0.7)
:   MDIIS step size.


[`mdiis_nvec`](#mmpbsa_ifv_mdiis_nvec){#mmpbsa_ifv_mdiis_nvec} (Default = 5)
:   Number of previous iterations MDIIS uses to predict a new solution.


[`mdiis_restart`](#mmpbsa_ifv_mdiis_restart){#mmpbsa_ifv_mdiis_restart} (Default = 10)
:   If the current residual is mdiis_restart times larger than the smallest residual in memory, then the MDIIS 
procedure is restarted using the lowest residual solution stored in memory. Increasing this number can sometimes 
help convergence.


[`maxstep`](#mmpbsa_ifv_maxstep){#mmpbsa_ifv_maxstep} (Default = 10000)
:   Maximum number of iterations allowed to converge on a solution.


[`npropagate`](#mmpbsa_ifv_npropagate){#mmpbsa_ifv_npropagate} (Default = 5)
:   Number of previous solutions propagated forward to create an initial guess for this solute atom configuration.

    * =0: Do not use any previous solutions
    * = 1..5: Values greater than 0 but less than 4 or 5 will use less system memory but may introduce artifacts to 
    the solution (_e.g._, energy drift).


### **Output**

[`polardecomp`](#mmpbsa_ifv_polardecomp){#mmpbsa_ifv_polardecomp} (Default = 0)
:   Decomposes solvation free energy into polar and non-polar components. Note that this typically requires 80% more 
computation time.

    * 0: Don’t decompose solvation free energy into polar and non-polar components. 
    * 1: Decompose solvation free energy into polar and non-polar components.

[`entropicdecomp`](#mmpbsa_ifv_entropicdecomp){#mmpbsa_ifv_entropicdecomp} (Default = 0)
:   Decomposes solvation free energy into energy and entropy components. Also performs temperature derivatives of other 
calculated quantities. Note that this typically requires 80% more computation time and requires a `.xvv` file version 
1.000 or higher (available within `GMXMMPBSA` data folder). 
See [§7.1.3](https://ambermd.org/doc12/Amber21.pdf#subsection.7.1.3) and 
[§7.3](https://ambermd.org/doc12/Amber21.pdf#section.7.3)

    * 0: No entropic decomposition
    * 1: Entropic decomposition

[`rism_verbose`](#mmpbsa_ifv_rism_verbose){#mmpbsa_ifv_rism_verbose} (Default = 0)
:   Level of output in temporary RISM output files. May be helpful for debugging or following convergence. 

    * 0: just print the final result
    * 1: additionally prints the total number of iterations for each solution
    * 2: additionally prints the residual for each iteration and details of the MDIIS solver (useful for debugging 
    and convergence analyses)
