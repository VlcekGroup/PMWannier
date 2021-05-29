# Manual for PMwannier package

This code performs unitary transformations on canonical orbitals to obtain a set of localized states. It adopts the Pipek-Mezey principle in the evaluation of the localization and employs steepest descent algorithm for the maximization of the objective function. It is primarily designed for the localization of occupied states but also available for virtual states if a specific virtual space is provided. For large complex ssystems, PWwannier code is able to search for the localized states on atoms (molecules) of interest, regardless of the localization of the rest. Further, the code can perform unitary transformations in a truncated subspace in order to lower the computational cost. In addition to the main program, tools in this package help determining the localized states on the selected atoms (molecules) and generate reconstructed orbitals. 

## Installation
In each folder of the source code, simply hit "make" to install the main code (wannier.x) and other tools. Note that for the main code, openmpi library is needed.

## Wannier.x
### Command to run a job
~/path_to_executable_file/wannier.x >& output_file_name &

### Files needed to run a job
(1) cnt.ini : provides the coordinates of atoms.\
(2) wf.bin : provides the canonical orbitals as the starting point.\
(3) input : specifies the parameters for the calculation, see Input Specifications. The name of this file must be "input".\
(4) wannier.bin : this file is needed if the job is restarted from the previous one.\

### Output files after running a job
(1) wannier.bin : contains all the unitary-transformed states in the occupied/specified subspace.\
(2) wannier_trun.bin : a truncated version of the wannier.bin, which contains wannier functions localized on the specified atoms only.\
(3) wannier_unocc.bin : contains all the unoccupied states used in the wannierization.\
(4) wannier_unocc_trun.bin : a truncated version of the wannier_unooc.bin with virtual wannierzation functions localized on the specified atoms only.\
(5) wannier_tmp.bin : a temporary wannier function bin file that is updated every iteration when the objective function change is positive.\
(6) wannier_sv.bin : this bin file will be generated when subspace_wannierization is true and restart_old is false. This will save the wannier functions before the subspace calculation is invoked. It helps optimizing the subspace parameters.\
(6) Output: by default the output information is written to screen. File with the name "Output" contains some of the input information, the objective function value at each iteration, the time spent at each step, and the delta_t value at each step.\

### Input Specification
#### Basic function
num_of_atoms\
(Integer. number of atoms. Default: calculated from cnt.ini file)\
\
num_of_occ\
(Integer. number of occupied states. Default: calculated from cnt.ini file)\
\
restart\
(Logical. restart the wannierization by reading the wannier.bin file. Default: F)\
\
max_iter\
(Integer. maximum iteration steps. Default: 50)\
\
 delta_t\
(Real. delta t parameter in the steepest descent process. Default: 1)\
\
delta_t_correction\
(Real. this will be triggered when the objective function change is negative. Delta_t = delta_t/corr_fac. Default: 1.1)\
\
num_of_threads\
(Integer. number of threads for MPI calculations. Default: 1)\
\
convergence_threshold\
(Real. This specifies the threshold to exit the iteration when the objective function change is less than the value. Default: 1E-8)\
\
#### Local wannierzation
local_wannierization\
(Logical. This is to perform wannierization on specified atoms. Default: F, If true, then in the following three lines you will need to specify: number of atoms for local wannierization, number of wannier functions of the specified atoms, the labels of the specified atoms)\
\
step_findwf\
(Integer. This defines the step interval of printing out the wannier functions localized on the specified atoms when using local wannierization. Default: 10)\

#### Localization of Virtual Orbitals
unocc_state\
(Logical. This is to localize the virtual orbitals. Default: F)\
\
num_unocc_lw\
(Integer. This is to specify the number of virtual wannier functions of interest on the specified atoms in local wannierzation)\
\
bulk\
(Logical. This is to obtain virtual wannier function together with local wannierzation, i.e., obtain virtual wannier functions on specified atoms. Default: F)\
\
num_of_unocc\
(Integer. If bulk is F, then you can choose to provide the number of virtual orbitals for the wannierization. Otherwise the program will take all the bound virtual states for wannierization)\

#### Subspace Wannerization
subspace_wannierzation\
(Logical. This is to perform wannierzation using an occupied subspace of KS states. Default: F)\
\
subspace_factor\
(Integer. This is to define the size of the occupied subspace. Suppose you are looking for N wannier functions localized on your specified atoms, then the size of the subspace will be N*sub_fac. Default: 3)\
\
subspace_threshold\
(Real. This defines when to start the subspace wannierzation. When the gain in the objective function is greater than this value, the subspace wannierzation will be triggered. Otherwise it will keep running using the full occupied space. Default: 1d0)\
\
restart_old\
(Logical. This is to connect the full space wannierization and the subspace wannierzation. Default: F. If this is true, it will read the wannier.bin file and continue with subspace wannierzation. Note that “restart” has to be true as well if this one is true.)\
\
specific_state\
(Logical. This is to perform wannierization using the KS states as specified in the following lines only. Default: F. If this is true, then in the next line provide the number of states used in the wannierzation and then followed by the state indices in one line.)\

## Post-PMwanner codes (tools)

### findWF.x
This is to generate wannier_trun.bin file from wannier.bin (or wannier_tmp.bin) if the job is terminated manually before all the iterations are finished.
#### Command to run a job
~/path_to_executable_file/findWF.x bin_file_name >& output_file_name &
#### Files needed to run a job
(1) findWF : this is an input file providing the number of atoms, number of wannier functions wanted, and the atom labels. See the specifications.\
(2) wannier.bin : defined as before.\
(3) cnt.ini : defined as before.\
#### Output files after running a job
(1) wannier_trun.bin : defined as before.\
(2) wannier_unocc_trun.bin : defined as before.\
(3) Output file : contains the information of the job.\
#### Input Specifications
num_occ_lw\
(Integer. This specifies the number of wannier functions to find on the specified atoms)\

num_atom_lw\
(Integer. This specifies the number of atoms of interest. In the next line provide the atom labels)\
\
num_of_threads\
(Integer. This specifies the number of cores for the MPI process. Default: 10)\
\
unocc_state\
(Logical. This is to get virtual wannier functions on the specified atoms.)\
\
num_unocc_lw\
(Integer. This specifies the number of virtual wannier functions of interest on the specified atoms.)\

### mkorb.x
This is to make reconstrcuted orbitals from the wannier functions of a composite system (atoms/molecules of interest plus the rest).\
#### Files needed to run a job
(1) wannier_trun.bin : defined as before.\
(1) wannier_mol.bin : contains the wannier functions of the isolated molecules.\
(2) wf_mol.bin : contains the canonical orbitals of the isolated molecules.\
#### Output files after running a job
(1) mdorb00X.bin : X represents the state index. This is a binary file with similar structure to the wf.bin.\
(2) orb00X.bin: X represent the state index. This is formatted file used for many-body calculation.\


