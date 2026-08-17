Program to calculate pitch angles for databases of collagen EL databases

1. Compiling software

Run compileProg in this directory, it requires a fortran compiler and cmake.
Run ./compileProg -h for more information.
An executble will be built and copied into this directory.


2. Running executable

Two arguments can be provided (but will need to be ordered as below):
./PitchAngle [nmin [file name]]

nmin - number of minima, if not provided this will be obtained from min.data.
file name - only used if nmin is provided, default is coords.prmtop (AMBER topology)

Additionally, two input files are required:
- one specifying the atoms that specify the helical axis
- a second specifying the terminal residues (no peptide bonds are used between terminal residues)

Examples for both files are provided in example_input
