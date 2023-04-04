MODULE COMMONS
   !> real, double precision 64-bit
   INTEGER, PARAMETER  :: REAL64 = SELECTED_REAL_KIND(15, 307)  
   ! number of residues
   INTEGER :: NRES = 0
   ! residue name
   CHARACTER(LEN=4), ALLOCATABLE :: RESNAMES(:)
   ! number of atoms
   INTEGER :: NATOMS = 0
   ! atom names
   CHARACTER(LEN=4), ALLOCATABLE :: ATOMNAMES(:)
   ! type for atoms in residue
   TYPE RESIDUE
      INTEGER :: NATS = 0
      INTEGER, ALLOCATABLE :: ATOMIDS(:)
      CHARACTER(LEN=4), ALLOCATABLE :: ATOMNAMES(:)
   END TYPE RESIDUE
   TYPE(RESIDUE), ALLOCATABLE :: ATOMDATA(:)
   INTEGER, ALLOCATABLE :: FIRSTAT(:)
   INTEGER, ALLOCATABLE :: LASTAT(:)

   !> atoms for the helical axis
   INTEGER :: NPOINT1, NPOINT2
   INTEGER, ALLOCATABLE :: ATSFOR1(:), ATSFOR2(:)

   !> number of peptide bonds
   INTEGER :: NPEPBONDS = 0
   !> identities of the atoms in these bonds
   INTEGER, ALLOCATABLE :: PEPBONDS(:,:)
   
   !> number of termini
   INTEGER :: NTERMINI = 0
   INTEGER, ALLOCATABLE :: TERMINI(:,:)

   !> number of methylene groups
   INTEGER :: NMETHYLENE = 0
   INTEGER, ALLOCATABLE :: METHYLENEGROUPS(:,:)

   !> pitch angle for peptides
   REAL(KIND=REAL64), ALLOCATABLE :: PEPANGLES(:,:)
   !> pitch angles for methylene groups
   REAL(KIND=REAL64), ALLOCATABLE :: METHYLENEANGLES(:,:)

END MODULE COMMONS