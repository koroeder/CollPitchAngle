PROGRAM CALCPITCH
   USE COMMONS
   USE UTILITIES, ONLY: GETUNIT
   USE CALCPROPERTIES, ONLY: GET_HELAX, DETERMINE_PEPTIDE, FIND_METHYLENE
   IMPLICIT NONE

   !> name of topology file 
   CHARACTER(LEN=25) :: TOPFILE="coords.prmtop"
   !> number of minima to be parsed
   INTEGER :: NMIN = 0

   !> current coordinates
   REAL(KIND=REAL64), ALLOCATABLE :: CURRX(:)
   !> unit for points.min
   INTEGER :: XUNIT
   INTEGER ::I, J

   ! initialise program
   CALL INIT_PROG()
   !initialise molecule numbers from files
   CALL GET_TERMINI()
   ! parse topology
   CALL READ_TOP(TOPFILE)
   !get atom ids for helical axis
   CALL ATOMSFORAXIS()
   ! get all peptide bonds
   CALL DETERMINE_PEPTIDE()
   ! get all methylene groups
   CALL FIND_METHYLENE()
   !allocate coordinate array
   ALLOCATE(CURRX(3*NATOMS))
   ! allocate arrays for storing data
   ALLOCATE(PEPANGLES(NMIN,NPEPBONDS))
   ALLOCATE(METHYLENEANGLES(NMIN,NMETHYLENE))

   ! iterate over all structures
   XUNIT = GETUNIT()
   OPEN(XUNIT,FILE='points.min',ACCESS='DIRECT',ACTION='READ',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)
   DO I = 1,NMIN
      READ(XUNIT,REC=I) (CURRX(J),J=1,3*NATOMS)
      CALL GET_HELAX(CURRX)
      CALL PEPTIDE_PITCH(I,CURRX)
      CALL METHYLENE_PITCH(I,CURRX)
   END DO
   CLOSE(XUNIT)

   CONTAINS
      SUBROUTINE INIT_PROG()
         USE UTILITIES, ONLY: GETUNIT
         IMPLICIT NONE

         INTEGER :: NARGS, TMPUNIT
         LOGICAL :: MINDFOUND
         CHARACTER(LEN=10) :: NMINDUMMY

         WRITE(*,*) "Program to calculate pitch angles in collagen for EL databases"
         WRITE(*,*) "Written by Dr K Roeder at KCL"
         WRITE(*,*) ""

         NARGS = COMMAND_ARGUMENT_COUNT()

         IF (NARGS.EQ.0) THEN
            WRITE(*,*) " parse_args> No command line arguments - determine number of minima from min.data"
            INQUIRE(FILE="min.data",EXIST=MINDFOUND)
            IF (.NOT.MINDFOUND) THEN
               WRITE(*,*) " parse_args> Cannot locate min.data - STOP"
               STOP
            END IF
            CALL EXECUTE_COMMAND_LINE("wc -l min.data > tmp_nmins")
            TMPUNIT = GETUNIT()
            OPEN(TMPUNIT, FILE="tmp_nmins", STATUS='OLD')
            READ(TMPUNIT,'(I8)') NMIN
            CALL EXECUTE_COMMAND_LINE("rm tmp_nmins")
            IF (NMIN.EQ.0) THEN
               WRITE(*,*) " parse_args> Number of minima is zero - something went wrong - STOP"
               STOP
            END IF
            WRITE(*,*) " parse_args> Number of minima to be parsed: ", NMIN
            WRITE(*,*) " parse_args> Using default name for topology"
         ELSE
            CALL GET_COMMAND_ARGUMENT(1,NMINDUMMY)
            READ(NMINDUMMY,'(I8)') NMIN
            WRITE(*,*) " parse_args> Number of minima to be parsed: ", NMIN
            IF (NARGS.GT.1) THEN
               CALL GET_COMMAND_ARGUMENT(2,TOPFILE)
               WRITE(*,*) " parse_args> Using topology ", TOPFILE
            ELSE
               WRITE(*,*) " parse_args> Using default name for topology"
            END IF
         END IF

      END SUBROUTINE INIT_PROG

      SUBROUTINE ATOMSFORAXIS()
         USE UTILITIES, ONLY: GETUNIT
         IMPLICIT NONE
         LOGICAL :: INPFOUND
         CHARACTER(LEN=25) :: INPUTFILE="atoms_helax.dat"
         INTEGER ::INUNIT, I

         INQUIRE(FILE=INPUTFILE,EXIST=INPFOUND)
         IF (.NOT.INPFOUND) THEN
            WRITE(*,*) " atomsforax> Cannot locate ", INPUTFILE, " - STOP"
            STOP
         END IF
         INUNIT = GETUNIT()
         OPEN(INUNIT, FILE=INPUTFILE, STATUS='OLD')
         READ(INUNIT,'(I6)') NPOINT1
         ALLOCATE(ATSFOR1(NPOINT1))
         READ(INUNIT,'(I6)') (ATSFOR1(I), I=1,NPOINT1)
         READ(INUNIT,'(I6)') NPOINT2       
         ALLOCATE(ATSFOR2(NPOINT2))
         READ(INUNIT,'(I6)') (ATSFOR2(I), I=1,NPOINT2)
         CLOSE(INUNIT)
      END SUBROUTINE ATOMSFORAXIS

      SUBROUTINE  GET_TERMINI()
         USE UTILITIES, ONLY: GETUNIT
         IMPLICIT NONE
         LOGICAL :: INPFOUND
         CHARACTER(LEN=25) :: INPUTFILE="termini.dat"
         INTEGER ::INUNIT, I

         INQUIRE(FILE=INPUTFILE,EXIST=INPFOUND)
         IF (.NOT.INPFOUND) THEN
            WRITE(*,*) " get_termini> Cannot locate ", INPUTFILE, " - will attmept to get termini from topology"
            RETURN
         END IF
         INUNIT = GETUNIT()
         OPEN(INUNIT, FILE=INPUTFILE, STATUS='OLD')
         READ(INUNIT,'(I6)') NTERMINI        
         ALLOCATE(TERMINI(NTERMINI,2))
         READ(INUNIT,'(I6)') (TERMINI(I,1), TERMINI(I,2), I=1,NTERMINI)         
         CLOSE(INUNIT)
      END SUBROUTINE GET_TERMINI

END PROGRAM CALCPITCH