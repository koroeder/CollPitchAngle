PROGRAM CALCPITCH
   USE PARSE_TOPOLOGY, ONLY: NATOMS, NRES, ATOMDATA, RESNAMES, READ_TOP
   USE UTILITIES, ONLY: GETUNIT
   IMPLICIT NONE
   !> real, double precision 64-bit
   INTEGER, PARAMETER  :: REAL64 = SELECTED_REAL_KIND(15, 307)

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
   ! parse topology
   CALL READ_TOP(TOPFILE)
   !allocate coordinate array
   ALLOCATE(CURRX(3*NATOMS))

   ! iterate over all structures
   XUNIT = GETUNIT()
   OPEN(XUNIT,FILE='points.min',ACCESS='DIRECT',ACTION='READ',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)
   DO I = 1,NMIN
      READ(XUNIT,REC=I) (CURRX(J),J=1,3*NATOMS)
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

END PROGRAM CALCPITCH