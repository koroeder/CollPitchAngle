MODULE PARSE_TOPOLOGY
   IMPLICIT NONE
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
    
   CONTAINS
      SUBROUTINE READ_TOP(TOPNAME)
         USE UTILITIES, ONLY: GETUNIT, READ_LINE
         IMPLICIT NONE
         CHARACTER(LEN=25), INTENT(IN) :: TOPNAME
         INTEGER :: TOPUNIT, J, IEND
         LOGICAL :: TOPEXISTS, CONTINUEPARSE, EOFT
         CHARACTER(LEN=100) :: LINE
         INTEGER, PARAMETER :: NWORDS=20
         CHARACTER(LEN=25) :: ENTRIES(NWORDS) = '', FLAG

         ! open file
         INQUIRE(FILE=TOPNAME, EXIST=TOPEXISTS)
         IF (.NOT.TOPEXISTS) THEN
            WRITE(*,*) " topology file ", TOPNAME, " not found - STOP"
            STOP
         END IF
         TOPUNIT = GETUNIT()
         OPEN(TOPUNIT, FILE=TOPNAME, STATUS='OLD')
         CONTINUEPARSE = .TRUE.

         DO WHILE (CONTINUEPARSE)
            !read the file line by line
            READ(TOPUNIT,'(A)', IOSTAT=IEND) LINE
            IF (IEND.LT.0) THEN
               CONTINUEPARSE = .FALSE.
               CYCLE
            ENDIF            
            CALL READ_LINE(LINE,NWORDS,ENTRIES)
            !use the %FLAG specifiers to find the correct input sections
            FLAG = ENTRIES(2)
            SELECT CASE (FLAG) 
               CASE("POINTERS")
                  !the next line is the format identifier - ignore
                  READ(TOPUNIT,*)
                  !first entry in the next line is the number of atoms
                  READ(TOPUNIT,'(A)') LINE
                  CALL READ_LINE(LINE,NWORDS,ENTRIES)
                  READ(ENTRIES(1),'(I8)') NATOMS
                  !the second entry in the next line is the number of residues
                  READ(TOPUNIT,'(A)') LINE
                  CALL READ_LINE(LINE,NWORDS,ENTRIES)               
                  READ(ENTRIES(2),'(I8)') NRES
                  ! call the allocator function
                  CALL ALLOC_TOP()
                  !we should now be ready to populate the remaining sections
               CASE("ATOM_NAME")
                  READ(TOPUNIT,'(20A4)') (ATOMNAMES(J), J=1,NATOMS)
               CASE("RESIDUE_LABEL")
                  READ(TOPUNIT,'(20A4)') (RESNAMES(J), J=1,NRES)
               CASE("RESIDUE_POINTER")
                  READ(TOPUNIT,'(12I6)') (FIRSTAT(J), J=1,NRES)
                  DO J=2,NRES
                     LASTAT(J-1) = FIRSTAT(J) - 1
                  END DO
                  LASTAT(NRES) = NATOMS
            END SELECT
         END DO
         CLOSE(TOPUNIT)

         ! create atomdata
         CALL POPULATE_ATOMDATA()
      END SUBROUTINE READ_TOP

      SUBROUTINE POPULATE_ATOMDATA()
         IMPLICIT NONE
         INTEGER :: I, J, FIRST, LAST, NATS

         DO I=1,NRES
            FIRST = FIRSTAT(I)
            LAST = LASTAT(I)
            NATS = LAST - FIRST + 1
            ATOMDATA(I)%NATS = NATS
            ALLOCATE(ATOMDATA(I)%ATOMIDS(NATS))
            ALLOCATE(ATOMDATA(I)%ATOMNAMES(NATS))
            DO J=1,NATS
               ATOMDATA(I)%ATOMIDS(J) = FIRST + J - 1
               ATOMDATA(I)%ATOMNAMES(J) = ATOMNAMES(FIRST + J - 1)
            END DO
         END DO

      END SUBROUTINE POPULATE_ATOMDATA

      SUBROUTINE ALLOC_TOP()
         CALL DEALLOC_TOP()
         ALLOCATE(RESNAMES(NRES))
         RESNAMES(1:NRES) = ''
         ALLOCATE(ATOMNAMES(NATOMS))
         ATOMNAMES(1:NATOMS) = ''
         ALLOCATE(ATOMDATA(NRES))
         ALLOCATE(FIRSTAT(NRES))
         FIRSTAT(1:NRES) = -1
         ALLOCATE(LASTAT(NRES))
         LASTAT(1:NRES) = -1
      END SUBROUTINE ALLOC_TOP

      SUBROUTINE DEALLOC_TOP()
         IF (ALLOCATED(RESNAMES)) DEALLOCATE(RESNAMES)
         IF (ALLOCATED(ATOMNAMES)) DEALLOCATE(ATOMNAMES)
         IF (ALLOCATED(ATOMDATA)) DEALLOCATE(ATOMDATA)
         IF (ALLOCATED(FIRSTAT)) DEALLOCATE(FIRSTAT)
         IF (ALLOCATED(LASTAT)) DEALLOCATE(LASTAT)
      END SUBROUTINE DEALLOC_TOP
         
END MODULE PARSE_TOPOLOGY