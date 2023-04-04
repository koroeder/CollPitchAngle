MODULE PARSE_TOPOLOGY
   USE COMMONS
   IMPLICIT NONE
   CONTAINS
      SUBROUTINE READ_TOP(TOPNAME)
         USE UTILITIES, ONLY: GETUNIT, READ_LINE
         IMPLICIT NONE
         CHARACTER(LEN=25), INTENT(IN) :: TOPNAME
         INTEGER :: TOPUNIT, J, IEND
         LOGICAL :: TOPEXISTS, CONTINUEPARSE
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
                  READ(TOPUNIT,*)
                  READ(TOPUNIT,'(20A4)') (ATOMNAMES(J), J=1,NATOMS)
               CASE("RESIDUE_LABEL")
                  READ(TOPUNIT,*)
                  READ(TOPUNIT,'(20A4)') (RESNAMES(J), J=1,NRES)
               CASE("RESIDUE_POINTER")
                  READ(TOPUNIT,*)
                  READ(TOPUNIT,'(10I8)') (FIRSTAT(J), J=1,NRES)
                  DO J=2,NRES
                     LASTAT(J-1) = FIRSTAT(J) - 1
                  END DO
                  LASTAT(NRES) = NATOMS
            END SELECT
         END DO
         CLOSE(TOPUNIT)

         ! create atomdata
         CALL POPULATE_ATOMDATA()
         ! find all molecules
         IF (.NOT.ALLOCATED(TERMINI)) CALL FIND_TERMINI()
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

      SUBROUTINE FIND_TERMINI()
         IMPLICIT NONE
         INTEGER :: TMPTER(NRES)
         INTEGER :: I, J

         NTERMINI = 0
         TMPTER(1:NRES) = -1
         DO I=1,NATOMS
            IF (ATOMNAMES(I).EQ."OXT") THEN
               DO J=1,NRES
                  IF ((FIRSTAT(J).LE.I).AND.(LASTAT(J).GE.I)) THEN
                     NTERMINI = NTERMINI + 1
                     TMPTER(NTERMINI) = J
                     EXIT
                  END IF
               END DO
            END IF
         END DO
         ALLOCATE(TERMINI(NTERMINI,2))
         DO I=1,NTERMINI
            TERMINI(I,2) = TMPTER(I)
            IF (I.EQ.1) THEN
               TERMINI(I,1) = 1
            ELSE
               TERMINI(I,1) = TERMINI(I-1,2) + 1
            END IF
         END DO
      END SUBROUTINE

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