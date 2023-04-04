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
       SUBROUTINE READ_TOP()
          ! read the first few lines to get the number of residues and atoms

          ! call the allocator function
          CALL ALLOC_TOP()

          ! read relevant entries


          ! create atomdata

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
         
        SUBROUTINE READ_LINE(LINE,NWORDS,WORDSOUT)
            CHARACTER(*), INTENT(IN) :: LINE
            INTEGER, INTENT(IN) :: NWORDS
            CHARACTER(*), DIMENSION(NWORDS), INTENT(OUT) :: WORDSOUT
            INTEGER:: J1,START_IND,END_IND,J2
            CHARACTER(25) :: WORD
            START_IND=0
            END_IND=0
            J1=1
            J2=0
            DO WHILE(J1.LE.LEN(LINE))
                IF ((START_IND.EQ.0).AND.(LINE(J1:J1).NE.' ')) THEN
                START_IND=J1
                ENDIF
                IF (START_IND.GT.0) THEN
                IF (LINE(J1:J1).EQ.' ') END_IND=J1-1
                IF (J1.EQ.LEN(LINE)) END_IND=J1
                IF (END_IND.GT.0) THEN
                    J2=J2+1
                    WORD=LINE(START_IND:END_IND)
                    WORDSOUT(J2)=TRIM(WORD)
                    START_IND=0
                    END_IND=0
                ENDIF
                ENDIF
                J1=J1+1
            ENDDO
        END SUBROUTINE READ_LINE
END MODULE PARSE_TOPOLOGY