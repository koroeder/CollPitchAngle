MODULE UTILITIES
   
   CONTAINS

      SUBROUTINE LOOKUP_ATOM(RESID,ATOM,ATOMID)
         USE COMMONS, ONLY: ATOMDATA
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: RESID
         CHARACTER(*), INTENT(IN) :: ATOM
         INTEGER, INTENT(OUT) :: ATOMID
         INTEGER :: NATS, I

         ATOMID = -1
         NATS = ATOMDATA(RESID)%NATS
         DO I=1,NATS
            IF (ATOM.EQ.ATOMDATA(RESID)%ATOMNAMES(I)) THEN
               ATOMID = ATOMDATA(RESID)%ATOMIDS(I)
               RETURN
            END IF
         ENDDO

      END SUBROUTINE LOOKUP_ATOM

      !> Function to find empty unit for output/input file
      INTEGER FUNCTION GETUNIT()
         IMPLICIT NONE
         LOGICAL :: INUSE
         INTEGER :: UNITNUM
      
         ! start checking for available units > 103, to avoid system default units
         ! 100, 101 and 102 are stdin, stdout and stderr respectively.

         INUSE=.TRUE.
         UNITNUM=103

         DO WHILE (INUSE)
            INQUIRE(UNIT=UNITNUM,OPENED=INUSE)
            IF (.NOT.INUSE) THEN
               GETUNIT=UNITNUM 
            ELSE     
               UNITNUM=UNITNUM+1
            ENDIF
         ENDDO
      END FUNCTION GETUNIT


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

END MODULE UTILITIES