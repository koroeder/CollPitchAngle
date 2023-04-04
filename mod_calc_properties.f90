MODULE CALCPROPERTIES
   USE VECTORS, ONLY: euc_norm, crossproduct, normed_cp2, vec_diff, normed_vec
   USE COMMONS
   IMPLICIT NONE
   REAL(KIND=REAL64), SAVE :: HELAX(3)

   CONTAINS
      SUBROUTINE GET_HELAX(X)
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: X(3*NATOMS)
         REAL(KIND=REAL64) :: ATS1(3), ATS2(3), NORM
         INTEGER :: I, J, IDX

         ATS1(1:3) = 0.0D0
         ATS2(1:3) = 0.0D0
         DO I=1,NPOINT1
            IDX = ATSFOR1(I)
            DO J=1,3
               ATS1(J) = ATS1(J) + X(3*(IDX-1)+J)
            END DO
         END DO
         DO I=1,NPOINT2
            IDX = ATSFOR2(I)
            DO J=1,3
               ATS2(J) = ATS2(J) + X(3*(IDX-1)+J)
            END DO
         END DO 
         ATS1 = ATS1/NPOINT1
         ATS2 = ATS2/NPOINT2
         
         CALL NORMED_VEC(VEC_DIFF(ATS1,ATS2),HELAX,NORM)
      END SUBROUTINE GET_HELAX


      SUBROUTINE DETERMINE_PEPTIDE()
         USE UTILITIES, ONLY: LOOKUP_ATOM
         IMPLICIT NONE
         
         INTEGER :: TMPBONDS(NRES,2)
         LOGICAL :: ISTERMINUS(NRES)
         INTEGER :: RES1, RES2, I, AT1, AT2
         
         ISTERMINUS(1:NRES) = .FALSE.

         DO I=1,NTERMINI
            ISTERMINUS(TERMINI(I,1)) = .TRUE.
            ISTERMINUS(TERMINI(I,2)) = .TRUE.
         END DO

         DO I=1,NRES-1
            RES1 = I
            RES2 = I + 1
            IF (ISTERMINUS(RES1).AND.ISTERMINUS(RES2)) CYCLE
            CALL LOOKUP_ATOM(RES1,"C",AT1)
            CALL LOOKUP_ATOM(RES2,"N",AT2)
            IF ((AT1.EQ.-1).OR.(AT2.EQ.-1)) CYCLE
            NPEPBONDS = NPEPBONDS + 1
            TMPBONDS(NPEPBONDS,1) = AT1
            TMPBONDS(NPEPBONDS,1) = AT2           
         END DO
         ALLOCATE(PEPBONDS(NPEPBONDS,2))
         DO I=1,NPEPBONDS
            PEPBONDS(I,1) = TMPBONDS(I,1)
            PEPBONDS(I,2) = TMPBONDS(I,2)
         END DO
      END SUBROUTINE DETERMINE_PEPTIDE
END MODULE CALCPROPERTIES