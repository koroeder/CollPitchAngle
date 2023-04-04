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

END MODULE CALCPROPERTIES