MODULE CALCPROPERTIES
   USE VECTORS, ONLY: euc_norm, crossproduct, normed_cp2, vec_diff, normed_vec
   USE COMMONS
   IMPLICIT NONE
   REAL(KIND=REAL64), SAVE :: HELAX(3)

   CONTAINS

      SUBROUTINE PEPTIDE_PITCH(IMIN, X)
         USE VECTORS, ONLY: ANGLE, VEC_DIFF
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: IMIN
         REAL(KIND=REAL64), INTENT(IN) :: X(3*NATOMS)
         REAL(KIND=REAL64) :: XC(3), XN(3), PEP(3)
         INTEGER :: I, AT1, AT2, IDX

         DO I=1,NPEPBONDS
            AT1=PEPBONDS(I,1)
            AT2=PEPBONDS(I,2)
            IDX = 3*(AT1-1)
            XC(1:3) = X(IDX+1:IDX+3)
            IDX = 3*(AT2-1)
            XN(1:3) = X(IDX+1:IDX+3)
            PEP(1:3) = VEC_DIFF(XC, XN)
            PEPANGLES(IMIN,I) = ANGLE(HELAX,PEP)
         END DO

      END SUBROUTINE PEPTIDE_PITCH

      SUBROUTINE METHYLENE_PITCH(IMIN,X)
         USE VECTORS, ONLY: ANGLE, VEC_DIFF
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: IMIN
         REAL(KIND=REAL64), INTENT(IN) :: X(3*NATOMS)
         REAL(KIND=REAL64) :: XC(3), XH1(3), XH2(3), XH12(3), VEC(3)
         INTEGER :: I, AT1, AT2, AT3, IDX

         DO I=1,NMETHYLENE
            AT1=METHYLENEGROUPS(I,1)
            AT2=METHYLENEGROUPS(I,2)
            AT3=METHYLENEGROUPS(I,3)
            IDX = 3*(AT1-1)
            XC(1:3) = X(IDX+1:IDX+3)
            IDX = 3*(AT2-1)
            XH1(1:3) = X(IDX+1:IDX+3)
            IDX = 3*(AT3-1)
            XH2(1:3) = X(IDX+1:IDX+3) 
            XH12(1:3) = (XH1+XH2)/2.0D0          
            VEC(1:3) = VEC_DIFF(XC, XH12)
            METHYLENEANGLES(IMIN,I) = ANGLE(HELAX,VEC)
         END DO
      END SUBROUTINE METHYLENE_PITCH

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
            TMPBONDS(NPEPBONDS,2) = AT2           
         END DO
         ALLOCATE(PEPBONDS(NPEPBONDS,2))
         DO I=1,NPEPBONDS
            PEPBONDS(I,1) = TMPBONDS(I,1)
            PEPBONDS(I,2) = TMPBONDS(I,2)
         END DO
      END SUBROUTINE DETERMINE_PEPTIDE

      SUBROUTINE FIND_METHYLENE()
         USE UTILITIES, ONLY: LOOKUP_ATOM
         IMPLICIT NONE
         
         INTEGER :: TMPMETH(3*NRES,3)
         INTEGER :: I, AT1, AT2, AT3
         
         DO I=1,NRES
            IF (RESNAMES(I).EQ."PRO") THEN
               CALL LOOKUP_ATOM(I,"CD",AT1)
               CALL LOOKUP_ATOM(I,"HD2",AT2)  
               CALL LOOKUP_ATOM(I,"HD3",AT3)  
               NMETHYLENE = NMETHYLENE + 1
               TMPMETH(NMETHYLENE,1) = AT1
               TMPMETH(NMETHYLENE,2) = AT2
               TMPMETH(NMETHYLENE,3) = AT3

               CALL LOOKUP_ATOM(I,"CG",AT1)
               CALL LOOKUP_ATOM(I,"HG2",AT2)  
               CALL LOOKUP_ATOM(I,"HG3",AT3)  
               NMETHYLENE = NMETHYLENE + 1
               TMPMETH(NMETHYLENE,1) = AT1
               TMPMETH(NMETHYLENE,2) = AT2
               TMPMETH(NMETHYLENE,3) = AT3

               CALL LOOKUP_ATOM(I,"CB",AT1)
               CALL LOOKUP_ATOM(I,"HB2",AT2)  
               CALL LOOKUP_ATOM(I,"HB3",AT3)  
               NMETHYLENE = NMETHYLENE + 1
               TMPMETH(NMETHYLENE,1) = AT1
               TMPMETH(NMETHYLENE,2) = AT2
               TMPMETH(NMETHYLENE,3) = AT3

            ELSE IF (RESNAMES(I).EQ."HYP") THEN
               CALL LOOKUP_ATOM(I,"CD",AT1)
               CALL LOOKUP_ATOM(I,"HD2",AT2)  
               CALL LOOKUP_ATOM(I,"HD3",AT3)  
               NMETHYLENE = NMETHYLENE + 1
               TMPMETH(NMETHYLENE,1) = AT1
               TMPMETH(NMETHYLENE,2) = AT2
               TMPMETH(NMETHYLENE,3) = AT3

               CALL LOOKUP_ATOM(I,"CB",AT1)
               CALL LOOKUP_ATOM(I,"HB2",AT2)  
               CALL LOOKUP_ATOM(I,"HB3",AT3)  
               NMETHYLENE = NMETHYLENE + 1
               TMPMETH(NMETHYLENE,1) = AT1
               TMPMETH(NMETHYLENE,2) = AT2
               TMPMETH(NMETHYLENE,3) = AT3
            END IF
         ENDDO
         ALLOCATE(METHYLENEGROUPS(NMETHYLENE,3))
         DO I=1,NMETHYLENE
            METHYLENEGROUPS(I,1) = TMPMETH(NMETHYLENE,1)
            METHYLENEGROUPS(I,2) = TMPMETH(NMETHYLENE,2)
            METHYLENEGROUPS(I,3) = TMPMETH(NMETHYLENE,3)
         END DO
      END SUBROUTINE FIND_METHYLENE
END MODULE CALCPROPERTIES