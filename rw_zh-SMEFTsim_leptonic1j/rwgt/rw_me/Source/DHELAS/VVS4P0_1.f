C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(1,2)*P(2,1) - P(-1,1)*P(-1,2)*Metric(1,2)
C     
      SUBROUTINE VVS4P0_1(V2, S3, COUP, M1, W1,V1)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 S3(*)
      REAL*8 P1(0:3)
      REAL*8 M1
      REAL*8 P2(0:3)
      REAL*8 W1
      COMPLEX*16 TMP20
      COMPLEX*16 DENOM
      COMPLEX*16 TMP28
      COMPLEX*16 COUP
      COMPLEX*16 V1(6)
      P2(0) = DBLE(V2(1))
      P2(1) = DBLE(V2(2))
      P2(2) = DIMAG(V2(2))
      P2(3) = DIMAG(V2(1))
      V1(1) = +V2(1)+S3(1)
      V1(2) = +V2(2)+S3(2)
      P1(0) = -DBLE(V1(1))
      P1(1) = -DBLE(V1(2))
      P1(2) = -DIMAG(V1(2))
      P1(3) = -DIMAG(V1(1))
      TMP20 = (P1(0)*V2(3)-P1(1)*V2(4)-P1(2)*V2(5)-P1(3)*V2(6))
      TMP28 = (P1(0)*P2(0)-P1(1)*P2(1)-P1(2)*P2(2)-P1(3)*P2(3))
      DENOM = COUP/(P1(0)**2-P1(1)**2-P1(2)**2-P1(3)**2 - M1 * (M1 -CI
     $ * W1))
      V1(3)= DENOM*S3(3)*(-CI*(P2(0)*TMP20)+CI*(V2(3)*TMP28))
      V1(4)= DENOM*S3(3)*(-CI*(P2(1)*TMP20)+CI*(V2(4)*TMP28))
      V1(5)= DENOM*S3(3)*(-CI*(P2(2)*TMP20)+CI*(V2(5)*TMP28))
      V1(6)= DENOM*S3(3)*(-CI*(P2(3)*TMP20)+CI*(V2(6)*TMP28))
      END


C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(1,2)*P(2,1) - P(-1,1)*P(-1,2)*Metric(1,2)
C     
      SUBROUTINE VVS4_5_6_7P0_1(V2, S3, COUP1, COUP2, COUP3, COUP4, M1
     $ , W1,V1)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 COUP2
      COMPLEX*16 S3(*)
      REAL*8 P1(0:3)
      REAL*8 M1
      REAL*8 P2(0:3)
      REAL*8 W1
      COMPLEX*16 COUP1
      COMPLEX*16 DENOM
      COMPLEX*16 COUP3
      COMPLEX*16 VTMP(6)
      INTEGER*4 I
      COMPLEX*16 V1(6)
      COMPLEX*16 COUP4
      CALL VVS4P0_1(V2,S3,COUP1,M1,W1,V1)
      CALL VVS5P0_1(V2,S3,COUP2,M1,W1,VTMP)
      DO I = 3, 6
        V1(I) = V1(I) + VTMP(I)
      ENDDO
      CALL VVS6P0_1(V2,S3,COUP3,M1,W1,VTMP)
      DO I = 3, 6
        V1(I) = V1(I) + VTMP(I)
      ENDDO
      CALL VVS7P0_1(V2,S3,COUP4,M1,W1,VTMP)
      DO I = 3, 6
        V1(I) = V1(I) + VTMP(I)
      ENDDO
      END


