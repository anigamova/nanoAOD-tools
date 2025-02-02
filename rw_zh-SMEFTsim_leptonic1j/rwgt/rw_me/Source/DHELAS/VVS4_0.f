C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(1,2)*P(2,1) - P(-1,1)*P(-1,2)*Metric(1,2)
C     
      SUBROUTINE VVS4_0(V1, V2, S3, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 S3(*)
      REAL*8 P1(0:3)
      REAL*8 P2(0:3)
      COMPLEX*16 TMP17
      COMPLEX*16 TMP20
      COMPLEX*16 TMP28
      COMPLEX*16 TMP24
      COMPLEX*16 VERTEX
      COMPLEX*16 COUP
      COMPLEX*16 V1(*)
      P1(0) = DBLE(V1(1))
      P1(1) = DBLE(V1(2))
      P1(2) = DIMAG(V1(2))
      P1(3) = DIMAG(V1(1))
      P2(0) = DBLE(V2(1))
      P2(1) = DBLE(V2(2))
      P2(2) = DIMAG(V2(2))
      P2(3) = DIMAG(V2(1))
      TMP24 = (P2(0)*V1(3)-P2(1)*V1(4)-P2(2)*V1(5)-P2(3)*V1(6))
      TMP20 = (P1(0)*V2(3)-P1(1)*V2(4)-P1(2)*V2(5)-P1(3)*V2(6))
      TMP17 = (V2(3)*V1(3)-V2(4)*V1(4)-V2(5)*V1(5)-V2(6)*V1(6))
      TMP28 = (P1(0)*P2(0)-P1(1)*P2(1)-P1(2)*P2(2)-P1(3)*P2(3))
      VERTEX = COUP*S3(3)*(-CI*(TMP20*TMP24)+CI*(TMP17*TMP28))
      END


C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(1,2)*P(2,1) - P(-1,1)*P(-1,2)*Metric(1,2)
C     
      SUBROUTINE VVS4_5_6_7_0(V1, V2, S3, COUP1, COUP2, COUP3, COUP4
     $ ,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 V2(*)
      COMPLEX*16 S3(*)
      REAL*8 P1(0:3)
      REAL*8 P2(0:3)
      COMPLEX*16 COUP1
      COMPLEX*16 COUP2
      COMPLEX*16 TMP
      COMPLEX*16 COUP3
      COMPLEX*16 VERTEX
      COMPLEX*16 V1(*)
      COMPLEX*16 COUP4
      CALL VVS4_0(V1,V2,S3,COUP1,VERTEX)
      CALL VVS5_0(V1,V2,S3,COUP2,TMP)
      VERTEX = VERTEX + TMP
      CALL VVS6_0(V1,V2,S3,COUP3,TMP)
      VERTEX = VERTEX + TMP
      CALL VVS7_0(V1,V2,S3,COUP4,TMP)
      VERTEX = VERTEX + TMP
      END


