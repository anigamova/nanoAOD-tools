C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     -2*Epsilon(1,2,3,-2)*P(-2,3)*P(-1,1)*P(-1,2) -
C      2*Epsilon(1,2,3,-2)*P(-2,2)*P(-1,1)*P(-1,3) -
C      2*Epsilon(1,2,3,-2)*P(-2,1)*P(-1,2)*P(-1,3) +
C      2*Epsilon(2,3,-1,-2)*P(-2,3)*P(-1,1)*P(1,2) +
C      2*Epsilon(2,3,-1,-2)*P(-2,2)*P(-1,1)*P(1,3) -
C      2*Epsilon(1,3,-1,-2)*P(-2,3)*P(-1,2)*P(2,1) -
C      2*Epsilon(1,3,-1,-2)*P(-2,1)*P(-1,2)*P(2,3) +
C      2*Epsilon(1,2,-1,-2)*P(-2,2)*P(-1,3)*P(3,1) +
C      2*Epsilon(1,2,-1,-2)*P(-2,1)*P(-1,3)*P(3,2) +
C      Epsilon(3,-1,-2,-3)*P(-3,2)*P(-2,3)*P(-1,1)*Metric(1,2) -
C      Epsilon(3,-1,-2,-3)*P(-3,1)*P(-2,3)*P(-1,2)*Metric(1,2) -
C      Epsilon(2,-1,-2,-3)*P(-3,3)*P(-2,2)*P(-1,1)*Metric(1,3) +
C      Epsilon(2,-1,-2,-3)*P(-3,1)*P(-2,2)*P(-1,3)*Metric(1,3) +
C      Epsilon(1,-1,-2,-3)*P(-3,3)*P(-2,1)*P(-1,2)*Metric(2,3) -
C      Epsilon(1,-1,-2,-3)*P(-3,2)*P(-2,1)*P(-1,3)*Metric(2,3)
C     
      SUBROUTINE VVV8_0(V1, V2, V3, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 TMP20
      REAL*8 P3(0:3)
      COMPLEX*16 TMP21
      COMPLEX*16 TMP32
      COMPLEX*16 TMP28
      COMPLEX*16 TMP52
      COMPLEX*16 TMP25
      COMPLEX*16 TMP18
      COMPLEX*16 V3(*)
      REAL*8 P2(0:3)
      COMPLEX*16 TMP46
      COMPLEX*16 TMP33
      COMPLEX*16 TMP53
      COMPLEX*16 TMP23
      COMPLEX*16 TMP24
      COMPLEX*16 TMP19
      COMPLEX*16 TMP49
      COMPLEX*16 V2(*)
      REAL*8 P1(0:3)
      COMPLEX*16 TMP45
      COMPLEX*16 TMP43
      COMPLEX*16 TMP54
      COMPLEX*16 TMP51
      COMPLEX*16 VERTEX
      COMPLEX*16 TMP41
      COMPLEX*16 TMP42
      COMPLEX*16 TMP48
      COMPLEX*16 TMP50
      COMPLEX*16 V1(*)
      COMPLEX*16 TMP22
      COMPLEX*16 TMP44
      COMPLEX*16 TMP17
      COMPLEX*16 COUP
      COMPLEX*16 TMP55
      COMPLEX*16 TMP47
      P1(0) = DBLE(V1(1))
      P1(1) = DBLE(V1(2))
      P1(2) = DIMAG(V1(2))
      P1(3) = DIMAG(V1(1))
      P2(0) = DBLE(V2(1))
      P2(1) = DBLE(V2(2))
      P2(2) = DIMAG(V2(2))
      P2(3) = DIMAG(V2(1))
      P3(0) = DBLE(V3(1))
      P3(1) = DBLE(V3(2))
      P3(2) = DIMAG(V3(2))
      P3(3) = DIMAG(V3(1))
      TMP24 = (P2(0)*V1(3)-P2(1)*V1(4)-P2(2)*V1(5)-P2(3)*V1(6))
      TMP25 = (P3(0)*V1(3)-P3(1)*V1(4)-P3(2)*V1(5)-P3(3)*V1(6))
      TMP41 = (-1D0)*(P3(0)*(V1(4)*(V3(6)*V2(5)-V3(5)*V2(6))+(V1(5)
     $ *(V3(4)*V2(6)-V3(6)*V2(4))+V1(6)*(V3(5)*V2(4)-V3(4)*V2(5))))
     $ +(P3(1)*(V1(3)*(V3(5)*V2(6)-V3(6)*V2(5))+(V1(5)*(V3(6)*V2(3)
     $ -V3(3)*V2(6))+V1(6)*(V3(3)*V2(5)-V3(5)*V2(3))))+(P3(2)*(V1(3)
     $ *(V3(6)*V2(4)-V3(4)*V2(6))+(V1(4)*(V3(3)*V2(6)-V3(6)*V2(3))
     $ +V1(6)*(V3(4)*V2(3)-V3(3)*V2(4))))+P3(3)*(V1(3)*(V3(4)*V2(5)
     $ -V3(5)*V2(4))+(V1(4)*(V3(5)*V2(3)-V3(3)*V2(5))+V1(5)*(V3(3)
     $ *V2(4)-V3(4)*V2(3)))))))
      TMP46 = (-1D0)*(P2(0)*(P3(1)*(V3(5)*V1(6)-V3(6)*V1(5))+(P3(2)
     $ *(V3(6)*V1(4)-V3(4)*V1(6))+P3(3)*(V3(4)*V1(5)-V3(5)*V1(4))))
     $ +(P2(1)*(P3(0)*(V3(6)*V1(5)-V3(5)*V1(6))+(P3(2)*(V3(3)*V1(6)
     $ -V3(6)*V1(3))+P3(3)*(V3(5)*V1(3)-V3(3)*V1(5))))+(P2(2)*(P3(0)
     $ *(V3(4)*V1(6)-V3(6)*V1(4))+(P3(1)*(V3(6)*V1(3)-V3(3)*V1(6))
     $ +P3(3)*(V3(3)*V1(4)-V3(4)*V1(3))))+P2(3)*(P3(0)*(V3(5)*V1(4)
     $ -V3(4)*V1(5))+(P3(1)*(V3(3)*V1(5)-V3(5)*V1(3))+P3(2)*(V3(4)
     $ *V1(3)-V3(3)*V1(4)))))))
      TMP47 = (-1D0)*(P1(0)*(P2(1)*(V3(6)*V1(5)-V3(5)*V1(6))+(P2(2)
     $ *(V3(4)*V1(6)-V3(6)*V1(4))+P2(3)*(V3(5)*V1(4)-V3(4)*V1(5))))
     $ +(P1(1)*(P2(0)*(V3(5)*V1(6)-V3(6)*V1(5))+(P2(2)*(V3(6)*V1(3)
     $ -V3(3)*V1(6))+P2(3)*(V3(3)*V1(5)-V3(5)*V1(3))))+(P1(2)*(P2(0)
     $ *(V3(6)*V1(4)-V3(4)*V1(6))+(P2(1)*(V3(3)*V1(6)-V3(6)*V1(3))
     $ +P2(3)*(V3(4)*V1(3)-V3(3)*V1(4))))+P1(3)*(P2(0)*(V3(4)*V1(5)
     $ -V3(5)*V1(4))+(P2(1)*(V3(5)*V1(3)-V3(3)*V1(5))+P2(2)*(V3(3)
     $ *V1(4)-V3(4)*V1(3)))))))
      TMP44 = (-1D0)*(P1(0)*(P3(1)*(V3(5)*V2(6)-V3(6)*V2(5))+(P3(2)
     $ *(V3(6)*V2(4)-V3(4)*V2(6))+P3(3)*(V3(4)*V2(5)-V3(5)*V2(4))))
     $ +(P1(1)*(P3(0)*(V3(6)*V2(5)-V3(5)*V2(6))+(P3(2)*(V3(3)*V2(6)
     $ -V3(6)*V2(3))+P3(3)*(V3(5)*V2(3)-V3(3)*V2(5))))+(P1(2)*(P3(0)
     $ *(V3(4)*V2(6)-V3(6)*V2(4))+(P3(1)*(V3(6)*V2(3)-V3(3)*V2(6))
     $ +P3(3)*(V3(3)*V2(4)-V3(4)*V2(3))))+P1(3)*(P3(0)*(V3(5)*V2(4)
     $ -V3(4)*V2(5))+(P3(1)*(V3(3)*V2(5)-V3(5)*V2(3))+P3(2)*(V3(4)
     $ *V2(3)-V3(3)*V2(4)))))))
      TMP45 = (-1D0)*(P1(0)*(P2(1)*(V3(5)*V2(6)-V3(6)*V2(5))+(P2(2)
     $ *(V3(6)*V2(4)-V3(4)*V2(6))+P2(3)*(V3(4)*V2(5)-V3(5)*V2(4))))
     $ +(P1(1)*(P2(0)*(V3(6)*V2(5)-V3(5)*V2(6))+(P2(2)*(V3(3)*V2(6)
     $ -V3(6)*V2(3))+P2(3)*(V3(5)*V2(3)-V3(3)*V2(5))))+(P1(2)*(P2(0)
     $ *(V3(4)*V2(6)-V3(6)*V2(4))+(P2(1)*(V3(6)*V2(3)-V3(3)*V2(6))
     $ +P2(3)*(V3(3)*V2(4)-V3(4)*V2(3))))+P1(3)*(P2(0)*(V3(5)*V2(4)
     $ -V3(4)*V2(5))+(P2(1)*(V3(3)*V2(5)-V3(5)*V2(3))+P2(2)*(V3(4)
     $ *V2(3)-V3(3)*V2(4)))))))
      TMP48 = (-1D0)*(P2(0)*(P3(1)*(V2(6)*V1(5)-V2(5)*V1(6))+(P3(2)
     $ *(V2(4)*V1(6)-V2(6)*V1(4))+P3(3)*(V2(5)*V1(4)-V2(4)*V1(5))))
     $ +(P2(1)*(P3(0)*(V2(5)*V1(6)-V2(6)*V1(5))+(P3(2)*(V2(6)*V1(3)
     $ -V2(3)*V1(6))+P3(3)*(V2(3)*V1(5)-V2(5)*V1(3))))+(P2(2)*(P3(0)
     $ *(V2(6)*V1(4)-V2(4)*V1(6))+(P3(1)*(V2(3)*V1(6)-V2(6)*V1(3))
     $ +P3(3)*(V2(4)*V1(3)-V2(3)*V1(4))))+P2(3)*(P3(0)*(V2(4)*V1(5)
     $ -V2(5)*V1(4))+(P3(1)*(V2(5)*V1(3)-V2(3)*V1(5))+P3(2)*(V2(3)
     $ *V1(4)-V2(4)*V1(3)))))))
      TMP49 = (-1D0)*(P1(0)*(P3(1)*(V2(6)*V1(5)-V2(5)*V1(6))+(P3(2)
     $ *(V2(4)*V1(6)-V2(6)*V1(4))+P3(3)*(V2(5)*V1(4)-V2(4)*V1(5))))
     $ +(P1(1)*(P3(0)*(V2(5)*V1(6)-V2(6)*V1(5))+(P3(2)*(V2(6)*V1(3)
     $ -V2(3)*V1(6))+P3(3)*(V2(3)*V1(5)-V2(5)*V1(3))))+(P1(2)*(P3(0)
     $ *(V2(6)*V1(4)-V2(4)*V1(6))+(P3(1)*(V2(3)*V1(6)-V2(6)*V1(3))
     $ +P3(3)*(V2(4)*V1(3)-V2(3)*V1(4))))+P1(3)*(P3(0)*(V2(4)*V1(5)
     $ -V2(5)*V1(4))+(P3(1)*(V2(5)*V1(3)-V2(3)*V1(5))+P3(2)*(V2(3)
     $ *V1(4)-V2(4)*V1(3)))))))
      TMP28 = (P1(0)*P2(0)-P1(1)*P2(1)-P1(2)*P2(2)-P1(3)*P2(3))
      TMP20 = (P1(0)*V2(3)-P1(1)*V2(4)-P1(2)*V2(5)-P1(3)*V2(6))
      TMP21 = (V3(3)*V1(3)-V3(4)*V1(4)-V3(5)*V1(5)-V3(6)*V1(6))
      TMP22 = (P3(0)*V2(3)-P3(1)*V2(4)-P3(2)*V2(5)-P3(3)*V2(6))
      TMP23 = (V3(3)*V2(3)-V3(4)*V2(4)-V3(5)*V2(5)-V3(6)*V2(6))
      TMP51 = (-1D0)*(P1(0)*(P2(1)*(V3(6)*P3(2)-V3(5)*P3(3))+(P2(2)
     $ *(V3(4)*P3(3)-V3(6)*P3(1))+P2(3)*(V3(5)*P3(1)-V3(4)*P3(2))))
     $ +(P1(1)*(P2(0)*(V3(5)*P3(3)-V3(6)*P3(2))+(P2(2)*(V3(6)*P3(0)
     $ -V3(3)*P3(3))+P2(3)*(V3(3)*P3(2)-V3(5)*P3(0))))+(P1(2)*(P2(0)
     $ *(V3(6)*P3(1)-V3(4)*P3(3))+(P2(1)*(V3(3)*P3(3)-V3(6)*P3(0))
     $ +P2(3)*(V3(4)*P3(0)-V3(3)*P3(1))))+P1(3)*(P2(0)*(V3(4)*P3(2)
     $ -V3(5)*P3(1))+(P2(1)*(V3(5)*P3(0)-V3(3)*P3(2))+P2(2)*(V3(3)
     $ *P3(1)-V3(4)*P3(0)))))))
      TMP50 = (-1D0)*(P1(0)*(P2(1)*(V3(5)*P3(3)-V3(6)*P3(2))+(P2(2)
     $ *(V3(6)*P3(1)-V3(4)*P3(3))+P2(3)*(V3(4)*P3(2)-V3(5)*P3(1))))
     $ +(P1(1)*(P2(0)*(V3(6)*P3(2)-V3(5)*P3(3))+(P2(2)*(V3(3)*P3(3)
     $ -V3(6)*P3(0))+P2(3)*(V3(5)*P3(0)-V3(3)*P3(2))))+(P1(2)*(P2(0)
     $ *(V3(4)*P3(3)-V3(6)*P3(1))+(P2(1)*(V3(6)*P3(0)-V3(3)*P3(3))
     $ +P2(3)*(V3(3)*P3(1)-V3(4)*P3(0))))+P1(3)*(P2(0)*(V3(5)*P3(1)
     $ -V3(4)*P3(2))+(P2(1)*(V3(3)*P3(2)-V3(5)*P3(0))+P2(2)*(V3(4)
     $ *P3(0)-V3(3)*P3(1)))))))
      TMP53 = (-1D0)*(P1(0)*(P2(1)*(P3(3)*V2(5)-P3(2)*V2(6))+(P2(2)
     $ *(P3(1)*V2(6)-P3(3)*V2(4))+P2(3)*(P3(2)*V2(4)-P3(1)*V2(5))))
     $ +(P1(1)*(P2(0)*(P3(2)*V2(6)-P3(3)*V2(5))+(P2(2)*(P3(3)*V2(3)
     $ -P3(0)*V2(6))+P2(3)*(P3(0)*V2(5)-P3(2)*V2(3))))+(P1(2)*(P2(0)
     $ *(P3(3)*V2(4)-P3(1)*V2(6))+(P2(1)*(P3(0)*V2(6)-P3(3)*V2(3))
     $ +P2(3)*(P3(1)*V2(3)-P3(0)*V2(4))))+P1(3)*(P2(0)*(P3(1)*V2(5)
     $ -P3(2)*V2(4))+(P2(1)*(P3(2)*V2(3)-P3(0)*V2(5))+P2(2)*(P3(0)
     $ *V2(4)-P3(1)*V2(3)))))))
      TMP52 = (-1D0)*(P1(0)*(P2(1)*(P3(2)*V2(6)-P3(3)*V2(5))+(P2(2)
     $ *(P3(3)*V2(4)-P3(1)*V2(6))+P2(3)*(P3(1)*V2(5)-P3(2)*V2(4))))
     $ +(P1(1)*(P2(0)*(P3(3)*V2(5)-P3(2)*V2(6))+(P2(2)*(P3(0)*V2(6)
     $ -P3(3)*V2(3))+P2(3)*(P3(2)*V2(3)-P3(0)*V2(5))))+(P1(2)*(P2(0)
     $ *(P3(1)*V2(6)-P3(3)*V2(4))+(P2(1)*(P3(3)*V2(3)-P3(0)*V2(6))
     $ +P2(3)*(P3(0)*V2(4)-P3(1)*V2(3))))+P1(3)*(P2(0)*(P3(2)*V2(4)
     $ -P3(1)*V2(5))+(P2(1)*(P3(0)*V2(5)-P3(2)*V2(3))+P2(2)*(P3(1)
     $ *V2(3)-P3(0)*V2(4)))))))
      TMP55 = (-1D0)*(P1(0)*(P2(1)*(P3(2)*V1(6)-P3(3)*V1(5))+(P2(2)
     $ *(P3(3)*V1(4)-P3(1)*V1(6))+P2(3)*(P3(1)*V1(5)-P3(2)*V1(4))))
     $ +(P1(1)*(P2(0)*(P3(3)*V1(5)-P3(2)*V1(6))+(P2(2)*(P3(0)*V1(6)
     $ -P3(3)*V1(3))+P2(3)*(P3(2)*V1(3)-P3(0)*V1(5))))+(P1(2)*(P2(0)
     $ *(P3(1)*V1(6)-P3(3)*V1(4))+(P2(1)*(P3(3)*V1(3)-P3(0)*V1(6))
     $ +P2(3)*(P3(0)*V1(4)-P3(1)*V1(3))))+P1(3)*(P2(0)*(P3(2)*V1(4)
     $ -P3(1)*V1(5))+(P2(1)*(P3(0)*V1(5)-P3(2)*V1(3))+P2(2)*(P3(1)
     $ *V1(3)-P3(0)*V1(4)))))))
      TMP54 = (-1D0)*(P1(0)*(P2(1)*(P3(3)*V1(5)-P3(2)*V1(6))+(P2(2)
     $ *(P3(1)*V1(6)-P3(3)*V1(4))+P2(3)*(P3(2)*V1(4)-P3(1)*V1(5))))
     $ +(P1(1)*(P2(0)*(P3(2)*V1(6)-P3(3)*V1(5))+(P2(2)*(P3(3)*V1(3)
     $ -P3(0)*V1(6))+P2(3)*(P3(0)*V1(5)-P3(2)*V1(3))))+(P1(2)*(P2(0)
     $ *(P3(3)*V1(4)-P3(1)*V1(6))+(P2(1)*(P3(0)*V1(6)-P3(3)*V1(3))
     $ +P2(3)*(P3(1)*V1(3)-P3(0)*V1(4))))+P1(3)*(P2(0)*(P3(1)*V1(5)
     $ -P3(2)*V1(4))+(P2(1)*(P3(2)*V1(3)-P3(0)*V1(5))+P2(2)*(P3(0)
     $ *V1(4)-P3(1)*V1(3)))))))
      TMP19 = (V3(3)*P2(0)-V3(4)*P2(1)-V3(5)*P2(2)-V3(6)*P2(3))
      TMP18 = (V3(3)*P1(0)-V3(4)*P1(1)-V3(5)*P1(2)-V3(6)*P1(3))
      TMP33 = (P3(0)*P2(0)-P3(1)*P2(1)-P3(2)*P2(2)-P3(3)*P2(3))
      TMP32 = (P3(0)*P1(0)-P3(1)*P1(1)-P3(2)*P1(2)-P3(3)*P1(3))
      TMP17 = (V2(3)*V1(3)-V2(4)*V1(4)-V2(5)*V1(5)-V2(6)*V1(6))
      TMP42 = (-1D0)*(P2(0)*(V1(4)*(V3(6)*V2(5)-V3(5)*V2(6))+(V1(5)
     $ *(V3(4)*V2(6)-V3(6)*V2(4))+V1(6)*(V3(5)*V2(4)-V3(4)*V2(5))))
     $ +(P2(1)*(V1(3)*(V3(5)*V2(6)-V3(6)*V2(5))+(V1(5)*(V3(6)*V2(3)
     $ -V3(3)*V2(6))+V1(6)*(V3(3)*V2(5)-V3(5)*V2(3))))+(P2(2)*(V1(3)
     $ *(V3(6)*V2(4)-V3(4)*V2(6))+(V1(4)*(V3(3)*V2(6)-V3(6)*V2(3))
     $ +V1(6)*(V3(4)*V2(3)-V3(3)*V2(4))))+P2(3)*(V1(3)*(V3(4)*V2(5)
     $ -V3(5)*V2(4))+(V1(4)*(V3(5)*V2(3)-V3(3)*V2(5))+V1(5)*(V3(3)
     $ *V2(4)-V3(4)*V2(3)))))))
      TMP43 = (-1D0)*(P1(0)*(V1(4)*(V3(6)*V2(5)-V3(5)*V2(6))+(V1(5)
     $ *(V3(4)*V2(6)-V3(6)*V2(4))+V1(6)*(V3(5)*V2(4)-V3(4)*V2(5))))
     $ +(P1(1)*(V1(3)*(V3(5)*V2(6)-V3(6)*V2(5))+(V1(5)*(V3(6)*V2(3)
     $ -V3(3)*V2(6))+V1(6)*(V3(3)*V2(5)-V3(5)*V2(3))))+(P1(2)*(V1(3)
     $ *(V3(6)*V2(4)-V3(4)*V2(6))+(V1(4)*(V3(3)*V2(6)-V3(6)*V2(3))
     $ +V1(6)*(V3(4)*V2(3)-V3(3)*V2(4))))+P1(3)*(V1(3)*(V3(4)*V2(5)
     $ -V3(5)*V2(4))+(V1(4)*(V3(5)*V2(3)-V3(3)*V2(5))+V1(5)*(V3(3)
     $ *V2(4)-V3(4)*V2(3)))))))
      VERTEX = COUP*2D0*(TMP17*1D0/2D0*(-CI*(TMP50)+CI*(TMP51))+(TMP21
     $ *1D0/2D0*(-CI*(TMP53)+CI*(TMP52))+(TMP23*1D0/2D0*(-CI*(TMP54)
     $ +CI*(TMP55))+(-CI*(TMP24*TMP44+TMP25*TMP45+TMP18*TMP48+TMP19
     $ *TMP49)+CI*(TMP28*TMP41+TMP32*TMP42+TMP33*TMP43+TMP20*TMP46
     $ +TMP22*TMP47)))))
      END


