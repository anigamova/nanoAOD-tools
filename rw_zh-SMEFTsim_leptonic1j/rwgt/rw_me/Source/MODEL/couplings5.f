ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP5()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_434 = (6.000000D+00*MDL_CH*MDL_COMPLEXI*MDL_VEVHAT__EXP__3)
     $ /MDL_LAMBDASMEFT__EXP__2
      GC_435 = (-6.000000D+00*MDL_CHBOX*MDL_COMPLEXI*MDL_LAM
     $ *MDL_VEVHAT__EXP__3)/MDL_LAMBDASMEFT__EXP__2
      GC_436 = (3.000000D+00*MDL_CHDD*MDL_COMPLEXI*MDL_LAM
     $ *MDL_VEVHAT__EXP__3)/(2.000000D+00*MDL_LAMBDASMEFT__EXP__2)
      GC_437 = (6.000000D+00*MDL_CHL3*MDL_COMPLEXI*MDL_LAM
     $ *MDL_VEVHAT__EXP__3)/MDL_LAMBDASMEFT__EXP__2
      GC_438 = (-3.000000D+00*MDL_CLL1*MDL_COMPLEXI*MDL_LAM
     $ *MDL_VEVHAT__EXP__3)/MDL_LAMBDASMEFT__EXP__2
      GC_443 = (MDL_CHBOX*MDL_EE__EXP__2*MDL_COMPLEXI
     $ *MDL_VEVHAT__EXP__3)/(2.000000D+00*MDL_CTH__EXP__2
     $ *MDL_LAMBDASMEFT__EXP__2*MDL_STH__EXP__2)
      GC_444 = (MDL_CHDD*MDL_EE__EXP__2*MDL_COMPLEXI
     $ *MDL_VEVHAT__EXP__3)/(8.000000D+00*MDL_CTH__EXP__2
     $ *MDL_LAMBDASMEFT__EXP__2*MDL_STH__EXP__2)
      GC_445 = -(MDL_CHL3*MDL_EE__EXP__2*MDL_COMPLEXI
     $ *MDL_VEVHAT__EXP__3)/(2.000000D+00*MDL_CTH__EXP__2
     $ *MDL_LAMBDASMEFT__EXP__2*MDL_STH__EXP__2)
      GC_446 = (MDL_CLL1*MDL_EE__EXP__2*MDL_COMPLEXI
     $ *MDL_VEVHAT__EXP__3)/(4.000000D+00*MDL_CTH__EXP__2
     $ *MDL_LAMBDASMEFT__EXP__2*MDL_STH__EXP__2)
      GC_447 = (-2.000000D+00*MDL_CHWB*MDL_COMPLEXI*MDL_VEVHAT)
     $ /MDL_LAMBDASMEFT__EXP__2+(4.000000D+00*MDL_CHWB*MDL_COMPLEXI
     $ *MDL_STH__EXP__2*MDL_VEVHAT)/MDL_LAMBDASMEFT__EXP__2
      GC_448 = (-2.000000D+00*MDL_CHWBTIL*MDL_COMPLEXI*MDL_VEVHAT)
     $ /MDL_LAMBDASMEFT__EXP__2+(4.000000D+00*MDL_CHWBTIL*MDL_COMPLEXI
     $ *MDL_STH__EXP__2*MDL_VEVHAT)/MDL_LAMBDASMEFT__EXP__2
      GC_509 = -((MDL_COMPLEXI*MDL_YC)/MDL_SQRT__2)
      GC_510 = -((MDL_CTH*MDL_CUBIM*MDL_YC)/(MDL_LAMBDASMEFT__EXP__2
     $ *MDL_SQRT__2))
      GC_511 = (MDL_CTH*MDL_CUBRE*MDL_COMPLEXI*MDL_YC)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2)
      GC_512 = -((MDL_CUGIM*MDL_YC)/(MDL_LAMBDASMEFT__EXP__2
     $ *MDL_SQRT__2))
      GC_513 = (MDL_CUGRE*MDL_COMPLEXI*MDL_YC)/(MDL_LAMBDASMEFT__EXP__2
     $ *MDL_SQRT__2)
      GC_518 = -((MDL_CTH*MDL_CUWIM*MDL_YC)/(MDL_LAMBDASMEFT__EXP__2
     $ *MDL_SQRT__2))
      GC_520 = (MDL_CTH*MDL_CUWRE*MDL_COMPLEXI*MDL_YC)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2)
      GC_526 = -((MDL_COMPLEXI*MDL_PROPCORR*MDL_YC)/MDL_SQRT__2)
      GC_532 = (MDL_CUBIM*MDL_STH*MDL_YC)/(MDL_LAMBDASMEFT__EXP__2
     $ *MDL_SQRT__2)
      GC_533 = -((MDL_CUBRE*MDL_COMPLEXI*MDL_STH*MDL_YC)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2))
      GC_534 = -((MDL_CUWIM*MDL_STH*MDL_YC)/(MDL_LAMBDASMEFT__EXP__2
     $ *MDL_SQRT__2))
      GC_535 = (MDL_CUWRE*MDL_COMPLEXI*MDL_STH*MDL_YC)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2)
      GC_536 = -((MDL_CTH*MDL_CUBIM*MDL_VEVHAT*MDL_YC)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2))
      GC_537 = (MDL_CTH*MDL_CUBRE*MDL_COMPLEXI*MDL_VEVHAT*MDL_YC)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2)
      END
