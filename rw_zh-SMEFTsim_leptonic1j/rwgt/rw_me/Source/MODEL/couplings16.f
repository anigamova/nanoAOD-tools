ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP16()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_1419 = -(MDL_CLEQU3IM*MDL_YTAU*MDL_YUP)/(2.000000D+00
     $ *MDL_LAMBDASMEFT__EXP__2)
      GC_1420 = -(MDL_CLEQU3IM*MDL_YTAU*MDL_YUP)/(4.000000D+00
     $ *MDL_LAMBDASMEFT__EXP__2)
      GC_1421 = (MDL_CLEQU3IM*MDL_YTAU*MDL_YUP)/(4.000000D+00
     $ *MDL_LAMBDASMEFT__EXP__2)
      GC_1422 = (MDL_CLEQU3IM*MDL_YTAU*MDL_YUP)/(2.000000D+00
     $ *MDL_LAMBDASMEFT__EXP__2)
      GC_1429 = -(MDL_CLEQU3RE*MDL_COMPLEXI*MDL_YTAU*MDL_YUP)
     $ /(4.000000D+00*MDL_LAMBDASMEFT__EXP__2)
      GC_1430 = (MDL_CLEQU3RE*MDL_COMPLEXI*MDL_YTAU*MDL_YUP)
     $ /(2.000000D+00*MDL_LAMBDASMEFT__EXP__2)
      END
