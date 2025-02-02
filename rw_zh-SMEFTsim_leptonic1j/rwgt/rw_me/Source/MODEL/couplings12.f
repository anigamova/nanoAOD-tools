ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP12()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_907 = (MDL_CHL3*MDL_COMPLEXI*MDL_VEVHAT__EXP__2*MDL_YS)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2)
      GC_908 = -(MDL_CLL1*MDL_COMPLEXI*MDL_VEVHAT__EXP__2*MDL_YS)
     $ /(2.000000D+00*MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2)
      GC_927 = -((MDL_CLEDQIM*MDL_YE*MDL_YS)/MDL_LAMBDASMEFT__EXP__2)
      GC_928 = (MDL_CLEDQIM*MDL_YE*MDL_YS)/MDL_LAMBDASMEFT__EXP__2
      GC_932 = (MDL_CLEDQRE*MDL_COMPLEXI*MDL_YE*MDL_YS)
     $ /MDL_LAMBDASMEFT__EXP__2
      GC_936 = -((MDL_CLEDQIM*MDL_YM*MDL_YS)/MDL_LAMBDASMEFT__EXP__2)
      GC_937 = (MDL_CLEDQIM*MDL_YM*MDL_YS)/MDL_LAMBDASMEFT__EXP__2
      GC_941 = (MDL_CLEDQRE*MDL_COMPLEXI*MDL_YM*MDL_YS)
     $ /MDL_LAMBDASMEFT__EXP__2
      GC_1113 = -((MDL_COMPLEXI*MDL_YTAU)/MDL_SQRT__2)
      GC_1119 = -((MDL_CEBIM*MDL_CTH*MDL_YTAU)/(MDL_LAMBDASMEFT__EXP__2
     $ *MDL_SQRT__2))
      GC_1120 = (MDL_CEBRE*MDL_CTH*MDL_COMPLEXI*MDL_YTAU)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2)
      GC_1121 = (MDL_CEWIM*MDL_CTH*MDL_YTAU)/(MDL_LAMBDASMEFT__EXP__2
     $ *MDL_SQRT__2)
      GC_1122 = -((MDL_CEWRE*MDL_CTH*MDL_COMPLEXI*MDL_YTAU)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2))
      GC_1126 = -((MDL_COMPLEXI*MDL_PROPCORR*MDL_YTAU)/MDL_SQRT__2)
      GC_1132 = (MDL_CEBIM*MDL_STH*MDL_YTAU)/(MDL_LAMBDASMEFT__EXP__2
     $ *MDL_SQRT__2)
      GC_1133 = -((MDL_CEBRE*MDL_COMPLEXI*MDL_STH*MDL_YTAU)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2))
      GC_1134 = (MDL_CEWIM*MDL_STH*MDL_YTAU)/(MDL_LAMBDASMEFT__EXP__2
     $ *MDL_SQRT__2)
      GC_1135 = -((MDL_CEWRE*MDL_COMPLEXI*MDL_STH*MDL_YTAU)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2))
      GC_1136 = (-3.000000D+00*MDL_CEHIM*MDL_VEVHAT*MDL_YTAU)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2)
      GC_1137 = (3.000000D+00*MDL_CEHRE*MDL_COMPLEXI*MDL_VEVHAT
     $ *MDL_YTAU)/(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2)
      GC_1141 = -((MDL_CEBIM*MDL_CTH*MDL_VEVHAT*MDL_YTAU)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2))
      GC_1142 = (MDL_CEBRE*MDL_CTH*MDL_COMPLEXI*MDL_VEVHAT*MDL_YTAU)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2)
      GC_1143 = (MDL_CEWIM*MDL_CTH*MDL_VEVHAT*MDL_YTAU)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2)
      GC_1144 = -((MDL_CEWRE*MDL_CTH*MDL_COMPLEXI*MDL_VEVHAT*MDL_YTAU)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2))
      GC_1153 = (MDL_CEBIM*MDL_STH*MDL_VEVHAT*MDL_YTAU)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2)
      END
