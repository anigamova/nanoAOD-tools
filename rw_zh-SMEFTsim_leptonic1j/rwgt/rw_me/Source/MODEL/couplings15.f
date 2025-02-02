ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP15()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_1298 = (MDL_CHL3*MDL_COMPLEXI*MDL_VEVHAT__EXP__2*MDL_YUP)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2)
      GC_1299 = -(MDL_CLL1*MDL_COMPLEXI*MDL_VEVHAT__EXP__2*MDL_YUP)
     $ /(2.000000D+00*MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2)
      GC_1300 = -((MDL_CUHIM*MDL_VEVHAT__EXP__2*MDL_YUP)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2))
      GC_1301 = (MDL_CUHRE*MDL_COMPLEXI*MDL_VEVHAT__EXP__2*MDL_YUP)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2)
      GC_1338 = -((MDL_CLEQU1IM*MDL_YE*MDL_YUP)/MDL_LAMBDASMEFT__EXP__2)
     $ 
      GC_1339 = (MDL_CLEQU1IM*MDL_YE*MDL_YUP)/MDL_LAMBDASMEFT__EXP__2
      GC_1343 = -((MDL_CLEQU1RE*MDL_COMPLEXI*MDL_YE*MDL_YUP)
     $ /MDL_LAMBDASMEFT__EXP__2)
      GC_1347 = -(MDL_CLEQU3IM*MDL_YE*MDL_YUP)/(2.000000D+00
     $ *MDL_LAMBDASMEFT__EXP__2)
      GC_1348 = -(MDL_CLEQU3IM*MDL_YE*MDL_YUP)/(4.000000D+00
     $ *MDL_LAMBDASMEFT__EXP__2)
      GC_1349 = (MDL_CLEQU3IM*MDL_YE*MDL_YUP)/(4.000000D+00
     $ *MDL_LAMBDASMEFT__EXP__2)
      GC_1350 = (MDL_CLEQU3IM*MDL_YE*MDL_YUP)/(2.000000D+00
     $ *MDL_LAMBDASMEFT__EXP__2)
      GC_1357 = -(MDL_CLEQU3RE*MDL_COMPLEXI*MDL_YE*MDL_YUP)/(4.000000D
     $ +00*MDL_LAMBDASMEFT__EXP__2)
      GC_1358 = (MDL_CLEQU3RE*MDL_COMPLEXI*MDL_YE*MDL_YUP)/(2.000000D
     $ +00*MDL_LAMBDASMEFT__EXP__2)
      GC_1365 = -((MDL_CLEQU1IM*MDL_YM*MDL_YUP)/MDL_LAMBDASMEFT__EXP__2)
     $ 
      GC_1366 = (MDL_CLEQU1IM*MDL_YM*MDL_YUP)/MDL_LAMBDASMEFT__EXP__2
      GC_1370 = -((MDL_CLEQU1RE*MDL_COMPLEXI*MDL_YM*MDL_YUP)
     $ /MDL_LAMBDASMEFT__EXP__2)
      GC_1374 = -(MDL_CLEQU3IM*MDL_YM*MDL_YUP)/(2.000000D+00
     $ *MDL_LAMBDASMEFT__EXP__2)
      GC_1375 = -(MDL_CLEQU3IM*MDL_YM*MDL_YUP)/(4.000000D+00
     $ *MDL_LAMBDASMEFT__EXP__2)
      GC_1376 = (MDL_CLEQU3IM*MDL_YM*MDL_YUP)/(4.000000D+00
     $ *MDL_LAMBDASMEFT__EXP__2)
      GC_1377 = (MDL_CLEQU3IM*MDL_YM*MDL_YUP)/(2.000000D+00
     $ *MDL_LAMBDASMEFT__EXP__2)
      GC_1384 = -(MDL_CLEQU3RE*MDL_COMPLEXI*MDL_YM*MDL_YUP)/(4.000000D
     $ +00*MDL_LAMBDASMEFT__EXP__2)
      GC_1385 = (MDL_CLEQU3RE*MDL_COMPLEXI*MDL_YM*MDL_YUP)/(2.000000D
     $ +00*MDL_LAMBDASMEFT__EXP__2)
      GC_1410 = -((MDL_CLEQU1IM*MDL_YTAU*MDL_YUP)/MDL_LAMBDASMEFT__EXP_
     $_2)
      GC_1411 = (MDL_CLEQU1IM*MDL_YTAU*MDL_YUP)/MDL_LAMBDASMEFT__EXP__2
      GC_1415 = -((MDL_CLEQU1RE*MDL_COMPLEXI*MDL_YTAU*MDL_YUP)
     $ /MDL_LAMBDASMEFT__EXP__2)
      END
