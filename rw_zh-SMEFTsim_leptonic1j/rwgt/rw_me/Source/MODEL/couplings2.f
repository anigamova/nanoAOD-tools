ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP2()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_223 = (4.000000D+00*MDL_COMPLEXI*MDL_GHAA)/MDL_VEVHAT
      GC_227 = (2.000000D+00*MDL_COMPLEXI*MDL_GHZA)/MDL_VEVHAT
      GC_245 = (2.000000D+00*MDL_COMPLEXI*MDL_GHZA*MDL_PROPCORR)
     $ /MDL_VEVHAT
      GC_260 = -6.000000D+00*MDL_COMPLEXI*MDL_LAM*MDL_VEVHAT
      GC_262 = (-3.000000D+00*MDL_CHBOX*MDL_COMPLEXI*MDL_VEVHAT)
     $ /MDL_LAMBDASMEFT__EXP__2
      GC_263 = -((MDL_CHDD*MDL_COMPLEXI*MDL_VEVHAT)
     $ /MDL_LAMBDASMEFT__EXP__2)
      GC_264 = (4.000000D+00*MDL_CHG*MDL_COMPLEXI*MDL_VEVHAT)
     $ /MDL_LAMBDASMEFT__EXP__2
      GC_265 = (-2.000000D+00*MDL_CHGTIL*MDL_COMPLEXI*MDL_VEVHAT)
     $ /MDL_LAMBDASMEFT__EXP__2
      GC_268 = (4.000000D+00*MDL_CHB*MDL_CTH__EXP__2*MDL_COMPLEXI
     $ *MDL_VEVHAT)/MDL_LAMBDASMEFT__EXP__2
      GC_269 = (-2.000000D+00*MDL_CHBTIL*MDL_CTH__EXP__2*MDL_COMPLEXI
     $ *MDL_VEVHAT)/MDL_LAMBDASMEFT__EXP__2
      GC_270 = (4.000000D+00*MDL_CHW*MDL_CTH__EXP__2*MDL_COMPLEXI
     $ *MDL_VEVHAT)/MDL_LAMBDASMEFT__EXP__2
      GC_271 = (-2.000000D+00*MDL_CHWTIL*MDL_CTH__EXP__2*MDL_COMPLEXI
     $ *MDL_VEVHAT)/MDL_LAMBDASMEFT__EXP__2
      GC_280 = -6.000000D+00*MDL_COMPLEXI*MDL_LAM*MDL_PROPCORR
     $ *MDL_VEVHAT
      GC_281 = -6.000000D+00*MDL_COMPLEXI*MDL_LAM*MDL_PROPCORR__EXP__2
     $ *MDL_VEVHAT
      GC_284 = (MDL_EE__EXP__2*MDL_COMPLEXI*MDL_VEVHAT)/(2.000000D+00
     $ *MDL_CTH__EXP__2*MDL_STH__EXP__2)
      GC_289 = (MDL_EE__EXP__2*MDL_COMPLEXI*MDL_PROPCORR*MDL_VEVHAT)
     $ /(2.000000D+00*MDL_CTH__EXP__2*MDL_STH__EXP__2)
      GC_291 = (MDL_EE__EXP__2*MDL_COMPLEXI*MDL_PROPCORR__EXP__2
     $ *MDL_VEVHAT)/(2.000000D+00*MDL_CTH__EXP__2*MDL_STH__EXP__2)
      GC_304 = (MDL_CHD*MDL_EE*MDL_COMPLEXI*MDL_VEVHAT)/(MDL_CTH
     $ *MDL_LAMBDASMEFT__EXP__2*MDL_STH)
      GC_305 = (MDL_CHE*MDL_EE*MDL_COMPLEXI*MDL_VEVHAT)/(MDL_CTH
     $ *MDL_LAMBDASMEFT__EXP__2*MDL_STH)
      GC_306 = (MDL_CHL1*MDL_EE*MDL_COMPLEXI*MDL_VEVHAT)/(MDL_CTH
     $ *MDL_LAMBDASMEFT__EXP__2*MDL_STH)
      GC_308 = (MDL_CHL3*MDL_EE*MDL_COMPLEXI*MDL_VEVHAT)/(MDL_CTH
     $ *MDL_LAMBDASMEFT__EXP__2*MDL_STH)
      GC_309 = (MDL_CHQ1*MDL_EE*MDL_COMPLEXI*MDL_VEVHAT)/(MDL_CTH
     $ *MDL_LAMBDASMEFT__EXP__2*MDL_STH)
      GC_310 = -((MDL_CHQ3*MDL_EE*MDL_COMPLEXI*MDL_VEVHAT)/(MDL_CTH
     $ *MDL_LAMBDASMEFT__EXP__2*MDL_STH))
      GC_311 = (MDL_CHQ3*MDL_EE*MDL_COMPLEXI*MDL_VEVHAT)/(MDL_CTH
     $ *MDL_LAMBDASMEFT__EXP__2*MDL_STH)
      GC_312 = (MDL_CHU*MDL_EE*MDL_COMPLEXI*MDL_VEVHAT)/(MDL_CTH
     $ *MDL_LAMBDASMEFT__EXP__2*MDL_STH)
      END
