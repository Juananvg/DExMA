test_ES <- function() {
  library(DExMA)
  checkEqualsNumeric(metaAnalysisDE(maObject, typeMethod = "REM", 
                    missAllow = 0.3, proportionData = 1)[1,1],
                    2.359380876, tolerance=1.0e-6)}