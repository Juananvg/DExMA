test_Pvalue <- function() {
  library(DExMA)
  checkEqualsNumeric(metaAnalysisDE(maObject, typeMethod = "Fisher", 
                                    missAllow = 0.3, proportionData = 1)[1,1],
                     30.48075, tolerance=1.0e-6)}