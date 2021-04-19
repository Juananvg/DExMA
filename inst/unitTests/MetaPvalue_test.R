test_Pvalue <- function() {
  library(DExMA)
  checkEqualsNumeric(metaAnalysisDE(maObject, typeMethod = "Fisher", 
                                    missAllow = 0.3, proportionData = 1)[1,1],
                                    34.358809, tolerance=1.0e-6)}