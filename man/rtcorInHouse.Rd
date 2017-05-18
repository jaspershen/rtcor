\name{rtcor}
\docType{data}
\alias{rtcorInHouse}
\title{RT Correction In house}
\description{
 RT Correction based on QCRT
}
\usage{
library(rtcor)
rtCorrectionInHouse(inputFilePath = '~/DDA_lib/NEG/idresults.csv',  # See ./data/idresults_example.csv
                    outputFileName = 'correctedResult.csv', 
					polyDegreeHighest = 5, outPath = '~/DDA_lib/test/',  # the highest degree of polynomial
					QCFilePath = '~/DDA_lib/QCRT/', 
					pol = 'NEG',  # polarity, 'NEG' or 'POS'
					method = 'polyline')  #method for correction RT, 'polyline' or 'loess'
}