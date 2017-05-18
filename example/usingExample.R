library(rtcor)
rtCorrectionInHouse(inputFilePath = '~/DDA_lib/NEG/idresults6.csv', outputFileName = 'correctedResult.csv', polyDegreeHighest = 5, outPath = '~/DDA_lib/test/', QCFilePath = '~/DDA_lib/QCRT/', pol = 'NEG', method = 'polyline')
