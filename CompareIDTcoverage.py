'''
Created On : 02/27/2015
@author: Ronak H Shah
'''
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
from collections import defaultdict
import sys
from collections import Counter
from tabulate import tabulate

def main():
    parser = argparse.ArgumentParser(prog='CompareIDTcoverage.py', description='Compare the coverage per probe sequence', usage='%(prog)s [options]')
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=True, help="make lots of noise [default]")
    parser.add_argument("-i1", "--covgFile1", action="store", dest="covgFile1", required=True, metavar='DOP-PROBE-BASED.sample_interval_summary', help="Location of the coverage file to be used")
    parser.add_argument("-i2", "--covgFile2", action="store", dest="covgFile2", required=True, metavar='DOP-PROBE-BASED.sample_interval_summary', help="Location of the coverage file to be used")
    parser.add_argument("-r", "--wellRef", action="store", dest="wellReferenceFile", required=True, metavar='someFile.txt', help="Location of the coordinate to well reference")
    parser.add_argument("-o", "--outputFilePrefix", action="store", dest="outFilePrefix", required=True, metavar='AnnotatedSV', help="Full path with prefix name for the output file")
    
    args = parser.parse_args()
    (WellDict) = Read_WellInfo(args.wellReferenceFile)
    (CovgDict1) = Read_CovgFile(args.covgFile1)
    (CovgDict2) = Read_CovgFile(args.covgFile2)
    AverageTheData(WellDict,CovgDict1,CovgDict2);
  
def Read_WellInfo(file):  
    dataDF = pd.read_csv(file, sep='\t', header=0, keep_default_na='True')
    d = defaultdict(list)
    for count, row in dataDF.iterrows():
        wellName = str(row.loc['NameWell Position'])
        chr = str(row.loc['Chr'])
        start = str(row.loc['Start'])
        end = str(row.loc['Stop'])
        val = chr + ":" + start + "-" + end
        d[wellName].append(val)
    return(d)

def Read_CovgFile(file):
    dataDF = pd.read_csv(file, sep='\t', header=0, keep_default_na='True')
    d = {}
    for count, row in dataDF.iterrows():
        target = str(row.loc['Locus'])
        covg = row.loc['Average_Depth_sample']
        d[target] = covg
    return(d)

def AverageTheData(wellDict,covgDict1,covgDict2):
    resultDict = defaultdict(list)
    for i,j in wellDict.iteritems():
        #print i,"\n"
        coverageForWell = []
        #print i,len(j)
        for coords in j:
            (chr,stend) = coords.split(":",1)
            (start,end) = stend.split("-",1)
            #print chr,",",start,",",end
            coverageForARange1 = []
            coverageForARange2 = []
            for items in xrange(int(start),int(end)):
                key = chr + ":" + str(items)
                if covgDict1.has_key(key):
                    val = covgDict1.get(key)
                    coverageForARange1.append(val)
                #else:
                    #print key,"1", "=>", 0
                if covgDict2.has_key(key):
                    val = covgDict2.get(key)
                    coverageForARange2.append(val)
                #else:
                    #print key,"2", "=>", 0
            averageRangeCoverage1 = avg(i,2,coverageForARange1)
            averageRangeCoverage2 = avg(i,2,coverageForARange2)
            ratio1 = None
            ratio2 = None
            if(averageRangeCoverage1 > averageRangeCoverage2 and averageRangeCoverage1 != 0):
                ratio1 = averageRangeCoverage2/float(averageRangeCoverage1)
            else:
                if(averageRangeCoverage2 != 0):
                    ratio2 = averageRangeCoverage1/float(averageRangeCoverage2)
            print i , chr, start, end , averageRangeCoverage1, averageRangeCoverage2, ratio1, ratio2
            #coverageForWell.append(averageRangeCoverage)
        #averageWellCoverage = avg(i,1,coverageForWell)
        #print averageWellCoverage
        #resultDict[i].append(averageWellCoverage)
    #print resultDict
    #headers = ["WellPosition", "Coverage"]
    #print tabulate(resultDict, headers="keys",tablefmt="plain")
        
def avg(i,num,list):
    sum = 0
    for elm in list:
        sum += elm
    if(i.startswith("6") and (num == 1)):
        divlen = 7
    else:
        divlen = len(list)
    return (sum/(divlen)*1.0)
        
if __name__ == "__main__":
    start_time = time.time()  
    main()
    end_time = time.time()
    print("Elapsed time was %g seconds" % (end_time - start_time))