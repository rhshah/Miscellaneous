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
    parser = argparse.ArgumentParser(
        prog='AnnalyzeIDTcoverage.py',
        description='Calculate the coverage per well',
        usage='%(prog)s [options]')
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        default=True,
        help="make lots of noise [default]")
    parser.add_argument(
        "-i",
        "--covgFile",
        action="store",
        dest="covgFile",
        required=True,
        metavar='DOP-PROBE-BASED.sample_interval_summary',
        help="Location of the coverage file to be used")
    parser.add_argument(
        "-r",
        "--wellRef",
        action="store",
        dest="wellReferenceFile",
        required=True,
        metavar='someFile.txt',
        help="Location of the coordinate to well reference")
    parser.add_argument(
        "-o",
        "--outputFilePrefix",
        action="store",
        dest="outFilePrefix",
        required=True,
        metavar='AnnotatedSV',
        help="Full path with prefix name for the output file")

    args = parser.parse_args()
    (WellDict) = Read_WellInfo(args.wellReferenceFile)
    (CovgDict) = Read_CovgFile(args.covgFile)
    AverageTheData(WellDict, CovgDict)


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


def AverageTheData(wellDict, covgDict):
    resultDict = defaultdict(list)
    for i, j in wellDict.iteritems():
        # print i,"\n"
        coverageForWell = []
        print i, len(j)
        for coords in j:
            (chr, stend) = coords.split(":", 1)
            (start, end) = stend.split("-", 1)
            # print chr,",",start,",",end
            coverageForARange = []
            for items in xrange(int(start), int(end)):
                key = chr + ":" + str(items)
                if key in covgDict:
                    val = covgDict.get(key)
                    coverageForARange.append(val)
                else:
                    print key, "=>", 0
            averageRangeCoverage = avg(i, 2, coverageForARange)
            coverageForWell.append(averageRangeCoverage)
        averageWellCoverage = avg(i, 1, coverageForWell)
        # print averageWellCoverage
        resultDict[i].append(averageWellCoverage)
    # print resultDict
    headers = ["WellPosition", "Coverage"]
    print tabulate(resultDict, headers="keys", tablefmt="plain")


def avg(i, num, list):
    sum = 0
    for elm in list:
        sum += elm
    if(i.startswith("6") and (num == 1)):
        divlen = 7
    else:
        divlen = len(list)
    return (sum / (divlen) * 1.0)

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    print("Elapsed time was %g seconds" % (end_time - start_time))
