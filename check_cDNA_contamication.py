"""
Created on 25/11/2014.

@Ronak Shah

"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
from collections import defaultdict
import sys
import math
from collections import Counter


def main():
    parser = argparse.ArgumentParser(
        prog='AnnalyzeSVs.py',
        description='Calculate the frequency of the Structural Variants',
        usage='%(prog)s [options]')
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        default=True,
        help="make lots of noise [default]")
    parser.add_argument(
        "-s",
        "--svFile",
        action="store",
        dest="svFilename",
        required=True,
        metavar='SVfile.txt',
        help="Location of the structural variant file to be used")
    parser.add_argument(
        "-o",
        "--outputFilePrefix",
        action="store",
        dest="outFilePrefix",
        required=True,
        metavar='AnnotatedSV',
        help="Full path with prefix name for the output file")
    parser.add_argument(
        "-r",
        "--ranges",
        action="store",
        nargs='*',
        dest="ranges",
        required=True,
        metavar='2,3,50,100',
        help="Number of ranges to be used to run Analysis; Multiple ranges are separated by space")
    args = parser.parse_args()
    (dataDF) = Read_SVfile(args.svFilename, args)


def Read_SVfile(file, args):
    # print file
    txt_fh = open(args.outFilePrefix, "wb")
    dataDF = pd.read_csv(file, sep='\t', header=0, keep_default_na='True')
    countData = defaultdict(list)
    # FirstPass
    txt_fh.write("chr1\tgene1\tpos1\trange\tchr2\tgene2\tpos2\trange\tcount\n")
    for count, row in dataDF.iterrows():
        #print row,"\n"
        chr1 = str(row.loc['site1_chrom'])
        chr2 = str(row.loc['site2_chrom'])
        pos1 = (row.loc['site1_pos'])
        pos2 = (row.loc['site2_pos'])
        gene1 = str(row.loc['site1_gene'])
        gene2 = str(row.loc['site2_gene'])
        if((math.isnan(pos1)) or (math.isnan(pos2))):
            continue
        #newEntry1 = chr1 + ":" + pos1 + ":" + gene1 + ":" + gene2
        #newEntry2 = chr2 + ":" + pos2 + ":" + gene2 + ":" + gene1
        #newEntry = chr1 + "\t" + str(pos1) + "\t" + str(pos1+1) + "\t" + chr2 + "\t" + str(pos2) + "\t" + str(pos2+1) + "\t" + gene1 + ":" + gene2
        #txt_fh.write(newEntry + "\n")
        # txt_fh.close()
        key = chr1 + ":" + gene1 + ":" + chr2 + ":" + gene2
        val = str(pos1) + ":" + str(pos2)
        countData[key].append(val)
        # countData[chr2].append(pos2)
    counts = Counter(countData)
    newData = list(set(countData))
    # print counts
    for item in newData:
        (chr1, gene1, chr2, gene2) = item.split(":")
        positions = counts[item]
        if(len(counts[item]) >= 2):
            # print '%s : %d' % (item, len(counts[item]))
            chr1pos = []
            chr2pos = []
            for pos in positions:
                print pos,"\n"
                (pos1, pos2) = pos.split(":")
                chr1pos.append(int(float(pos1)))
                chr2pos.append(int(float(pos2)))
            chr1avg = average(chr1pos)
            frac, whole = math.modf(chr1avg)
            if(frac >= 0.5):
                chr1avg = int(math.ceil(chr1avg))
            else:
                chr1avg = int(math.floor(chr1avg))
            chr2avg = average(chr2pos)
            frac, whole = math.modf(chr2avg)
            if(frac >= 0.5):
                chr2avg = int(math.ceil(chr2avg))
            else:
                chr2avg = int(math.floor(chr2avg))

            chr1variance = map(lambda x: (x - chr1avg) ** 2, chr1pos)
            chr2variance = map(lambda x: (x - chr2avg) ** 2, chr2pos)

            chr1standarddeviation = math.sqrt(average(chr1variance))
            frac, whole = math.modf(chr1standarddeviation)
            if(frac >= 0.5):
                chr1standarddeviation = int(math.ceil(chr1standarddeviation))
            else:
                chr1standarddeviation = int(math.floor(chr1standarddeviation))

            chr2standarddeviation = math.sqrt(average(chr2variance))
            frac, whole = math.modf(chr2standarddeviation)
            if(frac >= 0.5):
                chr2standarddeviation = int(math.ceil(chr2standarddeviation))
            else:
                chr2standarddeviation = int(math.floor(chr2standarddeviation))

            txt_fh.write(chr1 +
                         "\t" +
                         gene1 +
                         "\t" +
                         str(chr1avg) +
                         "\t" +
                         str(chr1standarddeviation) +
                         "\t" +
                         chr2 +
                         "\t" +
                         gene2 +
                         "\t" +
                         str(chr2avg) +
                         "\t" +
                         str(chr2standarddeviation) +
                         "\t" +
                         str(len(counts[item])) +
                         "\n")
        else:
            (pos1, pos2) = positions[0].split(":")
            # print pos1,pos2
            txt_fh.write(str(chr1) +
                         "\t" +
                         gene1 +
                         "\t" +
                         str(pos1) +
                         "\t0\t" +
                         str(chr2) +
                         "\t" +
                         gene2 +
                         "\t" +
                         str(pos2) +
                         "\t0\t" +
                         str(len(counts[item])) +
                         "\n")
    txt_fh.close()

'''
    keylist = countData.keys()
    keylist.sort()
    for k in keylist:
        sortedVal = sorted(countData[k],key=int)
        for val in sortedVal:
            newval = int(val)+1
            txt_fh.write(k + "\t" + val + "\t" + str(newval) + "\n")
    txt_fh.close()

    #Second Pass
    keysTodelete = []
    for count, row in dataDF.iterrows():
        #print row
        chr1 = str(row.loc['Chr1'])
        chr2 = str(row.loc['Chr2'])
        pos1 = (row.loc['Pos1'])
        pos2 = (row.loc['Pos2'])
        gene1 = row.loc['Gene1']
        gene2 = row.loc['Gene2']
        for cord,val in countData.items():
            (cChr,cPos,cGene1,cGene2) = cord.split(":")
            if (cChr == chr1):
                if((int(pos1)-int(args.range))>int(cPos)>(int(pos1)+int(args.range))):
                    newPos = int((pos1+cPos)/2)
                    newval = val + 1
                    newEntry = cChr + ":" + cPos + ":" + cGene1 + ":" + cGene2
                    counData[str(newEntry)] = newval
                    keysTodelete.append(cord)
            if (cChr == chr2):
                if((int(pos2)-int(args.range))>int(cPos)>(int(pos2)+int(args.range))):
                    newPos = int((pos2+cPos)/2)
                    newval = val + 1
                    newEntry = cChr + ":" + cPos + ":" + cGene1 + ":" + cGene2
                    counData[str(newEntry)] = newval
                    keysTodelete.append(cord)
    map(countData.__delitem__, keysTodelete)

    for cord,val in countData.items():
        print cord," = ", val



    counts = Counter(countData)
    newData = list(set(countData))
    #print counts
    for item in newData:
        if(counts[item]> 3):
            print '%s : %d' % (item, counts[item])
'''


def average(s):
    return sum(s) * 1.0 / len(s)

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    #print("Elapsed time was %g seconds" % (end_time - start_time))
