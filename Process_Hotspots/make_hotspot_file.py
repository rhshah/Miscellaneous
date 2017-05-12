#!/usr/bin/python
'''
@Description : This tool helps to modify the hotspot24k excel file to make a hotspot coordinate bed file 
@Created :  05/10/2017
@Updated : 05/12/2017
@author : Ronak H Shah

'''
from __future__ import division
import argparse
import sys
import os
import time
import logging
import tempfile

logging.basicConfig(
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.DEBUG)
logger = logging.getLogger('make_hotspot_file')
try:
    import coloredlogs
    coloredlogs.install(level='DEBUG')
except ImportError:
    logger.warning("make_hotspot_file: coloredlogs is not installed, please install it if you wish to see color in logs on standard out.")
    pass
try:
    import pandas as pd
except ImportError:
    logger.fatal("make_hotspot_file: pandas is not installed, please install pandas as it is required to run the removing process.")
    sys.exit(1)
try:
    import pybedtools
except ImportError:
    logger.fatal("make_hotspot_file: pybedtools is not installed, please install pybedtools as it is required to run the removing process.")
    sys.exit(1)

def main():
   parser = argparse.ArgumentParser(prog='make_hotspot_file.py', description='This tool helps to modify the hotspot24k excel file to make a hotspot coordinate bed file', usage='%(prog)s [options]')
   parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", help="make lots of noise")
   parser.add_argument("-ixlsx", "--input-excel", action="store", dest="inputXL", required=True, type=str, metavar='SomeID.xlsx', help="Input excel file which needs to be fixed")
   parser.add_argument("-obed","--output-bed", action="store", dest="outputBed", required=True, type=str, metavar='SomeID.bed', help="Output bed file name")
   parser.add_argument("-m", "--merge_intervals", action="store_true", dest="merge_intervals", help="Flag if given will output merged interval file")
   parser.add_argument("-o", "--outDir", action="store", dest="outdir", required=False, type=str, metavar='/somepath/output', help="Full Path to the output dir.")
   
   args = parser.parse_args()
   if(args.verbose):
       logger.info("make_hotspot_file: Started the run for replacing allele counts.")
   (snvDF,indelDF) = read_excel_file(args)
   (snvBed) = process_DF(snvDF,"SNVs")
   (indelBed) = process_DF(indelDF, "InDels")
   bedList = snvBed + indelBed
   makebed(args, bedList)
   if(args.verbose):
       logger.info("make_hotspot_file: Finished the run for replacing allele counts.")

def read_excel_file(args):
    snvDF = pd.read_excel(args.inputXL, sheetname="single_codon")
    indelDF = pd.read_excel(args.inputXL, sheetname="indels")
    return(snvDF,indelDF)

def process_DF(dataDF,tag):
    bedList = list()
    for index, row in dataDF.iterrows():
        gene = (str(row.loc['Hugo_Symbol'])).rstrip()
        coordinates = (str(row.loc['Genomic_Position'])).rstrip()
        aa = (str(row.loc['Variant_Amino_Acid'])).rstrip()
        coord_dict=dict([x.split("_") for x in coordinates.split("|")])
        for coord in coord_dict.keys():
            (chr,pos) = coord.split(':')
            #print chr,"\t",pos,"\t",pos,"\t",gene,":",aa,"\t","0","\t","+"
            bed = [chr,pos,pos,tag+";"+gene+";"+aa,".","+"]
            bedList.append(bed)
    return(bedList)

def makebed(args,bedList):
    temp = tempfile.NamedTemporaryFile()
    outfile = None
    if(args.outdir):
        outfile = os.path.join(args.outdir,args.outputBed)
    else:
        outfile = os.path.join(os.getcwd(),args.outputBed)
    for i in bedList:
        temp.write("\t".join(i) + "\n")
    temp.seek(0)
    newbed = pybedtools.BedTool(temp.name)
    srtbed = newbed.sort()
    if(args.merge_intervals):
        mrgbed = srtbed.merge(d=10)
        mrgbed.saveas(outfile)
    else:
        srtbed.saveas(outfile)
    # Automatically cleans up the file
    temp.close()
    
if __name__ == "__main__":
    start_time = time.time()  
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("make_hotspot_file: Elapsed time was %g seconds", totaltime)
    sys.exit(0)
