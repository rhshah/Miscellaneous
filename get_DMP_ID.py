"""
get_DMP_ID
~~~~~~~~~~

:Description: This module gets the information about Deidentified DMP ID associated with Identified DMP ID

"""

'''
Created on December 13, 2016
Description: This module gets the information about Deidentified DMP ID associated with Identified DMP ID
@author: Ronak H Shah
::Input::
listofid : List of identified ids
fileofmapping: csv file with identified and deidentified mapping
locationofbam: where the bam file will be found
::Output::
a tsv file with following header
SAMPLE_ID    D_SAMPLE_ID    SAMPLE_TYPE    GROUP_ID    BAM_LOCATION
'''

import argparse
import pandas as pd
import time
import os,sys
import coloredlogs
import logging
import glob
from collections import Counter


def main():
    parser = argparse.ArgumentParser(
        prog='get_DMP_ID.py',
        description='This module gets the information about Deidentified DMP ID associated with Identified DMP ID',
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
        "--ids",
        action="store",
        dest="idFilename",
        required=True,
        metavar='ids.txt',
        help="Location of the file containing identified sample ids")
    parser.add_argument(
        "-o",
        "--outputFile",
        action="store",
        dest="outFile",
        required=True,
        metavar='out.txt',
        help="Full path with name for the output file")
    parser.add_argument(
        "-r",
        "--refFile",
        action="store",
        dest="refFilename",
        required=True,
        metavar='refFile.txt',
        help="Location of the file containing identified and de-identified sample mapping")
    parser.add_argument(
        "-b",
        "--bamLocation",
        action="store",
        dest="bamLocation",
        required=True,
        metavar='/prefix/to/bam/file/location',
        help="Location where the deidentified bam files are")
    args = parser.parse_args()

    #Get the process is verbose or not
    verbose = args.verbose
    #Get where are we
    here = os.path.realpath('.')

    logger = logging.getLogger('get_DMP_ID')
    coloredlogs.install(level='DEBUG')
    
    if(os.path.isfile(args.idFilename)):
        if(verbose):
            logger.info("get_DMP_ID: %s is being Read...", args.idFilename)
        iDF = pd.read_table(args.idFilename,sep="\t",header=None)
        iDF.columns = ["SAMPLE_ID"]
    else:
        logging.error("get_DMP_ID: %s is not a file. We will exit. Please make sure you provide a valid sample id file before rerun.", args.idFilename)
        sys.exit(1)
        
    if(os.path.isfile(args.idFilename)):
        if(verbose):
            logger.info("get_DMP_ID: %s is being Read...", args.refFilename)
        dDF = pd.read_table(args.refFilename,sep=",",header=None,usecols=[1,2,3])
        dDF.columns = ["SAMPLE_ID","D_SAMPLE_ID","GROUP_ID"]
    else:
        logging.error("get_DMP_ID:%s is not a file. We will exit. Please make sure you provide a valid sample mapping file before rerun.", args.refFilename)
        sys.exit(1)
    
    if(os.path.isdir(args.bamLocation)):
        if(verbose):
            logger.info("get_DMP_ID: %s is a directory", args.bamLocation)
        else:
            pass
    else:
        logging.error("get_DMP_ID:%s is not a directory. We will exit. Please make sure you provide a valid bam location directory before rerun.", args.bamLocation)
       
    outDF = pd.DataFrame(
        columns=[
            "SAMPLE_ID",
            "D_SAMPLE_ID",
            "SAMPLE_TYPE",
            "GROUP_ID",
            "BAM_LOCATION"])
    
    for count,row in iDF.iterrows():
        iID = row.loc["SAMPLE_ID"]
        dDF_idx = dDF[dDF["SAMPLE_ID"]==iID].index().tolist()[0]
        dID = dDF.iloc[dDF_idx,"D_SAMPLE_ID"]
        gID = dDF.iloc[dDF_idx,"GROUP_ID"]
        bamFile = glob.glob(args.bamLocation +"/" + dID +"*.bam")
        outDF.loc[count,["SAMPLE_ID",
            "D_SAMPLE_ID",
            "SAMPLE_TYPE",
            "GROUP_ID",
            "BAM_LOCATION"]] = [iID,dID,gID,bamFile] 
    
    outDF.to_csv(args.outFile, sep='\t', index=False)
    if(verbose):
        logger.info("get_DMP_ID: Finished Mapping, Final data written in %s", args.outFile)
        
        
if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("iCallSV:Elapsed time was %g seconds", totaltime)
    
