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
import time
import os
import sys
import logging
import glob


logger = logging.getLogger('get_DMP_ID')
try:
    import coloredlogs
    coloredlogs.install(level='DEBUG')
except ImportError:
    logger.warning("get_DMP_ID: coloredlogs is not installed, please install it if you wish to see color in logs on standard out.")
    pass

try:
    import pandas as pd
except ImportError, e:
    logger.warning("get_DMP_ID: pandas is not installed, please install pandas as it is required to run the mapping.")
    sys.exit(1)


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
    
    if(os.path.isfile(args.idFilename)):
        if(verbose):
            logger.info("get_DMP_ID: %s is being Read...", args.idFilename)
        iDF = pd.read_table(args.idFilename,sep="\t",header=None,usecols=[0])
        iDF.columns = ["SAMPLE_ID"]
        if(verbose):
            logger.info("Identification DataFrame:\n%s",iDF.head(n=5))
    else:
        logger.error("get_DMP_ID: %s is not a file. We will exit. Please make sure you provide a valid sample id file before rerun.", args.idFilename)
        sys.exit(1)
        
    if(os.path.isfile(args.idFilename)):
        if(verbose):
            logger.info("get_DMP_ID: %s is being Read...", args.refFilename)
        dDF = pd.read_table(args.refFilename,sep=",",header=None,usecols=[0,1,2])
        dDF.columns = ["SAMPLE_ID","D_SAMPLE_ID","GROUP_ID"]
        if(verbose):
            logger.info("Mapping DataFrame:\n%s",dDF.head(n=5))
    else:
        logging.error("get_DMP_ID: %s is not a file. We will exit. Please make sure you provide a valid sample mapping file before rerun.", args.refFilename)
        sys.exit(1)
    
    if(os.path.isdir(args.bamLocation)):
        if(verbose):
            logger.info("get_DMP_ID: %s is a directory", args.bamLocation)
        else:
            pass
    else:
        logger.error("get_DMP_ID: %s is not a directory. We will exit. Please make sure you provide a valid bam location directory before rerun.", args.bamLocation)
        sys.exit(1)
    outDF = pd.DataFrame(
        columns=[
            "SAMPLE_ID",
            "D_SAMPLE_ID",
            "SAMPLE_TYPE",
            "GROUP_ID",
            "BAM_LOCATION"])
    i = 0
    for count,row in iDF.iterrows():
        iID = str(row.loc["SAMPLE_ID"])
        try:
            dDF_idx = dDF[dDF["SAMPLE_ID"] == iID].index.tolist()
        except IndexError:
            dDF_idx = None
        if(dDF_idx):
            dID = dDF.loc[dDF_idx[0],"D_SAMPLE_ID"]
            gID = dDF.loc[dDF_idx[0],"GROUP_ID"]
            bamFile = ",".join(glob.glob(args.bamLocation +"/" + str(dID) +"*.bam"))
            outDF.loc[i,["SAMPLE_ID",
                             "D_SAMPLE_ID",
                             "SAMPLE_TYPE",
                             "GROUP_ID",
                             "BAM_LOCATION"]] = [str(iID),str(dID),"TUMOR",str(gID),str(bamFile)]
            i = i+1
            #Get Normals
            gID_idx = dDF[dDF["GROUP_ID"] == gID].index.tolist()
            for gidx in gID_idx:
                dDF_iID = dDF.loc[gidx,"SAMPLE_ID"]
                if "-N" in dDF_iID:
                    dID = dDF.loc[gidx,"D_SAMPLE_ID"]
                    gID = dDF.loc[gidx,"GROUP_ID"]
                    bamFile = ",".join(glob.glob(args.bamLocation +"/" + str(dID) +"*.bam"))
                    outDF.loc[i,["SAMPLE_ID",
                                     "D_SAMPLE_ID",
                                     "SAMPLE_TYPE",
                                     "GROUP_ID",
                                     "BAM_LOCATION"]] = [str(dDF_iID),str(dID),"NORMAL",str(gID),str(bamFile)]
                    i = i + 1
                else:
                    continue
            
        else:
            logger.critical("Sample ID: %s is not present in mapping file we will print an empty entry.", iID)
            outDF.loc[count,["SAMPLE_ID",
                             "D_SAMPLE_ID",
                             "SAMPLE_TYPE",
                             "GROUP_ID",
                             "BAM_LOCATION"]] = [str(iID),None,"TUMOR",None,None] 
            
    
    outDF.sort_values(["GROUP_ID"], inplace=True, ascending=True)  
    outDF.drop_duplicates(inplace=True)
    outDF.to_csv(args.outFile, sep='\t', index=False)
    if(verbose):
        logger.info("get_DMP_ID: Finished Mapping, Final data written in %s", args.outFile)
        
        
if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("get_DMP_ID: Elapsed time was %g seconds", totaltime)
    
