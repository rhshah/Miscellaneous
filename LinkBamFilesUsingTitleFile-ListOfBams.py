'''
Created On : 05/27/2015
@author: Ronak H Shah
'''
import argparse
import pandas as pd
import os,sys
import time
import copy
def main():
    parser = argparse.ArgumentParser(
        prog='LinkBamFilesUsingTitleFile-ListOfBams.py',
        description='Link Files using a list of bams and there metainfo.',
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
        "--titleFile",
        action="store",
        dest="titleFile",
        required=True,
        metavar='title_file.txt',
        help="Location of the title file with meta info.")
    parser.add_argument(
        "-l",
        "--listOfBams",
        action="store",
        dest="listOfBams",
        required=True,
        metavar='somelist.txt',
        help="Location of the file of files with bams.")
    parser.add_argument(
        "-o",
        "--outputDir",
        action="store",
        dest="outDir",
        required=True,
        metavar='/some/path',
        help="Full path to output directory")

    args = parser.parse_args()
    (tDF) = Read_titleFile(args.titleFile)
    (bList) = Read_BamList(args.listOfBams)
    LinkTheData(tDF, bList,args)
    
def Read_titleFile(file):
    dataDF = pd.read_csv(file, sep='\t', header=0, keep_default_na='True')
    return (dataDF)
def Read_BamList(file):
    allFiles = []
    with open(file, 'r') as f:
        allFiles = [line.strip() for line in f]
    return(allFiles)

def LinkTheData(titleDF,bamList,args):
    for count,bamfile in enumerate(bamList):
        this_dir, bamFileName = os.path.split(bamfile)
        #print this_dir,this_filename
        baiFile = copy.copy(bamFileName)
        baiFile = this_dir + "/" + baiFile[:-1] + 'i'
        (sampleIdFromList,rest) = bamFileName.split("_",1)
        idx = int(titleDF[titleDF['Sample_ID'] == sampleIdFromList].index)
        barcode = str(titleDF.iloc[idx]['Barcode'])
        project = str(titleDF.iloc[idx]['Pool'])
        newSortBamName = args.outDir + "/" + sampleIdFromList + "_" + barcode + "_" + project + "_L000_mrg_cl_aln_srt.bam"
        newSortBaiName = args.outDir + "/" + sampleIdFromList + "_" + barcode + "_" + project + "_L000_mrg_cl_aln_srt.bai"
        newMdBamName =  args.outDir + "/" + sampleIdFromList + "_" + barcode + "_" + project + "_L000_mrg_cl_aln_srt_MD.bam"
        newMdBaiName =  args.outDir + "/" + sampleIdFromList + "_" + barcode + "_" + project + "_L000_mrg_cl_aln_srt_MD.bai"
        os.symlink(bamfile, newSortBamName)
        os.symlink(baiFile, newSortBaiName)
        os.symlink(bamfile, newMdBamName)
        os.symlink(baiFile, newMdBaiName)
    return
        
if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    print("Elapsed time was %g seconds" % (end_time - start_time))
