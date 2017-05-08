'''
dmpGermline2vcf
~~~~~~~~~~

:Description: This module converts multisample dmp maf like format to tumor/normal vcf

'''

'''
Created on March 21, 2017
Description: This module converts multisample dmp maf like format to tumor/normal vcf
@author: Ronak H Shah
'''

import argparse
import time
import os
import sys
import logging
import glob

logging.basicConfig(
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.DEBUG)
logger = logging.getLogger('dmpGermline2vcf')
try:
    import coloredlogs
    coloredlogs.install(level='DEBUG')
except ImportError:
    logger.warning("dmpGermline2vcf: coloredlogs is not installed, please install it if you wish to see color in logs on standard out.")
    pass

try:
    import pandas as pd
except ImportError, e:
    logger.warning("dmpGermline2vcf: pandas is not installed, please install pandas as it is required to run the mapping.")
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(prog='dmpGermline2vcf.py', description='This module converts multi-sample dmp maf like format to tumor/normal vcf', usage='%(prog)s [options]')
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", help="make lots of noise [default]")
    parser.add_argument("-dmp", "--dmp_input_list", action="store", dest="mInput", required=True, metavar='file.txt or listoffiles.list', help="Location of the mutation file or list of files to be used")
    parser.add_argument("-vcf", "--vcf_output", action="store", dest="vcf", required=True, metavar='AnnotatedSV', help="Full path with prefix name for the output file")
    args = parser.parse_args()
    (mDF) = ReadMfile(args.mInput)
    (header) = makeVCFheader()
    makeVCF(mDF,header,args.vcf)
   
def ReadMfile(input):
    dataDF = None
    if(input.endswith(".list")):
        allFiles = []
        listdf = []
        header = []
        with open(input, 'r') as f:
            allFiles = [line.strip() for line in f]
        for count,file in enumerate(allFiles):
            with open(file,'r') as filecontent:
                header = filecontent.readline().strip('\n').split('\t')
                #header = header[0:35]
                for line in filecontent:
                    data = line.strip('\n').split('\t')
                    #data = data[0:35]
                    listdf.append(data)
        dataDF = pd.DataFrame(listdf,columns=header)
    else:
        dataDF = pd.read_csv(input, sep='\t', header=0, keep_default_na='True')
    return(dataDF)

def makeVCFheader():
    header = ["##fileformat=VCFv4.3",
              "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">",
              "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Depth matching reference/alternate (REF/ALT) allele\">",
              "##FORMAT=<ID=ADF,Number=1,Type=Integer,Description=\"Depth matching forward reference/alternate (REF/ALT) allele\">",
              "##FORMAT=<ID=ADR,Number=1,Type=Integer,Description=\"Depth matching reverse reference/alternate (REF/ALT) allele\">",
              "##INFO=<ID=Sample,Number=1,Type=String,Description=\"Name of the Sample where the variant is called\">",
              "##INFO=<ID=CallMethod,Number=1,Type=String,Description=\"Name of the tool used to make the call\">",
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL"]
    header = "\n".join(header)
    return(header)
#`T_TotalDepth    T_RefCount    T_AltCount    T_AltFreq    T_Ref+    T_Ref-    T_Alt+    T_Alt-`
def makeVCF( dataDF, header, vcf):
    vcf_fh = open(vcf, "wb")
    vcf_fh.write(header + "\n")
    allelic_format = "DP:AD:ADF:ADR"
    for idx,row in dataDF.iterrows():
        chrom = row.loc['Chrom']
        pos = row.loc['Start']
        id = "."
        ref = row.loc['Ref']
        alt = row.loc['Alt']
        qual = "."
        filter = "PASS"
        sampleId = row.loc['Sample']
        callmethod = row.loc['CallMethod']
        info = "Sample=" + sampleId + ";" + "CallMethod=" + callmethod
        tdp = row.loc['T_TotalDepth']
        trd = row.loc['T_RefCount']
        tad = row.loc['T_AltCount']
        tvf = row.loc['T_AltFreq']
        trdp = row.loc['T_Ref+'] 
        trdn = row.loc['T_Ref-'] 
        tadp = row.loc['T_Alt+'] 
        tadn = row.loc['T_Alt-'] 
        tdpp = int(trdp) + int(tadp)
        tdpn = int(trdn) + int(tadn)
        tumor = str(tdp) + ":" + str(trd) + "," + str(tad) + ":" + str(trdp) + "," + str(tadp) + ":" + str(trdn) + "," + str(tadn)
        normal = "0:0,0:0,0:0,0"
        vcf_fh.write("\t".join([chrom,pos,id,ref,alt,qual,filter,info,allelic_format,tumor,normal]) + "\n")
        
    vcf_fh.close()
if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("dmpGermline2vcf: Elapsed time was %g seconds", totaltime)
    sys.exit(0)
