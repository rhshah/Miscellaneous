#!/usr/bin/python
'''
@Description : This tool helps to modify the hotspot24k excel file to make a hotspot coordinate bed file 
@Created :  05/10/2017
@Updated : 05/10/2017
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
   parser.add_argument("-o", "--outDir", action="store", dest="outdir", required=False, type=str, metavar='/somepath/output', help="Full Path to the output dir.")
   
   args = parser.parse_args()
   if(args.verbose):
       logger.info("make_hotspot_file: Started the run for replacing allele counts.")
   (snvDF,indelDF) = read_excel_file(args)
   (snvBed) = process_DF(snvDF)
   (indelBed) = process_DF(indelDF)
   bedList = snvBed + indelBed
   makebed(args, bedList)
   #write_output(args,cleanDF)
   if(args.verbose):
       logger.info("make_hotspot_file: Finished the run for replacing allele counts.")

def read_excel_file(args):
    snvDF = pd.read_excel(args.inputXL, sheetname="single_codon")
    indelDF = pd.read_excel(args.inputXL, sheetname="indels")
    return(snvDF,indelDF)

def process_DF(dataDF):
    bedList = list()
    for index, row in dataDF.iterrows():
        gene = (str(row.loc['Hugo_Symbol'])).rstrip()
        coordinates = (str(row.loc['Genomic_Position'])).rstrip()
        aa = (str(row.loc['Variant_Amino_Acid'])).rstrip()
        coord_dict=dict([x.split("_") for x in coordinates.split("|")])
        for coord in coord_dict.keys():
            (chr,pos) = coord.split(':')
            #print chr,"\t",pos,"\t",pos,"\t",gene,":",aa,"\t","0","\t","+"
            bed = [chr,pos,pos,gene+";"+aa,".","+"]
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
        #bed = pybedtools.create_interval_from_list(i)
        temp.write("\t".join(i) + "\n")
    temp.seek(0)
    newbed = pybedtools.BedTool(temp.name)
    srtbed = newbed.sort()
    #mrgbed = srtbed.merge(d=10)
    srtbed.saveas(outfile)
    # Automatically cleans up the file
    temp.close()

def get_index_for_fillout(filloutDF,m_chr,m_start,m_end,m_ref,m_alt):
    dataDF = filloutDF.loc[filloutDF['Chrom'] == m_chr]
    count = 0
    for index, row in dataDF.iterrows():
        
        vcf_chr = row.loc['Chrom']
        vcf_start = row.loc['Start']
        vcf_ref = (str(row.loc['Ref'])).rstrip()
        vcf_alt = (str(row.loc['Alt'])).rstrip()
        
        #Skip Complex
        if(len(vcf_ref) >= 2 and len(vcf_alt) >= 2 and (len(vcf_ref) != len(vcf_alt))):
            #print "Skipping",vcf_chr,vcf_start,vcf_ref,vcf_alt,m_chr,m_start,m_end,m_ref,m_alt,"\n"
            continue
        start = None
        end = None
        ref = None
        alt = None
        #SNV,DNP,ONP
        if(len(vcf_ref) == len(vcf_alt)):
            start = int(vcf_start)
            end = int(vcf_start) + len(vcf_ref) - 1
            ref = vcf_ref
            alt = vcf_alt
        else:
            #DELETIONS
            if(len(vcf_ref) > len(vcf_alt)):
                ref = vcf_ref[1:]
                if(len(vcf_alt) == 1):
                    alt = "-"
                if(len(vcf_alt) > 1):
                    alt = vcf_alt[1:]
                start = int(vcf_start) + 1
                end = int(vcf_start) + len(ref)
            #INSERTIONS
            if(len(vcf_ref) < len(vcf_alt)):
                if(len(vcf_ref) == 1):
                   ref = "-"
                if(len(vcf_ref) > 1):
                    ref = vcf_ref[1:]
                alt = vcf_alt[1:]
                start = int(vcf_start)
                end = int(vcf_start + 1)
        #print "IN:",vcf_chr,vcf_start,vcf_ref,vcf_alt,"Mod:",start,ref,alt,"MAF:",m_chr,m_start,m_end,m_ref,m_alt,"\n"
        if(m_start == start and m_end == end and m_ref == ref and m_alt == alt):
            #print index
            return(index)
        else:
            continue

def replace_allele_count(args,mafDF,filloutDF):
    mafDF_copy = mafDF.copy()
    for i_index, i_row in mafDF.iterrows():
        m_chr = i_row.loc['Chromosome']
        m_start = i_row.loc['Start_Position']
        m_end = i_row.loc['End_Position']
        m_type = i_row.loc['TYPE']
        m_vt = i_row.loc['Variant_Type']
        m_ref = (str(i_row.loc['Reference_Allele'])).rstrip()
        m_alt = (str(i_row.loc['Tumor_Seq_Allele2'])).rstrip()
        m_tumor_sample_name = i_row.loc['Tumor_Sample_Barcode']
        m_normal_sample_name = i_row.loc['Matched_Norm_Sample_Barcode']
        if(len(m_ref) >= 2 and len(m_alt) >= 2 and (len(m_ref) != len(m_alt))):
            logging.warn("make_hotspot_file: Skipping as complex event: %s, %d, %d, %s, %s", m_chr, m_start, m_end, m_ref, m_alt)
            continue
        else:
            filloutIndex = get_index_for_fillout(filloutDF, m_chr, m_start, m_end, m_ref, m_alt)
            if(filloutIndex != None):
                recordToReplace = filloutDF.iloc[[filloutIndex]]
                tcols = [col for col in recordToReplace.columns if m_tumor_sample_name in col]
                tcolName = tcols[0]
                ncols = [col for col in recordToReplace.columns if m_normal_sample_name in col]
                ncolName = ncols[0]
                t_record = recordToReplace.iloc[0][tcolName]
                tsampleDat=dict([x.split("=") for x in t_record.split(";")])
                n_record = recordToReplace.iloc[0][ncolName]
                nsampleDat=dict([x.split("=") for x in n_record.split(";")])
                mafDF_copy.set_value(i_index,"t_depth",tsampleDat.get("DP"))
                mafDF_copy.set_value(i_index,"t_ref_count",tsampleDat.get("RD"))
                mafDF_copy.set_value(i_index,"t_alt_count",tsampleDat.get("AD"))
                mafDF_copy.set_value(i_index,"n_depth",nsampleDat.get("DP"))
                mafDF_copy.set_value(i_index,"n_ref_count",nsampleDat.get("RD"))
                mafDF_copy.set_value(i_index,"n_alt_count",nsampleDat.get("AD"))
            else:
                logging.warn("make_hotspot_file: Skipping as could not find index in fillout as conversion of vcf to maf failed: %s, %d, %d, %s, %s", m_chr, m_start, m_end, m_ref, m_alt)
                continue
    return(mafDF_copy)

def write_output(args,output_DF):
    if(args.outdir):
        outFile = os.path.join(args.outdir,args.outputMaf)
    else:
        outFile = os.path.join(os.getcwd(),args.outputMaf)
    output_DF.to_csv(outFile, sep='\t', index=False)
    
if __name__ == "__main__":
    start_time = time.time()  
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("make_hotspot_file: Elapsed time was %g seconds", totaltime)
    sys.exit(0)
