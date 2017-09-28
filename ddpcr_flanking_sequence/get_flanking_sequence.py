"""
get_flanking_sequence
~~~~~~~~~~~~~~~~~~~~~

:Description: This module gets the information about 200bp flanking sequence given chromosome, start, end, reference allele and alternate allele based on maf

"""

'''
Created on February 1, 2017
Description: This module gets the information about 200bp flanking sequence given chromosome, start, end, reference allele and alternate allele based on maf
@author: Ronak H Shah
::Input::
chr : chromosome location of the alteration
position: position of the alteration
reference_allele: reference allele information
variant_allele: variant allele information
protein_ann: variant annotation for protein change
cDNA_ann: variant annotation for cDNA change
vcf: location of gzipped population frequency vcf
::Output::
a fasta file with following header
>protein_ann-cDNA_ann

::Example Run::
```
python get_flanking_sequence.py -chr 11 -s 64577305 -e 64577307 -ref AGC -alt GCTT -p_ann R92Qfs*25 -c_ann 120_125delins -r /Users/shahr2/Documents/PubData/Genome/hg19/Homo_sapiens_assembly19.fasta -vcf /Users/shahr2/Documents/PubData/Genome/hg19/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz -o /Users/shahr2/Desktop/ -of test_gfs
```
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
logger = logging.getLogger('get_flanking_sequence')
try:
    import coloredlogs
    coloredlogs.install(level='DEBUG')
except ImportError:
    logger.warning("get_flanking_sequence: coloredlogs is not installed, please install it if you wish to see color in logs on standard out.")
    pass

try:
    import pysam
except ImportError:
    logger.fatal("get_flanking_sequence: pysam is not installed, please install pysam as it is required to run the mapping. pysam version > 0.10.0")
    sys.exit(1)

try:
    import vcf as vcfpackage
except ImportError:
    logger.fatal("get_flanking_sequence: pyvcf is not installed, please install pyvcf as it is required to run the mapping. pysam version > 0.10.0")
    sys.exit(1)


#Run all sub function    
def main():
    parser = argparse.ArgumentParser(prog='get_flanking_sequence.py', description='This module gets the information about 200bp flanking sequence given chromosome, start, end, reference allele and alternate allele based on maf specification', usage='%(prog)s [options]')
    parser.add_argument("-chr", "--chromosome", action="store", dest="chr", type=str, required=True, metavar='1', help="Enter the information for what chromosome the event is on",choices=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]) 
    #parser.add_argument("-of", "--outputfileprefix", action="store", dest="outFilePrefix", required=True, metavar='OutFilePrefix', help="Output File Prefix for the flanking fasta file.")
    parser.add_argument("-v", "--verbose", action="store", required=True, dest="verbose",default="TRUE",metavar="TRUE", help="make lots of noise",choices=["TRUE","FALSE"])
    parser.add_argument("-s", "--start", action="store", dest="start", type=int, required=True, metavar=5, help="Start coordinate for the event")
    parser.add_argument("-e", "--end", action="store", dest="end", type=int, required=True, metavar=20, help="End coordinate for the event")
    parser.add_argument("-f", "--flank", action="store", dest="flank", type=int, required=False, default = 200, metavar=200, help="Number of bases to flank for the output [default=200]")
    parser.add_argument("-r", "--referenceFile", action="store", dest="refgenome", required=True, default="hg19",metavar='hg19', help="Full Path to the reference file with the fasta index.", choices=["hg19"])
    parser.add_argument("-ref", "--reference_allele", action="store", dest="ref", type=str, required=True, metavar='someref', help="reference allele for the event")
    parser.add_argument("-alt", "--alternate_allele", action="store", dest="alt", type=str, required=True, metavar='somealt', help="alternate allele for the event")
    parser.add_argument("-p_ann", "--protein_annotation", action="store", dest="pann", type=str, required=True, metavar='p.ann', help="Protein annotation as string")
    parser.add_argument("-c_ann", "--cDNA_annotations", action="store", dest="cann", type=str, required=True, metavar='c.ann', help="cDNA annotation as string")
    parser.add_argument("-gene", "--gene_name", action="store", dest="gene", type=str, required=True, metavar='GENE', help="GENE name as a string")
    #parser.add_argument("-o", "--outDir", action="store", dest="outdir", required=True, type=str, metavar='/somepath/output', help="Full Path to the output dir.")
    parser.add_argument("-vcf", "--pop_vcf", action="store", dest="vcf", type=str, required=False, default="/Users/shahr2/Documents/PubData/Genome/hg19/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz",metavar='/somepath/pop_vcf', help="Full Path to gzipped population frequency vcf")
    args = parser.parse_args()
    if(args.verbose == "TRUE"):
    	args.verbose = True
    else:
    	args.verbose = False
    #Get proper paths
    if(args.refgenome):
        if(args.refgenome == "hg19"):
            if(os.path.isdir("/ifs/")):
                args.refgenome = "/ifs/depot/resources/dmp/data/pubdata/hg-fasta/VERSIONS/hg19/Homo_sapiens_assembly19.fasta"
            else:
                args.refgenome = "/Users/shahr2/Documents/PubData/Genome/hg19/Homo_sapiens_assembly19.fasta"
        if(args.refgenome == "hg38"):
            if(os.path.isdir("/ifs/")):
                args.refgenome = "/ifs/depot/assemblies/H.sapiens/hg38/hg38.fasta"
            else:
                args.refgenome = "/Users/shahr2/Documents/PubData/Genome/hg19/Homo_sapiens_assembly19.fasta"
    else:
        if(args.verbose):
            logger.fatal("get_flanking_sequence: The reference genome input is not correct. Please input proper reference genome information")
        sys.exit(1)
    if(args.vcf):
        if(os.path.isdir("/opt/common/CentOS_6-dev")):
            args.vcf = "/opt/common/CentOS_6-dev/vep/cache/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz"
        else:
            pass
    else:
        if(args.verbose):
            logger.fatal("get_flanking_sequence: The population vcf input is not correct. Please input proper population vcf information")
        sys.exit(1)
    #Run Modules
    validate_inputs(args.chr, args.start, args.end, args.ref, args.alt, args.flank, args.refgenome, args.vcf, args.verbose)
    (mflank_left,mflank_right,flank_left,flank_right) = get_masked_flank_seq(args.chr, args.start, args.end, args.flank,args.refgenome,args.vcf,args.verbose)
    print_output(args.verbose,args.gene,args.pann,args.cann,mflank_left,mflank_right,flank_left,flank_right,args.ref,args.alt)
    #(outfile) = write_final_output(args.verbose,args.gene,args.pann,args.cann,mflank_left,mflank_right,flank_left,flank_right,args.ref,args.alt)
    if(args.verbose):
        logger.info("get_flanking_sequence: Finished generating flanking sequence")
        #logger.info("get_flanking_sequence: Flanking sequence is written in %s", outfile)

#Print output
def print_output(verbose,gene,pann,cann,mflank_left,mflank_right,flank_left,flank_right,ref,alt):
    header1 = ">" + gene + "-"+ pann + "-" + cann + "\n"
    record1 = str(flank_left) + "[" + ref + "/" + alt + "]" + str(flank_right) + "\n"
    header2 = ">Masked-" + gene + "-"+ pann + "-" + cann + "\n"
    record2 = str(mflank_left) + "[" + ref + "/" + alt + "]" + str(mflank_right) + "\n"
    print(header1 + record1 + header2 + record2)
    return()
#Write the final output       
def write_final_output(verbose,gene,pann,cann,mflank_left,mflank_right,flank_left,flank_right,ref,alt,prefix,outdir):
    outfile = os.path.join(outdir, prefix + ".fasta")
    fh = open(outfile,"w+")
    header1 = ">" + gene + "-"+ pann + "-" + cann + "\n"
    fh.write(header1)
    record1 = str(flank_left) + "[" + ref + "/" + alt + "]" + str(flank_right) + "\n"
    fh.write(record1)
    header2 = ">Masked-" + gene + "-"+ pann + "-" + cann + "\n"
    fh.write(header2)
    record2 = str(mflank_left) + "[" + ref + "/" + alt + "]" + str(mflank_right) + "\n"
    fh.write(record2)
    fh.close()
    return(outfile)

#Get the coordinates to mask
def get_masking_coord(vcf,chr,start,end):
    vcf_reader = vcfpackage.Reader(open(vcf, 'r'))
    positionToMask = []
    for record in vcf_reader.fetch(chr, start, end):
        if(record.is_snp):
            positionToMask.append(record.POS)
    return(positionToMask)


#Get the flanking sequence and mask for snps
def get_masked_flank_seq(chr,start,end,flank,refgenome,vcf,verbose):
    #Make pysam object of fasta
    fa = pysam.FastaFile(refgenome)
    #leftside
    #figure out start and end on left
    p_start_left = start - flank - 1
    if(p_start_left < 0):
        p_start = 1
    p_end_left = start-1
    flank_left = fa.fetch(chr, p_start_left, p_end_left)
    if(verbose):
        logging.info("get_flanking_sequence: Flanking sequence to left:%s", flank_left)
    left_pos_to_mask = get_masking_coord(vcf, chr, p_start_left, p_end_left)
    m_flank_left = flank_left
    if(len(left_pos_to_mask) > 1):
        for coord in left_pos_to_mask:
            index = abs(p_start_left-coord)
            m_flank_left = m_flank_left[:index-1] + "N" + m_flank_left[index:]
    if(verbose):
        logging.info("get_flanking_sequence: Masked sequence to left:%s", m_flank_left)
    #Righside 
    #figure out start and end on right
    p_start_right = end  
    p_end_right = end + flank 
    flank_right = fa.fetch(chr, p_start_right, p_end_right)
    if(verbose):
        logging.info("get_flanking_sequence: Flanking sequence to right:%s", flank_right)
    right_pos_to_mask = get_masking_coord(vcf, chr, p_start_right, p_end_right)
    m_flank_right = flank_right
    if(len(right_pos_to_mask) > 1):
        for coord in right_pos_to_mask:
            index = abs(p_start_right-coord)
            #print index,p_start_right,coord
            m_flank_right = m_flank_right[:index-1] + "N" + m_flank_right[index:]
            #print m_flank_right
    if(verbose):
        logging.info("get_flanking_sequence: Masked sequence to right:%s", m_flank_right)
    return(m_flank_left,m_flank_right,flank_left,flank_right)
    
#Check the inputs
def validate_inputs(chr,start,end,ref,alt,flank,refgenome,vcf,verbose):
    if(chr is 'X' or chr is 'Y'):
        pass
    else:
        chr = int(chr)
        if(1 <= chr <= 22):
            pass
        else:
            if(verbose):
                logger.fatal("get_flanking_sequence: The input for chromosome is not valid. Please input value between 1-22 or X or Y")
            sys.exit(1)
    if(isinstance( start, (int) )):
        pass
    else:
        if(verbose):
            logger.fatal("get_flanking_sequence: The input for start is not valid. Please input proper co-ordinate")
        sys.exit(1)
    if(isinstance( end, (int) )):
        pass
    else:
        if(verbose):
            logger.fatal("get_flanking_sequence: The input for end is not valid. Please input proper co-ordinate")
        sys.exit(1)
    if(isinstance( flank, (int) )):
        pass
    else:
        if(verbose):
            logger.fatal("get_flanking_sequence: The input for flank is not valid. Please input proper integer value as an input")
        sys.exit(1)
    if(len(ref) > 1):
        reflist = list(ref)
        for sref in reflist:
            if(sref in 'ATGC-'):
                pass
            else:
                if(verbose):
                    logger.fatal("get_flanking_sequence: The input for ref is not valid. Please input proper reference allele")
                sys.exit(1)
    else:
        if(ref in 'ATGC'):
            pass
        else:
            if(verbose):
                logger.fatal("get_flanking_sequence: The input for ref is not valid. Please input proper reference allele")
            sys.exit(1)
    if(len(alt) > 1):
        altlist = list(alt)
        for salt in altlist:
            if(salt in 'ATGC-'):
                pass
            else:
                if(verbose):
                    logger.fatal("get_flanking_sequence: The input for alt is not valid. Please input proper alternate allele")
                sys.exit(1)
    else:
        if(alt in 'ATGC-'):
            pass
        else:
            if(verbose):
                logger.fatal("get_flanking_sequence: The input for alt is not valid. Please input proper alternate allele")
            sys.exit(1)
    
    if(os.path.isfile(refgenome)):
        pass
    else:
        if(verbose):
            logger.fatal("get_flanking_sequence: The reference genome input location is not a file. Please input proper reference genome file")
        sys.exit(1)
    if(os.path.isfile(vcf)):
        pass
    else:       
        if(verbose):
            logger.fatal("get_flanking_sequence: The vcf input location is not a file. Please input proper vcf file")
        sys.exit(1)
        
'''            
    if(isinstance( outFilePrefix, (str) )):
        pass
    else:
        if(verbose):
            logger.fatal("get_flanking_sequence: The outputfileprefix is not a string . Please input proper ExAC vcf file")
        sys.exit(1)
    if(os.path.isdir(outdir)):
        pass
    else:
        if(verbose):
            logger.fatal("get_flanking_sequence: The outdir location is not a valid directory. Please input proper output directory")
        sys.exit(1)
'''
               
# Run the whole script
if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("get_flanking_sequence: Elapsed time was %g seconds", totaltime)
    sys.exit(0)

