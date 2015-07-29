'''
Created On : 06/01/2015
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
from cyordereddict import OrderedDict
from numpy import nan
import matplotlib.gridspec as gridspec
import scipy.spatial.distance as distance
import scipy.cluster.hierarchy as sch
import copy
import brewer2mpl
from gdata.data import GDFeed


def main():
    parser = argparse.ArgumentParser(
        prog='makeComprehensiveMutationList.py',
        description='Write Per Patient Mutation in all samples and write heatmap files',
        usage='%(prog)s [opFtions]')
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        default=True,
        help="make lots of noise [default]")
    parser.add_argument(
        "-m",
        "--mutationFile",
        action="store",
        dest="mutationFile",
        required=True,
        metavar='ImpactVaraints.txt',
        help="Location of mutation File")
    parser.add_argument(
        "-t",
        "--titleFile",
        action="store",
        dest="infoFile",
        required=True,
        metavar='title_file.txt',
        help="Location of the title File")
    parser.add_argument(
        "-o",
        "--outputFilePrefix",
        action="store",
        dest="outFilePrefix",
        required=True,
        metavar='AnnotatedSV',
        help="Full path with prefix name for the output file")
    args = parser.parse_args()
    (mutationsDF) = readExlFile(args.mutationFile)
    (sampleinfoDF) = readTxtFile(args.infoFile)
    (patientTosampleDF, sampleDF) = processSampleDF(sampleinfoDF)
    (newMutationDF) = modifyMutationDF(mutationsDF.copy(), sampleinfoDF, patientTosampleDF, args)
    WriteOutput(newMutationDF, "newMutations", args)
    (sampleinfoDF) = readTxtFile(args.infoFile)
    (patientTosampleDF, sampleDF) = processSampleDF(sampleinfoDF)
    (geneLevelDF) = makeGeneLevelDF(mutationsDF.copy(), sampleinfoDF, patientTosampleDF, args)
    WriteOutput(geneLevelDF, "GeneLevel", args)
    (sampleinfoDF) = readTxtFile(args.infoFile)
    (patientTosampleDF, sampleDF) = processSampleDF(sampleinfoDF)
    makePateintLevelDF(mutationsDF.copy(), sampleinfoDF, patientTosampleDF, args)
    (sampleinfoDF) = readTxtFile(args.infoFile)
    (patientTosampleDF, sampleDF) = processSampleDF(sampleinfoDF)
    #makePateintLevelDFForLichee(mutationsDF.copy(), sampleinfoDF, patientTosampleDF, args)
    #geneLevelDFcopy = geneLevelDF.copy()
    # makeHeatMaps(geneLevelDFcopy)
    


def readTxtFile(inFile):
    df = pd.read_csv(inFile, sep='\t', header=0, keep_default_na='True')
    print "Finished Reading Text File\n"
    return(df)


def readExlFile(inFile):
    df = pd.io.excel.read_excel(inFile, sheetname=0, keep_default_na='True')
    print "Finished Reading Excel File\n"
    return(df)


def processSampleDF(df):
    groupsOfPatients = df.groupby('Patient_ID')
    # Remove Normal Records
    gDF = groupsOfPatients.apply(lambda g: g[g['Class'] == 'Tumor'])
    # Remove Patients with only 1 sample
    gDF = gDF.groupby('Patient_ID').filter(lambda g: len(g) > 1)
    groupOfMultiSamplePatientsDF = gDF.groupby('Patient_ID')
    groupOfSamplesDF = df.groupby('Sample_ID')
    return(groupOfMultiSamplePatientsDF, groupOfSamplesDF)


def modifyMutationDF(mutationsDF, sampleDF, patientTosampleDF, args):
    sampleuniqdict = OrderedDict()
    names = mutationsDF.columns.values
    newMutationDF = pd.DataFrame(columns=names)
    indexForNewDF = 0
    for count, row in mutationsDF.iterrows():
        sampleName = str(row.loc['Sample'])
        sampleIndex = int(sampleDF[sampleDF['Sample_ID'] == sampleName].index)
        patientID = (sampleDF.iloc[sampleIndex]['Patient_ID'])
        # print sampleIndex, type(sampleDF), str(patientID)
        try:
            patientTosampleGrp = patientTosampleDF.get_group(patientID)
        except KeyError:
            if(args.verbose):
                print "The patient ID for the sample is not in the group:", patientID, "\n"
            continue
        NormalName = str(row.loc['NormalUsed'])
        chr = str(row.loc['Chrom'])
        start = str(row.loc['Start'])
        ref = str(row.loc['Ref'])
        alt = str(row.loc['Alt'])
        key = patientID + ";" + chr + ";" + start + ";" + ref + ";" + alt
        # print key,"\n"
        if(key in sampleuniqdict):
            if(args.verbose):
                print "Data Already Processed for the mutation:", key, "\n"
            continue
        else:
            for icount, irow in patientTosampleGrp.iterrows():
                newrecord = row.copy()
                gSampleID = str(irow.loc['Sample_ID'])
                (dp, rd, ad, vf) = (str(row.loc[gSampleID])).split(";")
                (n, dp) = dp.split("=")
                (n, rd) = rd.split("=")
                (n, ad) = ad.split("=")
                (n, vf) = vf.split("=")
                # print gSampleID,dp,rd,ad,vf,"\n"
                newrecord.loc['Sample'] = gSampleID
                newrecord.loc['T_TotalDepth'] = dp
                newrecord.loc['T_RefCount'] = rd
                newrecord.loc['T_AltCount'] = ad
                newrecord.loc['T_AltFreq'] = vf
                # print
                # newrecord.loc['Sample'],newrecord.loc['T_TotalDepth'],newrecord.loc['T_AltFreq']
                newMutationDF.loc[indexForNewDF] = newrecord
                indexForNewDF = indexForNewDF + 1
                # print newMutationDF
        sampleuniqdict[key] = NormalName
    return(newMutationDF)


def makeGeneLevelDF(mutationsDF, sampleDF, patientTosampleDF, args):
    if(args.verbose):
        print "Now running gene-level analysis\n"
    sampleuniqdict = OrderedDict()
    samplesToTraverse = sorted(sampleDF['Sample_ID'].tolist())
    
    colsForDF = copy.copy(samplesToTraverse)
    colsForDF.insert(0, "Gene-AminoAcid")
    geneLevelDF = pd.DataFrame(columns=colsForDF)
    gene_aa_dict = OrderedDict()
    # print geneLevelDF,"\n"
    geneIndex = 0
    for count, row in mutationsDF.iterrows():
        sampleName = str(row.loc['Sample'])
        #if((sampleName == "s-FFPE-Pooled-Tumor") or (sampleName == "s-EV-crc-043-P2") or (sampleName == "s-EV-crc-070-M3") or (sampleName == "s-EV-crc-039-P3") or (sampleName == "s-EV-crc-036-P3") or (sampleName == "s-EV-crc-058-M6") or (sampleName == "s-EV-crc-034-P3")):
         #   print "Skipping:", sampleName,"\n"
          #  continue
        geneName = str(row.loc['Gene'])
        if((geneName == "PDGFRA") or (geneName == "MEN1")):
            print "Skipping:", sampleName,"\n"
            continue
        callConfidence = str(row.loc['Call_Confidence'])
        if(callConfidence != 'HIGH'):
             print "Skipping:", sampleName,"\n"
             continue
        aachange = (str(row.loc['AAchange']))
        cdnaChange = (str(row.loc['cDNAchange']))
        NormalName = str(row.loc['NormalUsed'])
        chr = str(row.loc['Chrom'])
        start = str(row.loc['Start'])
        ref = str(row.loc['Ref'])
        alt = str(row.loc['Alt'])
        comment = str(row.loc['Comments'])
        if(aachange != "nan"):
                (p, aachange) = (str(row.loc['AAchange'])).split(".")
        else:
            aachange = "NA"
        if(cdnaChange != "nan"):
                (p, cdnaChange) = (str(row.loc['cDNAchange'])).split(".")
        else:
            cdnaChange = "NA"
        if(aachange == "NA"):
            gene_aa = geneName + "-" + cdnaChange 
        else:
            gene_aa = geneName + "-" + aachange 
        gene_aa_dict[gene_aa] = 0
        key =  gene_aa + ";" + chr + ";" + start + ";" + ref + ";" + alt
        # print key,"\n"
        if((key in sampleuniqdict)or ("Germline" in comment)):
            if(args.verbose):
                print "Data Already Processed for the mutation:", key, "\n"
            continue
        else:
            valList = []
            if(gene_aa in gene_aa_dict):
                    new_gene_aa = gene_aa + "-" + str(geneIndex)
                    gene_aa_dict[new_gene_aa] = gene_aa_dict[gene_aa]
            valList.append(new_gene_aa)
            for sampleName in samplesToTraverse:
                (dp, rd, ad, vf) = (str(row.loc[sampleName])).split(";")
                (n, dp) = dp.split("=")
                (n, rd) = rd.split("=")
                (n, ad) = ad.split("=")
                (n, vf) = vf.split("=")
                valList.append(vf)
            sampleuniqdict[key] = NormalName
            geneLevelDF.loc[geneIndex] = valList
            geneIndex = geneIndex + 1
    return(geneLevelDF)

def makePateintLevelDF(mutationsDF, sampleDF, patientTosampleDF, args):
    if(args.verbose):
        print "Making per patient heatmap Files\n"
    #sampleuniqdict = OrderedDict()
    for index,items in enumerate(patientTosampleDF):
        pid = items[0]
        df = items[1]
        samplesToTraverse = sorted(df['Sample_ID'].tolist())
        outsuffix = pid + "-" + "heatmapData"
        colsForDF = copy.copy(samplesToTraverse)
        colsForDF.insert(0, "Gene-AminoAcid")
        geneLevelDF = pd.DataFrame(columns=colsForDF)
        gene_aa_dict = OrderedDict()
        geneIndex = 0
        sampleuniqdict = OrderedDict()
        for count, row in mutationsDF.iterrows():
            sampleName = str(row.loc['Sample'])
            if(sampleName in samplesToTraverse):
                geneName = str(row.loc['Gene'])
                if((geneName == "PDGFRA") or (geneName == "MEN1")):
                    print "Skipping:", sampleName,"\n"
                    continue
                callConfidence = str(row.loc['Call_Confidence'])
                if(callConfidence != 'HIGH'):
                    print "Skipping:", sampleName,"\n"
                    continue
                aachange = (str(row.loc['AAchange']))
                cdnaChange = (str(row.loc['cDNAchange']))
                NormalName = str(row.loc['NormalUsed'])
                comment = str(row.loc['Comments'])
                chr = str(row.loc['Chrom'])
                start = str(row.loc['Start'])
                ref = str(row.loc['Ref'])
                alt = str(row.loc['Alt'])
                if(aachange != "nan"):
                        (p, aachange) = (str(row.loc['AAchange'])).split(".")
                else:
                    aachange = "NA"
                if(cdnaChange != "nan"):
                        (p, cdnaChange) = (str(row.loc['cDNAchange'])).split(".")
                else:
                    cdnaChange = "NA"
                if(aachange == "NA"):
                    gene_aa = geneName + "-" + cdnaChange 
                else:
                    gene_aa = geneName + "-" + aachange 
                gene_aa_dict[gene_aa] = 0
                key =  gene_aa + ";" + chr + ";" + start + ";" + ref + ";" + alt
                # print key,"\n"
                if((key in sampleuniqdict) or ("Germline" in comment)):
                    if(args.verbose):
                        print "Data Already Processed for the mutation:", key, "\n"
                    continue
                else:
                    valList = []
                    if(gene_aa in gene_aa_dict):
                        new_gene_aa = gene_aa + "-" + str(geneIndex)
                        gene_aa_dict[new_gene_aa] = gene_aa_dict[gene_aa]
                    valList.append(new_gene_aa)
                    for sampleName in samplesToTraverse:
                        (dp, rd, ad, vf) = (str(row.loc[sampleName])).split(";")
                        (n, dp) = dp.split("=")
                        (n, rd) = rd.split("=")
                        (n, ad) = ad.split("=")
                        (n, vf) = vf.split("=")
                        valList.append(vf)
                    sampleuniqdict[key] = NormalName
                    geneLevelDF.loc[geneIndex] = valList
                    geneIndex = geneIndex + 1
        WriteOutput(geneLevelDF, outsuffix, args)
    return
   
def makePateintLevelDFForLichee(mutationsDF, sampleDF, patientTosampleDF, args):
    if(args.verbose):
        print "Making per patient Lichee Files\n"
    #sampleuniqdict = OrderedDict()
    for index,items in enumerate(patientTosampleDF):
        pid = items[0]
        df = items[1]
        samplesToTraverse = sorted(df['Sample_ID'].tolist())
        outsuffix = pid + "-" + "Lichee"
        colsForDF = copy.copy(samplesToTraverse)
        colsForDF.insert(0, "Normal")
        colsForDF.insert(0, "EOA")
        colsForDF.insert(0, "alt")
        colsForDF.insert(0, "ref")
        colsForDF.insert(0, "pos")
        colsForDF.insert(0, "#chr")
        geneLevelDF = pd.DataFrame(columns=colsForDF)
        gene_aa_dict = OrderedDict()
        geneIndex = 0
        sampleuniqdict = OrderedDict()
        for count, row in mutationsDF.iterrows():
            sampleName = str(row.loc['Sample'])
            if(sampleName in samplesToTraverse):
                valList = []
                geneName = str(row.loc['Gene'])
                aachange = (str(row.loc['AAchange']))
                cdnaChange = (str(row.loc['cDNAchange']))
                NormalName = str(row.loc['NormalUsed'])
                comment = str(row.loc['Comments'])
                chr = str(row.loc['Chrom'])
                start = str(row.loc['Start'])
                ref = str(row.loc['Ref'])
                alt = str(row.loc['Alt'])
                valList.append(chr)
                valList.append(start)
                valList.append(ref)
                valList.append(alt)
                valList.append(geneName)
                valList.append("0")
                '''
                if(aachange != "nan"):
                        (p, aachange) = (str(row.loc['AAchange'])).split(".")
                else:
                    aachange = "NA"
                if(cdnaChange != "nan"):
                        (p, cdnaChange) = (str(row.loc['cDNAchange'])).split(".")
                else:
                    cdnaChange = "NA"
                if(aachange == "NA"):
                    gene_aa = geneName + "-" + cdnaChange 
                else:
                    gene_aa = geneName + "-" + aachange 
                gene_aa_dict[gene_aa] = 0
                '''
                key =  geneName + ";" + chr + ";" + start + ";" + ref + ";" + alt
                # print key,"\n"
                if((key in sampleuniqdict) or ("Germline" in comment)):
                    if(args.verbose):
                        print "Data Already Processed for the mutation:", key, "\n"
                    continue
                else:
                    
                    '''
                    if(gene_aa in gene_aa_dict):
                        new_gene_aa = gene_aa + "-" + str(geneIndex)
                        gene_aa_dict[new_gene_aa] = gene_aa_dict[gene_aa]
                    '''
                    
                    for sampleName in samplesToTraverse:
                        (dp, rd, ad, vf) = (str(row.loc[sampleName])).split(";")
                        (n, dp) = dp.split("=")
                        (n, rd) = rd.split("=")
                        (n, ad) = ad.split("=")
                        (n, vf) = vf.split("=")
                        valList.append(vf)
                    sampleuniqdict[key] = NormalName
                    geneLevelDF.loc[geneIndex] = valList
                    geneIndex = geneIndex + 1
        WriteOutput(geneLevelDF, outsuffix, args)
    return
     
def makeHeatMaps(geneLevelDF):
    axi = plt.imshow(geneLevelDF, interpolation='nearest', cmap=plt.cm.RdBu)
    ax = axi.get_axes()
    clean_axis(ax)


def WriteOutput(df, type, args):
    outputTxt = args.outFilePrefix + "-" + type + ".txt"
    outputExl = args.outFilePrefix + "-" + type + ".xlsx"
    # Print to TSV file
    df.to_csv(outputTxt, sep='\t', index=False)
    # Print to Excel
    df.to_excel(outputExl, sheet_name='Annotated_SVs', index=False)
    return
# helper for cleaning up axes by removing ticks, tick labels, frame, etc.


def clean_axis(ax):
    """Remove ticks, tick labels, and frame from axis."""
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    print("Elapsed time was %g seconds" % (end_time - start_time))
