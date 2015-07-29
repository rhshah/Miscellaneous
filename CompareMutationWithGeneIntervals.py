import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import math

def main():
    parser = argparse.ArgumentParser(prog='AnnalyzeSVs.py', description='Calculate the stats for the Mutation for Intervals', usage='%(prog)s [options]')
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=True, help="make lots of noise [default]")
    parser.add_argument("-m", "--mInput", action="store", dest="mInput", required=True, metavar='file.txt or listoffiles.list', help="Location of the mutation file or list of files to be used")
    parser.add_argument("-o", "--outputFilePrefix", action="store", dest="outFilePrefix", required=True, metavar='AnnotatedSV', help="Full path with prefix name for the output file")
    parser.add_argument("-g", "--annGeneIntervals", action="store", dest="annGeneIntervals", required=False, default='/dmp/data/mskdata/interval-lists/production/gene_intervals.list.annotated', metavar='Interval.txt', help="List of Annotated Gene Intervals.")
    parser.add_argument("-r", "--ranges", action="store", nargs='*', dest="ranges", required=True, metavar='2,3,50,100', help="Number of ranges to be used to run Analysis; Multiple ranges are separated by space")
    args = parser.parse_args()
    (mDF) = ReadMfile(args.mInput)
    ranges = GetRanges(args.ranges)
    intervalsDF = ReadIntervals(args.annGeneIntervals)
    FilterRecords(intervalsDF,mDF, ranges, args.outFilePrefix)
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
                header = header[0:35]
                for line in filecontent:
                    data = line.strip('\n').split('\t')
                    data = data[0:35]
                    listdf.append(data)
        dataDF = pd.DataFrame(listdf,columns=header)
    else:
        dataDF = pd.read_csv(input, sep='\t', header=0, keep_default_na='True')
    return(dataDF)

def GetRanges(range):
    if (range.__len__() > 1):
        ranges = ",".join(range)
    else:
        ranges = "".join(range)
    return(ranges)

def ReadIntervals(annGeneIntervals):
    Intervals = pd.DataFrame(columns=['chr', 'start', 'end', 'arm', 'gene', 'transcript', 'exon_number','exon_start','exon_end'])
    #Intervals = pd.DataFrame()
    with open(annGeneIntervals,'r') as filecontent:
        count = 0
        for line in filecontent:
            processedData = []
            if(line):
                data = line.rstrip('\n').split("\t")
                (chr,startend) = data[0].split(":")
                (start,end) = startend.split("-")
                #print data[2]
                if(data[2]):
                    if(',' in data[2]):
                        (type1,type2) = data[2].split(',')
                        (gene,transcript,exonNumber,exChr,exon_start,exon_end) = type1.split(":")
                        processedData = [chr,start,end,data[1],gene,transcript,exonNumber,exon_start,exon_end]
                        #Intervals.loc[count,['chr', 'start', 'end', 'arm', 'gene', 'transcript', 'exon_number','exon_start','exon_end']] = processedData
                        Intervals.loc[count] = processedData
                        count = count +1
                        
                        processedData = []
                        (gene,transcript,exonNumber,exChr,exon_start,exon_end) = type2.split(":")
                        processedData = [chr,start,end,data[1],gene,transcript,exonNumber,exon_start,exon_end]
                        #Intervals.loc[count,['chr', 'start', 'end', 'arm', 'gene', 'transcript', 'exon_number','exon_start','exon_end']] = processedData
                        Intervals.loc[count] = processedData
                        count = count +1
                    else:
                        (gene,transcript,exonNumber,exChr,exon_start,exon_end) = data[2].split(":")
                        processedData = [chr,start,end,data[1],gene,transcript,exonNumber,exon_start,exon_end]
                        #Intervals.loc[count,['chr', 'start', 'end', 'arm', 'gene', 'transcript', 'exon_number','exon_start','exon_end']] = processedData
                        Intervals.loc[count] = processedData
                        count = count + 1
                else:
                    gene = None
                    transcript = None
                    exonNumber = None
                    exChr = None
                    exon_start = None
                    exon_end = None
                    processedData = [chr,start,end,data[1],gene,transcript,exonNumber,exon_start,exon_end]
                    #Intervals.loc[count,['chr', 'start', 'end', 'arm', 'gene', 'transcript', 'exon_number','exon_start','exon_end']] = processedData
                    Intervals.loc[count] = processedData
                    count = count + 1
                #count = count + 1   
                #Intervals = Intervals.append(processedData)
                #Intervals['chr', 'start', 'end', 'arm', 'gene', 'transcript', 'exon_number','exon_start','exon_end'] = processedData
    return(Intervals)
            
    
def FilterRecords(annGeneIntervalsDF,mDF,ranges,fileprefix):
    mDFheaders = list(mDF.columns.values)
    orgFilterDF = pd.DataFrame(columns=mDFheaders)
    outorgFilterFile = fileprefix + ".xlsx"
    xlswriter = pd.ExcelWriter(outorgFilterFile)
    iRange = ranges.split(",")
    offTargetDF=pd.DataFrame(columns=mDFheaders)
    orgFilterList = []
    for qcount,qrow in mDF.iterrows():
        sampleName = qrow.loc['Sample']
        if("Pool" in sampleName):
            continue
        chrQ = qrow.loc['Chrom']
        startQ = qrow.loc['Start']
        idxList = annGeneIntervalsDF[annGeneIntervalsDF['chr'] == chrQ].index.tolist()
        for dindex in (idxList):
            startD = annGeneIntervalsDF.iloc[dindex]['start']
            endD = annGeneIntervalsDF.iloc[dindex]['end']
            exonstartD = annGeneIntervalsDF.iloc[dindex]['exon_start']
            exonendD = annGeneIntervalsDF.iloc[dindex]['exon_end']
            if(int(startD) <= int(startQ) <= int(endD)):
                #print qrow
                orgFilterList.append(qcount)
                #orgFilterDF.append(mDF.index[qcount])
                break
    orgFilterDF = mDF.ix[orgFilterList]
    orgFilterDF.to_excel(xlswriter, sheet_name='OrgFilter', index=False)
    by_VC = orgFilterDF.groupby(['VariantClass'])
    by_VF = by_VC['T_AltFreq']
    by_DP = by_VC['T_TotalDepth']
    by_SM = by_VC['Sample']
    orgFilterStatDF = pd.concat([by_SM.count(),by_VF.apply(np.median), by_VF.apply(np.std),by_DP.apply(np.median), by_DP.apply(np.std)],axis=1,keys=['Variant Counts','VF_Median', 'VF_Std', 'DP_Median','DP_Std'])
    titsList = GetTsTi(orgFilterDF, orgFilterStatDF )
    orgFilterStatDF['Ti/Ts'] = pd.Series(titsList,index=orgFilterStatDF.index)
    orgFilterStatDF.to_excel(xlswriter, sheet_name='OrgFilterStats', index=True)
    last = len(iRange) - 1
    for i, val in enumerate(iRange): 
        FilterDF = pd.DataFrame()
        FilterList = []
        FilterDF = pd.DataFrame(columns=mDFheaders)
        for qcount,qrow in mDF.iterrows():
            sampleName = qrow.loc['Sample']
            if("Pool" in sampleName):
                continue
            chrQ = qrow.loc['Chrom']
            startQ = qrow.loc['Start']
            idxList = annGeneIntervalsDF[annGeneIntervalsDF['chr'] == chrQ].index.tolist()
            for dindex in (idxList):
                startD = annGeneIntervalsDF.iloc[dindex]['start']
                endD = annGeneIntervalsDF.iloc[dindex]['end']
                exonstartD = annGeneIntervalsDF.iloc[dindex]['exon_start']
                exonendD = annGeneIntervalsDF.iloc[dindex]['exon_end']
                if((int(startD)-int(val)) <= int(startQ) <= (int(endD)+int(val))):
                    #Get the index that pass this condition
                    FilterList.append(qcount)
                    break
                
                #FilterDF = FilterDF.append(mDF.index[qcount])
        #Make a subset using the index
        #FilterDF = mDF.ix[FilterList]
        keepList = []
        for fIndex in FilterList:
            if(fIndex in orgFilterDF.index):
                continue
            else:
                keepList.append(fIndex)
        FilterDF = mDF.ix[keepList]
        datasheetName = "Filter_" + val
        FilterDF.to_excel(xlswriter, sheet_name=datasheetName, index=False)
        statsheetName = "Filter_" + val + "_Stats"
        by_VC = FilterDF.groupby(['VariantClass'])
        by_VF = by_VC['T_AltFreq']
        by_DP = by_VC['T_TotalDepth']
        by_SM = by_VC['Sample']
        FilterStatDF = pd.concat([by_SM.count(),by_VF.apply(np.median), by_VF.apply(np.std),by_DP.apply(np.median), by_DP.apply(np.std)],axis=1,keys=['Variant Counts','VF_Median', 'VF_Std', 'DP_Median','DP_Std'])
        
        titsList = GetTsTi(FilterDF, FilterStatDF)
        FilterStatDF['Ti/Ts'] = pd.Series(titsList,index=FilterStatDF.index)
        FilterStatDF.to_excel(xlswriter, sheet_name=statsheetName, index=True)
        offTargetList = []
        if(i == last):
            for fcount,frow in mDF.iterrows():
                sampleName = frow.loc['Sample']
                if("Pool" in sampleName):
                    continue
                else:
                    if(fcount in FilterList):
                        continue
                    else:
                        offTargetList.append(fcount)
            offTargetDF = mDF.ix[offTargetList]
    offTargetDF.to_excel(xlswriter, sheet_name='NoInterval', index=False)          
    by_VC = offTargetDF.groupby(['VariantClass'])
    by_VF = by_VC['T_AltFreq']
    by_DP = by_VC['T_TotalDepth']
    by_SM = by_VC['Sample']
    offTargetStatDF = pd.concat([by_SM.count(),by_VF.apply(np.median), by_VF.apply(np.std),by_DP.apply(np.median), by_DP.apply(np.std)],axis=1,keys=['Variant Counts','VF_Median', 'VF_Std', 'DP_Median','DP_Std'])
    titsList = GetTsTi(offTargetDF, offTargetStatDF )
    offTargetStatDF['Ti/Ts'] = pd.Series(titsList,index=offTargetStatDF.index)
    offTargetStatDF.to_excel(xlswriter, sheet_name='NoIntervalStats', index=True)
    xlswriter.save()
    return
    
def GetTsTi(rDF,sDF):
    ticount = 0
    tscount = 0
    titsList = []
    for count,row in rDF.iterrows():
        ref=row.loc['Ref']
        alt = row.loc['Alt']
        if(len(ref)==1 and len(alt)==1):
            if((ref=='C' and alt=='T') or (ref=='T' and alt=='C') or (ref=='A' and ref =='G') or (ref=='G' and ref =='A')):
                ticount = ticount + 1
            else:
                tscount = tscount+1
    titsRatio = 0
    if(tscount>0):
        titsRatio = float(ticount)/float(tscount)    
    else:
        titsRatio = 0
    for i in range(len(sDF.index)):
            titsList.append(titsRatio)
    return(titsList)





if __name__ == "__main__":
    start_time = time.time()  
    main()
    end_time = time.time()
    print("Elapsed time was %g seconds" % (end_time - start_time))