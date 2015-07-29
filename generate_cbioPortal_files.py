# This script will combine all the IMPACT clinical data based on the sample sheet and generate necessary cBio files

import sys, csv, subprocess, os.path, re, argparse
from collections import defaultdict
parser = argparse.ArgumentParser(description="Prepare .MAF file containing MSK-IMPACT sample list")
#parser.add_argument('-s', '--sample_list', dest='SAMPLE_LIST')
parser.add_argument('-v', '--variant_file_list', dest='VARIANT_FILE_LIST')

args = parser.parse_args()

#SAMPLE_LIST = args.SAMPLE_LIST
VARIANT_FILE_LIST = args.VARIANT_FILE_LIST

def parse(filename, delimiter = "\t"):
    file = open(filename, 'rU')
    tabbed_data = csv.reader(file, delimiter = delimiter)
    fields = tabbed_data.next()
    parsed_data = []
    for row in tabbed_data:
        parsed_data.append(dict(zip(fields, row)))
    file.close()
    return parsed_data

def ls(path = '.', args = ''):
     """ Invoke the shell ls command with args on path""" 
     a = subprocess.Popen(' ls' + ' ' + args + ' ' + path, stdout = subprocess.PIPE, shell = True)
     (output, err) = a.communicate()
     return output

def main():
    #sample_list = parse(SAMPLE_LIST)
    #print "Sample list is parsed..."
    #clinical_data = open("data_clinical.txt", "w")
    #clinical_data_writer = csv.writer(clinical_data, delimiter = "\t", quoting = csv.QUOTE_NONE)
    #clinical_data_writer.writerow(["SAMPLE_ID", "PATIENT_ID", "CANCER_TYPE", "SAMPLE_TYPE", "SAMPLE_CLASS", "METASTATIC_SITE", "PRIMARY_SITE", "CANCER_TYPE_DETAILED", "KNOWN_MOLECULAR_CLASSIFIER"])  # header for clinical_data file]
    
    analysis_output = open("allvariants.txt", "w")
    analysis_output_writer = csv.writer(analysis_output, delimiter = "\t", quoting = csv.QUOTE_NONE)
    analysis_output_writer.writerow(["Sample_ID", "Chr", "Pos", "Ref", "Alt", "Gene", "cDNAchange", "AAchange", "VF", "Exon", "VariantClass", "GeneralTumorType", "DetailedTumorType", "M-number", "SurgPathNumber", "Run"])

    variant_data = open("data_mutations_extended.txt", "w")
    variant_data_writer = csv.writer(variant_data, delimiter = "\t", quoting = csv.QUOTE_NONE)
    
    #all_sample_list = ["#sequenced_samples:"]
    #cbio_id_lookup = {}
    #for y in sample_list:
        #cbio_id_lookup[y['SampleSheetSampleID']] = y['SAMPLE_ID']
        #clinical_data_writer.writerow([y['SAMPLE_ID'], y['PATIENT_ID'], y['CancerType'], y['SampleType'], y['SampleClass'], y['MetastasisSite'], y['PrimarySite'], y['DetailedCancerType'], y['KnownMolecularClassifier']])  # This prints out the clinical info
        #all_sample_list.append(y['SAMPLE_ID'])
    #variant_data_writer.writerow([" ".join(all_sample_list)])
    variant_data_writer.writerow([ "Hugo_Symbol"  , "Entrez_Gene_Id"   , "Center"  , "NCBI_Build"  , "Chromosome"  , "Start_Position"   ,
                                   "End_Position"   , "Strand" , "Variant_Classification"   , "Variant_Type"  , "Reference_Allele" , "Tumor_Seq_Allele1"  ,
                                   "Tumor_Seq_Allele2" , "dbSNP_RS"  , "dbSNP_Val_Status" , "Tumor_Sample_Barcode"   , "Matched_Norm_Sample_Barcode" ,
                                   "Match_Norm_Seq_Allele1"   , "Match_Norm_Seq_Allele2"  , "Tumor_Validation_Allele1"  , "Tumor_Validation_Allele2"  ,
                                   "Match_Norm_Validation_Allele1"  , "Match_Norm_Validation_Allele2"   , "Verification_Status"   , "Validation_Status"   ,
                                   "Mutation_Status"  , "Sequencing_Phase"   , "Sequence_Source" , "Validation_Method"  , "Score"   , "BAM_File" ,
                                   "Sequencer"   , "t_ref_count"  , "t_alt_count"   , "n_ref_count"  , "n_alt_count"  , "DMP_cDNA_annotation" , "DMP_AA_annotation" ]) 
    # This part prints out the variant list
    with open(VARIANT_FILE_LIST, "r") as file:
        for line in file:
            file_path = line.strip()
            file_name = os.path.basename(file_path)
            #IMPACTv3-CLIN-20140001_AllSomaticMutInde
            variant_data = parse(file_path)
            for x in variant_data:
                #analysis_output_writer.writerow([x['Sample'], x['Chrom'], x['Start'], x['Ref'], x['Alt'], x['Gene'], x['cDNAchange'], x['AAchange'], x['T_AltFreq'], x['Exon'], x['VariantClass'], y['CancerType'], y['DetailedCancerType'], y['M-Number'], y['SurgPathNumber'], run_number])
                
                cbio_id = x['Sample']
                sampleComment = x['Comments']
                if('Likely Germline' in sampleComment):
                    sampleComment = 'Germline'
                if('Somatic' in sampleComment):
                    sampleComment = 'Somatic'
                # for SNVs
                if len(x['Ref']) == len(x['Alt']):
                    if(len(x['Ref']) == 1):
                        start_position = x['Start']
                        end_position = x['Start']
                        var_type = "SNP"
                        ref = x['Ref']
                        alt = x['Alt']
                    elif(len(x['Ref']) == 2):
                        start_position = x['Start']
                        end_position = int(x['Start']) + len(x['Ref']) -1
                        var_type = "DNP"
                        ref = x['Ref']
                        alt = x['Alt']
                    elif(len(x['Ref']) == 3):
                        start_position = x['Start']
                        end_position = int(x['Start']) + len(x['Ref']) -1
                        var_type = "TNP"
                        ref = x['Ref']
                        alt = x['Alt']
                elif len(x['Ref']) > len(x['Alt']):
                    var_type = "DEL"
                    ref = x['Ref'][1:]
                    if(len(x['Alt']) == 1):
                        alt = "-"
                    elif(len(x['Alt']) > 1):
                        alt = x['Alt'][1:]
                    start_position = x['Start']
                    end_position = int(x['Start']) + len(ref) - 1
                elif len(x['Ref']) < len(x['Alt']):
                    var_type = "INS"
                    if(len(x['Ref']) == 1):
                        ref = "-"
                    elif(len(x['Ref']) > 1):
                        ref = x['Ref'][1:]
                    alt = x['Alt'][1:]
                    start_position = x['Start']
                    end_position = int(x['Start']) + 1
                variant_data_writer.writerow([x['Gene'], "", "MSK-IMPACT", "hg19", x['Chrom'], start_position, end_position, "+" , x['VariantClass'], var_type, ref, alt, "", "", "", cbio_id,
                                              "", "", "", "", "", "", "", "", "", sampleComment, "", "", "", "", "", "HiSeq-2500", x['T_RefCount'], x['T_AltCount'], x['N_RefCount'], x['N_AltCount'], x['cDNAchange'], x['AAchange']])
                
    # This part combines the .seg files together
    print "Variant calls are processed...\n"
'''    
    print "Segmentation files are being processed..."
    seg_files = ls("/dmp/hot/zehira/IMPACTv3-CLIN-20140*/AllSampleResults/*.seg").split("\n")
    segmentation_data = open("mixed_dmp_MSK-IMPACT_2014_data_cna_hg19.seg", "w")
    segmentation_data_writer = csv.writer(segmentation_data, delimiter = "\t", quoting = csv.QUOTE_NONE)
    segmentation_data_writer.writerow(["ID", "chrom", "loc.start" , "loc.end" , "num.mark", "seg.mean"])
    # print seg_files
    samples_written_out = []
    for file in seg_files:
        if file:
            #print file + " is being processed..."
            file_name = os.path.basename(file)
            #M13-100-T_bc12_IMPACTv3-CLIN-0001_copynumber.seg
            if "68R"  in file_name:
                match_obj2 = re.match(r'.*IMPACTv3-CLIN-2*0*1*4*0+([1-9]-*[0-9]*R*)-*R*.*', file_name)
            elif "95R" in file_name:
                match_obj2 = re.match(r'.*IMPACTv3-CLIN-2*0*1*4*0+([1-9]-*[0-9]*R*)-*R*.*', file_name)
            elif "132R" in file_name:
                match_obj2 = re.match(r'.*IMPACTv3-CLIN-2*0*1*4*0+([1-9]-*[0-9]*R*)-*R*.*', file_name)
            elif "131R" in file_name:
                match_obj2 = re.match(r'.*IMPACTv3-CLIN-2*0*1*4*0+([1-9]-*[0-9]*R*)-*R*.*', file_name)
            else:
                match_obj2 = re.match(r'.*IMPACTv3-CLIN-2*0*1*4*0+([1-9]-*[0-9]*)-*R*.*', file_name)
            if not match_obj2:
                print file_name
            run_number = match_obj2.group(1)
            if run_number == "6-":
                run_number = "6"
            
            match_obj = re.match(r'(\d+-*[T,N][a-z]*[0-9]*R*r*p*t*-*[A-Z]*)_.*', file_name)
            if match_obj:
                sample_id = match_obj.group(1)
                
                for y in sample_list:
                    if sample_id in y['SampleSheetSampleID'] and y['Run'] == run_number:
                        samples_written_out.append(sample_id)
                        seg_file = parse(file, " ")
                        for x in seg_file:
                            try:
                                cbio_id = cbio_id_lookup[sample_id]
                            except Exception as e:
                                print "In Run: %s  -> %s is not a proper key: %s"%(file_name, sample_id, e)
                                
                            segmentation_data_writer.writerow([cbio_id, x['chrom'], x['loc.start'], x['loc.end'], x['num.mark'], x['seg.mean']])
                #else:
                    #print sample_id  # this prints out samples that failed - most of the time.                   
    
    for x in sample_list:
        if x['SampleSheetSampleID'] not in samples_written_out:
            print "WARNING: Sample :" + str(x['SampleSheetSampleID']) + " in run: " + str(x['Run']) + " is not processed (.seg file). Check! "  
'''
    
if __name__ == '__main__':
    main()

































