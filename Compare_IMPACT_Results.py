__author__ = 'zehira'

import csv, os, argparse
from collections import defaultdict

parser = argparse.ArgumentParser(prog = "Compare_IMPACT_Results.py", description="Compare results of two different IMPACT runs and create tables showing differences", usage = "%(prog)s [options]")
parser.add_argument('-f1', '--IMPACT-result-1', dest='IMPACT_1', help = "IMPACT results file #1 you want to compare")
parser.add_argument('-f2', '--IMPACT-result-2', dest='IMPACT_2', help = "IMPACT results file #2 you want to compare")
parser.add_argument('-n1', '--name-1', dest='NAME1', help = "Name for the first file")
parser.add_argument('-n2', '--name-2', dest='NAME2', help = "Name for the second file")
parser.add_argument('-p', '--pool', dest='POOL', help = "Name for the pool being compared")
parser.add_argument('-t', '--type', dest='TYPE', help = "The type of variants you are comparing. Eg: exonic, silent, non-panel, all...")
parser.add_argument('-k', '--kind', dest='KIND', help = "The kind of results you are comparing. Acceptable values: Mutation, SCNA, SV", choices = ['Mutation', 'SCNA', 'SV'])
parser.add_argument('-c', '--coverage', dest='COV', help = "A table showing the per sample mean coverage across canonical exons")

args = parser.parse_args()


def get_validation_exons():
    result = defaultdict()
    for line in open("/dmp/data/mskdata/two-tier-filters/production/validation_exons.txt", "r"):
        gene, exon = line.strip("\n").split("\t")
        key = "%s:%s" % (gene, exon)
        result[key] = 1
    return result



def parse(filename, delimiter = "\t"):
    file = open(filename, 'rU')
    tabbed_data = csv.reader(file, delimiter = delimiter)
    fields = tabbed_data.next()
    parsed_data = []
    for row in tabbed_data:
        parsed_data.append(dict(zip(fields, row)))
    file.close()
    return parsed_data


def make_data_dict(file):
    data = {}
    for line in file:

        key = "%s:%s:%s:%s:%s" % (line["Sample"], line["Chrom"], line["Start"], line["Ref"], line["Alt"])
        data[key] = line
    return data


def make_scna_data_dict(file):
    data = {}
    for line in file:
        if line['sig'] == "1":
            if line['fc'] > 0:
                key = "%s:%s:Amplification" % (line['sample'], line['region'])
                data[key] = line
            if line['fc'] < 0:
                key = "%s:%s:Deletion" % (line['sample'], line['region'])
                data[key] = line
    return data

def get_header(file):
    with open(file, 'r') as f:
        header = f.readline()
    return header


def get_sample_coverage():
    results_cv3 = defaultdict()
    results_cv5 = defaultdict()
    for line in open(args.COV, 'r'):
        row = line.strip("\n").split("\t")
        if "SampleName" in row:
            continue
        results_cv3[row[3]] = "%.0f" % float(row[1])
        results_cv5[row[3]] = "%.0f" % float(row[2])


    return results_cv3, results_cv5


def parse_scna_files(args):
    data1 = make_scna_data_dict(parse(args.IMPACT_1))
    data2 = make_scna_data_dict(parse(args.IMPACT_2))
    result = defaultdict()
    for key1 in data1:
        sample, gene, tmp = key1.split(":")
        if key1 in data2:
            if sample in result:
                result[sample]["common"] += 1
            else:
                result[sample] = defaultdict()
                result[sample]["common"] = 0
                result[sample]["unique" + str(args.NAME1)] = 0
                result[sample]["unique" + str(args.NAME2)] = 0
                result[sample]["common"] += 1
        else:
            result[sample] = defaultdict()
            result[sample]["common"] = 0
            result[sample]["unique" + str(args.NAME1)] = 0
            result[sample]["unique" + str(args.NAME2)] = 0
            result[sample]["unique" + str(args.NAME1)] += 1

    for key2 in data2:
        sample, gene, tmp = key2.split(":")
        if key2 not in data1:
            sample, gene, tmp = key1.split(":")
            if sample in result:
                result[sample]["unique" + str(args.NAME2)] += 1
            else:
                result[sample]["unique" + str(args.NAME2)] = 0
                result[sample]["common"] = 0
                result[sample]["unique" + str(args.NAME1)] = 0
                result[sample]["unique" + str(args.NAME2)] += 1
    for sample in result:
        print sample, result[sample]["common"], result[sample]["unique" + str(args.NAME1)], result[sample]["unique" + str(args.NAME2)]


def parse_files(args):
    #coverage_cv3, coverage_cv5 = get_sample_coverage()
    data1 = make_data_dict(parse(args.IMPACT_1))
    data2 = make_data_dict(parse(args.IMPACT_2))
    common_var_file1 = []
    common_var_file2 = []
    unique_to_file1 = []
    unique_to_file2 = []
    total_call_file1 = 0
    total_call_file2 = 0
    result = defaultdict()
    for key1 in data1:
        total_call_file1 += 1
        if key1 in data2:
            sample, chr, start, ref, alt = key1.split(":")
            if sample in result:
                result[sample]["common"] += 1
            else:
                result[sample] = defaultdict()
                result[sample]["common"] = 0
                result[sample]["unique" + str(args.NAME1)] = 0
                result[sample]["unique" + str(args.NAME2)] = 0
                result[sample]["common"] += 1

            common_var_file1.append(data1[key1])
            common_var_file2.append(data2[key1])
        else:
            sample, chr, start, ref, alt = key1.split(":")
            if sample in result:
                result[sample]["unique" + str(args.NAME1)] = 0
                result[sample]["unique" + str(args.NAME1)] += 1


            unique_to_file1.append(data1[key1])
    for key2 in data2:
        total_call_file2 += 1
        if key2 not in data1:
            sample, chr, start, ref, alt = key2.split(":")
            if sample in result:

                result[sample]["unique" + str(args.NAME2)] += 1
            else:
                result[sample] = defaultdict()
                result[sample]["common"] = 0
                result[sample]["unique" + str(args.NAME1)] = 0
                result[sample]["unique" + str(args.NAME2)] = 0
                result[sample]["unique" + str(args.NAME2)] += 1

            unique_to_file2.append(data2[key2])

    common_var_count = len(common_var_file1)
    unique_to_file1_count = len(unique_to_file1)
    unique_to_file2_count = len(unique_to_file2)


    print "Overall Results for pool %s:" % args.POOL
    print "-----------------"
    print "File\tTotal Calls\tCommon Calls\tUnique Calls"
    print "%s\t%s\t%s\t%s" % (args.NAME1, total_call_file1, common_var_count, unique_to_file1_count)
    print "%s\t%s\t%s\t%s" % (args.NAME2, total_call_file2, common_var_count, unique_to_file2_count)
    print "\n---------------------\n"
    print "Sample\tCommonMutationCount\t%sMutationCount\t%sMutationCount" % (args.NAME1,args.NAME2 )
    for sample in result:
        print "%s\t%s\t%s\t%s" % (sample, result[sample]["common"], result[sample]["unique" + str(args.NAME1)], result[sample]["unique" + str(args.NAME2)])

    print "____________________"
    '''
    print "Table2.1A"
    for sample in result:
        total_cv3_calls = result[sample]["common"] + result[sample]["unique" + str(args.NAME1)]
        total_cv5_calls = result[sample]["common"] + result[sample]["unique" + str(args.NAME2)]

        print "%s\t%s\t%s\t%s\t%s\t%s" % (sample, sample.replace("T", "N"), coverage_cv3[sample], total_cv3_calls,
                                          coverage_cv5[sample], total_cv5_calls)

    print "____________________"
    '''
    print "Table 2.1B"


    validation_exons = get_validation_exons()
    file_name = args.POOL + "MutationComparison.txt"
    t = open(file_name, "w")
    t.write("Sample\tGene\tExon\tcDNA\tAA\tDP_CV3\tAD_CV3\tVF_CV3\tDP_CV5\tAD_CV5\tVF_CV5\tValidationStatus\n")
    for key1 in data1:
        if key1 in data2:
            try:
                sample, chr, start, ref, alt = key1.split(":")
                cDNAchange = data1[key1]["cDNAchange"]
                aaCahnge = data1[key1]["AAchange"]
                gene = data1[key1]["Gene"]
                exon = data1[key1]['Exon']
                val = "NotValidated"
                key = "%s:%s" % (gene, exon)
                if key in validation_exons:
                    val = "Validated"
                #cv3_coverage = coverage_cv3[sample]
                cv3_DP = data1[key1]["T_TotalDepth"]
                cv3_AD = data1[key1]["T_AltCount"]
                cv3_VF = data1[key1]["T_AltFreq"]
                #cv3_NormCov = float(cv3_DP)/float(cv3_coverage)
                #cv5_coverage = coverage_cv5[sample]
                cv5_DP = data2[key1]["T_TotalDepth"]
                cv5_AD = data2[key1]["T_AltCount"]
                cv5_VF = data2[key1]["T_AltFreq"]
               # cv5_NormCov = float(cv5_DP)/float(cv5_coverage)
                line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (sample, gene, exon, cDNAchange, aaCahnge, cv3_DP, cv3_AD, cv3_VF, cv5_DP, cv5_AD, cv5_VF, val)
                t.write(line)
            except Exception as e:
                print "Encountered an error with this line. Error: %s. Common mutation" % e
                print "Line: %s" % data1[key1]
        else:
            try:
                print key1
                sample, chr, start, ref, alt = key1.split(":")
                cDNAchange = data1[key1]["cDNAchange"]
                aaCahnge = data1[key1]["AAchange"]
                gene = data1[key1]["Gene"]
                exon = data1[key1]['Exon']
                val = "NotValidated"
                key = "%s:%s" % (gene, exon)
                if key in validation_exons:
                    val = "Validated"
                #cv3_coverage = coverage_cv3[sample]
                cv3_DP = data1[key1]["T_TotalDepth"]
                cv3_AD = data1[key1]["T_AltCount"]
                cv3_VF = data1[key1]["T_AltFreq"]
                #cv3_NormCov = float(cv3_DP)/float(cv3_coverage)
                #cv5_coverage = ""
                cv5_DP = ""
                cv5_AD = ""
                cv5_VF = ""
                #cv5_NormCov = ""
                print key
                line =  "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (sample, gene, exon, cDNAchange, aaCahnge, cv3_DP, cv3_AD, cv3_VF, cv5_DP, cv5_AD, cv5_VF, val)
                t.write(line)
            except Exception as e:
                print "Encountered an error with this line. Error: %s. Uniq to CV3" % e
                print "Line: %s" % data1[key1]
    for key2 in data2:
        if key2 not in data1:
            try:
                sample, chr, start, ref, alt = key2.split(":")

                cDNAchange = data2[key2]["cDNAchange"]
                aaCahnge = data2[key2]["AAchange"]
                gene = data2[key2]["Gene"]
                exon = data2[key2]['Exon']
                val = "NotValidated"
                key = "%s:%s" % (gene, exon)
                if key in validation_exons:
                    val = "Validated"
                #cv3_coverage = coverage_cv3[sample]
                cv3_DP = ""
                cv3_AD = ""
                cv3_VF = ""
                #cv3_NormCov = ""
                #cv5_coverage = coverage_cv5[sample]
                cv5_DP = data2[key2]["T_TotalDepth"]
                cv5_AD = data2[key2]["T_AltCount"]
                cv5_VF = data2[key2]["T_AltFreq"]
                #cv5_NormCov = float(cv5_DP)/float(cv5_coverage)
                line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (sample, gene, exon, cDNAchange, aaCahnge, cv3_DP, cv3_AD, cv3_VF, cv5_DP, cv5_AD, cv5_VF, val)
                t.write(line)
            except Exception as e:
                print "Encountered an error with this line. Error: %s. Uniq to CV5" % e
                print "Line: %s" % data2[key2]
    t.close()
    return data1, data2

def write_out_comparisons(data1, data2, args):
    file_name_common_vars = "Common_variants_for_pool_%s_%s_vs_%s.txt" % (args.POOL, args.NAME1, args.NAME2)
    file_name_unique_to_file1 = "Unique_variants_for_pool_%s_%s.txt" % (args.POOL, args.NAME1)
    file_name_unique_to_file2 = "Unique_variants_for_pool_%s_%s.txt" % (args.POOL, args.NAME2)
    header = get_header(args.IMPACT_1)
    with open(file_name_common_vars, "w") as f:
        f.write(header)
        with open(args.IMPACT_1, 'r') as f1:
            for x in f1:
                x = x.strip()
                line = x.split("\t")
                key = "%s:%s:%s:%s:%s" % (line[0], line[2], line[3], line[4], line[5])
                if key in data2:
                    f.write("\t".join(line))
                    f.write("\n")
    f.close()
    with open(file_name_unique_to_file1, "w") as f:
        f.write(header)
        with open(args.IMPACT_1, 'r') as f1:
            for x in f1:
                x = x.strip()
                line = x.split("\t")
                key = "%s:%s:%s:%s:%s" % (line[0], line[2], line[3], line[4], line[5])
                if key not in data2:
                    f.write("\t".join(line))
                    f.write("\n")
    f.close()
    with open(file_name_unique_to_file2, "w") as f:
        f.write(header)
        with open(args.IMPACT_2, 'r') as f1:
            for x in f1:
                x = x.strip()
                line = x.split("\t")
                key = "%s:%s:%s:%s:%s" % (line[0], line[2], line[3], line[4], line[5])
                if key not in data1:
                    f.write("\t".join(line))
                    f.write("\n")
    f.close()


if __name__ == '__main__':
    if args.KIND == "Mutation":
        data1, data2 = parse_files(args)
        write_out_comparisons(data1, data2, args)
    elif args.KIND == "SCNA":
        parse_scna_files(args)
    else:
        parser.print_help()

        # python ../Scripts/GeneralScripts/Compare_IMPACT_Results.py -f1 /dmp/dms/qc-data/IMPACT/2014/IMPACTv3-CLIN-20140177/DEFAULT/IMPACTv3-CLIN-20140177_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedExonic.Filtered.txt -f2 /dmp/hot/zehira/CV5/IMPACTv5-VAL-20140177/IMPACTv5-VAL-20140177_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedExonic.Filtered.txt -n1 CV3 -n2 CV5 -p Clin177 -t exonic -k Mutation -c Clin177_SampleCoverage.txt
        # python ../Scripts/GeneralScripts/Compare_IMPACT_Results.py -f1 /dmp/dms/qc-data/IMPACT/2014/IMPACTv3-CLIN-20140167/20141123_21.42/IMPACTv3-CLIN-20140167_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedExonic.Filtered.txt -f2  /dmp/hot/zehira/CV5/IMPACTv5-VAL-20140167/IMPACTv5-VAL-20140167_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedExonic.Filtered.txt -n1 CV3 -n2 CV5 -p Clin167 -t exonic -k Mutation -c Clin167.SampleCoverage.txt
        # python ../Scripts/GeneralScripts/Compare_IMPACT_Results.py -f1 /dmp/dms/qc-data/IMPACT/2014/IMPACTv3-CLIN-20140163/20141119_21.28/IMPACTv3-CLIN-20140163_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedExonic.Filtered.txt -f2 /dmp/hot/zehira/CV5/IMPACTv5-VAL-20140163/IMPACTv5-VAL-20140163_AllSomaticMutIndel_withAlleleDepth_annovarAnnotatedExonic.Filtered.txt -n1 CV3 -n2 CV5 -p Clin163 -t exonic -k Mutation -c Clin163_SampleCoverage.txt
