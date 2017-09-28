# ddpcr\_flanking\_sequence
Scripts that help in getting flanking sequence for creating primer for ddPCR assay

## get\_flanking\_sequence.py

This module gets the information about 200bp flanking sequence given chromosome, start, end, reference allele and alternate allele based on **maf**

Source: http://github.com/mskcc/ddpcr_flanking_sequence

License: [Apache License 2.0](http://www.apache.org/licenses/LICENSE-2.0)

### Usage

```

python get_flanking_sequence.py --help
usage: get_flanking_sequence.py [options]

This module gets the information about 200bp flanking sequence given
chromosome, start, end, reference allele and alternate allele based on maf
specification

optional arguments:
  -h, --help            show this help message and exit
  -chr 1, --chromosome 1
                        Enter the information for what chromosome the event is
                        on
  -of OutFilePrefix, --outputfileprefix OutFilePrefix
                        Output File Prefix for the flanking fasta file.
  -v, --verbose         make lots of noise [default]
  -s 5, --start 5       Start coordinate for the event
  -e 20, --end 20       End coordinate for the event
  -f 200, --flank 200   Number of bases to flank for the output [default=200]
  -r /somepath/Homo_Sapeins_hg19.fasta, --referenceFile /somepath/Homo_Sapeins_hg19.fasta
                        Full Path to the reference file with the fasta index.
  -ref someref, --reference_allele someref
                        reference allele for the event
  -alt somealt, --alternate_allele somealt
                        alternate allele for the event
  -p_ann p.ann, --protein_annotation p.ann
                        Protein annotation as string
  -c_ann c.ann, --cDNA_annotations c.ann
                        cDNA annotation as string
  -gene GENE, --gene_name GENE
                        Gene name as string
  -o /somepath/output, --outDir /somepath/output
                        Full Path to the output dir.
  -vcf /somepath/pop_vcf, --pop_vcf /somepath/pop_vcf
                        Full Path to gzipped the population frequency vcf

```

### Inputs 
#### Based on Mutation Annotation Format ([MAF](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification))

chr : chromosome location of the alteration

start: start position of the alteration

end: end position of the alteration

reference_allele: reference allele information

variant_allele: variant allele information

protein_ann: variant annotation for protein change

cDNA_ann: variant annotation for cDNA change

gene: Gene name where the variant occurs

#### Poputation Frequecny VCF for masking

vcf: location of gzipped population frequency vcf

### Output
a fasta file with following header

```
>protein_ann-cDNA_ann

```

### Example Usage:

```
python get_flanking_sequence.py -chr 11 -s 64577305 -e 64577307 -ref AGC -alt GCTT -p_ann R92Qfs*25 -c_ann 120_125delins -r /Users/shahr2/Documents/PubData/Genome/hg19/Homo_sapiens_assembly19.fasta -vcf /Users/shahr2/Documents/PubData/Genome/hg19/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz -o /Users/shahr2/Desktop/ -of test_gfs -gene IDH2
```

### Requirements:
pysam : [v0.10.0](https://pysam.readthedocs.io/en/latest/)
pyvcf : [v0.6.8](http://pyvcf.readthedocs.io/en/latest/INTRO.html)
