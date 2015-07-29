#!/usr/bin/perl -w
# run_doc_on_samples.pl <bam.fof> <interval file> <working dir>---
# Author: Zehir <zehira@phoenix-h2.mskcc.org>
# Modified By Ronak H Shah
# Created: 13 May 2015
# Last Modified : 24 June 2015
# Version: 0.02
# Given a file of bam files, interval file and the working dir, this will run Depth of Coverage on the bam file using the interval file in the working dir
use warnings;
use strict;
use File::Basename;

die "Usage run_doc_on_samples.pl <bam.fof> <interval file> <working dir>\n" unless @ARGV == 3;
my $bams     = $ARGV[0];
my $interval = $ARGV[1];
my $cwd      = $ARGV[2];

open( IN, "<", $bams ) || die "Cannot open $bams, $!\n";
while (<IN>)
{
	chomp;
	my $bam = $_;
	if($bam =~ /\//){
		$bam = basename ($bam);
	}
	my ($sample) = $bam =~ /(.*)\.bam/;
	my $output = $sample . "_DOC_mapq";
	next if ( $bam =~ /Pool/ );
	next if ( $bam =~ /NTC/ );

	#if (!defined($sample)) {
	print "$bam\n";

	#}
	my $cmd =
"/common/sge/bin/lx-amd64/qsub -q test.q -V -wd $cwd -N DOC.${sample} -o ${sample}.DOC.stdout -e ${sample}.DOC.stderr -l h_vmem=8g,virtual_free=8g -pe smp 2 -b y '/dmp/resources/prod/tools/system/java/java-1.7.0-openjdk-1.7.0.9.x86_64/bin/java -Xmx4g -jar /dmp/resources/prod/tools/bio/gatk/production/GenomeAnalysisTK.jar -T DepthOfCoverage -L $interval -I $bam -R /dmp/data/pubdata/hg-fasta/production/Homo_sapiens_assembly19.fasta -o $output -rf BadCigar -mmq 20 -mbq 20'";
	print $cmd,"\n";
	`$cmd`;
}
close(IN)
__END__

=head1 NAME

run_doc_on_samples.pl <bam.fof> <interval file> <working dir>

=head1 SYNOPSIS

run_doc_on_samples.pl [options] args

      -opt --long      Option description

=head1 DESCRIPTION

Stub documentation for run_doc_on_samples.pl, 

=head1 AUTHOR

Zehir, E<lt>zehira@phoenix-h2.mskcc.orgE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2015 by Zehir

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
