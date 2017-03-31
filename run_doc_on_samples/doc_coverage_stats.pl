# Ahmet Zehir
# Modified By Ronak H Shah
# Jun 12, 2013
# given a depth of coverage file, a chromosome number; generates coverage statistics for the probes residing on that chromosome 
# Last Modified on: Jun 25 2015
use strict;
use warnings;
use Cwd;
use File::Basename;
use POSIX qw/ceil/;
#perl doc_coverage_stats.pl 35437182_DOC picard_targets.intervals

die "Usage doc_coverage_stats.pl <doc file> <interval list>\n" unless @ARGV == 2;

my $file = $ARGV[0];
my $intervalList = $ARGV[1];
if($file =~ /\//){
		$file = basename ($file);
}

my (%intervalHash, %docHash);
open IN, "<", $intervalList || die $!;
while(<IN>){
    chomp;
    next if($_=~/^@/);
    my @line = split("\t");
    my $chr = $line[0];
    my $start = $line[1];
    my $end = $line[2];
    my $strand = $line[3];
    my $name = $line[4];
    $intervalHash{$name} = $chr.":".$start.":".$end.":".$strand;
}
close IN;

print "interval file is parsed\n";

# foreach my $name (sort keys %intervalHash){
#     print "$name\t$intervalHash{$name}\n";

# }
# exit;

my $output = $file."_processed.txt";
my $wd = cwd();
my $outfile = $wd ."/".$output;

if (-e $outfile) {
    print STDERR "$output exists and will not be produced again\n";
    exit;
}   

my %coverage;

open IN, "<", "$file";
while(<IN>){
    chomp;
    next if($_=~/^Locus/);
    my @line = split("\t");
    my ($chr, $locus) = split(":", $line[0]);
    my $cov = $line[-1];
    $coverage{$chr}{$locus} = $cov;
}

close IN;

print "Coverge file is parsed\n";

open OUT, ">", $output;
    
foreach my $name (sort keys %intervalHash){
    my ($chr, $start, $end, $strand) = split(":", $intervalHash{$name});
    print STDERR "interval $name , chr $chr is being processed\n";
    #foreach my $locus (sort keys %{$coverage{$chr}}){
    #    if ($locus >= $start && $locus <= $end){
    for my $locus ($start..$end){
        my $cov = 0;
        if (exists($coverage{$chr}{$locus})){
            $cov = $coverage{$chr}{$locus};
        }
        push @{$docHash{$name}},$cov;   
    }
}

print STDERR "All data is parsed. Stats are being calculated\n";

    
print OUT "target\tchr\tinterval\tmin\tq10\tq25\tq50\tq75\tq90\tmax\tmean\tsd\tcv\tzscore_mean\tzscore_sd\tzscore_cv\n";
my @data =();
my @meanData = ();
my @sdData = ();
my @cvData = ();
my @chrData = ();
foreach my $name (sort keys %docHash){
    print STDERR "Stats for interval: $name are being calculated\n";
    my @sorted_doc = sort {$a <=> $b} (@{$docHash{$name}});
    my ($chr, $start, $end, $strand) = split(":", $intervalHash{$name});
    my $length = $#sorted_doc+1;
    
    my $index_q10 = int(0.1*$length);
    my $index_q25 = int(0.25*$length);
    my $index_q50 = int(0.5*$length);
    my $index_q75 = int(0.75*$length);
    my $index_q90 = int(0.90*$length);
    
    my $min = $sorted_doc[0];
    my $q10 = $sorted_doc[$index_q10];
    my $q25 = $sorted_doc[$index_q25];
    my $q50 = $sorted_doc[$index_q50];
    my $q75 = $sorted_doc[$index_q75];
    my $q90 = $sorted_doc[$index_q90];
    my $max = $sorted_doc[-1];
    my ($sum, $sumOfSq) = 0;
    
    foreach my $i (@sorted_doc){
        $sum += $i;
        $sumOfSq += $i*$i;
    }
    my $mean = $sum/$length;
    my $sqtotal =0;
    foreach my $i (@sorted_doc){
        $sqtotal += ($mean-$i)**2;
    }
    my $sd = ($sqtotal/$length)**0.5;
    my $cv;
    my $zscore;
    if ($mean == 0) {
        $cv =0;
    }else{
        $cv = $sd/$mean;
    }
    if ($sd == 0) {
        $zscore = 0;
    }else{
        $zscore= (($sum-$mean)/$sd);
    }
    push(@meanData,ceil($mean));
	push(@sdData,ceil($sd));
	push(@cvData,ceil($cv));
	push(@chrData,$chr);
	push(@data,"$name\t$chr:$start-$end\t$min\t$q10\t$q25\t$q50\t$q75\t$q90\t$max\t$mean\t$sd\t$cv"); 
	
    #print OUT "$name\t$chr:$start-$end\t$min\t$q10\t$q25\t$q50\t$q75\t$q90\t$max\t$mean\t$sd\t$cv\t$zscore\n";
}

my $mean_avg = average(\@meanData);
my $mean_sd = stdev(\@meanData);
my $sd_avg = average(\@sdData);
my $sd_sd = stdev(\@sdData);
my $cv_avg = average(\@cvData);
my $cv_sd = stdev(\@cvData);
my $count = 0;
for my $i (@data){
	my($name,$chrstr,$min,$q10,$q25,$q50,$q75,$q90,$max,$mean,$sd,$cv) = split("\t",$i);
	my $zscore_mean = ($mean - $mean_avg)/$mean_sd;
	my $zscore_sd = 0;
	if ($sd_sd != 0) {
		$zscore_sd = ($sd - $sd_avg)/$sd_sd;
	}
	my $zscore_cv = ($cv - $cv_avg)/$cv_sd;
	print OUT "$name\t$chrData[$count]\t$chrstr\t$min\t$q10\t$q25\t$q50\t$q75\t$q90\t$max\t$mean\t$sd\t$cv\t$zscore_mean\t$zscore_sd\t$zscore_cv\n";
	$count++
}
close(OUT);
#calculate Avg
sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty arrayn");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}
#calculate Std
sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}



