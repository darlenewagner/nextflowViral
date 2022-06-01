#!/apps/x86_64/perl/perl-5.16.1-MT/bin/perl
## a wrapper script for output of bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%AD]\t[%DP]\n' sample.vcf.gz
use strict;

my $count = 0;
my $sum = 0;

print "Ref_ID\tPosition\tRef\tAlt\tQual\tCoverage\tFraction\n"; 

while(<STDIN>)
{
    my @F = split(/\s+/, $_);
    my @half = split(/,/, $F[5]);
    printf "%\s\t%\d\t%\s\t%\s\t%\d\t%\d\t%0.4f\n", $F[0], $F[1], $F[2], $F[3], $F[4], $F[6], $half[1]/$F[6];
}


