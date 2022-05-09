#!/apps/x86_64/perl/perl-5.16.1-MT/bin/perl
use strict;
my $sam = $ARGV[0];

open(SAM, $sam) || die "Can't find .sam file, $sam $!";

my $count = 0;
my $sum = 0;

while(<SAM>)
{
    my @F = split(/\s+/, $_);
    if($F[8] > 0)
    {
	$sum = $sum + $F[8];
	$count = $count + 1;
    }
}

printf("%.2f\n", $sum/$count)
