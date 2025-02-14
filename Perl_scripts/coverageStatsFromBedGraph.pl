#!/apps/x86_64/perl/perl-5.16.1-MT/bin/perl
use strict;

my $bedGraph = $ARGV[0];

open(COVERAGE, $bedGraph) || die "Can't find .bedGraph file, $bedGraph $!";

my $count = 0;
my $sum = 0;
my $name = '';
my $found_name = 0;
my $start = 0;
my $end = 0;
my $prevCoverage = 0;
my @F = ();
my @Regions = ();


while(<COVERAGE>)
 {
    @F = split(/\s+/, $_);

    if(($F[0] =~ /^[A-Z].+/) && ($found_name = 0))
       {
	 $name = $F[0];
	 $found_name = 1;
       }

    my $coverage = $F[2];
    chomp $coverage;
    
    if(($coverage > 0) && ($prevCoverage == 0))
      {
	  push @Regions, $F[1];
	  $prevCoverage = $coverage;
	  $sum = $sum + $coverage;
	  $count++;
      }
    elsif($coverage > 0)
      {
	  $prevCoverage = $coverage;
	  $sum = $sum + $coverage;
	  $count++;
      }
    elsif(($coverage == 0) && ($prevCoverage > 0))
      {
	  push @Regions, $F[1];
	  $prevCoverage = $coverage;
      }
    else
      {
          $prevCoverage = $coverage;
      }
 }

 my $str = sprintf("%0.2f", $sum/$count);
 print "Coverage, ", $str, "x\n";

 print "\nSegments:\n";

 for( my $r = 1; $r < scalar @Regions; $r = $r + 2)
   {
      print $Regions[$r] - $Regions[$r - 1], "\t";
   }


print "\n\n";
