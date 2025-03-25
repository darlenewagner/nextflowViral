#!/apps/x86_64/perl/perl-5.16.1-MT/bin/perl
use strict;
use Getopt::Long;
use List::Util qw(max);

my $bedGraph = $ARGV[0];

open(COVERAGE, $bedGraph) || die "Can't find .bedGraph file, $bedGraph $!";

my $first = rindex($bedGraph, '/');
my $second = index($bedGraph, '.');
my $reads = substr($bedGraph, $first + 1, $second - $first - 1);

my $verbose = 0;

## --verbose flag sets $verbose to 1
GetOptions('verbose' => \$verbose);

my $count = 0;
my $sum = 0;
my $name = '';
my $found_name = 0;
my $start = 0;
my $end = 0;
my $prevCoverage = 0;
my $total = 0;
my @F = ();
my @Regions = ();


while(<COVERAGE>)
 {
    @F = split(/\s+/, $_);

    if($F[0] =~ /^\w/)
    {
	
	$name = $F[0];
	#$found_name = 1;
	$total++;
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
	  #$sum = $sum + $coverage;
	  $count++;
      }
    else
      {
	  #$sum = $sum + $coverage;
          $prevCoverage = $coverage;
      }
 }

push @Regions, $count;

if($verbose)
  {
   my $str = sprintf("%0.1f", $sum/$count);
   print $reads, " mapped to ", $name, ", average depth of coverage, ", $str, "x\n";

   print "\nSegments for breadth of coverage (percent):\n";
   
   for( my $r = 1; $r < scalar @Regions; $r = $r + 2)
     {
         my $percentages = sprintf("%0.1f", 100*($Regions[$r] - $Regions[$r - 1])/$total);
         print $percentages, "\t";
     }
  }
else
  {
      my $str = sprintf("%0.1f", $sum/$count);

     my @percentage = ();
      
   for( my $r = 1; $r < scalar @Regions; $r = $r + 2)
     {
         push @percentage, 100*($Regions[$r] - $Regions[$r - 1])/$total;
        # print $percentages, "\t";
     }

      my $maxBreadth = sprintf("%0.1f", max(@percentage));
      print $name, "\tAvg. Depth: ", $str, "\tMax. Breadth: ", $maxBreadth, "%\n";
      
  }
print "\n";
