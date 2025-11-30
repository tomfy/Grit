#!/usr/bin/perl -w
use strict;
use List::Util qw (min max sum);
use Getopt::Long;

# given a permutation (specified as e.g.: -target_perm '1,-3,-4,-5,6,7,2' )
my $signed = 0;
my $target_perm_string = '1,-5,-4,-3,-6,2,7'; # comma separated
my $lambda_max = 35;
my $max_length = 8;
my $count_1marker_invs = 0;
my $count_wholechrom_invs = 0;
my $identify_perm_and_reversed_perm = 1;

GetOptions(
           'target_permution=s' => \$target_perm_string,
           'signed=s' => \$signed,
           'lambda_max=f' => \$lambda_max,
           'max_length|max_inversions=i' => \$max_length,
	   'single_marker_inversions=i' => \$count_1marker_invs,
	   'whole_chrom_inversions=i' => \$count_wholechrom_invs,
	   'identify_permrevperm=i' => \$identify_perm_and_reversed_perm,
          );

$target_perm_string =~ s/\s//g; # remove whitespace
my @target_perm = split(",", $target_perm_string);
if (!$signed) {
  @target_perm = map(abs($_), @target_perm);
 
}

print "# Signed? ", ($signed)? 'signed' : 'unsigned', "\n";
print "# Target permutation: $target_perm_string\n";
print "# Max inversions: $max_length\n";
print "# lambda_max: $lambda_max\n";
print "# Count single-marker inversions: ", ($count_1marker_invs)? "Yes" : "No", "\n";
my $min_rev_size  = ($signed  or  $count_1marker_invs)? 1 : 2;
print "# Min inversion length: $min_rev_size markers.\n";
print "# Count whole chromosome inversion: ", ($count_wholechrom_invs)? "Yes" : "No", "\n";
print "# Identify permutation and same permutation reversed? ", ($identify_perm_and_reversed_perm)? "Yes" : "No", "\n";
my $Nmrk = scalar @target_perm; # number of markers
my $Nmrk_factorial = 1;
for my $i (1..$Nmrk) {
  $Nmrk_factorial *= $i;
}

my $Nrevs = ($signed or $count_1marker_invs)?  $Nmrk*($Nmrk+1)/2  : $Nmrk*($Nmrk-1)/2;
if (!$count_wholechrom_invs) {
  $Nrevs--;
}	    # don't count inversion of all markers (i.e. swap of ends)
print "# N inversions: $Nrevs\n";

# set initial perm to 1,2,3,4,5,...,$Nmrk
my @init_perm = (); for my $i (1..$Nmrk) { push @init_perm, $i; }

my $sum = 1.0;
my $sum_arg = 1.0;
my $Linv_factorial = 1.0;
my @relprobs = (0) x $max_length;


my @xarray = @init_perm;
@init_perm = @target_perm;
@target_perm = @xarray;
$target_perm_string = join(",", @target_perm);
print "# init perm: ", join(', ', @init_perm), "\n";
print "# target perm: ", join(', ', @target_perm), "\n";

print STDERR "# 0 0 0 0 0\n";
my $perm_pathcount = {join(",", @init_perm) => 1}; # after 0 inversions this (init perm) is the only perm reached
for my $Linv (1..$max_length) { # get the perms reached after 1,2, ... 5 inversions
  my $newperm_pathcount = {}; # {0,0,0,0,0,0,0}; # keys: permutations (string), values: number of inversion sequences leading to the perm.
  
  while (my($perm, $pathcount) = each %$perm_pathcount) {
    my @parray = split(",", $perm);
    for my $left (0..$Nmrk-1) {
      my $upper_index = $Nmrk;
      if ($left == 0  and !$count_wholechrom_invs) {
	$upper_index = $Nmrk-1;
      } 
      for my $right ($left+$min_rev_size..$upper_index) {
	
	my $new_perm = ($signed)?
	  invert_signed(\@parray, $left, $right) :
	  invert_unsigned(\@parray, $left, $right);

	#print STDERR "# GGG   $Linv  $perm  $new_perm  $pathcount\n" if($new_perm eq $target_perm_string);
	$newperm_pathcount->{$new_perm} += $pathcount;
	# print "$Linv   $new_perm   $target_perm_string\n" if(
	if ($Linv <=4  and  $new_perm eq $target_perm_string) {
	  #print "!!!! $Linv inversions;  \n$perm\n $new_perm    $pathcount\n";
	}
	#print "$Linv   $new_perm  $pathcount  ", $newperm_pathcount->{$new_perm}, "\n";
      }
    }
  } # after loop to add all possible inversions to previous set of paths.
  $Linv_factorial *= $Linv;
  $sum_arg *= $lambda_max/$Linv;
  $sum += $sum_arg;		# $lambda_max**$Linv/$Linv_factorial;
 
  #print STDERR "$Linv  $XXX  ", $lambda_max**$Linv, "  $Linv_factorial ", $lambda_max**$Linv/$Linv_factorial, "\n";  
   
  # while(my($newperm, $pathcount) = each %$newperm_pathcount){ print "$inv   $newperm  $pathcount\n"; }
  my $IoLf = 1.0 - exp(-$lambda_max)*$sum; # the integral divided by L!
  #print STDERR "# AAA: $Linv  $sum_arg   ", $lambda_max**$Linv/$Linv_factorial, "  $sum   $IoLf\n";
  #print STDERR "ZZZ: $ZZZ\n";
  $perm_pathcount = $newperm_pathcount;
  if (exists $perm_pathcount->{$target_perm_string}) {
    
    my $totalNpaths = $Nrevs**$Linv;

    my $K = $Nmrk_factorial * (($signed)? 2**$Linv : 1);
    my $Mapprox = $totalNpaths/$Nmrk_factorial;
    $Mapprox /= 2**$Linv if($signed);
    my $Npaths = $perm_pathcount->{$target_perm_string};
    
    # if($Linv == 3 ){ print "xxx $Npaths \n"; }
    my $relprob = $IoLf*$Npaths/$totalNpaths;
    #$cumeprob += $relprob;
    print STDERR "# $Linv   $Npaths  $totalNpaths   $IoLf  $relprob   $K $Mapprox  ", $Npaths/$Mapprox, "\n";
    $relprobs[$Linv] = $relprob;
  } else {
    print STDERR "# $Linv  0 - - 0 - \n";
    $relprobs[$Linv] = 0;
  }

}
my $cumeprob = sum(@relprobs);
print STDERR join(", ", @relprobs), "\n";
print STDERR "# $cumeprob \n";
print "0  0 0 \n";
#my $K = 
for my $i (1..$max_length) {
  print "$i  $relprobs[$i] ", $relprobs[$i]/$cumeprob, "\n";
}


sub invert_signed{
  my $parray = shift;		# permutation as array ref
  my $l = shift;		# left break-point
  my $r = shift;		# right break-point
  # invert elements $l through $r-1
  my @revparray = @$parray;
  for my $i ($l..$r-1) {
    $revparray[$i] = -1*$parray->[$r-1 - ($i-$l)];
  }
  if ($identify_perm_and_reversed_perm) {
    if (abs($revparray[0]) > abs($revparray[-1])) {
      @revparray = reverse @revparray;
      @revparray =  map(-1*$_, @revparray);
    }
  }
  my $rev_perm = join(",", @revparray);
  # print "$p  $l $r  $rev_perm\n";
  return $rev_perm;
}

sub invert_unsigned{
  my $p = shift;		# permutation as array ref
  my $l = shift;		# left break-point
  my $r = shift;		# right break-point
  # invert elements $l through $r-1
  my @revparray = @$p;
  for my $i ($l..$r-1) {
    $revparray[$i] = $p->[$r-1 - ($i-$l)];
  }
  if ($identify_perm_and_reversed_perm) {

    if (abs($revparray[0]) > abs($revparray[-1])) {
      @revparray = reverse @revparray;
    }
  }
  my $rev_perm = join(",", @revparray);
  # print "$p  $l $r  $rev_perm\n";
  return $rev_perm;
}
