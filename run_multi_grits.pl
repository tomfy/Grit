#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $nprocs = 1;
my $input_data_file = undef; # file with the marker orders or  2 genomes
my $control_file = undef; # control file
my $nchains = undef; # number chains per process
my $seed = time();
my $signed = undef;

GetOptions(
	   'input_file|data_file=s' => \$data_file,
	   'control_file=s' => \$control_file,
	   'nprocs|processes=i' => \$nprocs,
	   'nchains=i' => \$nchains,  # number of chains per process
	   'seed|rng_seed=i' => \$seed,
	   'signed!' => $signed, # -nosigned -> unsigned analysis
	  );
print STDERR "nprocs: $nprocs\n";
my $pid = -1;
my $proc_index = 0;
my $prefix;
my $short_long = '--long';

while (1) {
  if ($pid == 0) { # child process
    $prefix = $proc_index . "_"; # $proc_index 1..$nproc-1
    $short_long = ($proc_index % 2  ==  0)? '--short' : '--long';
    print STDERR "Child process $proc_index. Running grit with prefix:  $prefix, seed: $seed\n";
    my $grit_command = "~/Grit/grit --prefix $prefix  --rng $seed  $short_long ";

    if(defined $data_file){ $grit_command .= " --data $data_file "; }	
    if(defined $control_file){ $grit_command .= " --control $control_file "; }
    if(defined $signed){ $grit_command .= " $signed "; }
    
    system"$grit_command";
    print STDERR "Child process $proc_index. grit completed.\n";
    last;
  } elsif (defined $pid) { # parent process
    # $proc_index++;
    $proc_index++;
    $seed += 10000;
    #sleep($sleeptime);
    print STDERR "parent process with pid: $pid  proc_index:  $proc_index  seed: $seed\n";
    if($proc_index < $nprocs){
      print STDERR "parent process with pid: $pid  proc_index:  $proc_index  forking.\n";
      $pid = fork();
    }else{ # the last one, just run in parent process
      # $proc_index++;
      $prefix = $proc_index . "_"; #$proc_index == $nproc
      #sleep($sleeptime);
      $short_long = ($proc_index % 2  ==  0)? '--short' : '--long';
      print STDERR "running grit from parent process, with prefix: $prefix  seed: $seed\n";
      my $grit_command = "~/Grit/grit --prefix $prefix  --rng $seed  $short_long";
      if(defined $data_file){ $grit_command .= " --data $data_file "; }	
      if(defined $control_file){ $grit_command .= " --control $control_file "; }
      if(defined $signed){ $grit_command .= " $signed "; }
      system "$grit_command";
      print STDERR "parent process $proc_index grit completed.\n";
      last;
    }
  } else { 
    die "fork failed.\n";
    last;
  }
 
}
