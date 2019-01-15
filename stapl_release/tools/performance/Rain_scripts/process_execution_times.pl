#!/usr/bin/perl
=head1 SYNOPSIS

  process_execution_times.pl
    This script processes .results files generated by the -generate flag of the 
    frontend script test_driver.pl. It scans the local directory for all
    .results files, extracts the timinng data, and stores it into the database.
    
  Usage:
    process_execution_times.pl [options]

  Options:
    -debug
      Shows debug info

    -magic
      Doesn't do anything, but it's still magic.

    -num_cores
      Tells us how many cores the resultswere run on.

    -help
      Shows this help information.
=cut

use Pod::Usage;
use Getopt::Long;
use DBI;
use Net::Domain qw (hostname hostfqdn hostdomain);

use warnings;
use strict;

sub trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };

# Global Variables
my $opt_help;
my $opt_magic;
my $opt_num_cores;
my $opt_debug;
my $opt_data_size;
my @result_files; # list of all .results files.
my $function_path; # Path to our script utilities.
my @timings;

# Get command line options.
my $options_result =
GetOptions ("magic"       => \$opt_magic,
            "num_cores=s" => \$opt_num_cores,
            "data_size=s"  => \$opt_data_size,
            "debug"        => \$opt_debug,
            "help"        => \$opt_help);

die "Invalid option specification" if (!$options_result);

# Display help is -help option was provided.
if ($opt_help)
{ pod2usage(); exit; }

if (not defined ($opt_num_cores))
{
  die "You need to define the number of cores with -num_cores=# \nSee -help for more help.\n";
}

# Find all .results files in local directory.
for my $file (glob("*.results"))
{
  push(@result_files, $file) if (-f $file);
  if ($opt_debug) {print "Push count: @result_files\n";}
}

# Database Setup
my $dsn = 'DBI:mysql:stapl_new_perf:pdbsrv.cs.tamu.edu';
my $db_user_name = 'nthomas';
my $db_password = 'howdy';
my $dbh =  DBI->connect($dsn, $db_user_name, $db_password)
  || die "connect did not work: $DBI::errstr";

# Load .results timings into hash table
foreach my $file (@result_files)
{
  open(TIMES, '<', $file);
  if ($opt_debug) { print "$file opened.\n"; }
  
  # Set up the data structure for this file in the timings list.
  my @tests = ();
  my $file_record = {'file' => $file, 'tests' => \@tests};
  push(@timings, $file_record);
  
  # Process each line in the file and put the timings in the hash table.
  my $test_held = 0; # 0 is empty, 1 is building, 2 is completed/commitable
  my $timing_held = 0; # 0 is empty, 1 is building, 2 is completed/commitable  
  while (my $line = <TIMES>)
  {
    # A line will either begin with "Test", "Status", "Version", "Time", or be junk data.
    chomp $line;

    #Valid lines are of the form "Keyword : Value"
    my $value_offset = 2; #3 offset is 2 for "test : value" lines and 1 for "test: value" lines.
    my @line_entries;
    
    if(index($line, " : ") != -1)
    {
      @line_entries = split(" ", $line);
    }
    else 
    {
      @line_entries = split(":", $line);
      $value_offset = 1;
    }
    my $entry_count = $#line_entries + 1;
    
    if( (($entry_count > 2) && ($line_entries[1] eq ":")) || ($entry_count == 2) ) # If not, this line is probably junk
    {
      # Process a Test line
      if(($entry_count > $value_offset) && ($line_entries[0] eq "Test")) 
      {
        push(@{$timings[-1]{tests}}, {}); # Pushes a reference to an empty hash into the tests array
        
        # Populate test structure. 
        ${$timings[-1]->{tests}[-1]}{'test'} = trim($line_entries[$value_offset]);
        ${$timings[-1]->{tests}[-1]}{'result'} = "Undefined";
        ${$timings[-1]->{tests}[-1]}{'timings'} = [];
        $test_held = 1;
      }

      # Process a Status line
      if(($test_held == 1) && ($entry_count > $value_offset) && ($line_entries[0] eq "Status")) 
      {
        ${$timings[-1]->{tests}[-1]}{'result'} = trim($line_entries[$value_offset]);
      }
      
      # Process a Version line
      if(($test_held > 0) && ($entry_count > $value_offset) && ($line_entries[0] eq "Version")) 
      {      
        push (@{${$timings[-1]->{tests}[-1]}{'timings'}}, {});
        ${${$timings[-1]->{tests}[-1]}{'timings'}[-1]}{'version'} =trim($line_entries[$value_offset]);
        ${${$timings[-1]->{tests}[-1]}{'timings'}[-1]}{'time'}="Pending";
        $timing_held = 1;
      }

      # Process a Time line
      if(($timing_held == 1) && ($entry_count > $value_offset) && ($line_entries[0] eq "Time")) 
      {      
        ${${$timings[-1]->{tests}[-1]}{'timings'}[-1]}{'time'}= trim($line_entries[$value_offset]);
        $test_held = 2;
        $timing_held = 0;
      }
    }
  }
  close(TIMES);
}

if($opt_debug) 
{
  foreach my $this_file (@timings) 
  {
    print "\n\nFILE: $this_file->{file}\n";
    foreach my $this_test ( @{$this_file->{tests}} ) 
    {
      print "  $this_test->{test}\n";
      foreach my $this_timing (@{$this_test->{timings}})
      {
        print "    $this_timing->{version}";
        print "  $this_timing->{time} \n";
      }
    }
  }
}

# Get machine id. 
my $machine_name = hostfqdn();
my $machine_id = -1;

my $query = $dbh->prepare(qq{ SELECT machine_id, name FROM machines });
$query->execute();
while (my ($id, $name) = $query->fetchrow_array()) 
{
  if ($name eq $machine_name) {$machine_id = $id};
}
$query->finish();

#Add to machines table if not found
if ($machine_id == -1)
{
  $query = $dbh->prepare("INSERT INTO machines (`name`, `description`) VALUES (?, ?)");
  $query->bind_param( 1, $machine_name );
  $query->bind_param( 2, $machine_name );
  $query->execute();
  $machine_id = $query->{mysql_insertid};
  $query->finish();
}

# Get User ID
my $username = getpwuid( $< );
my $user_id = -1;

$query = $dbh->prepare(qq{ SELECT user_id, csuser FROM users });
$query->execute();
while (my ($id, $name) = $query->fetchrow_array()) 
{
  if ($name eq $username) {$user_id = $id};
}
$query->finish();

#Add to users table if not found
if ($user_id == -1)
{
  $query = $dbh->prepare("INSERT INTO users (`csuser`, `name`) VALUES (?, ?)");
  $query->bind_param( 1, $username );
  $query->bind_param( 2, $username );
  $query->execute();
  $user_id = $query->{mysql_insertid};
  $query->finish();
}

#
# Insert times into the database
my $run_id = -1;
my $test_id = -1;
my $num_results = 0;
foreach my $this_file (@timings) 
{
  #find application_id
  my $app_name = $this_file->{file};
  my $app_id = -1;

  $query = $dbh->prepare(qq{ SELECT app_id, name FROM applications });
  $query->execute();
  while (my ($id, $name) = $query->fetchrow_array()) 
  {
    if ($name eq $app_name) {$app_id = $id};
  }
  $query->finish();

  #Add to applications table if not found
  if ($app_id == -1)
  {
    $query = $dbh->prepare("INSERT INTO applications (`name`, `description`) VALUES (?, ?)");
    $query->bind_param( 1, $app_name );
    $query->bind_param( 2, $app_name );
    $query->execute();
    $app_id = $query->{mysql_insertid};
    $query->finish();
  } 
  if ($opt_debug) { print "Inserting into runs (machine_id, app_id, user_id, opt_num_cores) VALUES $machine_id, $app_id, $user_id, $opt_num_cores \n"; }
  #insert into runs table
  $query = $dbh->prepare("INSERT INTO runs (`machine_id`, `application_id`, `user_id`, `description`, `processor_count`, `data_size`) VALUES (?, ?, ?, 'Nightly Rel_alpha Timing', ?, ?)");
  $query->bind_param( 1, $machine_id );
  $query->bind_param( 2, $app_id );
  $query->bind_param( 3, $user_id );
  $query->bind_param( 4, $opt_num_cores );
  $query->bind_param( 5, $opt_data_size ); 
  $query->execute();
  $run_id = $query->{mysql_insertid};
  $query->finish();

  foreach my $this_test ( @{$this_file->{tests}} ) 
  {
    foreach my $this_timing (@{$this_test->{timings}})
    {
      if ($opt_debug) { print "Inserting into tests (run_id, this_test->{test}, this_timing->{version}) VALUES $run_id, $this_test->{test}, $this_timing->{version} \n"; }
      #insert into tests table
      $query = $dbh->prepare("INSERT INTO tests (`run_id`, `name`, `version`) VALUES (?, ?, ?)");
      $query->bind_param( 1, $run_id );
      $query->bind_param( 2, $this_test->{test} );
      $query->bind_param( 3, $this_timing->{version} );
      $query->execute();
      $test_id = $query->{mysql_insertid};
      $query->finish();
      
      if ($opt_debug) { print "Inserting into runs (test_id, this_timing->{time}) VALUES $test_id, $this_timing->{time} \n"; }
      #insert into test_results table
      $query = $dbh->prepare("INSERT INTO test_results (`test_id`, `stat_id`, `value`) VALUES (?, ?, ?)");
      $query->bind_param( 1, $test_id );
      $query->bind_param( 2, 5 );
      $query->bind_param( 3, $this_timing->{time} );
      $query->execute();
      $query->finish(); 
      
      $num_results++;  
    }
  }
}

if ($num_results == 0)
{
  print "\n\nFAILURE TO INSERT DATA RESULTS!! VERIFY SUCCESSFUL EXDECUTION OF TESTS!!!\n\n";
}
else
{
  print " $num_results runs inserted in results table. ";
}

#Close down DB connection gracefully.
$dbh->disconnect();
