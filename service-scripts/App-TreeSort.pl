#
# The TreeSort application
#

use Bio::KBase::AppService::AppScript;
use Bio::KBase::AppService::AppConfig;

use strict;
use Data::Dumper;
use File::Basename;
use File::Slurp;
use File::Temp;
use LWP::UserAgent;
use JSON::XS;
use IPC::Run qw(run);
use Cwd;
use Clone;

my $script = Bio::KBase::AppService::AppScript->new(\&process_treesort, \&preflight);

my $rc = $script->run(\@ARGV);

exit $rc;


sub preflight
{
   # Declare and assign local variables for the parameters passed to the preflight function.
   my($app, $app_def, $raw_params, $params) = @_;

   print STDERR "Pre-flight TreeSort ", Dumper($params, $app);

   return {
      cpu => 2,
      memory => "64G",
      runtime => 18000,
      storage => 0,
   };
}

sub process_treesort
{
   # Declare and assign local variables for the parameters passed to the process_treesort function.
   my($app, $app_def, $raw_params, $params) = @_;

   # Uncomment to troubleshoot the parameters.
   # warn Dumper($app_def, $raw_params, $params);

   my $token = $app->token();
   my $ws = $app->workspace();

   # Uncomment to print user info.
   #my @cmd = ("p3-whoami");
   #IPC::Run::run(\@cmd);

   # Create a temp directory for intermediate calculations/results.
   # DMD TEST: Don't clean up the temp directory when we're testing.
   my $cwd = File::Temp->newdir(CLEANUP => 0); # CLEANUP => 1 );

   # Create an "input" subdirectory for the input FASTA file, etc.
   my $input_dir = "$cwd/input";
   -d $input_dir or mkdir $input_dir or die "Cannot mkdir $input_dir: $!";

   # Create a "staging" subdirectory for output files that will be copied to the workspace.
   my $stage_dir = "$cwd/stage";
   -d $stage_dir or mkdir $stage_dir or die "Cannot mkdir $stage_dir: $!";

   # Create a "working" subdirectory for the input and intermediate file(s).
   my $work_dir = "$cwd/work";
   -d $work_dir or mkdir $work_dir or die "Cannot mkdir $work_dir: $!";
   
   # TODO: Are these needed?
   my $data_api = Bio::KBase::AppService::AppConfig->data_api_url;
   my $dat = { data_api => $data_api };
   my $sstring = encode_json($dat);

   # Clone the input parameters.
   my $params_to_app = Clone::clone($params);

   # Encode the input parameters as a job description JSON file.
   my $job_desc = "$cwd/jobdesc.json";
   open(JOB_DESC, ">", $job_desc) or die "Cannot write $job_desc: $!";
   print JOB_DESC JSON::XS->new->pretty(1)->encode($params_to_app);
   close(JOB_DESC);

   my $parallel = $ENV{P3_ALLOCATED_CPU};

   # Run the Python script that runs TreeSort.
   my @cmd = ("run_treesort.py", "-i", $input_dir, "-j", $job_desc, "-s", $stage_dir, "-w", $work_dir);
   my $ok = run(\@cmd);

   # TODO: This is for testing purposes and can be deleted.
   print "Contents of the staging directory after running TreeSort:\n";
   print `ls -l $stage_dir`."\n";

   # Was the command successful?
   if (!$ok)
   {
      die "Command failed: @cmd\n";
      exit 1;
   }

   # A modifiable version of the workspace's result folder name.
   my $result_folder = $app->result_folder;

   # Make sure the result folder name doesn't end with "/.".
   if ($result_folder =~ m{/\.$})
   {
      $result_folder = substr $result_folder, 0, -2;
   }

   print "Result folder is now $result_folder";

   # Map file extensions to BV-BRC file types.
   my %suffix_map = (aln => 'aligned_dna_fasta',
                     csv => 'csv',
                     pdf => 'pdf',
                     tre => 'nwk',
                     tsv => 'tsv');

   my @suffix_map = map { ("--map-suffix", "$_=$suffix_map{$_}") } keys %suffix_map;


   # Make sure the staging directory exists.
   if (! -d $stage_dir) {
      die "Staging directory $stage_dir does not exist\n";
      exit 1;
   }

   # If the result directory doesn't exist, create it.
   if (! -d "ws:$result_folder") {

      # Make sure the result folder (output path) exists.
      my @cmd = ("p3-mkdir", "ws:$result_folder");
      print "@cmd\n";

      my $ok = IPC::Run::run(\@cmd);
      if (!$ok)
      {
         die "Error $? creating directory ws:$result_folder\n";
         exit 1;
      }
   }

   # Use the p3 utility to copy the staged files to the user's workspace.
   my @cmd = ("p3-cp", "-r", "-f", @suffix_map, "$stage_dir/", "ws:$result_folder");
   print "@cmd\n";
   my $ok = IPC::Run::run(\@cmd);
   if (!$ok)
   {
      warn "Error $? copying output with @cmd\n";
      exit 1;
   }
   
   exit 0;
}