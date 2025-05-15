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

   warn Dumper($app_def, $raw_params, $params);

   my $token = $app->token();
   my $ws = $app->workspace();

   # Uncomment to print user info.
   my @cmd = ("p3-whoami");
   IPC::Run::run(\@cmd);

   # Create a temp directory for intermediate calculations/results.
   my $cwd = File::Temp->newdir( CLEANUP => 1 );

   # Create a "working" subdirectory for the input and intermediate file(s).
   my $work_dir = "$cwd/work";
   
   # Create a "staging" subdirectory for output files that will be copied to the workspace.
   my $stage_dir = "$cwd/stage";

   -d $work_dir or mkdir $work_dir or die "Cannot mkdir $work_dir: $!";
   -d $stage_dir or mkdir $stage_dir or die "Cannot mkdir $stage_dir: $!";

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
   my @cmd = ("run_treesort.py", "-j", $job_desc, "-s", $stage_dir, "-w", $work_dir);
   my $ok = run(\@cmd);
   if (!$ok)
   {
      die "Command failed: @cmd\n";
   }

   # TODO: This is for testing purposes and can be deleted.
   print `ls -l $stage_dir`."\n";

   # A modifiable version of the workspace's result folder name.
   my $result_folder = $app->result_folder

   # Make sure the result folder name doesn't end with a period.
   if (substr $result_folder, -1 eq ".")
   {
      $result_folder = substr $result_folder, 0, -1;
   }

   # Map file extensions to BV-BRC file types.
   my %suffix_map = (aln => 'aligned_dna_fasta',
                     csv => 'csv',
                     pdf => 'pdf',
                     tre => 'nwk',
                     tsv => 'tsv');

   my @suffix_map = map { ("--map-suffix", "$_=$suffix_map{$_}") } keys %suffix_map;

   if (opendir(my $dh, $stage_dir))
   {
      # Iterate over the files in the staging directory.
      while (my $p = readdir($dh))
      {
         next if $p =~ /^\./;

         # Use the p3 utility to copy the output files to the workspace.
         my @cmd = ("p3-cp", "-r", "-f", @suffix_map, "$stage_dir/$p", "ws:" . $result_folder);
         print "@cmd\n";
         my $ok = IPC::Run::run(\@cmd);
         if (!$ok)
         {
               warn "Error $? copying output with @cmd\n";
         }
      } 
      closedir($dh);
   }
   else
   {
      warn "Staging directory $stage_dir does not exist\n";
   }
}