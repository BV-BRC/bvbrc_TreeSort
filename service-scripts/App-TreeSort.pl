#
# The TreeSort application
#

# DMD: How do I run this in my local environment? What needs to be set up?
use Bio::KBase::AppService::AppScript;
use Bio::KBase::AppService::AppConfig;

# TODO: Are all of these needed?
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

   #warn Dumper($app_def, $raw_params, $params);

   my $token = $app->token();
   my $ws = $app->workspace();

   # Uncomment to print user info.
   #my @cmd = ("p3-whoami");
   #IPC::Run::run(\@cmd);

   #
   # Create an output directory under the current dir. The app service is meant to invoke
   # the app script in a working directory; we create a folder here to encapsulate
   # the job output.
   #
   # We also create a staging directory for the input files from the workspace.
   #

   # TODO: may not need a staging directory

   # Create a temp directory for intermediate calculations/results.
   my $cwd = File::Temp->newdir( CLEANUP => 1 );
   my $work_dir = "$cwd/work";
   my $stage_dir = "$cwd/stage";

   -d $work_dir or mkdir $work_dir or die "Cannot mkdir $work_dir: $!";
   -d $stage_dir or mkdir $stage_dir or die "Cannot mkdir $stage_dir: $!";

   # DMD: What is this for? $sstring isn't used!
   my $data_api = Bio::KBase::AppService::AppConfig->data_api_url;
   my $dat = { data_api => $data_api };
   my $sstring = encode_json($dat);

   # Clone the input parameters.
   my $params_to_app = Clone::clone($params);

   #
   # Encode the input parameters as a job description JSON file.
   #
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

   print `ls -l $stage_dir`;

   my %suffix_map = (aln => 'aligned_dna_fasta',
                     csv => 'csv',
                     pdf => 'pdf',
                     tsv => 'tsv');

   my @suffix_map = map { ("--map-suffix", "$_=$suffix_map{$_}") } keys %suffix_map;

   if (opendir(my $dh, $stage_dir))
   {
      while (my $p = readdir($dh))
      {
         next if $p =~ /^\./;

         my @cmd = ("p3-cp", "-r", "-f", @suffix_map, "$stage_dir/$p", "ws:" . $app->result_folder);
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
      warn "Output directory $stage_dir does not exist\n";
   }
}