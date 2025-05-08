
from job_data import JobData
from treesortrunner.common import DEFAULT_REF_SEGMENT, InferenceType, InputSource, InputParameter, INVALID_FASTA_CHARS, \
   Method, safeTrim, ScriptOption, VALID_SEGMENTS
import re
import os
import subprocess
import sys

# TODO: use https://github.com/BV-BRC/bvbrc_subspecies_classification/blob/master/scripts/run_subspecies_classification.py as an example.


# A class that is responsible for 
class TreeSortRunner:

   # The base URL
   base_url: str = None

   # The default name of the input FASTA file.
   default_input_filename: str = "input.fasta"

   # The name of the FASTA file to use as input. This is specified in job_data.
   input_filename: str = None

   # The JSON data for the job.
   job_data: JobData = None
   
   # The directory where the scripts will be run.
   work_directory = None


   # C-tor
   def __init__(self, job_data: JobData, work_directory: str):
         
      # Set and validate the work directory.
      self.work_directory = work_directory
      if not self.work_directory or len(self.work_directory) < 1:
         raise ValueError("The work directory parameter is invalid")
      
      # Determine the base URL.
      if "P3_BASE_URL" in os.environ:
         self.base_url = os.environ["P3_BASE_URL"]
      else:
         self.base_url = "https://www.bv-brc.org"

      # Validate the job data.
      if not job_data:
         raise ValueError("The job data parameter is invalid")
      
      self.job_data = job_data

      if not TreeSortRunner.is_job_data_valid():
         raise ValueError("Job data in the constructor is invalid")
      
      # TODO: Set the output path's absolute path?
      #self.job_data.output_path = os.path.abspath(self.job_data.output_path)
      
      
   # Is the JobData instance valid?
   def is_job_data_valid(self) -> bool:

      try:
         if not self.job_data or not isinstance(self.job_data, JobData):
            raise Exception("job_data is empty")
         
         # Validate the input_source.
         if self.job_data.input_source not in [i.value for i in InputSource]:
            raise ValueError("job_data.input_source is not a valid input source") 

         # If input_source is fasta_data, make sure an input_fasta_data value was provided.
         if self.job_data.input_source == InputSource.FastaData.value and \
         not self.job_data.input_fasta_data:
            raise ValueError("The input FASTA data is invalid")
         
         # If input_source is fasta_file, make sure an input_fasta_file value was provided.
         if self.job_data.input_source == InputSource.FastaFile.value and \
         not self.job_data.input_fasta_file:
            raise ValueError("The input FASTA file is invalid")
         
         # If input_source is fasta_file_id, make sure an input_fasta_file_id value was provided.
         if self.job_data.input_source == InputSource.FastaFileID.value and \
         not self.job_data.input_fasta_file_id:
            raise ValueError("The input FASTA file ID is invalid")

         # Validate the method
         if self.job_data.method not in [m.value for m in Method]:
            raise ValueError("job_data.method is not a valid method") 

         # Validate the output path.
         if not self.job_data.output_path:
            raise ValueError("The output path is invalid")

         # Validate the reference segment and provide a default if not provided.
         refSegment = safeTrim(self.job_data.ref_segment)
         if not refSegment:
            refSegment = DEFAULT_REF_SEGMENT
         elif not refSegment in VALID_SEGMENTS:
            raise ValueError(f"Invalid reference segment: {refSegment}")

         # Validate the segments
         segments = safeTrim(self.job_data.segments)
         if len(segments) > 0:
            for segment in segments.split(","):
               if not segment in VALID_SEGMENTS:
                  raise ValueError(f"Invalid segment: {segment}")

      except Exception as e:
         sys.stderr.write(f"Invalid job data:\n {e}\n")
         return False

      return True
   

   # Build alignments and trees and compile a descriptor file.
   def prepare_dataset(self) -> bool:

      sys.stdout.write("In prepare_dataset\n\n")
      
      try:
         # Example usage: ./prepare_dataset.sh --segments "HA,NA" segments.fasta HA myoutdir
         cmd = ["prepare_dataset.sh"]

         # Should --fast be added?
         if (self.job_data.inference_type and self.job_data.inference_type == InferenceType.FastTree.value):
            cmd.append(ScriptOption.FastTree.value)

         # The segments (optional)
         if self.job_data.segments:
            cmd.append(ScriptOption.Segments.value)
            cmd.append(self.job_data.segments)

         # The input FASTA file.
         cmd.append(self.input_filename)

         # The reference segment
         refSegment = self.job_data.ref_segment
         if refSegment:
            cmd.append(refSegment)

         # The output path
         cmd.append(self.job_data.output_path)

         # TEST
         print(f"{' '.join(cmd)}\n\n")

         result = subprocess.check_call(cmd, shell=False)
         print(f"Result = {str(result)}\n\n")

      except ValueError as e:
         sys.stderr.write(f"Error preparing dataset:\n {e}\n")
         return False
         
      return True
     

   # Determine the source of the FASTA input file and prepare it for use.
   def prepare_input_file(self) -> bool:

      sys.stdout.write("In prepare_input_file\n\n")

      try:
         # The input source determines how the input file is provided.
         input_source = safeTrim(self.job_data.input_source)

         if input_source == InputSource.FastaData.value:
         
            # Copy user data to the input file.
            try:
               with open(self.default_input_filename, "w+") as input_file:
                  input_file.write(self.job_data.input_fasta_data)

            except Exception as e:
               raise IOError(f"Error copying FASTA data to input file:\n {e}\n")

         elif input_source == InputSource.FastaFile.value:

            # The input file will be in the local filesystem.
            self.input_filename = safeTrim(self.job_data.input_fasta_file)
            if len(self.input_filename) == 0:
               raise ValueError("Invalid input filename")
            
         elif input_source == InputSource.FastaFileID.value:

            # Fetch the input file from the workspace.
            try:
               # TODO: Is this syntax correct?
               fetch_fasta_cmd = ["p3-cp", f"ws:{self.job_data.input_fasta_file_id}", self.input_filename]
               subprocess.check_call(fetch_fasta_cmd, shell=False)

            except Exception as e:
               raise IOError("Error copying FASTA file from workspace:\n %s" % str(e))

         else:
            raise ValueError(f"Invalid input source: {input_source}")

         # Set the input filename and include its full path.
         self.input_filename = os.path.abspath(self.input_filename)

         # Validate the input FASTA file.
         if not os.path.exists(self.input_filename) or os.path.getsize(self.input_filename) == 0:
            raise IOError("Input FASTA file is invalid or empty")

         # Remove invalid characters from FASTA headers.
         with open(self.input_filename, 'r+') as f:
            data = ""
            for line in f:
               if line.startswith(">"):
                  line = line.strip()
                  line = re.sub(INVALID_FASTA_CHARS, "", line)
                  line = re.sub(" ", "_", line)
                  line += "\n"

               data += line

            f.seek(0)
            f.write(data)
            f.truncate()

      except Exception as e:
         sys.stderr.write(f"Error processing input file:\n {e}\n")
         return False
      
      return True


   # Run TreeSort on the command line.
   def tree_sort(self) -> bool:

      sys.stdout.write("In tree_sort\n\n")
      
      try:
         cmd = ["treesort"]

         # The clades output path
         clades_path = safeTrim(self.job_data.clades_path)
         if len(clades_path) > 0:
            cmd.append(ScriptOption.CladesPath.value)
            cmd.append(clades_path)

         # The descriptor path.
         descriptor_path = safeTrim(self.job_data.descriptor_path)
         if len(descriptor_path) > 0:
            cmd.append(ScriptOption.DescriptorPath.value)
            cmd.append(descriptor_path)

         # The "match on" options are mutually exclusive.
         if self.job_data.match_on_strain:
            cmd.append(ScriptOption.MatchOnStrain.value)

         elif self.job_data.match_on_epi:
            cmd.append(ScriptOption.MatchOnEPI.value)

         elif self.job_data.match_on_regex:
            cmd.append(ScriptOption.MatchOnRegex.value)
            cmd.append(self.job_data.match_on_regex)

         # No collapse
         if self.job_data.no_collapse:
            cmd.append(ScriptOption.NoCollapse.value)

          # Always add the output path.
         cmd.append(ScriptOption.OutputPath.value)
         cmd.append(self.job_data.output_path)

         # Equal rates
         if self.job_data.equal_rates:
            cmd.append(ScriptOption.EqualRates.value)

         # Is time scaled (timetree)
         if self.job_data.is_time_scaled:
            cmd.append(ScriptOption.IsTimeScaled.value)

         # TODO: uncomment
         # result = subprocess.check_call(cmd, shell=False)
         print(f"{' '.join(cmd)}\n\n")

      except ValueError as e:
         sys.stderr.write(f"Error preparing dataset:\n {e}\n")
         return False
         
      return True