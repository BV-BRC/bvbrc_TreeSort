
from .common import DEFAULT_REF_SEGMENT, InferenceType, InputSource, InputParameter, INVALID_FASTA_CHARS, \
   safeTrim, ScriptOption, VALID_SEGMENTS
import re
import os
import subprocess
import sys
import traceback

# TODO: use https://github.com/BV-BRC/bvbrc_subspecies_classification/blob/master/scripts/run_subspecies_classification.py as an example.


# A class that is responsible for 
class TreeSortRunner:

   # The base URL
   base_url: str = None

   # The input filename with full path.
   input_filename: str = "input.fasta"

   # The JSON data for the job.
   job_data: dict = None

   # The output filename with full path.
   output_filename: str = None
   

   # C-tor
   def __init__(self, job_data, output_filename):
         
      # Set and validate the job data.
      self.job_data = job_data
      if not self.job_data:
         raise ValueError("Invalid job data")
      
      # Set and validate the output filename.
      self.output_filename = safeTrim(output_filename)
      if self.output_filename == "":
         raise ValueError("Invalid output filename")

      # Determine the base URL.
      if "P3_BASE_URL" in os.environ:
         self.base_url = os.environ["P3_BASE_URL"]
      else:
         self.base_url = "https://www.bv-brc.org"
         
      
   @staticmethod
   def is_job_data_valid(job_data: dict) -> bool:

      try:
         if not job_data or not isinstance(job_data, dict) or len(job_data) == 0:
            raise Exception("job_data is empty")
         
         # Validate the reference segment and provide a default if not provided.
         refSegment = safeTrim(job_data[InputParameter.RefSegment.value])
         if not refSegment:
            refSegment = DEFAULT_REF_SEGMENT

         if not refSegment in VALID_SEGMENTS:
            raise ValueError(f"Invalid reference segment: {refSegment}")

         # Validate the segments
         segments = safeTrim(job_data[InputParameter.Segments.value])
         if len(segments) > 0:
            for segment in segments.split(","):
               if not segment in VALID_SEGMENTS:
                  raise ValueError(f"Invalid segment: {segment}")

      except Exception as e:
         sys.stderr.write(f"Job data is invalid:\n {e}\n")
         return False

      return True
   

   # Build alignments and trees and compile a descriptor file.
   def prepare_dataset(self) -> bool:

      try:
         # Example usage: ./prepare_dataset.sh --segments "HA,NA" segments.fasta HA myoutdir
         cmd = ["prepare_dataset.sh"]

         inference_type = self.job_data[InputParameter.InferenceType.value]
         if (inference_type and inference_type == InferenceType.FastTree.value):
            cmd.append(ScriptOption.FastTree.value)

         segments = self.job_data[InputParameter.Segments.value]
         if segments:
            # TODO: validate segments
            cmd.append(ScriptOption.Segments.value)
            cmd.append(segments)

         cmd.append(self.input_filename)

         refSegment = self.job_data[InputParameter.RefSegment.value]
         if refSegment:
            # TODO: validate refSegment
            cmd.append(refSegment)

         cmd.append(self.output_filename)

         # TODO: uncomment
         # result = subprocess.check_call(cmd, shell=False)
         print(cmd)

      except ValueError as e:
         sys.stderr.write(f"Error preparing dataset:\n {e}\n")
         return False
         
      return True
     

   # Determine the source of the FASTA input file and prepare it for use.
   def prepare_input_file(self) -> bool:

      try:
         # The input source determines how the input file is provided.
         input_source = self.job_data[InputParameter.InputSource.value]

         if input_source == InputSource.FastaData.value:
         
            # Copy user data to the input file.
            try:
               with open(self.input_filename, "w+") as input_file:
                  input_file.write(self.job_data[InputParameter.InputFastaData.value])

            except Exception as e:
               raise IOError("Error copying FASTA data to input file:\n {e}\n")

         elif input_source == InputSource.FastaFile.value:

            # The input file will be in the local filesystem.
            try:
               self.input_filename = safeTrim(self.job_data[InputParameter.InputFastaFile])
               if len(self.input_filename) == 0:
                  raise IOError(f"Invalid input file")

            except Exception as e:
               raise RuntimeError("Error copying FASTA file from workspace:\n {e}\n")
            
         elif input_source == InputSource.FastaFileID.value:

            # Fetch the input file from the workspace.
            try:
               fetch_fasta_cmd = ["p3-cp", "ws:%s" %(self.job_data[InputParameter.InputFastaFileID.value]), self.input_filename]
               subprocess.check_call(fetch_fasta_cmd, shell=False)

            except Exception as e:
               raise RuntimeError("Error copying FASTA file from workspace:\n %s" % str(e))

         else:
            raise ValueError(f"Invalid input source: {input_source}")

         # TEST
         sys.stdout.write(f"self.input_filename = {self.input_filename}, cwd = {os.getcwd()}\n")

         # Validate the input file.
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

      try:
         cmd = ["treesort"]

         # Always add the descriptor path.
         cmd.append(ScriptOption.DescriptorPath)
         cmd.append(self.job_data[InputParameter.DescriptorPath])

         # The "match on" options are mutually exclusive.
         if self.job_data[InputParameter.MatchOnStrain]:
            cmd.append(ScriptOption.MatchOnStrain)

         elif self.job_data[InputParameter.MatchOnEPI]:
            cmd.append(ScriptOption.MatchOnEPI)

         elif self.job_data[InputParameter.MatchOnRegex]:
            cmd.append(ScriptOption.MatchOnRegex)
            cmd.append(self.job_data[InputParameter.MatchOnRegex])

         # No collapse
         if self.job_data[InputParameter.NoCollapse]:
            cmd.append(ScriptOption.NoCollapse)

          # Always add the output path.
         cmd.append(ScriptOption.OutputPath)
         cmd.append(self.job_data[InputParameter.OutputPath])

         # Equal rates
         if self.job_data[InputParameter.EqualRates]:
            cmd.append(ScriptOption.EqualRates)

         # Is time scaled (timetree)
         if self.job_data[InputParameter.IsTimeScaled]:
            cmd.append(ScriptOption.IsTimeScaled)

         # TODO: uncomment
         # result = subprocess.check_call(cmd, shell=False)
         print(cmd)

      except ValueError as e:
         sys.stderr.write(f"Error preparing dataset:\n {e}\n")
         return False
         
      return True