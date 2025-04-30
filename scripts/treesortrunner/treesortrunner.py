
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

   # The full path to the input directory.
   #input_path: str = None

   # The JSON data for the job.
   job_data: dict = None

   # The full path to the output directory.
   #output_path: str = None
   

   # C-tor
   def __init__(self, job_data):
         
      # Determine the base URL.
      if "P3_BASE_URL" in os.environ:
         self.base_url = os.environ["P3_BASE_URL"]
      else:
         self.base_url = "https://www.bv-brc.org"

      # Set and validate the job data.
      self.job_data = job_data
      if not self.job_data:
         raise ValueError("No job data was provided")
      elif not TreeSortRunner.is_job_data_valid(self.job_data):
         raise ValueError("The job data is invalid")
      
      
   @staticmethod
   def is_job_data_valid(job_data: dict) -> bool:

      try:
         if not job_data or not isinstance(job_data, dict) or len(job_data) == 0:
            raise Exception("job_data is empty")
         
         # Validate the input_source.
         if job_data[InputParameter.InputSource.value] not in InputSource:
            raise ValueError("job_data.input_source is not a valid input source") 

         # If input_source is fasta_data, make sure an input_fasta_data value was provided.
         if job_data[InputParameter.InputSource.value] == InputSource.FastaData.value and \
         not not job_data[InputParameter.InputFastaData]:
            raise ValueError("The input FASTA data is invalid")
         
         # If input_source is fasta_file, make sure an input_fasta_file value was provided.
         if job_data[InputParameter.InputSource.value] == InputSource.FastaFile.value and \
         not not job_data[InputParameter.InputFastaFile]:
            raise ValueError("The input FASTA file is invalid")
         
         # If input_source is fasta_file_id, make sure an input_fasta_file_id value was provided.
         if job_data[InputParameter.InputSource.value] == InputSource.FastaFileID.value and \
         not not job_data[InputParameter.InputFastaFileID]:
            raise ValueError("The input FASTA file ID is invalid")

         # Validate the method
         if job_data[InputParameter.Method] not in Method:
            raise ValueError("job_data.method is not a valid method") 

         # Validate the output path.
         if not job_data[InputParameter.OutputPath]:
            raise ValueError("The output path is invalid")

         # Validate the reference segment and provide a default if not provided.
         refSegment = safeTrim(job_data[InputParameter.RefSegment.value])
         if not refSegment:
            refSegment = DEFAULT_REF_SEGMENT
         elif not refSegment in VALID_SEGMENTS:
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

      sys.stdout.write("In prepare_dataset\n\n")
      
      try:
         # Example usage: ./prepare_dataset.sh --segments "HA,NA" segments.fasta HA myoutdir
         cmd = ["prepare_dataset.sh"]

         inference_type = self.job_data[InputParameter.InferenceType.value]
         if (inference_type and inference_type == InferenceType.FastTree.value):
            cmd.append(ScriptOption.FastTree.value)

         segments = self.job_data[InputParameter.Segments.value]
         if segments:
            cmd.append(ScriptOption.Segments.value)
            cmd.append(segments)

         cmd.append(self.input_path)

         refSegment = self.job_data[InputParameter.RefSegment.value]
         if refSegment:
            # TODO: validate refSegment
            cmd.append(refSegment)

         cmd.append(self.output_path)

         # TODO: uncomment
         # result = subprocess.check_call(cmd, shell=False)
         print(f"{"".join(cmd)}\n\n")

      except ValueError as e:
         sys.stderr.write(f"Error preparing dataset:\n {e}\n")
         return False
         
      return True
     

   # Determine the source of the FASTA input file and prepare it for use.
   def prepare_input_file(self) -> bool:

      sys.stdout.write("In prepare_input_file\n\n")

      try:
         # The input source determines how the input file is provided.
         input_source = safeTrim(self.job_data[InputParameter.InputSource.value])

         if input_source == InputSource.FastaData.value:
         
            # Copy user data to the input file.
            try:
               with open(self.default_input_filename, "w+") as input_file:
                  input_file.write(self.job_data[InputParameter.InputFastaData.value])

            except Exception as e:
               raise IOError(f"Error copying FASTA data to input file:\n {e}\n")

         elif input_source == InputSource.FastaFile.value:

            # The input file will be in the local filesystem.
            self.input_filename = safeTrim(self.job_data[InputParameter.InputFastaFile])
            if len(self.input_filename) == 0:
               raise ValueError("Invalid input filename")
            
         elif input_source == InputSource.FastaFileID.value:

            # Fetch the input file from the workspace.
            try:
               # TODO: Is this syntax correct?
               fetch_fasta_cmd = ["p3-cp", f"ws:{self.job_data[InputParameter.InputFastaFileID.value]}", self.input_filename]
               subprocess.check_call(fetch_fasta_cmd, shell=False)

            except Exception as e:
               raise IOError("Error copying FASTA file from workspace:\n %s" % str(e))

         else:
            raise ValueError(f"Invalid input source: {input_source}")

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
         print(f"{"".join(cmd)}\n\n")

      except ValueError as e:
         sys.stderr.write(f"Error preparing dataset:\n {e}\n")
         return False
         
      return True