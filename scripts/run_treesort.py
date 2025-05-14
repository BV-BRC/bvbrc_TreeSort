#!/usr/bin/env python

import argparse
from dataclasses import dataclass
from enum import Enum
import json
import os
import re
import subprocess
import sys
import traceback
from typing import Optional

#-----------------------------------------------------------------------------------------------------------------------------
# Define constants
#-----------------------------------------------------------------------------------------------------------------------------

# The default reference segment.
DEFAULT_REF_SEGMENT = "HA"

# Characters to remove from FASTA headers.
INVALID_FASTA_CHARS = "[\\['\"(),;|:\\]]"

# A list of valid segments for the influenza virus.
VALID_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]


#-----------------------------------------------------------------------------------------------------------------------------
# Define enums
#-----------------------------------------------------------------------------------------------------------------------------

class Filename(str, Enum):
   Descriptor = "descriptor.csv"

class InferenceType(str, Enum):
   FastTree = "FastTree"
   IQTree = "IQ-Tree"
 
class InputSource(str, Enum):
   FastaData = "fasta_data"
   FastaFile = "fasta_file"
   FastaFileID = "fasta_file_id"
   
class Method(str, Enum):
   Local = "local"
   MinCut = "mincut"

# Command-line script flags / optional arguments.
class ScriptOption(str, Enum):
   CladesPath = "--clades"
   DescriptorPath = "-i"
   Deviation = "--dev"
   EqualRates = "--equal-rates"
   FastTree = "--fast"
   IsTimeScaled = "--timetree"
   MatchOnEPI = "--match-on-epi"
   MatchOnRegex = "--match-on-regex"
   MatchOnStrain = "--match-on-strain"
   Method = "-m"
   NoCollapse = "--no-collapse"
   OutputPath = "-o"
   PValue = "--pvalue"
   Segments = "--segments"


#-----------------------------------------------------------------------------------------------------------------------------
# Define utility functions
#-----------------------------------------------------------------------------------------------------------------------------

# Trim a string that's possibly null and always return a trimmed, non-null value.
def safeTrim(text: str):
   if not text:
      return ""
   else:
      return text.strip()


#-----------------------------------------------------------------------------------------------------------------------------
# Define classes
#-----------------------------------------------------------------------------------------------------------------------------

# The contents of the jobdesc.json file.
@dataclass
class JobData:
   clades_path: str
   descriptor_path: str
   deviation: float
   equal_rates: bool
   inference_type: InferenceType
   input_fasta_data: str
   input_fasta_file: str
   input_fasta_file_id: str
   input_source: InputSource
   is_time_scaled: bool
   match_on_epi: bool
   match_on_regex: Optional[str]
   match_on_strain: bool
   method: Method
   no_collapse: bool
   output_path: str
   p_value: float
   ref_segment: str
   segments: str


# The class responsible for processing input parameters, preparing the input file, and 
# running the TreeSort application.
class TreeSortRunner:

   # The base URL
   base_url: str = None

   # The default name of the input FASTA file.
   default_input_filename: str = "input.fasta"

   # The name of the FASTA file to use as input. This is specified in job_data.
   input_filename: str = None

   # The JSON data for the job.
   job_data: JobData = None
   
   # The directory where the output files will be created.
   staging_directory = None

   # The directory where the scripts will be run.
   work_directory = None


   # C-tor
   def __init__(self, job_data: JobData, staging_directory: str, work_directory: str):
         
      # Set and validate the staging directory.
      self.staging_directory = staging_directory
      if not self.staging_directory or len(self.staging_directory) < 1:
         raise ValueError("The staging directory parameter is invalid")
      
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

      # Validate the job data.
      if not self.is_job_data_valid():
         raise ValueError("Job data in the constructor is invalid")
      
      
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
         if self.job_data.input_source == InputSource.FastaFileID.value:
            if not self.job_data.input_fasta_file_id:
               raise ValueError("The input FASTA file ID is invalid")
            elif self.job_data.input_fasta_file_id.startswith("ws:"):
               self.job_data.input_fasta_file_id = self.job_data.input_fasta_file_id[3:]

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
         cmd.append(self.staging_directory)

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

            # Use the default filename.
            self.input_filename = self.default_input_filename

            # Copy the input file from the workspace to the working directory.
            try:
               fetch_fasta_cmd = ["p3-cp", f"ws:{self.job_data.input_fasta_file_id}", f"{self.work_directory}/{self.input_filename}"]
               subprocess.check_call(fetch_fasta_cmd, shell=False)

            except Exception as e:
               raise IOError("Error copying FASTA file from workspace:\n %s" % str(e))

         else:
            raise ValueError(f"Invalid input source: {input_source}")

         # Set the input filename and include its full path.
         self.input_filename = f"{self.work_directory}/{self.input_filename}"
         
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


   # Run TreeSort in a sub-process.
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




def main(argv=None):
    
   if argv is None:
      argv = sys.argv[1:]  # Exclude the script name

   # Create an argument parser.
   parser = argparse.ArgumentParser(description="A script to run TreeSort")
   parser.add_argument("-j", "--job-filename", dest="job_filename", help="A JSON file with the job description", required=True)
   parser.add_argument("-s", "--staging-directory", dest="staging_directory", help="The directory where output files will be created", required=True)
   parser.add_argument("-w", "--work-directory", dest="work_directory", help="The directory where the scripts will be run", required=True)
   
   args = parser.parse_args()

   # Validate the job filename parameter.
   job_filename = safeTrim(args.job_filename)
   if len(job_filename) == 0:
      traceback.print_exc(file=sys.stderr)
      sys.stderr.write("Invalid job filename parameter\n")
      sys.exit(-1)

   # Validate the staging directory parameter.
   staging_directory = safeTrim(args.staging_directory)
   if len(staging_directory) == 0:
      traceback.print_exc(file=sys.stderr)
      sys.stderr.write("Invalid staging directory parameter\n")
      sys.exit(-1)

   # Validate the work directory parameter.
   work_directory = safeTrim(args.work_directory)
   if len(work_directory) == 0:
      traceback.print_exc(file=sys.stderr)
      sys.stderr.write("Invalid work directory parameter\n")
      sys.exit(-1)

   # Load job data
   job_data = None
   try:
      with open(job_filename, "r", encoding="utf-8") as job_file:
         job_dict = json.load(job_file)
         job_data = JobData(**job_dict)

   except Exception as e:
      traceback.print_exc(file=sys.stderr)
      sys.stderr.write(f"Invalid job file:\n{e}\n")
      sys.exit(-1)

   try:
      # Create a TreeSortRunner instance
      runner = TreeSortRunner(job_data, staging_directory, work_directory)

   except Exception as e:
      traceback.print_exc(file=sys.stderr)
      sys.stderr.write(f"Unable to create an instance of TreeSortRunner:\n{e}\n")
      sys.exit(-1)
   
   # Prepare the input file
   if not runner.prepare_input_file():
      traceback.print_exc(file=sys.stderr)
      sys.stderr.write("An error occurred in prepare_input_file\n")
      sys.exit(-1)

   # Prepare the dataset
   if not runner.prepare_dataset():
      traceback.print_exc(file=sys.stderr)
      sys.stderr.write("An error occurred in prepare_dataset\n")
      sys.exit(-1)

   # Run TreeSort
   """if not runner.tree_sort():
      traceback.print_exc(file=sys.stderr)
      sys.stderr.write("An error occurred in tree_sort\n")
      sys.exit(-1)"""


if __name__ == "__main__" :
   main()
   
