#!/usr/bin/env python

import argparse
from dataclasses import dataclass
from enum import Enum
import json
import os
import re
#import shutil # For deleting files and directories
import subprocess
import sys
import time
import traceback
from typing import Optional

#-----------------------------------------------------------------------------------------------------------------------------
# Define constants
#-----------------------------------------------------------------------------------------------------------------------------

# The default reference segment.
DEFAULT_REF_SEGMENT = "HA"

# The name of the descriptor file.
DESCRIPTOR_FILE_NAME = "descriptor.csv"

# The default name of the input FASTA file.
INPUT_FASTA_FILE_NAME = "input.fasta"

# Characters to remove from FASTA headers.
INVALID_FASTA_CHARS = "[\\['\"(),;:\\]]"

# A list of valid segments for the influenza virus.
VALID_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]


#-----------------------------------------------------------------------------------------------------------------------------
# Define enums
#-----------------------------------------------------------------------------------------------------------------------------

class InputSource(str, Enum):
   FastaData = "fasta_data"
   FastaFileID = "fasta_file_id"
   PreparedFiles = "prepared_files"
   
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

class TreeInference(str, Enum):
   FastTree = "FastTree"
   IQTree = "IQ-Tree"
   # RAxML = "RAxML" # Not implemented yet
 

#-----------------------------------------------------------------------------------------------------------------------------
# Define utility functions
#-----------------------------------------------------------------------------------------------------------------------------

# Trim a string that's possibly null and always return a trimmed, non-null value.
def safeTrim(text: str|None):
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
   deviation: float
   equal_rates: bool
   input_existing_directory: Optional[str]
   input_fasta_data: Optional[str]
   input_fasta_file_id: Optional[str]
   input_source: InputSource
   is_time_scaled: bool
   match_on_epi: bool
   match_on_regex: Optional[str]
   match_on_strain: bool
   method: Method
   no_collapse: bool
   output_path: str
   prepare_dataset: bool
   p_value: float
   ref_segment: str
   ref_tree_inference: TreeInference
   segments: Optional[str]
   
   

# The class responsible for processing input parameters, preparing the input file, and 
# running the TreeSort application.
class TreeSortRunner:

   # The base URL
   base_url: str

   # Segments found in the input FASTA file.
   fasta_segments: list[str]

   # the name of the directory containing the input FASTA file.
   input_directory: str

   # The name of the FASTA file to use as input. This is specified in job_data.
   input_filename: str

   # The JSON data for the job.
   job_data: JobData
   
   # The directory where the output files will be created.
   staging_directory: str

   # The directory where the scripts will be run.
   work_directory: str


   # C-tor
   def __init__(self, input_directory: str, job_data: JobData, staging_directory: str, work_directory: str):
         
      # Set and validate the input directory.
      self.input_directory = input_directory
      if not self.input_directory or len(self.input_directory) < 1:
         raise ValueError("The input directory parameter is invalid")
      
      # Set and validate the staging directory.
      self.staging_directory = staging_directory
      if not self.staging_directory or len(self.staging_directory) < 1:
         raise ValueError("The staging directory parameter is invalid")
      
      # Set and validate the work directory.
      self.work_directory = work_directory
      if not self.work_directory or len(self.work_directory) < 1:
         raise ValueError("The work directory parameter is invalid")
      
      # Initialize member variables.
      self.fasta_segments = []

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
      
      
   # Find the segment name in a FASTA header.
   def get_segment_from_header(self, header: str) -> str | None:

      result = None

      for segment_name in VALID_SEGMENTS:
         if re.search(rf"\|{segment_name}\|", header, re.IGNORECASE):
            result = segment_name
            break

      return result
   

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
         
         if self.job_data.prepare_dataset:

            # The dataset will be prepared before running TreeSort.

            if self.job_data.input_source == InputSource.FastaData.value:

               # Make sure an input_fasta_data value was provided.
               if not self.job_data.input_fasta_data:
                  raise ValueError("The input FASTA data is invalid")

            elif self.job_data.input_source == InputSource.FastaFileID.value:

               # Make sure an input_fasta_file_id value was provided.
               if not self.job_data.input_fasta_file_id:
                  raise ValueError("The input FASTA file ID is invalid")
               
               elif self.job_data.input_fasta_file_id.startswith("ws:"):

                  # Remove the "ws:" prefix from the directory name.
                  self.job_data.input_fasta_file_id = self.job_data.input_fasta_file_id[3:]
         else:

            # If input_source is prepared_files, make sure an input_existing_directory value was provided.
            if self.job_data.input_source != InputSource.PreparedFiles.value:
               raise ValueError(f"The input source {self.job_data.input_source} is invalid")
            
            if not self.job_data.input_existing_directory:
               raise ValueError("An existing directory of prepared files is required")
            
            if self.job_data.input_existing_directory.startswith("ws:"):

               # Remove the "ws:" prefix from the directory name.
               self.job_data.input_existing_directory = self.job_data.input_existing_directory[3:]

         # Validate the inference method
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
   

   def prepare_dataset(self) -> bool:

      sys.stdout.write("In prepare_dataset\n\n")
      
      # The result status defaults to false.
      result_status = False

      # NOTE: prepare_dataset.sh deleted and recreated the working directory
      # at this point, but we don't need to do that here as the working directory
      # was just created and is empty.
      
      try:

         # Split the input FASTA file into multiple files by segment.
         self.split_fasta_by_segment()

         # Run mafft to align the segments.


         # Remove the segment-specific FASTA files.


         # Iterate over all valid segments and:
         # 1. Build reference trees for each segment. Note that the reference segment
         #    might be handled differently.
         #
         # 2. Root trees with treetime


         # Create the descriptor file.

         result_status = True

      except Exception as e:
         sys.stderr.write(f"Error preparing dataset:\n {e}\n")
         return False

      return result_status


   # Determine the source of the FASTA input file and prepare it for use.
   def prepare_input_file(self) -> bool:

      sys.stdout.write("In prepare_input_file\n\n")

      try:
         # The input source determines how the input file is provided.
         input_source = safeTrim(self.job_data.input_source)

         if input_source == InputSource.FastaFileID.value:

            # Populate the input filename, including its full path.
            self.input_filename = f"{self.input_directory}/{INPUT_FASTA_FILE_NAME}"

            try:
               # Copy the input file from the workspace to the working directory.
               fetch_fasta_cmd = ["p3-cp", f"ws:{self.job_data.input_fasta_file_id}", self.input_filename]
               subprocess.call(fetch_fasta_cmd, shell=False)

            except Exception as e:
               raise IOError("Error copying FASTA file from workspace:\n %s" % str(e))

         elif input_source == InputSource.FastaData.value:
      
            # Create the input filename, including its full path.
            self.input_filename = f"{self.input_directory}/{INPUT_FASTA_FILE_NAME}"

            try:
               # Create a file that contains the input FASTA data.
               with open(self.input_filename, "w+") as input_file:
                  input_file.write(str(self.job_data.input_fasta_data))

            except Exception as e:
               raise IOError(f"Error copying FASTA data to input file:\n {e}\n")

         else:
            raise ValueError(f"Invalid input source: {input_source}")

         # Validate the input filename.
         if not self.input_filename:
            raise ValueError("Invalid input filename (empty)")
         
         # Validate the input FASTA file.
         if not os.path.exists(self.input_filename) or os.path.getsize(self.input_filename) == 0:
            raise IOError("Input FASTA file is invalid or empty")

         # Remove invalid characters from FASTA headers.
         # TODO: Since we're already iterating over the file contents, this might be a good place
         # to split the file into multiple files by segment.
         """with open(self.input_filename, 'r+') as f:
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
            f.truncate()"""

      except Exception as e:
         sys.stderr.write(f"Error processing input file:\n {e}\n")
         return False
      
      return True


   # Run prepare_dataset.sh to build alignments and trees and compile a descriptor file.
   def run_prepare_dataset(self) -> bool:

      sys.stdout.write("In run_prepare_dataset\n\n")
      
      # The result status defaults to false.
      result_status = False

      try:
         # Example usage: ./prepare_dataset.sh segments.fasta HA myoutdir
         # TODO: Replace this with the actual script name!
         cmd = ["prepare_dataset_051625.sh"]

         # Should --fast be added?
         # TODO: Support FastTree, IQ-Tree, and RAxML
         if (self.job_data.ref_tree_inference and self.job_data.ref_tree_inference == TreeInference.FastTree.value):
            cmd.append(ScriptOption.FastTree.value)

         # The segments (optional)
         if self.job_data.segments:
            cmd.append(ScriptOption.Segments.value)
            cmd.append(self.job_data.segments)

         # The input FASTA file
         cmd.append(self.input_filename)

         # The reference segment
         refSegment = self.job_data.ref_segment
         if refSegment:
            cmd.append(refSegment)

         # The output directory
         cmd.append(self.work_directory)

         # TEST
         print(f"{' '.join(cmd)}\n\n")

         result = subprocess.call(cmd, shell=False)
         if result == 0:
            result_status = True

      except ValueError as e:
         sys.stderr.write(f"Error preparing dataset:\n {e}\n")
         return False
         
      return result_status
     

   # Split a FASTA file containing sequences for multiple segments into individual 
   # files for each segment found in the FASTA headers.
   def split_fasta_by_segment(self) -> bool:

      # A dictionary of lines of FASTA keyed by segment name.
      fasta_by_segment = {}

      with open(self.input_filename, 'r') as f:
         
         current_segment = None

         # Iterate over every line in the file.
         for line in f:

            #print(f"{line}\n")

            if line.startswith(">"):
               header = line.strip()
               # header = re.sub(INVALID_FASTA_CHARS, "", header)
               # header = re.sub(" ", "_", header)
               current_segment = self.get_segment_from_header(header)

               #print(f"In the header, current segment = {current_segment}")
               
               line = header

            if not current_segment:
               continue
            
            # Append a new line character to the end of the line.
            line = f"{line}\n"

            # Update the list of FASTA segments we have encountered.
            if not current_segment in self.fasta_segments:
               self.fasta_segments.append(current_segment)

               print(f"Just added segment {current_segment}")

            # Found the segment in the header.
            if current_segment in fasta_by_segment:
               fasta_by_segment[current_segment] += line
            else:
               fasta_by_segment[current_segment] = line


         # Write the contents of fasta_by_segment to new FASTA files.
         for segment in fasta_by_segment.keys():

            print(f"segment key = {segment}\n")

            fasta = fasta_by_segment.get(segment)
            if not fasta or len(fasta) < 1:
               print("no fasta for this key")
               continue 
            
            try:
               with open(f"{self.work_directory}/{segment}-{INPUT_FASTA_FILE_NAME}", "w") as fasta_file:
                  fasta_file.write(fasta)

            except Exception as e:
               sys.stderr.write(f"Error creating FASTA file for {segment}:\n {e}\n")
               return False
         
      return True



   # Run TreeSort in a sub-process.
   def tree_sort(self) -> bool:

      sys.stdout.write("In tree_sort\n\n")
      
      # The result status defaults to false.
      result_status = False

      try:
         cmd = ["treesort"]

         # The clades output path
         clades_path = safeTrim(self.job_data.clades_path)
         if len(clades_path) > 0:
            cmd.append(ScriptOption.CladesPath.value)
            cmd.append(clades_path)

         # The descriptor file is in the working directory.
         cmd.append(ScriptOption.DescriptorPath.value)
         cmd.append(f"{self.work_directory}/{DESCRIPTOR_FILE_NAME}")

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

          # Always add the output path (staging directory)
         cmd.append(ScriptOption.OutputPath.value)
         cmd.append(self.staging_directory)
         
         # Equal rates
         if self.job_data.equal_rates:
            cmd.append(ScriptOption.EqualRates.value)

         # Is time scaled (timetree)
         if self.job_data.is_time_scaled:
            cmd.append(ScriptOption.IsTimeScaled.value)

         # TEST
         print(f"{' '.join(cmd)}\n\n")

         # Run the command
         result = subprocess.call(cmd, shell=False)
         if result == 0:
            result_status = True

      except ValueError as e:
         sys.stderr.write(f"Error preparing dataset:\n {e}\n")
         return False
         
      return result_status



def main(argv=None):
    
   if argv is None:
      argv = sys.argv[1:]  # Exclude the script name

   # Create an argument parser.
   parser = argparse.ArgumentParser(description="A script to run TreeSort")
   parser.add_argument("-i", "--input-directory", dest="input_directory", help="The directory that will contain the FASTA input file(s)", required=True)
   parser.add_argument("-j", "--job-filename", dest="job_filename", help="A JSON file with the job description", required=True)
   parser.add_argument("-s", "--staging-directory", dest="staging_directory", help="The directory where output files will be created", required=True)
   parser.add_argument("-w", "--work-directory", dest="work_directory", help="The directory that will contain generated intermediate files", required=True)
   
   args = parser.parse_args()

   # Validate the input directory parameter.
   input_directory = safeTrim(args.input_directory)
   if len(input_directory) == 0:
      traceback.print_exc(file=sys.stderr)
      sys.stderr.write("Invalid input directory parameter\n")
      sys.exit(-1)

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
      runner = TreeSortRunner(input_directory, job_data, staging_directory, work_directory)

   except Exception as e:
      traceback.print_exc(file=sys.stderr)
      sys.stderr.write(f"Unable to create an instance of TreeSortRunner:\n{e}\n")
      sys.exit(-1)
   
   # Should we prepare the dataset?
   if runner.job_data.prepare_dataset:
      
      # Prepare the input file
      if not runner.prepare_input_file():
         traceback.print_exc(file=sys.stderr)
         sys.stderr.write("An error occurred in prepare_input_file\n")
         sys.exit(-1)

      # Prepare the dataset
      if not runner.run_prepare_dataset():
         traceback.print_exc(file=sys.stderr)
         sys.stderr.write("An error occurred in prepare_dataset\n")
         sys.exit(-1)

   # Run TreeSort
   if not runner.tree_sort():
      traceback.print_exc(file=sys.stderr)
      sys.stderr.write("An error occurred in tree_sort\n")
      sys.exit(-1)


if __name__ == "__main__" :
   main()
   
