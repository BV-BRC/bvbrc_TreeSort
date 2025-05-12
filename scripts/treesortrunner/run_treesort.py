#!/usr/bin/env python

import argparse
from treesortrunner.common import safeTrim
from treesortrunner.treesortrunner import TreeSortRunner
from job_data import JobData
import json
import os
import sys
import traceback


def main(argv=None):
    
   if argv is None:
      argv = sys.argv[1:]  # Exclude the script name

   # Create an argument parser.
   parser = argparse.ArgumentParser(description="A script to run TreeSort")
   parser.add_argument("-j", "--job-filename", dest="job_filename", help="A JSON file for the job", required=True)
   parser.add_argument("-w", "--work-directory", dest="work_directory", help="The directory where the scripts will be run", required=True)
   
   args = parser.parse_args()

   # Validate the job filename parameter.
   job_filename = safeTrim(args.job_filename)
   if len(job_filename) == 0:
      traceback.print_exc(file=sys.stderr)
      sys.stderr.write("Invalid job filename parameter\n")
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
      runner = TreeSortRunner(job_data, work_directory)

   except Exception as e:
      traceback.print_exc(file=sys.stderr)
      sys.stderr.write(f"Unable to create an instance of TreeSortRunner:\n{e}\n")
      sys.exit(-1)
   
   # Prepare the input file
   if not runner.prepare_input_file():
      traceback.print_exc(file=sys.stderr)
      sys.stderr.write("An error occurred in prepare_input_file\n")
      sys.exit(-1)

   # Create the output directory if it doesn't exist.
   output_path = os.path.abspath(job_data.output_path)
   if not os.path.exists(output_path):
      os.mkdir(output_path)

   # Go to the output directory.
   #os.chdir(output_path)

   # Prepare the dataset
   if not runner.prepare_dataset():
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
   
