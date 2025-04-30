
import argparse
from .common import safeTrim
from .treesortrunner import TreeSortRunner
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
   parser.add_argument("-o", "--output-dir-name", dest="output_dir_name", help="The output directory name (defaults to current directory)", required=False, default=".")

   args = parser.parse_args()

   # Validate the job filename
   job_filename = safeTrim(args.job_filename)
   if len(job_filename) == 0:
      traceback.print_exc(file=sys.stderr)
      sys.stderr.write("Invalid job filename parameter\n")
      sys.exit(-1)

   # Load job data
   job_data = None
   try:
      with open(job_filename, "r", encoding="utf-8") as job_file:
         job_data = json.load(job_file)

   except Exception as e:
      traceback.print_exc(file=sys.stderr)
      sys.stderr.write(f"Invalid job file:\n {e}\n")
      sys.exit(-1)

   
   # For debugging   
   #print(job_data)

   # Validate the job data
   if not TreeSortRunner.is_job_data_valid(job_data):
      sys.exit(-1)

   # Validate the output directory name
   output_dir_name = safeTrim(args.output_dir_name)
   if len(output_dir_name) == 0:
      traceback.print_exc(file=sys.stderr)
      sys.stderr.write("Invalid output directory name parameter\n")
      sys.exit(-1)

   # Create the output directory if it doesn't exist.
   output_dir_name = os.path.abspath(output_dir_name)
   if not os.path.exists(output_dir_name):
      os.mkdir(output_dir_name)

   # Create a TreeSortRunner instance
   runner = TreeSortRunner(job_data, output_dir_name)

   # Process the input file
   if not runner.prepare_input_file():
      traceback.print_exc(file=sys.stderr)
      sys.stderr.write("An error occurred in prepare_input_file\n")
      sys.exit(-1)

   # Go to the output directory.
   os.chdir(output_dir_name)

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
   
