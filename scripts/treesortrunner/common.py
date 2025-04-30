
from enum import Enum

#-----------------------------------------------------------------------------------------------------------------------------
# Define constants
#-----------------------------------------------------------------------------------------------------------------------------

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

class InputParameter(str, Enum):
   CladesPath = "clades_path"
   DescriptorPath = "descriptor_path"
   Deviation = "deviation"
   EqualRates = "equal_rates"
   InferenceType = "inference_type"
   InputFastaData = "input_fasta_data"
   InputFastaFile = "input_fasta_file"
   InputFastaFileID = "input_fasta_file_id"
   InputSource = "input_source"
   IsTimeScaled = "is_time_scaled"
   MatchOnEPI = "match_on_epi"
   MatchOnRegex = "match_on_regex"
   MatchOnStrain = "match_on_strain"
   Method = "method"
   NoCollapse = "no_collapse"
   OutputPath = "output_path"
   PValue = "p_value"
   RefSegment = "ref_segment"
   Segments = "segments"
 
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
# Define functions
#-----------------------------------------------------------------------------------------------------------------------------

# Trim a string that's possibly null and always return a non-null value.
def safeTrim(text: str):
   if not text:
      return ""
   else:
      return text.strip()
