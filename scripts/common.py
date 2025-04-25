
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
   FastaFileID = "fasta_file_id"
   
class Method(str, Enum):
   Local = "local"
   MinCut = "mincut"




#-----------------------------------------------------------------------------------------------------------------------------
# Define functions
#-----------------------------------------------------------------------------------------------------------------------------

# Trim a string that's possibly null and always return a non-null value.
def safeTrim(text: str):
   if not text:
      return ""
   else:
      return text.strip()
