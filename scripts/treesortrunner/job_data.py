
from common import InferenceType, InputSource, Method
from dataclasses import dataclass
from typing import Optional

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