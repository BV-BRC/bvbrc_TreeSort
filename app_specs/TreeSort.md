# Application specification: TreeSort

This is the application specification for the service with identifier TreeSort.

The backend script implementing the application is [App-TreeSort.pl](../service-scripts/App-TreeSort.pl).

The raw JSON file for this specification is [TreeSort.json](TreeSort.json).

This service performs the following task: TreeSort infers both recent and ancestral reassortment events along the branches of a phylogenetic tree of a fixed genomic segment.

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------  |
| clades_path | Clades path | string | ??? | |
| deviation | Deviation | float? | | 2 |
| equal_rates | Equal rates | boolean | | |
| input_existing_directory | Existing directory | string | | |
| input_fasta_data | Input FASTA data | string | | |
| input_fasta_file_id | Input FASTA file ID | wsid | | |
| input_source | Input source | fasta_data, fasta_file_id | :heavy_check_mark: | fasta_file_id |
| is_time_scaled | Is time scaled? | boolean | | false |
| match_on_epi  | Match on "EPI_ISL_XXX" field | boolean |  |  |
| match_on_regex | Match on RegEx | string (regex) |  | ??? |
| match_on_strain | Match on strain? | boolean |  | |
| method | Method | local, mincut | :heavy_check_mark: | local |
| no_collapse | No collapse | boolean | | |
| output_path | Output path | string | :heavy_check_mark: | |
| prepare_dataset | Prepare a FASTA file first? | boolean |  |  |
| p_value | P-value | float | | 0.001 |
| ref_segment | Reference Segment | string | :heavy_check_mark: | HA |
| ref_tree_inference | Reference tree inference method | FastTree, IQTree | | IQTree |
| segments | Segments | string | :heavy_check_mark: | All segments |






