# `CGB_input`
A script to create a [`cgb`](https://github.com/ErillLab/cgb) ready input file
## About
`CGB_input` exists to create and populate a correctly formatted JSON file, complete with all necessary information required to successfully initiate an analysis of a transcription factor, or group of related transcripted factors, with `cgb`: *a python library for comparative genomics of transcriptional regulation in Bacteria*.

Upon completion of a successful input file generation, the resulting JSON file will be both `cgb` and human readable. Human readability facilitates the option to make last minute minor adjustments to the generated file with no need to repeat the entire process.
## Getting Started
### Using The Input File Generator:
* To use this software, you (the end user/researcher) should have a basic understanding of navigating your computer's directories from the command line/terminal.
* It's generally recommended you run this command:
  * `git clone https://github.com/shane806/CGB_input.git`
* Running this command will download the entire repository and create a directory called `CGB_input` in your current working directory.
* Then change your current working directory to this newly created directory:
  * `cd CGB_input`
* In the main `CGB_input` directory are various files and directories which vary in their purpose, such as:
	* `test_input.json` and `test_input.md`: 
		* The `test_input.json` file is an example of the format of the JSON file which you, the end user/researcher must provide. 
		* The `test_input.md` file is an overview/detailed explanation of the  JSON input file you will provide.
		* **MAKE SURE YOUR INPUT JSON FILE IS IN THE MAIN `CGB_input` DIRECTORY**
    * `cgb_input_run.py`: The `cgb_input_run.py` file is, other than the input JSON file you provide, is the only file you will (or should) modify when running the software. On the 3rd line of `cgb_input_run.py`, you will enter the name of the file, enclosed by quotation marks, to the right of the equal sign operator, where the left side is INPUT_FILE.
      * At the time of cloning the repository, this is the line in the file: `INPUT_FILE = "test_input.json"`
      * **THIS MUST BE A JSON FILE AND THE FILE MUST BE FORMATTED FOLLOWING THE PROTOCOL LAID OUT IN `test_input.md`**
	* `environment.yml`: This is the dependencies file for an anaconda environment.
    * `README.md`: The explanation of the software's purpose and basic instructions for using this software.
    * `cgb_input/`: The `cgb_input/` directory, which holds all the python files used to generate the `cgb` ready input file. 
	* `log_files/`: The `log_files/` directory contains, or will contain, files which are generated throughout various phases of the input file generation process, and will show you when, where, why, and how many hits are removed/filtered out.
    * `ouput_files/`: The `output_files/` directory simply contains the final product of the input file generation process. The newly created file placed in this directory after a successful input file generation can be used with no further editing in the `cgb` software.

### TASK 1:
* Using the provided Transcription Factor (TF) protein accession numbers , the CGB_input program performs a BLASTP search for each TF, limited by e-value and (optionally) by a taxonomic organism ID. (See input_test.md for details)

* If no taxonomic ID is specified, BLASTP hits will be obtained up to a maximum specified e-value (1e-10 by default). If, however,  a taxonomic ID is specified, the BLASTP search will be constrained to that taxon (e.g. class, order...) well as by the specified e-value.

* Identified HSPs from BLASTP searches are consolidated into a list of unique hits, relative to each other, and to the TFs under study.

### TASK 2:
* Once HSPs for all TFs have been returned, the next task is to identify and retrieve the 'best' genome record to which the HSP belongs, and collect the following coding sequence (CDS) data:
	1. genome acccession.version ID
	2. start position
	3. stop position
	4. strand directionality
	5. p_score (defined in genome_record_retrieval() function)
	* 'Best' is defined in detail in the genome_record_retrieval() function, but generally, complete RefSeq genomes are most desired.
* Once all pieces of information are collected, they are organized into a dictionary of dictionaries, where the main, 'outer', key is the protein accession ID of the HSP, and all the genome and CDS information is stored in a dictionary as the value for the HSP
  * ex.
```
{ 'ZZ_123456.7' :{
                 'genome_accession' : 'AA_987654.3',
                 'start_pos' : 1234,
                 'stop_pos' : 5678,
                 'strand' : '+',
                 'p_score' : 7
                 }
}
```
* Following genome retrieval, the Entrez database is then searched for any unidentified genome records such as plasmids for complete records and contig records for WGS records.

### TASK 3:
  * Following an exhaustive search for all relevant genome records for each HSP, the plasmid or contig record(s) are added to the corresponding HSP under a new key called, 'genomes'.

### TASK 4:
  * The last major task to complete is sampling HSPs
  * There are three sampling methods from which the user may choose, and will indicate which of the three method(s) in the input_parameters section of the input file. (see test_input.md)

  * #### METHODS
    1. Select only one genome record per taxonomic level (e.g. One genome record per genus)
	2. Select only records for HSPs with promoter region identities below a user-specified threshold, compared following a pair-wise protocol. (e.g. 90% identity)

	3. Select all HSPs identified in BLASTP search (No post-BLASTP search filtering)

### TASK 5:
  * Following the protocol indicated in TASK 4, all relevant information is compiled into a JSON file which will be 100% cgb ready, and will still be human readable to allow for post-hoc edits.
