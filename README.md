# ``cgb_input_generator``
A script to create a [`cgb`](https://github.com/ErillLab/cgb) ready input file
## About
``cgb_input_generator`` exists to create and populate a correctly formatted JSON file, complete with all necessary information required to successfully initiate an analysis of a transcription factor, or group of related transcripted factors, with ``cgb``: _a python library for comparative genomics of transcriptional regulation in Bacteria_.

Upon completion of a successful input file generation, the resulting JSON file will be both ``cgb`` and human readable. Human readability facilitates the option to make last minute minor adjustments to the generated file with no need to repeat the entire process.
## Getting Started
### TASK 1:
  * Using the provided TF protein accession numbers , cgb_input_generator 
    performs a BLASTP search for each TF, limited by e-value and (optionally) 
    by a taxonomic organism ID. (See input_test.md for details)
  
  * If no taxonomic ID is specified, BLAST hits will be obtained up to a
    maximum specified e-value (1e-10 by default). If, however,  a taxonomic ID 
    is specified, the BLASTP search will be constrained to that taxon (e.g.
    class, order...) well as by the specified e-value.
  
  * Identified HSPs from BLASTP searches are consolidated into a list of
    unique hits, relative to each other, and to the TFs under study 
### TASK 2:
  * Once HSPs for all TFs have been returned, the next task is to identify and retrieve the 'best' genome record to which the HSP belongs, and collect the following coding sequence (CDS) data: 
      1. genome acccession.version ID 
      2. start position
      3. stop position
      4. strand directionality
      5. p_score (defined in genome_record_retrieval() function)
      * 'Best' is defined in detail in the genome_record_retrieval() function,
       but generally, complete RefSeq genomes are most desired.     
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
  * Following genome retrieval, the Entrez database is then searched for any 
    unidentified genome records such as plasmids for complete records and contig
    records for WGS records
### TASK 3:
  * Following an exhaustive search for all relevant genome records for each HSP, the plasmid or contig record(s) are added to the corresponding HSP under a new key called, 'genomes'.
### TASK 4:
  * The last major task to complete is sampling HSPs
     * There are three sampling methods from which the user may choose, and will
       indicate which of the three method(s) in the input_parameters section
       of the input file. (see test_input.md)
       
       * **METHODS**
         1. Select only one genome record per taxonomic level (e.g. genus)
	
         2. Select only records for HSPs with pair-wise promoter region identities below a user-specified threshold (e.g. 90% identity)
         3. Select all HSPs identified

### TASK 5:
  * Following the protocol indicated in TASK 4, all relevant information is compiled into a JSON file which will be 100% cgb ready, and will still be human readable to allow for post-hoc edits.
