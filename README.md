# ``cgb_input_generator``
A script to create a [`cgb`](https://github.com/ErillLab/cgb) ready input file
## About
``cgb_input_generator`` exists to create and populate a correctly formatted JSON file, complete with all necessary information required to successfully initiate an analysis of a transcription factor, or group of related transcripted factors, with ``cgb``: _a python library for comparative genomics of transcriptional regulation in Bacteria_.

Upon completion of a successful input file generation, the resulting JSON file will be both ``cgb`` and human readable. Human readability facilitates the option to make last minute minor adjustments to the generated file with no need to repeat the entire process.

(Might change the Getting Started section below into more of an overview of the program's functionality rather than the current input file overview? And integrate current Getting Started section into test_input.md?)
## Getting Started
###  ``cgb_input_generator`` Input File

1. **A list of Transcription Factors**
   * Transcription Factor proteins (TFs) are represented by an NCBI protein accession.version identifier. 
     * The accession component of the identifier will be alphanumeric, at a minimum, and some will also have an underscore. 
     * The version component is much simpler. It will contain only single digit numeric values.
		
2. **A list of each Transcription Factor’s respective binding sites.**
   * Each list contains experimentally validated DNA binding sites for the corresponding TF protein, represented as strings. 
   * Each string is a DNA sequence (all uppercase) represented as they're commonly represented: 'ACTG'.
  
As long as these two fields are present and accurate, ``cgb_input_generator`` can run and return a meaningful result using default values for parameters.

3. **A list of ``cgb_input_generator`` parameters**
   * Actually a dictionary, where the key is the name of the parameter, and the value is the specified value for that parameter.
4. **A list of ``cgb`` configuration parameters**
   * Also a dictionary, the format of which follows what is outlined above in ``cgb_input_generator`` parameters.
   
See **test_input.md** for a more detailed explanation of the input file for ``cgb_input_generator`` and see **input_skeleton.json** for an example of what your input file should look like. In fact, the user is encouraged to download and use **input_skeleton.json** as a base when preparing their input file.
