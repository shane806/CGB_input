# Input File Overview
  * The ``cgb_input_generator`` input file is specified in JSON format, which defines elements as dictionaries and lists thereof. The JSON file contains four primary entries (`"TF"`, `"motifs"`, `"input_parameters"`, and `"cgb_parameters"`).
## TF
### The `TF` field  specifies the transcription factor family for which the analysis is being performed.
`"TF" : "transcription_factor_name_here"`

## Motifs
  * A list of motifs
  * Each motif is defined as a dictionary with two key:value pairs.

### 1st key:value pair
  * **_Key_** = `"protein_accession"`
  * **_Value_** = `"NP_000000.0"`  
    * A Transcription Factor (TF) protein represented as an NCBI accession.version identifer for proteins
    
    * Example: `[{"protein_accession : "NP.000000.0", ... , "protein_accession" : "Q01234.1"}]`
    		    
### 2nd key:value pair
  * **_Key_** = `"sites"`
  * **_Value_** =  `["ATGCATCG", "ATCGGCTA"]`
    * A list of aligned TF-binding sites reported (or inferred by the user) for the TF under analysis reported in the genome accessions specified for this organism
    
  * Example: `"sites" : ["ATGCATCG", "ATCGGCTA"]`
  
### 1st and 2nd key:value pair together in motifs

```{"motifs" : [
	{ "protein_accession" : "WP_012345.1",
	  "sites" : ["ATGCATCG", 
	  	    "ATCGGCTA"]
	},
	...
	,
	{ "protein_accession" : "NP_987654.9",
	  "sites" : ["ATCGGCTA",
	  	     "GCTACGTA"]
	}   
     ]
  }
```
## `"input_parameters"`
**A dictionary of ``input_parameters``**

1. `"BLAST_eval"`
   * The evalue cutoff for NCBI BLAST search
     * float value between 0 and 1
     * default value: `1e-10`
     * See https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp#expect for more details.
   
2. `"BLAST_dbase"`
   * The database used for NCBI BLAST search
     * database options: `nr`, `refseq_protein`, `landmark`, `swissprot`, `pataa`, `pdb`, `env_nr`, and `tsa_nr`
     * default option: `nr`
     
3. `"max_hits"`
   * The max number of BLAST hits (HSPs) returned in NCBI BLAST search for orthologs
     * integer value between 0 and 500
     * default value: 50
     
4. `"selected_taxon"`
   * The taxonomic level/rank used to sample genomes on a 'one per' basis (e.g. one genome per genus if `"selected_taxon"` is `genus`)
     * level/rank options: `species`, `genus`, `family`, `order`, `class`, `phylum`
     * default option: `null`
5. `"maxID"`
   * The maximum allowed identity between promoter sequences in terms of sequence similarity, expressed as a decimal (float)
     * A float value between 0 and 1
     * default value: `0.8`
     
6. `"up_region"`

   * The number of upstream bases relative to TLS, used in comparing and selectingpromoters below a maximum percent ID relative to other promoters in the sample population
     * A positive integer value between 50 and 300
     * default value: 250
      
7. `"dw_region"`

   * The number of downstream bases relative to TLS, used in comparing and selecting promoters below a maximum percent ID relative to other promoters in the sample population  
     * A positive integer value between 10 and 100
     * default value: 25

8. `"tax_ID"`
   * A taxonomic identifier which, like the aformentioned protein_accession identifiers, is an ID number. `"tax_ID"` is an ID number which represents an 'organism' to limit BLAST hits by that organisms clade
     * A positive integer value 

9. `"TF_family"` (**This parameter must have a value entered by the user! `cgb_input_generator` will return `None` otherwise**) 
   * A name for the TF under study expressed as as string
     * A string value
     * no default

10. `"min_cover"`
   * coverage is the ratio of total number of amino acid residues used in a pairwise alignment to the total number of amino acid residues in the original query sequence. A low coverage value for an HSP has a high probability of being a protein domain which is not the target for `cgb.`
   
   * Therefore, `"min_cover"` is used to minimize the risk of false positive BLAST hits in the form of protein domains.
     * float value between 0.25 and 1
     * default value: `0.75`
     
11. `"sleepy"`
   * Additional sleep time, in seconds, for queries if NCBI returns HTTP error for too many queries
      * A positive float or integer value less than or equal to 3.
      * default value: `0.5'
      
12. `"outputfile_name"`
   * A file path for `cgb_input_generator` to return the results of the generated input file, ending with a .json extension (or one will be applied in its absence), which ultimately will become the filepath for the `cgb_input_generator`
   
     * A string value 
### `cgb_parameters`
4. **A dictionary of ``cgb`` configuration parameters**

#### Site search and operon prediction methods
Several parameters control differenet aspects of CGB behavior regarding how it searches for sites and predicts operons.

- `prior_regulation_probability`
	- Type: `mixed`
	- Allowed values:
		- Positive real number in [0..1] interval
		- "No" string
    - Default value: `"No"`
	- Recommended value: `0.03`
		- Alternatively, any estimate by the user on a reference species, using the number of experimentally reported binding sites as the putative number of regulated operons, and dividing this by the number of known (or predicted) operons in that reference genome.
	- Rationale:
		- This parameter allows the user to define a the prior probability of regulation, which is used for the computation of the posterior probability of regulation in CGB.
		- When "No" is specified, CGB will estimate the prior probability of regulation independently in each genome, based on the information content of the TF-binding motif transferred to that genome as a proxy for the predicted number of regulated operons.

- `alpha`
	- Type: `float`
	- Allowed values:
		- Positive real number in [0..1] interval
    - Default value: `0.00285714285714` [1/300]
    - Rationale:
    	- CGB assumes that, in the promoter region of a regulated gene, the distribution of PSSM scores is a mixture of two distributions: the background score distribution and the distribution of scores in functional sites. The former can be inferred from the genome and the latter for the collection (or more generally, PSWM model) of known (or transferred) sites.
    	- For the mixture distribution `R`, the `alpha` parameter defines the mixing ratio between model (`M`) and background (`B`) distributions, so that `R ~ alpha*M + (1-alpha)*B`.
    	- For a *regular* bacterial transcription factor, we expect 1 site per regulated promoter, and an average promoter length of 300 bp, so the default prior is set to 1/300.

- `promoter_up_distance`
	- Type: `integer`
	- Allowed values:
		- Positive integer under 1,000
    - Default value: `300`
    - Rationale:
    	- This parameter controls how far upstream of a gene's translational start site (TLS) CGB will scan for putative TF-binding sites.

- `promoter_dw_distance`
	- Type: `integer`
	- Allowed values:
		- Positive integer under 1,000
    - Default value: `50`
    - Rationale:
    	- This parameter controls how far downstream of a gene's  translational start site (TLS)  CGB will scan for putative TF-binding sites.

- `phylogenetic_weighting`
	- Type: `boolean`
	- Allowed values:
		- `true / false`
    - Default value: `true`
    - Recommended value: `true`
    - Rationale:
    	- This paramater determines whether phylogenetic weighting is activated. If activated, the genetic distance between species (estimated from a CLUSTALO) alignment of the homologous TF protein sequences in each target and reference species) is used to weigh how much each reference TF-binding motif contributes to the inferred TF-binding model in each target species.

- `site_count_weighting`
	- Type: `boolean`
	- Allowed values:
		- `true / false`
    - Default value: `true`
    - Rationale:
    	- This paramater determines whether site-count weighting is activated. If activated, the number of sites reported in each reference species is used to weigh the contribution to the inferred TF-binding model in target genomes of different reference species (weighing them up contributions from species with many sites, and down for species with few sites). The aim of this weighting is to counteract artifactual changes in the TF-binding motif due to undersampling.

- `posterior_probability_threshold_for_reporting`
 	- Type: `float`
	- Allowed values:
		- Positive real number in [0..1] interval
    - Default value: `0.5`
    - Rationale:
    	- This parameter defines the minimum posterior probability value required to report an operon as regulated by the transcription factor.

- `operon_prediction_probability_threshold`
 	- Type: `float`
	- Allowed values:
		- Positive real number in [0..1] interval
    - Default value: `0.5`
    - Rationale:
    	- This parameter defines the minimum posterior probability value required to split an operon. Predicted operons can be split if one of their "internal" genes appears to be regulated. This parameter determines how high the posterior probability of regulation must be in order to split the operon.

- `operon_prediction_distance_tuning_parameter`
 	- Type: `float`
	- Allowed values:
		- Positive real number in [0.5..5] interval
    - Default value: `1.0`
    - Rationale:
    	- This parameter modulates the maximum intergenic distance used to define operons. CGB computes the mean intergenic distance between the first two genes of all opposing directons. This mean intergenic distance is then considered the maximum distance between the end of a gene and the start of the next in order to consider those genes as belonging to the same operon.
    	- Values for this tuning parameter larger than 1.0 will therefore increase the maximum intergenic distance within operons, allowing more distant genes to belong to the same operon and decreasing the overall number of predicted operons.
    	- Values lower than 1.0 will decrease the maximum intergenic distance, increasing the tendency to split genes into multiple operons and hence increasing the number of predicted operons.

#### Ancestral state reconstruction

Two parameters control the ancestral state reconstruction process.

- `ancestral_state_reconstruction`
	- Type: `boolean`
	- Allowed values:
		- `true / false`
    - Default value: `false`
    - Rationale:
    	- This paramater determines whether ancestral state reconstruction will be performed.

- `bootstrap_replicates`
	- Type: `integer`
	- Allowed values:
		- Positive integer under 10,000
    - Default value: `100`
    - Recommended value: `100` for fast, preliminary analysis; `1,000` for robust ancestral state reconstruction
    - Rationale:
    	- This parameter sets the number of boostrap replicates generated when performing ancestral state reconstruction. Since ancestral state reconstruction operates with discrete traits, for any given gene the inferred probability of regulation in each species is discretized into 1/0 values (or A for absence of ortholog). This is done `bootstrap_replicates` times, and ancestral reconstruction is run on each replicate. The resulting values for ancestral states in each non-terminal node are averaged.

#### Plots and output

Several parameters control the plots and output files generated by CGB.

- `heatmap_plot`
	- Type: `boolean`
	- Allowed values:
		- `true / false`
    - Default value: `true`
    - Rationale:
    	- Generates a heatmap plot, with a species phylogram on the x axis and the ortholog groups on the y axis. The color scheme indicates the average posterior probability of regulation for the gene in each species.

- `motif_plot`
	- Type: `boolean`
	- Allowed values:
		- `true / false`
    - Default value: `true`
    - Rationale:
    	- Plots the inferred TF-binding model in each target species using WebLogo.

- `gene_regulation_plot`
	- Type: `boolean`
	- Allowed values:
		- `true / false`
    - Default value: `false`
    - Rationale:
    	- Plots the reconstructed ancestral states as a phylogram for each orthologous group.

- `taxon_regulation_plot`
	- Type: `boolean`
	- Allowed values:
		- `true / false`
    - Default value: `false`
    - Rationale:
    	- Plots the reconstructed ancestral states for each taxon.

- `network_size_plot`
	- Type: `boolean`
	- Allowed values:
		- `true / false`
    - Default value: `false`
    - Rationale:
    	- Plots a phylogram depicting the size of the network in each taxon.

- `site_printout`
	- Type: `boolean`
	- Allowed values:
		- `true / false`
    - Default value: `true`
    - Rationale:
    	- Determines whether identified putative TF-binding sites for each species will be reported as a CSV file.

#### NCBI Entrez API

- `Entrez email`
	- Type: `string`
	- Allowed values:
		- `any valid email address`
    - Default value: `None` [program will halt if not provided]
    - Rationale:
    	- Entrez API requires that a valid email be provided. This is used to warn users of excessive usage previous to blacklisting. You want to provide an email address in case of overuse; otherwise, if you do not receive the notification and persist, your IP will be blacklisted.

- `Entrez API key`
	- Type: `string`
	- Allowed values:
		- `any valid NCBI Entrez API key`
    - Default value: `None`
    - Rationale:
    	- NCBI allows you to generate an Entrez API key. This allows your scripts to make more requests per second.

- `Entrez API key`
	- Type: `float`
	- Allowed values:
		- `0-10 seconds`
    - Default value: `0`
    - Rationale:
    	- In case of overload, NCBI servers will return a "HTTP Error 429: Too Many Requests" error message. This parameter enables you to add an additional "wait" time in between NCBI API requests (beyond the one already specified in the BioPython library, based on presence or absence of API key), so that you can adjust request intervals to meet server load.

#### NCBI BLAST

- `TF eval`
	- Type: `float`
	- Allowed values:
		- `positive values between 0 and 1`
    - Default value: `0.001`
    - Rationale:
    	- TF_eval allows the user to control the BLAST e-value used to identify instances of the TF factor via tBLASTN using the input reference TF accessions and the specified target genomes.

- `Homolog eval`
	- Type: `float`
	- Allowed values:
		- `positive values between 0 and 1`
    - Default value: `0.001`
    - Rationale:
    	- homolog_eval allows the user to control the BLAST e-value used to identify homologs as reciprocal best blast hits among the specified target genomes.

#### HMMER

- `HMMER eval`
	- Type: `float`
	- Allowed values:
		- `positive values between 0 and 1`
    - Default value: `0.00001`
    - Rationale:
    	- hmmer_eval allows the user to specify the maximum e-value used to identify hits for an orthologous group against the COG, eggNOG or PFAM databases.

- `COG_search`
	- Type: `boolean`
	- Allowed values:
		- `true / false`
    - Default value: `false`
    - Rationale:
    	- COG_search determines whether a HMMER search against a local copy of the COG database will be performed.

- `NOG_search`
	- Type: `boolean`
	- Allowed values:
		- `true / false`
    - Default value: `false`
    - Rationale:
    	- NOG_search determines whether a HMMER search against a local copy of the eggNOG database will be performed.

- `PFAM_search`
	- Type: `boolean`
	- Allowed values:
		- `true / false`
    - Default value: `false`
    - Rationale:
    	- PFAM_search determines whether a HMMER search against a local copy of the eggNOG database will be performed.

- `COG database name`
	- Type: `string`
	- Allowed values:
		- `a valid path including the file name and extension for the COG HMMER database`
    - Default value: `None`
    - Rationale:
    	- COG_dbname specifies the path (including file name) for the COG HMMER database. For instance, "/home/theuser/HHMERdbs/COG_database.hmm". Full path recommended (relative should also work).
    	- 
- `eggNOG database name`
	- Type: `string`
	- Allowed values:
		- `a valid path including the file name and extension for the eggNOG HMMER database`
    - Default value: `None`
    - Rationale:
    	- eggNOG_dbname specifies the path (including file name) for the eggNOG HMMER database. For instance, "/home/theuser/HHMERdbs/bact.hmmer". Full path recommended (relative should also work).

- `PFAM database name`
	- Type: `string`
	- Allowed values:
		- `a valid path including the file name and extension for the eggNOG PFAM database`
    - Default value: `None`
    - Rationale:
    	- PFAM_dbname specifies the path (including file name) for the PFAM HMMER database. For instance, "/home/theuser/HHMERdbs/Pfam-A.hmm". Full path recommended (relative should also work).

- `Orthologous group evalue jump`
	- Type: `float`
	- Allowed values:
		- `positive values`
    - Default value: `5`
    - Rationale:
    	- OGejump indicates how many logs can separate the first and second HMMER hits against the COG/eggNOG database. For instance, when set to 5, if the first hit has 1e-30 e-value, the second hit must have 1e-25 or smaller in order to be reported.

- `COG max hit report`
	- Type: `float`
	- Allowed values:
		- `positive values`
    - Default value: `2`
    - Rationale:
    	- maxCOG indicates how many hits against the COG database are reported, at the most, for any orthologous group.

- `NOG max hit report`
	- Type: `float`
	- Allowed values:
		- `positive values`
    - Default value: `2`
    - Rationale:
    	- maxNOG indicates how many hits against the eggNOG database are reported, at the most, for any orthologous group.

- `PFAM max hit report`
	- Type: `float`
	- Allowed values:
		- `positive values`
    - Default value: `2`
    - Rationale:
    	- maxNOG indicates how many hits against the PFAM database are reported, at the most, for any orthologous group.
