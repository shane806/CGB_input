#co-created by: Dr. Ivan Erill and Shane Humphrey
"""
Reads in a JSON file containing at least five required elements:
1. Transcription factor (TF)family name 
2. NCBI protein identifiers for transcription factors to be studied
3. Reported TF-binding sites for the proteins. 
4. input_generator parameters which govern the scope of the BLASTP
   search for orthologs and the sampling process
5. cgb_parameters which govern the scope of cgb
Basic operation:
TASK 1:
  a. Using the provided TF protein accession numbers , cgb_input_generator 
    performs a BLASTP search for each TF, limited by e-value and (optionally) 
    by a taxonomic organism ID. (See input_test.md for details)
  b. If no taxonomic ID is specified, BLAST hits will be obtained up to a
    maximum specified e-value (1e-10 by default). If, however,  a taxonomic ID 
    is specified, the BLASTP search will be constrained to that taxon (e.g.
    class, order...) well as by the specified e-value.
  c. Identified HSPs from BLASTP searches are consolidated into a list of
    unique hits, relative to each other, and to the TFs under study 
TASK 2:
  a. Once HSPs for all TFs have been returned, the next task is to identify 
    and retrieve the 'best' genome record to which the HSP belongs, and collect
    the following coding sequence (CDS) data: 
      1. genome acccession.version ID 
      2. start position
      3. stop position
      4. strand directionality
      5. p_score (defined in genome_record_retrieval() function)
      * 'Best' is defined in detail in the genome_record_retrieval() function,
       but generally, complete RefSeq genomes are most desired.
  b. Once all pieces of information are collected, they are organized into 
     a dictionary of dictionaries, where the main, 'outer', key is the protein
     accession ID of the HSP, and all the genome and CDS information is stored
     in a dictionary as the value for the HSP
     ex.
      { 'ZZ_123456.7' : { 'genome_accession' : 'AA_987654.3', 
                          'start_pos : 1234,
                          'stop_pos : 5678,
                          'strand : '+',
                          'p_score' : 7}}
  c. Following genome retrieval, the Entrez database is then searched for any 
    unidentified genome records such as plasmids for complete records and contig
    records for WGS records
TASK 3:
  a. Following an exhaustive search for all relevant genome records for each 
  HSP, the plasmid or contig record(s) are added to the corresponding HSP under
  a new key called, 'genomes'.
TASK 4:
  a. - The last major task to complete is sampling HSPs
     - There are three sampling methods from which the user may choose, and will
       indicate which of the three method(s) in the input_parameters section
       of the input file. (see test_input.md)
       METHODS
        1. Select only one genome record per taxonomic level (e.g. genus)
        2. Select only records for HSPs with pair-wise promoter region 
           identities below a user-specified threshold (e.g. 90% identity)
           
        3. Select all HSPs identified
TASK 5:
    a. Following the protocol indicated in TASK 4, all relevant information is
      compiled into a JSON file which will be 100% cgb ready, and will still be
       human readable to allow for post-hoc edits.
"""

# Make sure you have: biopython and urllib2 (or urllib3, but still import 
# urllib2 installed in your environment)

# All others are pre-installed if you're using a conda environment
from Bio import SeqIO, Entrez
import collections, json, time
from NCBI_funcs import blast_search
from GenomeFuncs import prettyeze_spname, genome_record_retrieval, contig_accessions
from TaxFilt import filter
from ParamValidation import input_check
from MaxIDFuncs import maxID_filt
from MkLogDirs import make_log_directories

def cgb_input_main(inputfile_name):
    """ This function takes 1 input: A file name of the JSON file as a string.
        (To reduce the liklihood of value or IO errors, it's highly recommended
        that you use the template provided in the github repository')

        Its goal is to create an input file for the cgb pipeline. These are the
        necessary elements which should be in your JSON input file being passed
        into the cgb pipeline input generator:

            1. TF protein accession.version identifiers for every TF under
               analysis.

            2. A list of putative binding sites for each TF protein accession

            3. outputfile (A file path, ending with a .json extension which can
               be directly used to populate CGB upon the  generating the
               relevant set of genomes

            4. user_email (required for NCBI Entrez API access)

            5. NCBI API key (for more queries per second)

            6. min_cover (alignment coverage for BLAST HSPs)

            7. TF_family (a name for the TF under study)

            8. evalue (evalue cutoff for NCBI BLASTP)

            9. txid (taxon 'organism' limit for NCBI BLASTP)

            10. nhits (max number of NCBI BLASTP results)

            11. selected_taxon (the highest taxonomic level to sample)

            12. upstream and downstream positions relative to TLS, for both
             comparing promoters and for CGB

            13. maxID to compare promoters pairwise and retrieve only those with
             less than maxID

            14. Additional sleep time for queries if NCBI returns HTTP error for
                too many queries
    """
    
    log_file = {}

    try:

        IGPs_checked = input_check(inputfile_name)

    except:

        print 'input check fail'

        return None

    log_file["input_parameters"] = IGPs_checked

    user_email = IGPs_checked['user_email']

    apikey = IGPs_checked['apikey']

    evalue = IGPs_checked['evalue']

    blast_db = IGPs_checked['blast_db']

    nhits = IGPs_checked['nhits']

    selected_taxon = IGPs_checked['selected_taxon']

    maxID = IGPs_checked['maxID']

    up_region = IGPs_checked['up_region']

    dw_region = IGPs_checked['dw_region']

    txID = IGPs_checked['tax_ID']

    min_cover = IGPs_checked['min_cover']

    sleepy = IGPs_checked['sleepy']

    TF_family = IGPs_checked['TF_family']

    outputfile = IGPs_checked['outputfile_name']
    
    cgb_pipeline_parameters = IGPs_checked['cgb_parameters']

    TF_accessions = IGPs_checked['TF_accessions']

    motifs = IGPs_checked['motif_object']

    Entrez.email = user_email

    Entrez.api_key = apikey
    
    log_file_directory, output_dir = make_log_directories(TF_family)
    
    outputfile = output_dir + outputfile

    # make a general 'log_files' directory for all log files for the TF

    # output dictionary with all the information included in JSON file
    # what will be changed herein is the motifs and genomes lists

    # we use OrderedDict so that the order in the JSON file is maintained for
    # human readibilty/editability when the JSON output is generated
    output = collections.OrderedDict([("TF", TF_family), ("motifs", []),
                                      ("genomes", [])])

    for key, val in cgb_pipeline_parameters.iteritems():
        
        output[key] = val

    # create dictionary that will store information about the orthologs,
    # which will be the target species for CGB
    orthologs = {}

    # GET TF ACCESSION NUMBERS AND GENOMES
    # first, create a list of TFs, download each IPG record, then look for the
    # specific TF provided by the user and extract its encoding genome
    # Access the reference TF accession numbers (protein), get their IPG records,
    # and create a list of them
    for item in motifs:

        skipprot = False

        # Get genome accessions for the TFs
        for cnt in range(5):

            try:

                records = Entrez.read(Entrez.efetch(db="protein",
                                                    id=item['protein_accession'],
                                                    rettype='ipg',
                                                    retmode='xml'))

                time.sleep(sleepy)

                break

            except:

                print 'NCBI exception raised.\n Reattempt iteration: ' + \
                      str(cnt + 1)

                if cnt == 4:
                    print '\n***Warning: Check reference TF accession: ' + \
                          item['protein_accession']

                    print '\nUnable to retrieve record through Entrez after' \
                          ' 5 attempts'

                    continue

        for prot_rec in records['IPGReport']['ProteinList']:

            if prot_rec.attributes['accver'] == item['protein_accession']:

                TF_genome = prot_rec['CDSList'][0].attributes['accver']

                break

            else:

                print '***Warning: Check the accession: ' + \
                      item['protein_accession']

                print '***Could not access its genome file!'

                print '*** might be missing the version number (e.g. xxxx.1)'

                skipprot = True

        if skipprot: continue

        # GET THE NAME OF THE TF

        # First, download each protein record and create a variable for their
        # respective names
        for cnt in range(5):

            try:

                protein_record = SeqIO.read(Entrez.efetch(db="protein",
                                                          id=item['protein_accession'],
                                                          rettype='gp',
                                                          retmode="text"), 
                                            'genbank')

                time.sleep(sleepy)

                break

            except:

                print 'NCBI exception raised. Reattempt iteration: ' + str(cnt + 1)

                if cnt == 4:
                    log_file[item['protein_accession']] = 'Reference TF fail on ' \
                                                          '"get name" efetch'

                    print 'Could not download record after 5 attempts'

                    skipprot = True
                    
                    continue

        if skipprot: continue

        # get the species name from the annotations "source", which identifies
        # the organism, and prettyeze it
        name = prettyeze_spname(protein_record.annotations['source'])

        # create a compound name using the TF family (e.g. LexA_) and the
        # organism name
        TF_name = TF_family + '_' + name

        print '\nProcessing: ' + TF_name

        # TF_accession (item['protein_accession]), TF_name, TF_genome, and
        # binding sites (item['sites']) are compiled into a dictionary, then
        # added to output for CGB's motif object.
        motif_object = {"name": TF_name, "genome_accession": TF_genome,
                        "protein_accession": item['protein_accession'],
                        "sites": item['sites']}

        output["motifs"].append(motif_object)

    # SECOND LOOP: SEARCH FOR ORTHOLOGS
    # BLAST search is performed using the blast_search function.
    # each hit stored in a dictionary called orthologs as keys with an empty
    # dictionary as its value
    for accession in TF_accessions:
        
        print '\nBLASTing: ' + accession

        blast_hits, blast_log = blast_search(accession, blast_db, 
                                                        evalue, nhits,
                                                        min_cover, txID, sleepy,
                                                        log_file_directory)
        # Ensure no repeating accession numbers in the dictionary
        for hit in blast_hits:
            
            if hit not in orthologs and hit not in TF_accessions:
                
                orthologs[str(hit)] = {}
            
            if hit not in orthologs and hit not in blast_log:
                
                blast_log[hit] = 'blast hit not unique or same as one of the ' \
									'reference TFs'
        
        for log in blast_log:
            
            if log not in log_file:
                
                log_file[log] = blast_log[log]

    print "\nAll reference protein records in input file pre-processed"\
        " and BLASTed."
    
    num_log = collections.OrderedDict()

    orth_start_num = len(orthologs)

    num_log['hits before processing'] = orth_start_num
    print orth_start_num, 'hits before processing'
    

    # THIRD LOOP: work with orthologs (targets)
    # For each ortholog, we want to pick one genome (from the IPG record of that
    # protein) We may also want to restrict things further by imposing that
    # only one ortholog per "clade" (as specified by the user) is pulled into
    # the cgb config file

    toberemoved = []

    for ortholog in orthologs:

        print '\nProt: ' + ortholog

        # use the function `genome_record_retrieval` to get the best selected
        # encoding for the TF ortholog identified by BLASTP
        # this function will implementprioritization (best=complete RefSeq)
        # and return the "best" matching genome with a CDS encoding the ortholog

        cds_accession = genome_record_retrieval(ortholog, sleepy,
                                                             log_file_directory)

        # make a list of the genome records for getting taxonomic info later

        # use the function `contig_accessions` to get the complete genome
        # (all chromids) or the complete set of contigs (for a WGS assembly)
        # for the nucleotide record mapping to the ortholog
        if cds_accession:
            
            try:

                print '     Genome: ' + cds_accession['acc'] + ' with score: ' + \
                      str(cds_accession['p_score'])
    
                target_name, target_range = contig_accessions(
                    cds_accession['acc'], cds_accession['p_score'], sleepy)
            
            except Exception as ex:
                
                log_file[ortholog] = 'Error message reads: ' + ex.message\
                    +'...\n' + 'check for issue in genome_record_retrieval'\
                        ' or contig accessions'
            
                toberemoved.append(ortholog)
                
                continue
            
            if not target_range:

                toberemoved.append(ortholog)

                log_file[ortholog] = 'from contig_accessions '\
                    '(function returned None for target range)'

            # store the genomic cds data, contigs, p_score, and name in
            # orthologs
            orthologs[ortholog]['nuc_record'] = cds_accession
            orthologs[ortholog]['genome'] = target_range
            orthologs[ortholog]['genome_score'] = cds_accession['p_score']
            orthologs[ortholog]['genome_name'] = target_name

        else:

            toberemoved.append(ortholog)

            log_file[ortholog] = 'from genome retrieval'

    # remove orthologs that did not pan out
    gen_ret_cnt, contig_cnt, other_cnt = 0, 0, 0
    for ortholog in toberemoved:
        
        if log_file[ortholog] == 'from genome retrieval':
            
            gen_ret_cnt += 1
        
        elif log_file[ortholog].startswith('from contig_accessions'):
            
            contig_cnt += 1
        
        else:
            
            other_cnt += 1
        
        del orthologs[ortholog]
            
    orth_dir_prefilt = log_file_directory + 'orthologs.json'
    
    post_contigs_cnt = len(orthologs)
    
    num_log['orthologs removed in genome retrieval'] = gen_ret_cnt
    
    num_log['orthologs removed in contig accessions'] = contig_cnt
    
    num_log['orthologs removed in genome retrieval or contig accessions'\
            ' (unclear)'] = other_cnt
    
    with open(orth_dir_prefilt, 'w') as f:
        
        json.dump(orthologs, f, indent=2, separators=(',', ': '))
    
    cnt_lost_pre_TaxFilt = orth_start_num - post_contigs_cnt
    
    print '\nTotal number of valid orthologs detected: ' + str(post_contigs_cnt)

    # This control flow (Beginning here and continuing through the next if/else
    # statement, is how we can sample by both taxonomic level AND max %ID, or
    # just sample by one of those two methods, or neither, and always perform
    # the sampling in the correct order

    # If the user wants to select the best ortholog from each group related to
    # a specific taxonomical level (e.g. get one record for each species,
    # genus, etc.)

    if selected_taxon:
        
        num_log['orthologs lost pre taxon filtering'] = \
            cnt_lost_pre_TaxFilt


        filtered_orthologs = filter(selected_taxon, orthologs, sleepy,
                                        log_file_directory)
        
        cnt_post_TaxFilt = len(filtered_orthologs)
        
        cnt_lost_post_TaxFilt = post_contigs_cnt - cnt_post_TaxFilt
        
        num_log['orthologs lost due to taxon filtering'] = cnt_lost_post_TaxFilt
        
        num_log['orthologs remaining after taxon filtering'] = cnt_post_TaxFilt

    # If not sampling by a tax. level, define the orthologs dictionary, as it
    # currently is, as a new dictionary, filtered_orthologs.
    else:

        filtered_orthologs = orthologs

    # if restricting by maxID
    log_maxID = {}
    if maxID:
        
        cnt_pre_maxID = len(filtered_orthologs)
        
        num_log['orthologs before maxID filtering'] = cnt_pre_maxID
        
        
        filtered_list = maxID_filt(maxID, filtered_orthologs, up_region, 
                                   dw_region, sleepy)

        filtered_orthologs_final = {}
        
        rem_by_maxID = 0
        for ortholog in filtered_orthologs:

            if ortholog in filtered_list:

                filtered_orthologs_final[ortholog] = filtered_orthologs[ortholog]

            else:
                
                rem_by_maxID += 1
                
                log_maxID[ortholog] = 'maxID filter'
                
        num_log['orthologs removed by maxID filter'] = rem_by_maxID

        num_log['orthologs remaining after maxID filter'] = \
            len(filtered_orthologs_final)
            
        max_ID_log = log_file_directory + 'maxID.json'

        with open(max_ID_log, 'w') as f:

            json.dump(log_maxID, f, indent=2)

    # if all orthologs are wanted, without further study:
    else:

        filtered_orthologs_final = filtered_orthologs
    
    other_cnt = 0
    for ortholog in orthologs:

        if ortholog not in filtered_orthologs_final:

            if ortholog not in log_file and ortholog not in log_maxID:
                
                
                
                other_cnt += 1
                
                log_file[ortholog] = 'other/reason unclear'
    
    num_log['other'] = other_cnt
    
    num_dir = log_file_directory + 'num_log.json'
    
    with open(num_dir, 'w') as f:
        
        json.dump(num_log, f, indent=2, separators=(',', ': '))
    
    for ortholog in filtered_orthologs_final:
    
        target = {"name": filtered_orthologs_final[ortholog]['genome_name'],
                  "accession_numbers": filtered_orthologs_final[ortholog]['genome']}

        output["genomes"].append(target)

    print '\nTotal number of valid orthologs after filtering', \
        len(output['genomes'])

    # write out the json file
    print '\nWriting out JSON input file for CGB'

    with open(outputfile, "w") as f:

        json.dump(output, f, indent=2, separators=(',', ': '))

    others_log = log_file_directory + 'others_and_parameters.json'

    with open(others_log, 'w') as f:

        json.dump(log_file, f, indent=2)
    return output

### ENTER YOUR INPUT FILE PATH HERE ###
# =====================================
inputfile_name = '/Users/shanehumphrey/anaconda3/envs/new/pseudomonadales_test.json'
# =====================================
# ### ENTER YOUR INPUT FILE PATH HERE ###

final = cgb_input_main(inputfile_name)
