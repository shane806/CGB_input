from Bio import SeqIO, Entrez
import collections, json, time
import NCBI_funcs, GenomeFuncs, ParamValidation, MkLogDirs, TaxFilt, MaxIDFuncs


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

        IGPs_checked = ParamValidation.input_check(inputfile_name)

    except:

        print 'input check fail'

        return None

    log_file["input_parameters"] = IGPs_checked

    user_email = IGPs_checked[0]

    apikey = IGPs_checked[1]

    evalue = IGPs_checked[2]

    blast_db = IGPs_checked[3]

    nhits = IGPs_checked[4]

    selected_taxon = IGPs_checked[5]

    maxID = IGPs_checked[6]

    up_region = IGPs_checked[7]

    dw_region = IGPs_checked[8]

    txID = IGPs_checked[9]

    min_cover = IGPs_checked[10]

    sleepy = IGPs_checked[11]

    TF_family = IGPs_checked[12]

    outputfile = IGPs_checked[13]
    
    cgb_pipeline_parameters = IGPs_checked[14]

    TF_accessions = IGPs_checked[15]

    motifs = IGPs_checked[16]

    Entrez.email = user_email

    Entrez.api_key = apikey
    
    log_file_directory, output_dir = MkLogDirs.make_log_directories(TF_family)
    
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
        name = GenomeFuncs.prettyeze_spname(protein_record.annotations['source'])

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

        blast_hits, blast_log = NCBI_funcs.blast_search(accession, blast_db, 
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

    print len(orthologs), 'hits before processing'

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
        # this function will implement prioritization (best=complete RefSeq)
        # and return the "best" matching genome with a CDS encoding the ortholog

        cds_accession = GenomeFuncs.genome_record_retrieval(ortholog, sleepy,
                                                             log_file_directory)

        # make a list of the genome records for getting taxonomic info later

        # use the function `contig_accessions` to get the complete genome
        # (all chromids) or the complete set of contigs (for a WGS assembly)
        # for the nucleotide record mapping to the ortholog
        if cds_accession:
            
            try:

                print '     Genome: ' + cds_accession['acc'] + ' with score: ' + \
                      str(cds_accession['p_score'])
    
                target_name, target_range = GenomeFuncs.contig_accessions(
                    cds_accession['acc'], cds_accession['p_score'], sleepy)
            
            except Exception as ex:
                
                log_file[ortholog] = 'Error message reads: ' + ex.message\
                    +'...\n' + 'check for issue in genome_record_retrieval'\
                        ' or contig accessions'
            
                toberemoved.append(ortholog)
                
                continue
            
            if target_range == None:

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
    for ortholog in toberemoved:
        
        if ortholog in orthologs:
            
            del orthologs[ortholog]
            
    orth_dir_prefilt = log_file_directory + 'orthologs.json'
    
    with open(orth_dir_prefilt, 'w') as f:
        
        json.dump(orthologs, f, indent=2, separators=(',', ': '))

    print '\nTotal number of valid orthologs detected: ' + str(len(orthologs))

    # This control flow (Beginning here and continuing through the next if/else
    # statement, is how we can sample by both taxonomic level AND max %ID, or
    # just sample by one of those two methods, or neither, and always perform
    # the sampling in the correct order

    # If the user wants to select the best ortholog from each group related to
    # a specific taxonomical level (e.g. get one record for each species,
    # genus, etc.)

    if selected_taxon:

        filtered_orthologs = TaxFilt.filter(selected_taxon, orthologs, sleepy,
                                        log_file_directory)

    # If not sampling by a tax. level, define the orthologs dictionary, as it
    # currently is, as a new dictionary, filtered_orthologs.
    else:

        filtered_orthologs = orthologs

    # if restricting by maxID
    if maxID:

        log_maxID = {}
        
        filtered_list = MaxIDFuncs.maxID_filt(maxID, filtered_orthologs, 
                                               up_region, dw_region, sleepy)

        filtered_orthologs_final = {}

        for ortholog in filtered_orthologs:

            if ortholog in filtered_list:

                filtered_orthologs_final[ortholog] = filtered_orthologs[ortholog]

            else:

                log_maxID[ortholog] = 'maxID filter'
                

        max_ID_log = log_file_directory + 'maxID.json'

        with open(max_ID_log, 'w') as f:

            json.dump(log_maxID, f, indent=2)

    # if all orthologs are wanted, without further study:
    else:

        filtered_orthologs_final = filtered_orthologs

    for ortholog in orthologs:

        if ortholog not in filtered_orthologs_final:

            if ortholog not in log_file:
                log_file[ortholog] = 'other/reason unclear'

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
inputfile_name = '/Users/shanehumphrey/anaconda3/envs/new/GcrA_Alphaproteobacteria_ctrl_input.json'
# =====================================
# ### ENTER YOUR INPUT FILE PATH HERE ###

final = cgb_input_main(inputfile_name)
