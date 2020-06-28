import NCBI_funcs
from Bio import Entrez
import time, json

def filter(selected_taxon, orthologs, sleepy, log):
    """
    This function takes in the current dictionary of valid orthologs
    (orthologs), the taxonomic level the user wants to sample by
    (selected_taxon), and the sleep cycle period (sleepy)

    The function loops through every ortholog, fetches its nucleotide record,
    and from that nucleotide record, accesses its taxonomic identifier.

    With the taxonomic identifier, the taxonomy database is queried and the
    lineage of the organism gathered.

    With the gathered lineage of the ortholog, the current dictionary of valid
    orthologs is updated to contain only those whose lineage contains the user-
    specified selected_taxon.

    Once the orthologs dictionary is updated, a new dicitonary of taxonomic
    groups is created which will sort the orthologs by unique taxons.

    (e.g. if selected_taxon is genus, the orthologs will be sorted by
    the different genera in the orthologs dictionary)

    Once sorting is complete, the best genome in each group is selected,
    complete genomes being selected first, or if no complete genomes available,
    the longest WGS record is chosen

    The selected genomes for each group are added to a list and this list is
    the object returned by the function
    """
    tax_log_file = {}
    
    for ortholog in orthologs:
    
        # download the nucleotide record with the genome record just obtained
        # from the genome record retrieval function
        record = NCBI_funcs.try_read_and_efetch('nucleotide',
                                      orthologs[ortholog]['nuc_record']['acc'],
                                      None, sleepy, 'xml')
    
        if not record:
    
            # update tax_tax_log_file and continue to next ortholog
            tax_log_file[ortholog] = 'efetch of nucleotide record failed'
            continue
    
        # Downloading the nucleotide record is done for the purpose of retrieving
        # the taxonomic identifier for the ortholog.
        # The'features' variable is very useful in that it allows for a quick and precise
        # extraction of the feature table object, which is the last object
        # to parse before arriving at the taxonomic identifier
        features = record[0]['GBSeq_feature-table'][0]
    
        # loop through features until a 'key: value' pair which value starts
        # with 'taxon'. Split the at the ':' , and the remaining element is the
        # taxonomic identifier
    
        for quals in features['GBFeature_quals']:
    
            if quals['GBQualifier_value'].split(':')[0] == 'taxon':
    
                tax_id = quals['GBQualifier_value'].split(':')[1]
    
                break
    
        # Download the taxonomy file with the tax_id
        taxon_record = NCBI_funcs.try_read_and_efetch(database='taxonomy',
                                            identifier='txid'+tax_id+'[Organism]',
                                            ret_type='xml', sleepy=sleepy,
                                            ret_mode='xml')
    
        if not taxon_record:
            
            # update tax_log_file and continue to next ortholog
            tax_log_file[ortholog] = "efetch of taxonomy database record" \
                                  " failed"
            continue
    
        # The 'taxon_record object has in it's main dictionary a key, 'Rank',
        # which shows it's most exclusive taxonomic level. We first check the
        # value for this key to see if it's the user's provided 'selected_taxon'.
        # If not, we loop through the 'LineageEx' dictionary, which is also 
        # found in the 'taxon_record' main dictionary. 
        # Each dictionary in 'LineageEx' has 3 keys: 'ScientificName', 'TaxId', and 'Rank'. We're interested in
        # values associated with 'Rank' and 'ScientificName'.
    
        # The 'Rank' key will be a taxonomic level (rank here is synonymous with
        # level) and its value will be one of the 8 putative taxonomic levels,
        # but we're only interested in the levels defined in the list defined
        # above, levels.
        # And if the value for the 'Rank' key is, for example, 'species', then
        # the key, 'ScientificName', in that dictionary, is the name of that
        # species
    
        orthologs[ortholog]['taxonomy'] = None
        
        # get the ScientificName for each taxonomic rank (level)
        if taxon_record[0]['Rank'] == selected_taxon:
    
            selected_name = taxon_record[0]['ScientificName']
            orthologs[ortholog]['taxonomy'] = selected_name
    
        else:
    
            for taxon_rank in taxon_record[0]['LineageEx']:
    
                if taxon_rank['Rank'] == selected_taxon:
    
                    selected_name = taxon_rank['ScientificName']
    
                    orthologs[ortholog]['taxonomy'] = selected_name
                    
                    print 'The ' + selected_taxon + ' for ' + ortholog\
                        + ' is ' + selected_name
        
        if not orthologs[ortholog]['taxonomy']:
            
            print ortholog + ' lacks necessary taxonomic data'
            
            tax_log_file[ortholog] = 'lacks necessary taxonomic data'
        
    taxon_groups = {}
    
    for ortholog in tax_log_file:
        
        if ortholog in orthologs:
            
            del orthologs[ortholog]
    
    # phase_3 starting here: The filtering process
    # organize the remaining orthologs by unique taxonomic groups
    # i.e. if selected_taxon = 'species', organize orthologs by different
    # species
    
    # update orthologs dictionary with the taxonomy
    
    # The following code is where the filtering process actually starts.
    # There will likely be records in orthologs which have the same
    # taxonomic classification, especially if selected_taxon is set to
    # phylum, order, or class.
    # So, we add to taxon_groups all instances of each different taxon at
    # the level of selected_taxon, and repeated instances of the SAME
    # taxonomic classification are added to a list which will be sorted
    # by best genome record score
    for ortholog in orthologs:
        
        if orthologs[ortholog]['taxonomy'] not in taxon_groups:
        
            taxon_groups[orthologs[ortholog]['taxonomy']]\
                = [ortholog]
    
        # these are the repeated instances being added to the list
        else:
    
            taxon_groups[orthologs[ortholog]['taxonomy']]\
                .append(ortholog)
    
    # final list of orthologs after taxonomic filtering
    selected_orthologs = []
    
    # select one genome from each group:
    # choose a complete one if available, or the longest WGS otherwise
    for group in taxon_groups:
        
        print '\nProcessing group: ' + str(group)
    
        potential_selected = []
    
        best_priority_score = 0
    
        for ortholog in taxon_groups[group]:
    
            # first choose best priorization score
            if orthologs[ortholog]['genome_score'] > best_priority_score:
                best_priority_score = orthologs[ortholog]['genome_score']
                
    
                potential_selected = [ortholog]
                
            elif orthologs[ortholog]['genome_score'] == best_priority_score:
                
                potential_selected.append(ortholog)
        
        if len(potential_selected) < 2 or best_priority_score >= 6:

            # if only one ortholog w/ score < 6 or multiple orthologs with 
            # complete genome records for a taxon_group, pick one at random 
            # (the first in the list)

            final_selected = potential_selected[0]
            print '\n-->Best score for: ' + group + ' = ' + \
                str(best_priority_score) + ', from record: ' +\
                    final_selected
                  
            selected_orthologs.append(final_selected)
            
            if best_priority_score >= 6:
                
                for ortholog in potential_selected[1:]:
                    
                    tax_log_file[ortholog] = 'multiple complete genomes for'\
                        + str(group) + '. First genome in list selected, this'\
                        'genome and others in ' + str(group) + ' discarded.'
            
            continue
                
        else:
            
            # if no complete genome available, select the longest WGS record
            # by summing the length, in base pairs, of all nucleotide records
            # associated with each ortholog
            print '\n-->Best scores for: ' + group + ' = ' + \
                str(best_priority_score) + ', from records: '
     
            for ortholog in potential_selected:
                
                print ortholog
    
            max_len = 0
            
            for ortholog in potential_selected:
    
                current_len = 0
    
                # for each genome associated with the ortholog
                print '\nDetermining longest genome record for', ortholog
                
                genome_fail_cnt = 0

                for target_genome in orthologs[ortholog]['genome']:
                    
                    if  genome_fail_cnt <= 10:
                                        
                        try:
                            
                            genome_rec = Entrez.read(Entrez.efetch(
                                db='nuccore',id = target_genome,rettype='docsum',
                                retmode='xml'))
        
                            time.sleep(sleepy)
        
                            record_len = int(genome_rec[0]['Length'])
            
                            current_len += record_len
            
                            if current_len > max_len:
                                
                                max_len = current_len
            
                                final_selected = ortholog
                        except:
                            
                            genome_fail_cnt += 1
                            
                            print '\n|--> Efetch failed; Skipping ' \
                                + target_genome
                            
                            continue  
                    else:
                        
                        print '\n|--> At least 10 WGS Efetch attempts failed'\
                            ' for this ortholog; Removing ortholog from'\
                                ' consideration'
                            
                        tax_log_file[ortholog] = 'WGS contigs Efetch failure'
                        
                        break
    
            print '----> Longest cumulative nucleotide record for ' + \
                  group + ' = ' + str(max_len) + ' from ' + final_selected
    
            selected_orthologs.append(final_selected)
            
            for ortholog in potential_selected:
            
                if ortholog != final_selected:
                    
                    if ortholog not in tax_log_file:
                        
                        tax_log_file[ortholog] = 'shorter genome than longest'\
                            ' genome for ' + group
            
    for ortholog in tax_log_file:
        
        if ortholog in orthologs:
            
            del orthologs[ortholog]
    
    taxon_level_filter_log = log + 'taxonomic_filter_log.json'
    
    with open(taxon_level_filter_log, 'w') as f:
        json.dump(tax_log_file, f, indent=2)

    return orthologs

