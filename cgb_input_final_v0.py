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

from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO, Entrez, pairwise2
import collections
import json
import time
import urllib2

### ENTER YOUR INPUT FILE PATH HERE ###
# =====================================
INPUT_FILEPATH = 'your_file_path.json'
# =====================================
### ENTER YOUR INPUT FILE PATH HERE ###

def try_efetch(database, identifier, ret_mode, sleepy, \
                        ret_type = 'text'):
    """
    Simple function to query entrez database, and safeguard against 
    premature errors returned by NCBI if their server gets overwhelmed.
    """
    for cnt in range(5):
        try:
            record = Entrez.efetch(db= database, id = identifier, \
                                   rettype = ret_type, retmode= ret_mode)
            time.sleep(sleepy)
            break
        
        except:
            print 'NCBI exception raised. Reattempt iteration: ' + str(cnt + 1)
    
    if cnt == 4:
        print 'Could not download record after 5 attempts...'
        print '\nreturning None'
        return None
    
    return record


def try_read_and_efetch(database, identifier, ret_mode, sleepy, \
                        ret_type = 'text'):
    """
    This function will also query the entrez database, and safeguard  
    against premature errors returned by NCBI if their server gets   
    overwhelmed. The difference here is the record will be read 
    within this function, with Entrez.read(), as opposed to being
    read by some other function outside this function.
    """
    for cnt in range(5):
        try:
            record = Entrez.read(Entrez.efetch(db= database, \
                                                id = identifier,\
                                                rettype = ret_type,\
                                                retmode= ret_mode))
            time.sleep(sleepy)
            break
        
        except:
            print 'NCBI exception raised. Reattempt iteration: ' + str(cnt+1)
    
    if cnt == 4:
        print 'Could not download record after 5 attempts...'
        print '\nreturning None'
        return None
    
    return record

def try_read_and_esearch(database, term_val, ID_type, sleepy):
    
    """
    Very similar to the try_read_and_efetch() function defined above,
    except this function calls Entrez with esearch and its parameters
    """
    for cnt in range(5):
        try:
            record = Entrez.read(Entrez.esearch(db= database, \
                                                term = term_val,\
                                                idtype = ID_type))
                                                
            time.sleep(sleepy)
            break
        
        except:
            print 'NCBI exception raised. Reattempt iteration: ' + str(cnt+1)
    
    if cnt == 4:
        print 'Could not download record after 5 attempts...'
        print '\nreturning None'
        return None
    
    return record


def blast_search(TF_accession, cutoff, nhits, min_cover, tax_id, sleepy,\
                 dbase):
    """
    - Remote BLASTp search to detect orthologs.
    Receives a list of protein accession numbers, limiting evalue,      
    max number of HSPs (hits) to be retrieved, a tax_id number that's used to
    constrain the BLASTP search to a database encompassing only the sequences annotated 
    to the taxon identifier via the entrez_query [organism] modifier.
    - Makes remote call to NCBI BLASTP API.
 	If only RefSeq is targeted, then the parameter dbase should be 
 	'refseq_protein' (see https://www.biostars.org/p/129932/)
    - Returns a list containing the protein accessions for the BLASTP hits.
    """

    # gets a protein/AA sequence using the given accessions
    # although this is not strictly necessary (NCBI BLAST can search with 
    # accession), this service often goes down, leading to BLAST returning no 
    # results

    #store all the hits in a list
    handle = try_efetch('protein', TF_accession, 'text', sleepy, \
                        ret_type = 'fasta')
    if handle == None:
        return None

    protrec = SeqIO.read(handle, "fasta")

    protseq = protrec.format('fasta')
 
    # if taxon filtering
    if tax_id!= None:
        taxon = "txid" + str(tax_id) + "[orgn]"
    
        #perform BLASTP search and parse results
        for cnt in range(5):
            try:
                handleresults = NCBIWWW.qblast(program = 'blastp', \
                                       database = dbase, \
                                       sequence = protseq, \
                                       entrez_query = taxon,\
                                       expect = cutoff, hitlist_size = nhits)
                time.sleep(sleepy)
                
                break
            except:
                print 'NCBI exception raised on attempt: ' + str(cnt+1) + \
                    '\nreattempting now...'
                
                if cnt == 4:
                    print 'Could not download record after 5 attempts'
                    return None
    else:
        for cnt in range(5):
            try:
        
                #perform BLASTP search and parse results
                handleresults = NCBIWWW.qblast(program = 'blastp', \
                                       database = dbase, \
                                       sequence = protseq, expect = cutoff,\
                                        hitlist_size = nhits)
                time.sleep(sleepy)
                
                break
            except:
                print 'NCBI exception raised on attempt: ' + str(cnt+1) + \
                    '\nreattempting now...'
                
                if cnt == 4:
                    print 'Could not download record after 5 attempts'
                    return None
    
    blast_records = list(NCBIXML.parse(handleresults))
    
    orthologs = []
    
    # for each blast hit
    for record in blast_records[0].alignments:
        for hsp in record.hsps:
                        
            # if min_cover is a parameter the user specified in the input JSON
            # file
            if min_cover:
                
                # calculate coverage to weed out false positives from protein domains
                query_len = float(record.length) # <- length of TF's complete AA seq.
                align_len = float(hsp.align_length) # <- length of alignment 
                coverage = float(align_len / query_len)
                
                # if min. coverage is less than coverage calculated:
                if coverage >= min_cover:
                    
                    # add the blast hit to orthologs list
                    orthologs.append(record.hit_id.split('|')[-2])
            else:
                orthologs.append(record.hit_id.split('|')[-2])
        
    return orthologs
######################################################################
def split_accession(acc):

    """Splits a genbank/refseq record into prefix (e.g. NZ_MDWE or MDWE) and
       numeric component (e.g. 0100092), and returns them as a list

    """
    #go through the record until we get the first numeric digit

    #store the prefix and use its length to get the numeric region

    prefix = ''
    for c in acc:
        if c.isdigit():
            break
        else:
            prefix=prefix + c
            
    suffix = acc[len(prefix):len(acc)]
    return [prefix, suffix]

######################################################################

def genome_record_to_seq(grecord, upstream, downstream, sleepy):

    """Gets a genome record consisting of accession, start position and strand.
       Queries NCBI to retrieve as many positions upstream and downstream as
       desired from TLS and returns a sequence object.

    """

    #assign positions in record, according to strand
    if  grecord['strand'] == '+':
        s_start = int(grecord['start'])-upstream
        s_stop = int(grecord['start'])+downstream
        s_strand = 1

    else:
        s_stop = int(grecord['stop'])+upstream
        s_start = int(grecord['stop'])-downstream
        s_strand = 2

    #Download FASTA record containing the desired upstream sequence
    for cnt in range(5):
        try:
            net_handle = Entrez.efetch(db="nuccore",id=grecord['acc'], \
                                   strand=s_strand, seq_start=s_start, \
                                   seq_stop=s_stop, rettype='fasta',\
                                   retmode="txt")
            time.sleep(sleepy)
            break
        
        except:
            print 'NCBI exception raised. Reattempt iteration: ' + str(cnt+1)    
    if cnt == 4:
        print 'Could not download record after 5 attempts'
        return None


    gnome_record= SeqIO.read(net_handle, "fasta")   

    time.sleep(sleepy)
    
    return gnome_record


######################################################################

def id_below_maxid_perc(el1, el2, max_percent_id):
    """Aligns two sequence elements and determines whether they are
       more than %ID identical (false) or not (true)
       Scoring: Match:+2, Mismatch:-1, GapO: -2, GapE: -0.2
    """
    
    al=pairwise2.align.globalms(el1.seq, el2.seq, 2, 0, -2, -.5,\
                                one_alignment_only=True, \
                                penalize_end_gaps=False)
    
    matches=0
    gapless=0
    
    #for each position in the alignment
    for ch_pair in zip(al[0][0],al[0][1]):
        
        #if this is a non-gapped position
        if '-' not in ch_pair:
            gapless=gapless+1
            
            #if it's a match, count it
            if ch_pair[0]==ch_pair[1]:
                matches=matches+1
        
    perID = (float(matches)/float(gapless))
    
    # return true or false depending on percent identity
    if perID <= float(max_percent_id):
        return(True)
    else:
        return(False)

def identity_filter_list(orthos, percent_id):
    """Gets a list of orthologs with upstream sequence and a max percent id.
       Goes through the list, removing any records with more than %ID.
       Returns the trimmed list.
    """
    filt_list=[]

    #create a list of protein accession + promoter pairs
    item_list=[]
    
    # orthos = orthologs
    
    # percent_id = 75
    for key, value in orthos.iteritems():

        item_list.append([key,value['promoter'], value['genome_score']])

    #sort item_list so that the first elements are the "best" records
    #(i.e. those with better genome score)
    #this  way, they are likely to be selected if they are not >maxID
    item_list.sort(key = lambda k: k[2], reverse=True)

    cnt=0
    while cnt < len(item_list):
        #get next first element in upstream seq list, removing it from list

        current_item=item_list.pop(0)

        #and add it to the filtered list

        filt_list.append(current_item)

        #check against all remaining elements in item_list; remove them from
        #the list if they are not below threshold of identity
        #at each pass revised item_list hence contains only elements less 
        #than %ID  identical to previously processed elements now stored in 
        #filt_list
        item_list[:]=[upel for upel in item_list if \
                    id_below_maxid_perc(upel[1],current_item[1], percent_id)]
    
    #generate a filtered dictionary
    filt_orthos={}
    
    for item in filt_list:

        filt_orthos[item[0]]=orthos[item[0]]

    return filt_orthos
  
######################################################################

def genome_record_retrieval(ortholog_acc, sleepy):
    """Takes a protein accession as an input. Retrieves its IPG record.
        The idea here is to obtain, prioritarily, data from complete genome 
        records if they exist, from RefSeq (AC_ and NC_ accessions) . If no 
        RefSeq is available, then select complete genome records from GenBank 
        (AE, CP, CY accessions). Otherwise, select contigs or WGS scaffolds from 
        RefSeq (NT_, NW_, NZ_). If that fails, get contigs or WGS scaffolds from 
        GenBank (AAAA-AZZZ). Only when nothing else is available, select direct 
        submissions from GenBank (U, AF, AY, DQ, and other 2 letter prefixes
        identified below).
    """
    # Attempt to download the ortholog's IPG record
    records = try_read_and_efetch('protein', ortholog_acc, 'xml', sleepy, \
                                  ret_type = 'ipg')
    if records == None:
        return None

    #create scoring for priorization for non-complete genomes:
    #genbank, EMBL, & DDBJ nucleotide accession
    # numbers indexed here:
    # 'https://www.ncbi.nlm.nih.gov/Sequin/acc.html'
    # 'http://www.nslc.wustl.edu/elgin/genomics/bio4342/1archives/2006/AccReference.pdf'
    
    # and the list of refseq . numbers:
    #  https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_numbers_and_mole/?report=objectonly
    priority = {'NC_': 7, 'AC_': 7, "AE": 6, "CP": 6, "CY": 6, "NZ_" : 5, \
                "NT_": 5, "NW_": 5, "U": 3, "AF": 3, "AY": 3,"DQ": 3, \
                'EF' : 3, ' EU' : 3, 'FJ' : 3, 'GQ' : 3, 'GU' : 3, \
                'HM' : 3, 'HQ' : 3, 'JF' : 3, 'JN' : 3, 'JQ' : 3, \
                'JX' : 3, 'KC' : 3, 'KF' : 3, ' KJ' : 3, 'KM' : 3, 'KP' : 3,\
                'KR' : 3, ' KT' : 3, 'JX' : 3, 'KC' : 3, 'KF' : 3, ' KJ' : 3, \
                'KM' : 3, 'KP' : 3, 'KR' : 3, ' KT' : 3, 'KU' : 3, 'KX' : 3, \
                'KY' : 3, ' MF' : 3, 'MG' : 3, 'MH' : 3, 'MK' : 3, ' MN' : 3, \
                'MT' : 3}

    # list to hold all genomes and dictionary that is eventually returned as
    # the 'best' option out of all genomes
    genomelist= []
    best_genome = {}
    
    #from the IPG report, retrieve all the genome records from all CDS  
    #listed, keeping only accession number, location of start and stop positions
    if 'ProteinList' in records['IPGReport'].keys():
        for protein_rec in records['IPGReport']['ProteinList']:
            if 'CDSList' in protein_rec.keys():
                for cds in protein_rec['CDSList']:

                    cds_dict = {'acc' : cds.attributes['accver'], 'start' : \
                                cds.attributes['start'], 
                                'stop' : cds.attributes['stop'], \
                                'strand' : cds.attributes['strand'] , \
                                'p_score' : 0}
                        
                    # assign score according to the above defined criteria
                    for key in priority:
                        if cds_dict['acc'].startswith(key):
                            cds_dict['p_score'] = priority[key]
                        if cds_dict['p_score'] > 6.5: 
                            best_genome = cds_dict # Complete refseq
                            return best_genome

                    #special case: NZ_ records that map to "complete" WGS
                    #NZ_ records with a shorter ( <7 ) segment of trailing digits 
                    # e.g. NZ_CP030158,
                    if cds_dict['acc'].startswith('NZ_'):
                        pr,sf = split_accession(cds_dict['acc'])

                        #complete WGS genomes seem to always have < 7 trailing 
                        #digits in their accession, so we pick as "complete" 
                        #anything with less than 7
                        if len(sf.split('.')[0]) < 7:
                            cds_dict['p_score'] = 6.5
                            print '     Atypical NZ_: ' + cds_dict['acc'] + \
                                ' upgraded to RefSeq complete genome.'
                                
                    #finally, score gb WGS ("AAAA-AZZZ" prefixes)
                    elif cds_dict['acc'][0:3].isalpha():
                        cds_dict['p_score'] = 4
                    # add each cds_dict to genomelist
                    genomelist.append(cds_dict)
                    
                # here, we choose a random genome to be our 'best_genome' and
                # compare all other genomes to it. If the p_score of the genome
                # being compared to best_genome is higher than the current 
                # best_genome, make it the best_genome
                best_genome = genomelist[0]
                for genome in genomelist:
                    
                    if genome['p_score'] > best_genome['p_score']:
                        best_genome = genome
                return best_genome
            
            else:
                # There are often many protein records in the ProteinList section
                # of the IPGReport. Most of them will have the cooresponding 
                # CDS data we need. However, more often than not, there's
                # a protein record which has >= 1 identical proteins, but one or
                # more of the cooresponding CDS regions will not exist. 
                #
                # So, to get all CDS regions available, simply continue into the
                # the next protein record. 
                # (EXAMPLE: WP_150851050.1)
                continue
    else:
        print str(ortholog_acc) + " has no  ProteinList feature in its "\
        "IPGReport"
        return None

######################################################################

def prettyeze_spname(spname):
    """Takes a nucleotide record, accesses its 'GBSeq_definition' tag to
       get the species name and then removes characters that might not work
       well with cgb (more appropriately, with programs that cgb calls, such
       as clustalw and BLAST)
    """

    mod_name = spname.split(',')[0].replace(' ','_').replace('.','_').replace('-','_')
    
    return(mod_name)

##############################################################################


def create_master(acc):

    """A protein's IPG record provides information about the contig where the 
    CDS that encodes the protein is located, not the master record for the
    corresponding genome project. 
    
    For WGS genome records, contigs have numerals within a range, and the
    master record is identified with the same prefix but a zero trail.
    
    For instance, the master record NZ_NFJD00000000 encompasses contigs
    NZ_NFJD01000001 to NZ_NFJD01000014.
    
    This function transforms a genome accession into its correspondent master 
    record accession, flipping all its digits into 0"""

    master = ''
    #e.g. "NZ_NFJD01000014" in "NZ_NFJD01000014.1"
    acc_preperiod = acc[:acc.find('.')]

    #get the part after the period in the accession
    #e.g. "NZ_NFJD01000014" in "NZ_NFJD01000014.1"
    acc_postperiod = acc[acc.find('.'):]


    #iterate through accession prefix (e.g. "NZ_NFJD01000014")
    for character in acc_preperiod:

        #if the character is a numeric digit, add a zero to master accession
        if character.isdigit():
            master += '0'

        #otherwise (e.g. NZ_NFJD) add the same character
        else:
            master += character

    #put together the master accession by combining new prefix and suffix
    master = master + acc_postperiod
    
    return master

##############################################################################

def contig_accessions(nuc_record_acc, nuc_record_score, sleepy):
    
    """
        This function receives a nucleotide record accession and its 
        prioritization score.
       
        It then acts differentially depending on the nature of the input genome:
        -   Given a nucleotide WGS master record accession, it identifies
            and returns as a list all the contigs associated with it.
            Details on the genbank nucleotide record structure can be found at:
            https://www.ncbi.nlm.nih.gov/data_specs/schema_alt/NCBI_GBSeq.mod.xsd
    
        -  If the input record is not a master record but a contig record, a 
            master will be created from itself (using create_master function)
    
        -   If the input record is a complete genome record, there is no need to 
            create a master record. However, many bacterial species have 
            multiple chromosomes (and plasmids), and these need to be obtained.
            See for instance: Agrobacterium fabrum str. C58, which has four
            complete sequences (NC_003062.2, NC_003063.2, , NC_003064.2 and 
            NC_003065.3).
    
        -   In such cases (in fact, in all complete sequence cases), the 
            function links out via the Assembly accesion associated to the 
            nucleotide accession obtained and uses Esearch to identify and 
            grab all complete records for that organism.
    """
    
    # determine whether this is a genbank or refseq accession
    #-> genbank accessions do NOT contain '_' characters, refseq accessions do
    # https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.
    #refseq_accession_numbers_and_mole/?report=objectonly
    # http://www.nslc.wustl.edu/elgin/genomics/bio4342/1archives/2006/AccReference.pdf
    # https://www.ncbi.nlm.nih.gov/genbank/wgs/

    refseq = 0
    if '_' in nuc_record_acc:
        refseq=1
    
    #for complete genome records, search other chromosomes or plasmids. 
    #This can be performed retrieving the Assembly database  accession and 
    #using it in ESearch
        
    #Download record information from Nucleotide
    nuc_rec = try_read_and_efetch(database = 'nucleotide', identifier = \
                                  nuc_record_acc, ret_mode = 'xml', \
                                  sleepy = sleepy, ret_type = None)
	
    if nuc_rec == None:
        return None

    sp_name = prettyeze_spname(nuc_rec[0]['GBSeq_definition'])
    
    # if nuc_rec is part of a complete genome record (p_score >= 6), we need
    # to see if there are any plasmids or other unique genome records associated
    # with the organism, so we first attempt to esearch with a BioSample ID,
    # if its in nuc_rec. If not, we search with the BioProject ID. If that's
    # not in nuc_rec, or the search yields zero results, we just use the 
    # record obtained in genome_record_retrieval
    if nuc_record_score >= 6:
        BioSample_found = False
        BioProject_found = False

        for element in nuc_rec[0]["GBSeq_xrefs"]:
            
            if "BioSample" in element["GBXref_dbname"]:
                
                BioSample = element["GBXref_id"] + "[BioSample]"
                BioSample_found = True
            
            else:
                
                if 'BioProject' in element['GBXref_dbname']:
                    BioProject = element['GBXref_id'] + '[BioProject]'
                    BioProject_found = True

        # in case more genome records exist. e.g. plasmids not associated with 
        # the complete circular chromosome
        if BioSample_found:
            Id_List = try_read_and_esearch(database = 'nucleotide', \
                                           term_val = BioSample, \
                                           ID_type = 'acc', \
                                           sleepy = sleepy)['IdList']
                    
        elif BioProject_found:
            
            bioP_id_list = try_read_and_esearch(database = 'BioSample', \
                                                term_val = BioProject, \
                                                ID_type = 'acc', \
                                               sleepy = sleepy)['IdList']
                

            if len(bioP_id_list) > 0:
                for Id in bioP_id_list:
                    
                    BioSample = Id + '[Biosample]'
                    
                    Id_List = try_read_and_esearch(database = 'nucleotide', \
                                           term_val = BioSample, \
                                           ID_type = 'acc', \
                                           sleepy = sleepy)['IdList']
                    			    	  
                    if nuc_record_acc in Id_List:
                        break
                    
                    else:
                        print '*** Warning: attempt to get all genomes from '\
                            'BioProject database for ' + nuc_record_acc + \
                                ' unsuccessful; using only the nucleotide accession'\
                                    ' identified through genome retrieval.'
                        return sp_name, nuc_record_acc
            else:
                print '*** Warning: attempt to get all genomes from '\
                    'BioProject database for ' + nuc_record_acc + \
                        ' unsuccessful; using only the nucleotide accession'\
                            ' identified through genome retrieval.'
                return sp_name, nuc_record_acc
        
        else:
            print '*** Warning: attempt to get all genomes from either '\
                'BioProject or BioSample Databasaes for ' + \
                    nuc_record_acc + ' proved unsuccessful; using only the '\
                        'nucleotide accession identified through genome '\
                            'retrieval.'
                            
            # return only the accession, as it was given
            return sp_name, nuc_record_acc
        
        # Check genome records returned in Id_list. If the original genome record
        # is a refseq record, grab only refseq records. If the original 
        # genome record is a genbank record, grab only genbank records.
    	contigs=[]
        for nuc_rec in Id_List:
            
            if refseq == True and '_' in nuc_rec:
                contigs.append(nuc_rec)
                return sp_name, contigs
                
            elif '_' not in nuc_rec:
                contigs.append(nuc_rec)
                return sp_name, contigs        
    
    
    #for WGS records:
    elif nuc_record_score <=5:
    #if the record is neither a primary record nor a master record, then 
    #it must be a contig record --> a master record will be created using 
    #the `create_master` function

        if "GBSeq_alt-seq" not in nuc_rec[0].keys():
            print '             ' + nuc_record_acc + ' is not a master record'
    
            #generate a XXX0000000-type master record accession from contig acc
            master_record_acc = create_master(nuc_record_acc)
    
            print '             --> Master record accession created: ' \
                    + master_record_acc
    
            #in principle, records missing "GBSeq_alt-seq" should be WGS 
            #contigs (with an associated master), but not all are.
            #so first we check whether the generated master record can be 
            #obtained; 
            #grab the nucleotide associated with the master record
            for cnt in range (5):
                try:
                    handle=Entrez.efetch(db = "nucleotide", \
                                         id = master_record_acc,\
                                     retmode = "xml")
                    time.sleep(sleepy)
                    break
                
                except:
                    if cnt == 4:
                        print '             --> Master record ' + \
                            master_record_acc + ' did not pan out'
                
                    #if we cannot get to the "master" we created
                    #get species name, replacing spaces by underscores
                    sp_name = prettyeze_spname(nuc_rec[0]['GBSeq_definition'])
            
                    #return only the accession, as None [to discard it]
                    return sp_name, None
    
            #if the master record accssion worked out, grab it
            master_rec = Entrez.read(handle)
    
        #if it was already a master record 
        #assign the previously downloaded nucleotide info    
        else:
            master_rec = nuc_rec
        
        #Double check that this is indeed a master record
        #("GBSeq_alt-seq" must be a key of the record)
        if ("GBSeq_alt-seq" in master_rec[0].keys()):

        #ranges deals with the fact that some master records do not hold 
        #a contiguous set of contigs
        #for instance: the WGS mater record NZ_AUIV00000000.1 encompasses:
        # WGS_SCAFLD  NZ_AUIV01000001-NZ_AUIV01000002

        # WGS_SCAFLD  NZ_AUIV01000005-NZ_AUIV01000015

        # WGS_SCAFLD  NZ_AUIV01000020-NZ_AUIV01000026

        # WGS_SCAFLD  NZ_AUIV01000029-NZ_AUIV01000057

        #each interval is kept as an element in ranges
            ranges = []

        #process the record depending on whether it is genbank or refseq
            if refseq:
    
                #get first and last WGS_SCAFLD contig records listed in master
                #For non-contiguous records (example: NZ_AUIV00000000.1): 
                #Each interval is presented as an element of the ranges list. 
                #Then, several elements will contain the WGS_SCAFLD tag. So,
                #we iterate over the list and add those intervals to ranges
                for item in master_rec[0]["GBSeq_alt-seq"]:
                    if item['GBAltSeqData_name']=='WGS_SCAFLD':
    
                        #if this range contains more than one item
                        if "GBAltSeqItem_last-accn" in \
                            item["GBAltSeqData_items"][0].keys():
                            first = item["GBAltSeqData_items"][0]\
                                        ["GBAltSeqItem_first-accn"]
    
                            last = item["GBAltSeqData_items"][0]\
                                        ["GBAltSeqItem_last-accn"]
                        else:
    
                            first = item["GBAltSeqData_items"][0]\
                                        ["GBAltSeqItem_first-accn"]
                            last = first
                        ranges.append([first, last])
    
            #if genbank
            else:
                #get first and last WGS contig records listed in master rec
                #the same procedure as for RefSeq, but looking for 'WGS' tags
                
                for item in master_rec[0]["GBSeq_alt-seq"]:
                    if item['GBAltSeqData_name']=='WGS':
    
                        #if this range contains more than one item
                        if "GBAltSeqItem_last-accn" in \
                            item["GBAltSeqData_items"][0].keys():
    
                            first = str(item["GBAltSeqData_items"][0]\
                                        ["GBAltSeqItem_first-accn"])
                                        
                            last = str(item["GBAltSeqData_items"][0]\
                                        ["GBAltSeqItem_last-accn"])
                        else:
                            first = str(item["GBAltSeqData_items"][0]\
                                        ["GBAltSeqItem_first-accn"])
    
                            last = first
    
                        ranges.append([first, last])
    
            #get the list of contigs associated to the WGS record
            contig_list = []
    
            #for each interval in range
            for interval in ranges:
                #get the first and last accessions
    
                first = interval[0]
                last = interval[1]
    
                #get the prefix and the numbers for start and end 
                prefix, first_num = split_accession(first)
                prefix, last_num = split_accession(last)
                suffix_len = len(first_num)
    
                #for each possible number between first and last, generate the
                #corresponding accession
                for number in range(int(first_num),int(last_num)+1):
    
                    #get number from first to last, transform number to string,
                    #look at size difference between number and the "normal"
                    #suffix (numeric) length for these records
                    #pad with zeroes on the left to get the completed string
                    #corresponding to the numeric part
                    num_str = str(number)
                    
                    padding = '0' * (suffix_len-len(num_str))
    
                    num_str=padding+num_str
    
                    #create the composite accession (prefix+suffix) and add 
                    #to contig list
                    contig_list.append(prefix + num_str)
    
            #get species name, replacing spaces with underscores
            sp_name=prettyeze_spname(master_rec[0]['GBSeq_definition'])
    
            return sp_name, contig_list

    	else:
            print '             --> Master record does not contain adequate fields.'

            #if the master record that we generated does not work well
            #get species name, replacing spaces by underscores

            sp_name = prettyeze_spname(nuc_rec[0]['GBSeq_definition'])

            #return only the accession, as None [to discard it]
            return sp_name, None

# contig_accessions('ADXF01000600.1', 4, .5)

##############################################################################


def tax_level_filter(selected_taxon, orthologs, sleepy):
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
    
    remove_list = []
    for ortholog in orthologs:
        
        # download the nucleotide record with the genome record just obtained
        # from the genome record retrieval function
        record = try_read_and_efetch(database = 'nucleotide', \
                                     identifier = \
                                     orthologs[ortholog]['nuc_record']['acc'],\
                                     ret_mode = 'xml', sleepy = sleepy, \
                                     ret_type = None)
            
        if record == None:
	
            # update remove_list and continue to next ortholog
            remove_list.append(ortholog)
            continue
        
        # Downloading the nucleotide record is done for the purpose of retrieving
        # the taxonomic identifier for the ortholog.         
        # features is very useful in that it allows for a quick and precise 
        # extraction of the feature table object, which is the last object
        # to parse before arriving at the taxonomic identifier
        features = record[0]['GBSeq_feature-table'][0]['GBFeature_quals']
        
        # loop through features until a 'key: value' pair which value starts 
        # with 'taxon'. Split the at the ':' , and the remaining element is the 
        # taxonomic identifier
        for quals in features:
            if quals['GBQualifier_value'].split(':')[0] == 'taxon':
                
                tax_id = quals['GBQualifier_value'].split(':')[1]
    
        # Download the taxonomy file with the tax_id
        taxonomy = try_read_and_efetch(database = 'taxonomy', \
                                       identifier = 'txid'+ tax_id + '[Organism]', \
                                       ret_mode = 'xml', sleepy = sleepy,\
                                       ret_type = 'xml')
        if taxonomy == None:
            
            # update remove_list and continue to next ortholog
            remove_list.append(ortholog)
            continue
        
        # LineageEx, defined below, is a list of dictionaries. Each dictionary 
        # has 3 keys: 'ScientificName', 'TaxId', and 'Rank'. We're interested in
        # values associated with 'Rank' and 'ScientificName'.
        
        # The 'Rank' key will be a taxonomic level (rank here is synonymous with
        # level) and its value will be one of the 8 putative taxonomic levels,
        # but we're only interested in the levels defined in the list defined
        # above, levels.
        # And if the value for the 'Rank' key is, for example, 'species', then
        # the key, 'ScientificName', in that dictionary, is the name of that 
        # species
        
        ranks_list = taxonomy[0]['LineageEx']
        
        # get the ScientificName for each taxonomic rank (level)
        tax_ranks = {}
        for ranks in ranks_list:
            
        # the taxonomic levels we're interested in are:
        # species, genus, family, order, class, and phylum
                if ranks['Rank'] == 'species':
                    tax_ranks['species'] = str(ranks['ScientificName'])
                    
                elif ranks['Rank'] == 'genus':
                    tax_ranks['genus'] = str(ranks['ScientificName'])
                    
                elif ranks['Rank'] == 'family':
                    tax_ranks['family'] = str(ranks['ScientificName'])
                    
                elif ranks['Rank'] == 'order':
                    tax_ranks['order'] = str(ranks['ScientificName'])
                    
                elif ranks['Rank'] == 'class':
                    tax_ranks['class'] = str(ranks['ScientificName'])
                    
                elif ranks['Rank'] == 'phylum':
                    tax_ranks['phylum'] = str(ranks['ScientificName'])
        
        
        print "\nLineage of " + orthologs[ortholog]['genome_name'] + ":"
        for level, name in tax_ranks.iteritems():
            print '\n ---->'+  level.upper() + ": " + name
        
        # add orthologs to remove_list which don't contain the required taxonomic
        # data
        if selected_taxon not in tax_ranks:
            print ortholog, 'lacks necessary taxonomic data and will be '\
                'removed from analysis.'
            remove_list.append(ortholog)
   
    # update orthologs accordingly
    for ortholog in remove_list:
        del orthologs[ortholog]
        
    # phase_3 starting here: The sampling
    # organize the remaining orthologs by unique taxonomic groups
    # i.e. if selected_taxon = 'species', organize orthologs by different
    # species
    taxon_groups = {}
    
    # update orthologs dictionary with the taxonomy
    for ortholog in orthologs:
        orthologs[ortholog]['taxonomy'] = tax_ranks
            
        # The following code is where sampling actually starts
        # there will likely be records in orthologs which have the same 
        # taxonomic rank, especially if selected_taxon is set to phylum,
        # order, or class.
        # so, we add to taxon_groups all instances of each different 
        # selected_taxon, and repeated instances of the SAME selected_taxon are
        # all added to a list which will be sorted by best genome record score
        if orthologs[ortholog]['taxonomy'][selected_taxon] not in taxon_groups:
            taxon_groups[orthologs[ortholog]['taxonomy'][selected_taxon]]\
                = [ortholog]
         
        # these are the repeated instances being added to the list 
        else:
            taxon_groups[orthologs[ortholog]['taxonomy']\
                         [selected_taxon]].append(ortholog)
    
    # final list of orthologs after taxonomic filtering         
    selected_orthologs = []
    
    #select one genome from each group:
    #choose a complete one if available, or the longest WGS otherwise
    for group in taxon_groups:
        print '\nProcessing group: ' + str(group)
        selected = None
        
        bestpriority = 0
        for ortholog in taxon_groups[group]:
            
            #first choose best priorization score
            if orthologs[ortholog]['genome_score'] > bestpriority:
                
                bestpriority = orthologs[ortholog]['genome_score']
                
                selected = ortholog
        
        print '\n-->Best score for: ' + group + ' = ' + str(bestpriority) + \
              ', from record: ' + selected
        
        #if no complete genome available, select the longest WGS record
        #by summing the length, in base pairs, of all nucleotide records 
        # associated with each ortholog
        if bestpriority < 6:            
            max_len = 0            
            for ortholog in taxon_groups[group]:
                
                current_len = 0
                
                #for each genome associated with the ortholog
                print '\nDetermining longest genome record for', ortholog
                for target_genome in orthologs[ortholog]['genome']:
                    
                    rec = try_read_and_efetch('nucleotide', target_genome, \
                                              'xml', sleepy, ret_type= None)
                        
                    if rec == None:
			
                        print '\n|--> Efetch failed; Skipping ' + target_genome
                        continue
                    
                    record_len = rec[0]['Length']
                    current_len += int(record_len)
                    time.sleep(sleepy)
                    
                    if current_len > max_len:
                        max_len  = current_len
                        selected = ortholog
            
            print '----> Longest cumulative nucleotide record for ' + \
                  group + ' = ' + str(max_len) + ' from ' + selected
                
        selected_orthologs.append(selected)
        
    return selected_orthologs

##############################################################################
def maxID_filt(maxID, prefilt_orthologs, up_region, dw_region, sleepy):

    print 'Obtaining upstream nucleotide sequences for homologs'

    for ortholog in prefilt_orthologs:
        print '  |-> Grabbing upstream nucleotide sequence for: ' + ortholog
        
        #grab the sequence corresponding to that record
        seq = genome_record_to_seq(prefilt_orthologs[ortholog]['nuc_record'], \
                                  up_region, dw_region, sleepy)
        
        prefilt_orthologs[ortholog]['promoter'] = seq
    
    # filter the list of upstream sequences
    filtered_orthologs = identity_filter_list(prefilt_orthologs, maxID)
    
    return filtered_orthologs

##############################################################################

def input_check(inputfile_name):
    """ This function takes 1 input: A file name of the JSON file as a string.
        It will be the first function called when the main cgb_input function
        is called.
        
        The goal is to validate each of the 16 user-provided parameters
        involved in generating the genome accession numbers of the orthologs 
        that ultimately are passed into the next phase of the CGB pipeline. 
        A few key points about this function:
        
            - If an improper/incompatible value(s) is(are) given. The function
              will prompt the user an adhoc update to that specific parameter's
              value.
              
            - If a proper/compatible value can't be obtained within 3 or less
              prompts, a default value is defined for each parameter and will be
              used in place of a user-provided value.
              
             - If null (null in JSON == None in python) is the user-provided 
             value for any given parameter, a predetermined default value will
             be set in it's place
              
        STRUCTURE:
        The JSON file contains 1 main JSON object (python dictionary)
        containing 4  K:V pairs):
                
            Main JSON object: 
                
            1st K:V pair - The name of the TF family being analyzed 
                             as "TF : "TF family name"
                             
            2nd K:V pair - "motifs" : JSON array (python list) of JSON objects
            
                * Each JSON object in the array has 2 K:V pairs
                            
                "motifs" : [{"protein_accession" : "accession_number_1",
                             "sites" : ["SITE ONE", "SITE TWO", ..., "SITE N"]},
            
                            {"protein_accession" : "accession_number_2", 
                             "sites" : ["SITE ONE", "SITE TWO", ..., "SITE N"]}, 
                        
                            ...,
                        
                            {"protein_accession" : "accession_number_N", 
                             "sites" : ["SITE_ONE", "SITE_TWO", ..., "SITE_N"]}
                            
            3rd K:V pair -
            
                "input_parameters" : JSON object with 14 unordered K:V pairs
                        
                    * Each K:V pair in "input_parameters" is a parameter 
                     used in generating the genome files which
                      the cgb pipeline ultimately takes as input.
                          
            4th K:V pair - 
                
                "cgb_parameters" : JSON object with 33 unordered K:V pairs
                
                    * Each K:V pair within "cgb_pipeline_parameters" represents 
                      a parameter used directly in the CGB pipeline
					  
    """        
    # Initialize list which will ultimately be the returned object containing
    # all validated parameters for the cgb_input function
    IGPs_checked = []
    
    ### *** JSON input file IO *** ###
    try:
        with open(inputfile_name, 'r') as path:
            in_file = json.load(path)

    # If unsuccessful, grab the raised exception and initialize in new variable  
    # to be called later
    except Exception as ex:
        except_type = type(ex)

        # Create a general template to show why an error was raised
        temp = "A(n) {} was raised. \nArguments: {!r}"
        print temp.format(type(ex).__name__, ex.args)
        
        
        if except_type == ValueError:
            print "\n Your JSON file isn't structured properly somewhere. This"\
                    " type of error can often be tough to pinpoint."\
                    "\n A potentially time saving solution to this issue"\
                    " is to copy/paste the contents of your JSON file into a web"\
                    " based JSON formatter and/or validator such as:"\
                    "\n https://jsonformatter.curiousconcept.com/ "\
                    "\n to quickly identify what syntax error is"\
                    " preventing a successful read in."
            return None
                    
        elif except_type == IOError:
            print "\n Remember, the file name passed into the script must include "\
            " the .json extension. For example, If your file's name is"\
            " 'input_f','input_f.json' would be the correct syntax used "\
            " (assuming input_f is a JSON file, a necessary condition of the input "\
            " file in order to be compatible with this program, and located in "\
            " your current working directory), including the quotation marks, "\
            " for the file to be properly recognized and read in by python." \
            
            print "\n If the file is not in your current directory, you " \
                " must use a relative or an absolute path." 
            
            print "\n e.g. ../currentworkingdirectory/myjsonfile.json or"\
                " home/User/otherdirectory/myjsonfile.json"
            return None
    
    CGB_parameters = in_file['cgb_parameters']
        
    print '\nWorking with parameters as follows:'
    
    ### *** internet connection test *** ###
    try:
        urllib2.urlopen("https://www.google.com/", timeout = 1)
    except:
        print "It is necessary to have a working internet connection to use"\
            " this program. Connect to an available active network and restart"\
                " the program."
        return None

    ### *** NCBI account email address check *** ###
    if 'entrez_email' in CGB_parameters:
        user_email = str(CGB_parameters['entrez_email'])
        if "@" in user_email and "." in user_email:
                user_email_checked = user_email
                print "\n -NCBI account email address:", user_email_checked
                IGPs_checked.append(user_email_checked)
            
        else:
            print "The email address provided under 'entrez_email' within"\
                " the 'cgb_parameters' object of your JSON file is invalid. Please"\
                    " update that value to the email address used to register "\
                        "for an NCBI account at www.ncbi.nlm.nih.gov"
            return None
    else:
        print "It's necessary for a key, 'entrez_email' to be in the "\
            "CGB_parameters dictionary of your JSON file, and for that email to"\
                " contain a valid NCBI account email address."
        print "\nUpdate your JSON file with this key and cooresponding value."
        return None
    

    ### *** NCBI account api key check *** ###
    if 'entrez_apikey' in CGB_parameters:
        apikey = str(CGB_parameters['entrez_apikey'])
    
        # use urllib2 module to access the entrez database with the given
        # apikey. If the api key is valid, the request will go through, otherwise,
        # an excepytion will be raised.
        api_url_prefix = "https://www.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?api_key="
        
        try:
            urllib2.urlopen(api_url_prefix + apikey)
        except Exception as badAPI:
           temp = "A(n) {} was raised. \nArguments: {!r}"
           print temp.format(type(badAPI).__name__, badAPI.args)
           print "Update the entrez_apikey parameter in your JSON file to a valid"\
               " api key."
           return None
    
        apikey_checked = apikey
        print "\n -NCBI account api key:", apikey_checked
    else:
        print "Non NCBI account api key specified. Setting parameter to none."
        apikey_checked = None
        
    IGPs_checked.append(apikey_checked)
    
    IGPs = in_file['input_parameters']

    ### *** BLAST_eval check *** ###
    
    # isinstance() checks the type and returns a boolean
    if 'BLAST_eval' in IGPs:
        if isinstance(IGPs['BLAST_eval'], float) and IGPs['BLAST_eval'] != None:
            BLAST_eval = IGPs['BLAST_eval']
            
            # if the e-value is pos. and less than 1, validate it
            if BLAST_eval >= 0 and BLAST_eval <= 1:
                BLAST_eval_checked = BLAST_eval
                print "\n -BLAST limiting e-value set to: ", BLAST_eval_checked
            
            # if not, return None
            else:
                print '\nThe "BLAST_eval" parameter reflects "the statistical '\
                    'significance threshold for reporting matches against database'\
                        ' sequences." For the purposes of this program, the '\
                            '"BLAST_eval" parameter should be a positive number '\
                                'less than 1.'
                                
                print "\nSee https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD="\
                    "Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp#expect for more "\
                        "details."
                return None
        
        # if user wants default
        elif IGPs['BLAST_eval'] == None:
            BLAST_eval_checked = 1e-10 # default
            print "\n -BLAST_eval parameter in your JSON set to null (default.)"
            print "\n -Setting parameter to default: 1e-10"
        
        # if not the right type, return None
        else:
            print '\nThe "BLAST_eval" parameter reflects "the statistical '\
                'significance threshold for reporting matches against database'\
                    ' sequences." For the purposes of this program, the '\
                        '"BLAST_eval" parameter should be a positive number '\
                            'less than 1.'
                            
            print "\nSee https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD="\
                "Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp#expect for more "\
                    "details."
            return None
    else:
        BLAST_eval_checked = 1e-10 # default
        print " -No BLAST_eval parameter in your JSON file. Setting parameter"\
            " to default: 1e-10"
    IGPs_checked.append(BLAST_eval_checked)
    
    ### *** BLAST_dbase check *** ###
    if 'BLAST_dbase' in IGPs:
        if isinstance(IGPs['BLAST_dbase'], unicode) and IGPs['BLAST_dbase'] != None:
            BLAST_dbase = str(IGPs['BLAST_dbase'])
            
            # database options
            databases = ['nr', 'refseq_protein', 'landmark', 'swissprot', \
                         'pataa', 'pdb', 'env_nr', 'tsa_nr']
                
            # if BLAST_dbase given is one of the databases available, validate it
            if BLAST_dbase in databases:
                BLAST_dbase_checked = BLAST_dbase
                print "\n -Database to be queried in search of orthologs:"\
                    , BLAST_dbase_checked
                    
            # if not, print the available options and return None
            else:
                print "The 'BLAST_dbase' parameter reflects the NCBI protein "\
                    "database which will be used to find orthologs. Here's a list"\
                        " of the available databases."
                for database in databases:
                    print "\n -" + database
                print "\nUpdate the 'BLAST_dbase' parameter in your JSON file with "\
                    "one of the above databases, or null to set to default."
                return None
            
        # if user wants default
        elif IGPs['BLAST_dbase'] == None:
            BLAST_dbase_checked = 'nr' # default
            print "\n -The default BLAST database, nr, will be used."
        
        # if parameter not the right type
        else: 
            print "The 'BLAST_dbase' parameter reflects the NCBI protein "\
                "database which will be used to find orthologs. Here's a list"\
                    " of the available databases."
                    
            for database in databases:
                print " \n" + database
                
            print "\nUpdate the 'BLAST_dbase' parameter in your JSON file with "\
                "one of the above databases."
            return None
    else:
        BLAST_dbase_checked = 'nr' # default
        print "No database specified. Setting BLAST_dbase parameter to "\
            "default: nr"
    
    IGPs_checked.append(BLAST_dbase_checked)
        
    ### *** BLAST max_hits check ***
    if 'max_hits' in IGPs:
        nhits = IGPs['max_hits']
        if nhits != None and isinstance(nhits, int):
            
            # if positive and less than 500, validate it
            if nhits > 0 and nhits < 500:
                nhits_checked = nhits
                print "\n -Max number of hits to be returned in BLAST search:"\
                    , nhits_checked
            else:
                print 'This parameter limits the max number of BLAST hits returned, so'\
                    ' a positive integer is the only compatible type for this parameter.'\
                        ' Update the value for max_hits to a whole number between 1 and 100.'
                return None
        
        elif nhits == None:
            nhits_checked = 50
            print "\N -Max number of hits to be returned in BLAST search set to "\
                "default: 50"
            
        else:
            print 'This parameter limits the max number of BLAST hits returned, so'\
                ' a positive integer is the only compatible type for this parameter.'\
                    ' Update the value for max_hits to a whole number between 1 and 100.'
            return None
    else:
        nhits_checked = 50
        print "no max_hits parameter in JSON file. Setting parameter to default"\
            , nhits_checked

    IGPs_checked.append(nhits_checked)
    
    ### *** selected_taxon check *** ###    
    if 'selected_taxon' in IGPs:

        if IGPs['selected_taxon'] == None:
            selected_taxon_checked = None # default
            print "\n -Null (default) value for selected_taxon parameter given"\
                " in JSON."
            print "\n  --Default setting will be used: None"
                    
        # Define the taxon options available
            
        #  if selected_taxon is one of the options, validate it. Or if 'null'
        # given as value in JSON, validate default
        else:
            taxon_levels = ['species', 'genus', 'family', 'order','class', 'phylum']

            if IGPs['selected_taxon'] in taxon_levels:
                selected_taxon = str(IGPs['selected_taxon'])
                selected_taxon_checked = selected_taxon

                print '\n -Sampling results by taxonomic level: '\
                    + selected_taxon_checked
                    
            
        
            
        # if not validated, print the available options and return None
            else:
                print "\nThe taxonomic level indicated in selected_taxon is "\
                    "incompatible with this program."
                print "\nPlease update your JSON with one of the "\
                    "following:"
                            
                for taxon in taxon_levels:
                    print " -" + taxon
                return None
            
    # if parameter not in IGPs, set to default
    else:
        print '\n -Orthologs will not be sampled according to a taxonomic '\
            'level.'
        selected_taxon_checked = None # default

    IGPs_checked.append(selected_taxon_checked)

	        
        # if ID_filter was set to True but None given for maxID, set to default
    if'maxID' in IGPs:
        if isinstance(IGPs['maxID'], (float)):
            maxID = float(IGPs['maxID'])

            # if maxID a pos. number between 0 and 1, validate it, otherwise,
            # return none
            if maxID > 0 and maxID < 1:            
                maxID_checked = maxID
                print "\n -Orthologs used in analysis will be less than "\
                    + str(maxID_checked*100) + "% identical to TFs"
                    
                IGPs_checked.append(maxID_checked)

            else:
                print 'The maxID parameter must be a float number between zero '\
                        'and one...'
                print "\nUpdate your JSON with a value for maxID within these"\
                    " constraints."
                return None

        
        elif IGPs['maxID'] == None:
            print '\n -maxID parameter in your JSON file set to null.'
            print "\n -Setting parameter to default: 0.8"
            maxID_checked = 0.8 # default

            IGPs_checked.append(maxID_checked)
            
        else:
            
                print 'The maxID parameter must be a float number between zero '\
                    'and one...'
                print "\nUpdate your JSON with a value for maxID within these"\
                    " constraints."
                return None
    else:
        print "\n -No maxID parameter in your JSON file."
        print "\n -Setting parameter to default: None"
        maxID_checked = 0.8 # default
        
        IGPs_checked.append(maxID_checked)
        
    ### *** up_region and dw_region checks *** ###
    
    if 'up_region' in IGPs:
            
        if isinstance(IGPs['up_region'], int):            
            up_region = IGPs['up_region']
            
            # if upregion between 25 and 1000, validate it, otherwise, return None
            if up_region >= 25 and up_region <= 1000:
                up_region_checked = up_region
                print "\n -Number of bases to check in the upstream region: ",\
                    up_region_checked
            else:
                print "the 'up_region' parameter takes only positive integer "\
                    "values from 25 to 1000. Update the value of 'up_region' "\
                        "in your JSON file such that it falls within these "\
                            "constraints"
                return None
        
        # if up_region set to None, validate it as default
        elif IGPs['up_region'] == None:
            print "\n -up_region parameter set to null in your JSON."
            print "\n -Setting parameter to default: 250"
            up_region_checked = 250
        
        # if up_region not the right type, return None
        else:
            print "the 'up_region' parameter takes only positive integer values from"\
                " 25 to 1000. Update the value of 'up_region' in your JSON file such"\
                    " that it falls within these constraints"
            return None
    else:
        print "\n -No up_region parameter in your JSON file."
        print "\n -Setting parameter to default: 250"
        up_region_checked = 250

    IGPs_checked.append(up_region_checked)

    if 'dw_region' in IGPs:
        if isinstance(IGPs['dw_region'], int):                   
            dw_region = IGPs['dw_region']
            
            # if dw_region between 10 and 100, validate it, otherwise, return None        
            if dw_region >= 10 and dw_region <= 100:
                dw_region_checked = dw_region
                print "\n -Number of bases to check in the downstream region: ",\
                    dw_region_checked
    
            else:
                print "The 'dw_region' parameter takes only positive integer "\
                    "values from 10 to 100. Update the value of 'dw_region' "\
                        "in your JSON file such that it falls within these "\
                            "constraints."
                return None
            
        # if dw_region set to None, validate as default        
        elif IGPs['dw_region'] == None:
            print "\n -dw_region parameter in your JSON file set to null."
            print "\n -Setting parameter to default: 25"
            dw_region_checked = 25
        
        # if not the right type, return None
        else:
            print "the 'dw_region' parameter takes only positive integer values from"\
                " 25 to 1000. Update the value of 'dw_region' in your JSON file such"\
                    " that it falls within these constraints"
            return None
    else:
        print "\n -No dw_region parameter in your JSON file."
        print "\n -Setting parameter to default: 25"
        dw_region_checked = 25
        
    IGPs_checked.append(dw_region_checked)
    
    ### *** tax_ID checks *** ###
    if 'tax_ID' in IGPs:
        if isinstance(IGPs['tax_ID'], int):
            tax_ID = str(IGPs['tax_ID'])

            # Taxonomic identifiers are numerical values ranging from 1 digit to
            # 7 digits, if tax_ID satisfies that req., validate it; otherwise,
            # return None        
            if len(tax_ID) < 1 or len(tax_ID)> 7:
                print "The 'tax_ID' parameter reflects an NCBI taxonomic "\
                    "identifier which consists of all digits from 1 to 7 "\
                        "digits in length. Update the 'tax_ID' parameter in "\
                        "your JSON file such that it falls within these "\
                            "constraints."
                return None
            
            # Call the Entrez taxonomy database with the tax_ID given, if no
            # results returned, i.e. an empty list, return None
            Entrez.email = user_email_checked
            
            Entrez.apikey = apikey_checked
            
            tax_ID_check = Entrez.read(Entrez.efetch(db = 'taxonomy', \
                                                     id = tax_ID))
            if len(tax_ID_check) == 0:
                print "\n -The tax_ID parameter, " + tax_ID + ", in your JSON"\
                    " file isn't a valid taxonomic identifier."
                print " \n -Update your JSON file with a valid taxonomic "\
                    "identifier, or null, to set the default of None"
                return None
            
            else:
                tax_ID_checked = tax_ID
                
                IGPs_checked.append(tax_ID_checked)

                print "\n -BLAST results will be restricted to those within the group"\
                    " specified by taxonomic identifier: ", tax_ID_checked
        
        # if tax_ID set to None, validate the default
        elif IGPs['tax_ID'] == None:
            tax_ID_checked = None
            print "\n -The tax_ID parameter set to null in your JSON file."
            print "\n -BLAST results will not be restricted to any taxonomic group."
            IGPs_checked.append(tax_ID_checked)
        
        # if not the right type, return None
        else:
            print "The 'tax_ID' parameter reflects an NCBI taxonomic identifier "\
                " which consists of all digits from 1 to 7 digits in length. Update"\
                    " the 'tax_ID' parameter in your JSON file such that it falls"\
                    " within these constraints."
    else:
        print "\n -No tax_ID parameter in your JSON file."
        print "\n -Setting parameter to default: None"
        tax_ID_checked = None
        IGPs_checked.append(tax_ID_checked)
    
    ### *** min_cover check *** ###
    if 'min_cover' in IGPs:
        if isinstance(IGPs['min_cover'], (float, int)):
            min_cover = IGPs['min_cover']
            
            # if min_cover between .25 and 1, validate it, otherwise, return None
            if min_cover >= 0.25 and min_cover <= 1:
                min_cover_checked = min_cover
                print "\n -Minimum amino acid sequence coverage set to: "\
                    + str(min_cover_checked * 100) + "%"
                
            else:
                print "\nThe 'min_cover' parameter reflects a ratio of the BLAST"\
                    " alignment length to the length of the query amino acid sequence"\
                        " expressed as a float value between 0 and 1. Update the"\
                            " 'min_cover' parameter in your JSON file such that it"\
                                " falls within these constraints."
                return None
        # if min_cover set to None, validate the default
        elif IGPs['min_cover'] == None:
            min_cover_checked = 0.75
            print "\n -The min_cover parameter in your JSON file set to null."
            print "\n -Setting parameter to default: 0.75"

        # if not the right type, return None
        else:
            print "\nThe 'min_cover' parameter reflects a ratio of the BLAST"\
                " alignment length to the length of the query amino acid sequence"\
                    " expressed as a float value between 0 and 1. Update the"\
                        " 'min_cover' parameter in your JSON file such that it"\
                            " falls within these constraints."
            return None
    else:
        print "\n -No min_cover parameter in your JSON file."
        print "\n -Setting parameter to default: 0.75"
        min_cover_checked = 0.75
    IGPs_checked.append(min_cover_checked)
    
    ### *** sleepy check *** ###
    if 'sleepy' in IGPs:
        if isinstance(IGPs['sleepy'], (float, int)):
            sleepy = IGPs['sleepy']
    
            # if sleepy between 0 and 5 seconds, validate it; otherwise, return None
            if sleepy > 0 and sleepy <= 5:
                sleepy_checked = sleepy
                print "\n -Additional sleep time for queries (in seconds): ", \
                    sleepy_checked
    
            else:
                print "\nThe 'sleepy' parameter reflects an amount of time, "\
                    "in seconds, that the program will wait to query the NCBI "\
                        "database again in the event an HTTP error is returned. "\
                            "Update the 'sleepy' parameter such that it falls within"\
                                " these constraints (positive integer or float "\
                                    "<= 5)."
                return None
        # if sleepy set to None, validate the default
        elif IGPs['sleepy'] == None:
            sleepy_checked = 0.5 # default
            print "\n -The sleepy parameter in your JSON file set to null."
            print "\n -Setting parameter to default: 0.5 seconds"
        
        # if not the right type, return None
        else:
            print "\nThe 'sleepy' parameter reflects an amount of time, "\
                "in seconds, that the program will wait to query the NCBI "\
                    "database again in the event an HTTP error is returned. "\
                        "Update the 'sleepy' parameter such that it falls within"\
                            " these constraints (positive integer or float "\
                                "<= 5)."
            return None
    else:
        sleepy_checked = 0.5 # default
        print "\n -No sleepy parameter in your JSON file."
        print "\n -Setting parameter to default: 0.5 seconds"
        
    IGPs_checked.append(sleepy_checked)
        
        ### *** TF_family checks *** ###
    if 'TF_family' in IGPs:
        if isinstance(IGPs['TF_family'], unicode):
            TF_family_checked = str(IGPs['TF_family'])
            print "\n -TF family: " + TF_family_checked

            IGPs_checked.append(TF_family_checked)
            
        else:
            print "\nThe 'TF_family' parameter reflects a name for the TF to be "\
                "analyzed, which should contain a string type. Update your JSON file "\
                    "such that it falls within this constraint."
            return None
    else:
        print "For this program to work effectively, their must be a name for "\
            "your TF to be analyzed, please update your input parameters to "\
                "include a key, 'TF_family', and its value should be the name."
        return None
    
    ### *** outputfile_name checks *** ###
    if 'outputfile_name' in IGPs:
        
        if isinstance(IGPs['outputfile_name'], unicode):
            outputfile_name = str(IGPs['outputfile_name'])
            

            # if outputfile_name ends with .json, validate it; otherwise, return
            # None
            if outputfile_name.endswith('.json'):
                outputfile_name_checked = outputfile_name
                print "\n -Name of output JSON file: "\
                        + outputfile_name_checked
            
            else:
                outputfile_name_checked = outputfile_name + '.json'
                print "\n -Name of output JSON file: "\
                        + outputfile_name_checked

        else:
            print "\n -Name for output file given in your JSON file is invalid."
            print "\n - Setting output file name to default: "\
                "\n -cgb_input_generated.json"
            
            outputfile_name_checked = "cgb_input_generated.json"
        
    # if not the right type, return None
    else:
        print "\n -No outputfile_name parameter in your JSON file."
        print "\n -Setting output file name to default:"\
            " \n -cgb_input_generated.json"
        outputfile_name_checked = "cgb_input_generated.json"
    
    IGPs_checked.append(outputfile_name_checked)
    

    # Finally, add the CGB_parameters dict, a list of reference protein 
    # accession numbers, and the whole 'motifs' object in the input file
    IGPs_checked.append(CGB_parameters)
    
    TF_accessions = []
    for acc in in_file['motifs']:
        TF_accessions.append(acc['protein_accession'])
        
    IGPs_checked.append(TF_accessions)
    
    motif_object = in_file['motifs']
    
    IGPs_checked.append(motif_object)
    
    return IGPs_checked
    
def cgb_input(inputfile_name):
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
    try:
        IGPs_checked = input_check(inputfile_name)
    except:
        return None
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
    outputfile =  IGPs_checked[13]
    cgb_pipeline_parameters = IGPs_checked[14] 
    TF_accessions = IGPs_checked[15]
    motifs = IGPs_checked[16]
    Entrez.email = user_email
    Entrez.api_key = apikey

    #output dictionary with all the information included in JSON file
    #what will be changed herein is the motifs and genomes lists

    #we use OrderedDict so that the order in the JSON file is maintained for 
    #human readibilty/editability when the JSON output is generated 
    output = collections.OrderedDict([("TF", TF_family), ("motifs",[]), \
                                      ("genomes",[])])

    for key, val in cgb_pipeline_parameters.iteritems():
        output[key] = val

    #create dictionary that will store information about the orthologs, 
    #wich will be the target species for CGB
    orthologs = {}
    
    # GET TF ACCESSION NUMBERS AND GENOMES        
    # first, create a list of TFs, download each IPG record, then look for the specific TF provided by
    # the user and extract its encoding genome
    # Access the reference TF accession numbers (protein), get their IPG records,
    # and create a list of them
    for item in motifs:
        skipprot=False
        
        # Get genome accessions for the TFs
        for cnt in range(5):
            try:
                records = Entrez.read(Entrez.efetch(db = "protein", \
                                            id = item['protein_accession'], \
                                            rettype = 'ipg', retmode = 'xml'))
                    
                time.sleep(sleepy)
                break

            except:
                print 'NCBI exception raised.\n Reattempt iteration: ' + \
                str(cnt+1)
                if cnt == 4:
                    print 'Could not download record after 5 attempts'
                    # return None
        
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

        #GET THE NAME OF THE TF

        # First, download each protein record and create a variable for their 
        # respective names
        for cnt in range(5):
            try:
                protein_record = SeqIO.read(Entrez.efetch(db="protein", \
                                            id = item['protein_accession'], \
                                            rettype = "gp", \
                                            retmode = "text"), 'genbank')
                time.sleep(sleepy)
                break
            
            except:
                print 'NCBI exception raised. Reattempt iteration: ' + str(cnt+1)    
           
            if cnt == 4:
                print 'Could not download record after 5 attempts'
                continue

        #get the species name from the annotations "source", which identifies 
        #the organism, and prettyeze it
        name = prettyeze_spname(protein_record.annotations['source'])

        #create a compound name using the TF family (e.g. LexA_) and the 
        #organism name
        TF_name = TF_family + '_' + name
        
        print '\nProcessing: ' + TF_name
        
        # TF_accession (item['protein_accession]), TF_name, TF_genome, and 
        # binding sites (item['sites']) are compiled into a dictionary, then 
        # added to output for CGB's motif object.
        motif_object = {"name" : TF_name, "genome_accession" : TF_genome, \
                        "protein_accession" : item['protein_accession'], "sites" : \
                        item['sites']}
        
        output["motifs"].append(motif_object)
        
    #SECOND LOOP: SEARCH FOR ORTHOLOGS
    # BLAST search is performed using the blast_search function.
    # each hit stored in a dictionary called orthologs as keys with an empty 
    # dictionary as its value
    for accession in TF_accessions:
        print '\nBLASTing: ' + accession
        
        blast_hits = blast_search(accession, evalue, nhits, min_cover, \
                                txID, sleepy, blast_db)
        
        # Ensure no repeating accession numbers in the dictionary
        for hit in blast_hits:
            hit = str(hit)
            if hit not in orthologs and hit not in TF_accessions:
                orthologs[hit] = {}
        
    print "\nAll reference protein records in input file pre-processed and BLASTed."
    
    print len(orthologs), 'hits before processing'
    
    #THIRD LOOP: work with orthologs (targets)
 	#For each ortholog, we want to pick one genome (from the IPG record of that
    # protein) We may also want to restrict things further by imposing that 
    # only one ortholog per "clade" (as specified by the user) is pulled into 
    # the cgb config file
    
    toberemoved=[]
        
    for ortholog in orthologs:
    
        print '\nProt: ' + ortholog
    
        #use the function `genome_record_retrieval` to get the best selected 
        #encoding for the TF ortholog identified by BLASTP
        #this function will implement prioritization (best=complete RefSeq) 
        #and return the "best" matching genome with a CDS encoding the ortholog
        
        cds_accession = genome_record_retrieval(ortholog, sleepy)
        
        # make a list of the genome records for getting taxonomic info later
        
        #use the function `contig_accessions` to get the complete genome 
        #(all chromids) or the complete set of contigs (for a WGS assembly)
        #for the nucleotide record mapping to the ortholog
        if cds_accession != None:
            # cds_accessions.append(cds_accession['acc'])
            
            print '     Genome: ' + cds_accession['acc'] + ' with score: ' + \
                  str(cds_accession['p_score'])
            
            target_name, target_range = \
                contig_accessions(cds_accession['acc'],\
                                  cds_accession['p_score'], sleepy)
            
            if target_range == None:
                
                toberemoved.append(ortholog)
                
            # store the genomic cds data, contigs, p_score, and name in 
            # orthologs
            orthologs[ortholog]['nuc_record'] = cds_accession
            orthologs[ortholog]['genome'] = target_range
            orthologs[ortholog]['genome_score'] = cds_accession['p_score']
            orthologs[ortholog]['genome_name'] = target_name

        else:
            toberemoved.append(ortholog)
            
            
    # remove orthologs that did not pan out
    for ortholog in toberemoved:
        del orthologs[ortholog]
    
    print '\nTotal number of valid orthologs detected: ' + str(len(orthologs))

    # This control flow (Beginning here and continuing through the next if/else 
    # statement, is how we can sample by both taxonomic level AND max %ID, or
    # just sample by one of those two methods, or neither, and always perform
    # the sampling in the correct order

    # If the user wants to select the best ortholog from each group related to 
    # a specific taxonomical level (e.g. get one record for each species, 
    # genus, etc.)
    if selected_taxon:
        
        filtered_list = tax_level_filter(selected_taxon, orthologs, sleepy)
        
        filtered_orthologs = {}
        for ortholog in orthologs:
            
            if ortholog in filtered_list:
                
                filtered_orthologs[ortholog] = orthologs[ortholog]
                
    # If not sampling by a tax. level, define the orthologs dictionary, as it
    # currently is, as a new dictionary, filtered_orthologs.
                
    else:
        
        filtered_orthologs = orthologs
        
    # if restricting by maxID
    if maxID:
        
        filtered_list = maxID_filt(maxID, filtered_orthologs, up_region, \
                                    dw_region, sleepy)
        
        filtered_orthologs_final = {}
        for ortholog in filtered_orthologs:
            if ortholog in filtered_list:
                filtered_orthologs_final[ortholog] = filtered_orthologs[ortholog]

    # if all orthologs are wanted, without further study:
    else:        
        filtered_orthologs_final = filtered_orthologs
        
    
    for ortholog in filtered_orthologs_final:
        target = {"name": filtered_orthologs_final[ortholog]['genome_name'], \
                  "accession_numbers": filtered_orthologs_final[ortholog]['genome']}
            
        output["genomes"].append(target)

    print '\nTotal number of valid orthologs after filtering', \
        len(output['genomes'])

    #write out the json file
    print '\nWriting out JSON input file for CGB'
    json.dump(output, open(outputfile, "w"), indent = 4, \
              separators = (' , ', ': '))

    return output

final = cgb_input(input_filepath)
