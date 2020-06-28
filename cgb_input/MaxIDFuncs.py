from Bio import Entrez, SeqIO, pairwise2
import time

def genome_record_to_seq(grecord, upstream, downstream, sleepy):
    """Gets a genome record consisting of accession, start position and strand.
       Queries NCBI to retrieve as many positions upstream and downstream as
       desired from TLS and returns a sequence object.

    """

    # assign positions in record, according to strand
    if grecord['strand'] == '+':
        s_start = int(grecord['start']) - upstream
        s_stop = int(grecord['start']) + downstream
        s_strand = 1

    else:
        s_stop = int(grecord['stop']) + upstream
        s_start = int(grecord['stop']) - downstream
        s_strand = 2

    # Download FASTA record containing the desired upstream sequence
    for cnt in range(5):

        try:

            net_handle = Entrez.efetch(db="nuccore", id=grecord['acc'], \
                                       strand=s_strand, seq_start=s_start, \
                                       seq_stop=s_stop, rettype='fasta', \
                                       retmode="txt")

            time.sleep(sleepy)

            break

        except:

            print 'NCBI exception raised. Reattempt iteration: ' + str(cnt + 1)

            if cnt == 4:
                print 'Could not download record after 5 attempts'
    
                return None
    
    gnome_record = SeqIO.read(net_handle, "fasta")

    time.sleep(sleepy)

    return gnome_record


######################################################################
def id_below_maxid_perc(el1, el2, max_percent_id):
    """Aligns two sequence elements and determines whether they are
       more than %ID identical (false) or not (true)
       Scoring: Match:+2, Mismatch:-1, GapO: -2, GapE: -0.2
    """

    al = pairwise2.align.globalms(el1.seq, el2.seq, 2, 0, -2, -.5, \
                                  one_alignment_only=True, \
                                  penalize_end_gaps=False)
    
    matches = 0
    gapless = 0
    
    # for each position in the alignment
    for ch_pair in zip(al[0][0], al[0][1]):

        # if this is a non-gapped position
        if '-' not in ch_pair:

            gapless = gapless + 1

            # if it's a match, count it
            if ch_pair[0] == ch_pair[1]:
                
                matches = matches + 1

    perID = (float(matches) / float(gapless))
    
    # return true or false depending on percent identity
    if perID <= float(max_percent_id):

        return (True)

    else:

        return (False)


def identity_filter_list(orthos, percent_id):
    """Gets a list of orthologs with upstream sequence and a max percent id.
       Goes through the list, removing any records with more than %ID.
       Returns the trimmed list.
    """
    filt_list = []

    # create a list of protein accession + promoter pairs
    item_list = []

    for key, value in orthos.iteritems():
        item_list.append([key, value['promoter'], value['genome_score']])

    # sort item_list so that the first elements are the "best" records
    # (i.e. those with better genome score)
    # this  way, they are likely to be selected if they are not >maxID
    item_list.sort(key=lambda k: k[2], reverse=True)

    cnt = 0

    while cnt < len(item_list):
        # get next first element in upstream seq list, removing it from list

        current_item = item_list.pop(0)

        # and add it to the filtered list

        filt_list.append(current_item)

        # check against all remaining elements in item_list; remove them from22
        # the list if they are not below threshold of identity
        # at each pass revised item_list hence contains only elements less
        # than %ID  identical to previously processed elements now stored in
        # filt_list
        item_list[:] = [upel for upel in item_list if \
                        id_below_maxid_perc(upel[1], current_item[1], percent_id)]

    # generate a filtered dictionary
    filt_orthos = {}

    for item in filt_list:
        filt_orthos[item[0]] = orthos[item[0]]

    return filt_orthos

########################################################################
def maxID_filt(maxID, prefilt_orthologs, up_region, dw_region, sleepy):
    print 'Obtaining upstream nucleotide sequences for homologs'

    for ortholog in prefilt_orthologs:
        print '  |-> Grabbing upstream nucleotide sequence for: ' + ortholog

        # grab the sequence corresponding to that record
        seq = genome_record_to_seq(prefilt_orthologs[ortholog]['nuc_record'], \
                                   up_region, dw_region, sleepy)

        prefilt_orthologs[ortholog]['promoter'] = seq

    # filter the list of upstream sequences
    filtered_orthologs = identity_filter_list(prefilt_orthologs, maxID)

    return filtered_orthologs

