from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import time, json


def try_efetch(database, identifier, ret_mode, sleepy, ret_type='text'):
    """
    Simple function to query entrez database, and safeguard against
    premature errors returned by NCBI if their server gets overwhelmed.
    """
    for cnt in range(5):

        try:

            record = Entrez.efetch(db=database, id=identifier,
                                   rettype=ret_type, retmode=ret_mode)
            time.sleep(sleepy)

            break

        except:

            print 'NCBI exception raised. Reattempt iteration: ' + str(cnt + 1)

            print 'Could not download record after 5 attempts...'

            print '\nreturning None'

            return None

    return record


def try_read_and_efetch(database, identifier, ret_type, sleepy,
                        ret_mode='text'):
    """
    This function will also query the entrez database, and safeguard
    against premature errors returned by NCBI if their server gets
    overwhelmed. The difference here is the record will be read
    within this function, with Entrez.read(), as opposed to being
    read by some other function outside this function.
    """
    for cnt in range(5):

        try:

            record = Entrez.read(Entrez.efetch(db=database,
                                               id=identifier,
                                               rettype=ret_type,
                                               retmode=ret_mode))
            time.sleep(sleepy)

            break

        except:

            print 'NCBI exception raised. Reattempt iteration: ' + str(cnt + 1)

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

            record = Entrez.read(Entrez.esearch(db=database,
                                                term=term_val,
                                                idtype=ID_type))

            time.sleep(sleepy)
            break

        except:

            print 'NCBI exception raised. Reattempt iteration: ' + str(cnt + 1)

            if cnt == 4:
                print 'Could not download record after 5 attempts...'

                print '\nreturning None'

                return None

    return record


def try_qblastp(db, seq, e_val, nhits, sleepy, tax_id=None):
    """
    Implements the NCBIWWW.qblast function according to the specific needs
    of the cgb_input program
    """
    # if taxon filtering
    if tax_id:

        taxon = "txid" + str(tax_id) + "[orgn]"

        # perform BLASTP search and parse results
        for cnt in range(5):

            try:

                handleresults = NCBIWWW.qblast(program='blastp',
                                               database=db,
                                               sequence=seq,
                                               entrez_query=taxon,
                                               expect=e_val, hitlist_size=nhits)

                time.sleep(sleepy)
                break

            except:
                print 'NCBI exception raised on attempt: ' + str(cnt + 1) + \
                      '\nreattempting now...'

                if cnt == 4:
                    print 'Could not download record after 5 attempts'

                    return None
    else:
        
        for cnt in range(5):

            try:
                # perform BLASTP search and parse results
                handleresults = NCBIWWW.qblast(program='blastp',
                                               database=db,
                                               sequence=seq, expect=e_val,
                                               hitlist_size=nhits)

                time.sleep(sleepy)
                break

            except:

                print 'NCBI exception raised on attempt: ' + str(cnt + 1) + \
                      '\nreattempting now...'

                if cnt == 4:
                    print 'Could not download record after 5 attempts'

                    return None

    blast_records = list(NCBIXML.parse(handleresults))

    return blast_records


def blast_search(TF_accession, dbase, cutoff, nhits, min_cover, tax_id, sleepy):

    log_file = {}
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

    # store all the hits in a list
    handle = try_efetch('protein', TF_accession, 'text', sleepy,
                        ret_type='fasta')
    if not handle:

        return None

    protrec = SeqIO.read(handle, "fasta")

    protseq = protrec.format('fasta')

    blast_records = try_qblastp(dbase, protseq, cutoff, nhits, sleepy, tax_id)

    orthologs = []

    # for each blast hit
    for record in blast_records[0].alignments:

        for hsp in record.hsps:

            # calculate coverage to weed out false positives from protein domains
            if min_cover:

                query_len = float(record.length)  # <-len. of TF's complete AA seq.

                align_len = float(hsp.align_length)  # <- len. of alignment

                coverage = float(align_len / query_len)

                # if min. coverage is less than coverage calculated:
                if coverage >= min_cover:

                    # add the blast hit to orthologs list
                    orthologs.append(record.hit_id.split('|')[-2])

                else:

                    log_file[record.hit_id.split('|')[-2]] = 'hit failed' \
                                                             ' coverage test'

            else:

                orthologs.append(record.hit_id.split('|')[-2])

    return orthologs, log_file
