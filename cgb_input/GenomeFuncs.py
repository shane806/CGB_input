from Bio import Entrez
from NCBI_funcs import try_read_and_efetch, try_read_and_esearch
import json


def genome_record_retrieval(ortholog_acc, sleepy, log_dir):
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
    log_file = {}

    # genome_ret_dir = log_dir + 'genome_record_retrieval.json'

    # Attempt to download the ortholog's IPG record
    records = try_read_and_efetch('protein', ortholog_acc, 'ipg', sleepy, 'xml')

    if not records:
        # log_file[ortholog_acc] = 'genome retrieval: efetch of IPG record failed'

        return None

    # create scoring for priorization for non-complete genomes:
    # genbank, EMBL, & DDBJ nucleotide accession
    # numbers indexed here:
    # 'https://www.ncbi.nlm.nih.gov/Sequin/acc.html'
    # 'http://www.nslc.wustl.edu/elgin/genomics/bio4342/1archives/2006/AccReference.pdf'

    # and the list of refseq . numbers:
    # https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_numbers_and_mole/?report=objectonly
    priority = {'NC_': 7, 'AC_': 7, "AE": 6, "CP": 6, "CY": 6, "NZ_": 5,
                "NT_": 5, "NW_": 5, "U": 3, "AF": 3, "AY": 3, "DQ": 3,
                'EF': 3, ' EU': 3, 'FJ': 3, 'GQ': 3, 'GU': 3,
                'HM': 3, 'HQ': 3, 'JF': 3, 'JN': 3, 'JQ': 3,
                'JX': 3, 'KC': 3, 'KF': 3, ' KJ': 3, 'KM': 3, 'KP': 3,
                'KR': 3, ' KT': 3, 'KU': 3, 'KX': 3, 'KY': 3, ' MF': 3,
                'MG': 3, 'MH': 3, 'MK': 3, ' MN': 3, 'MT': 3}

    # list to hold all genomes and dictionary that is eventually returned as
    # the 'best' option out of all genomes
    genomelist = []

    best_genome = {}

    # from the IPG report, retrieve all the genome records from all CDS
    # listed, keeping only accession number, location of start and stop positions
    if 'ProteinList' in records['IPGReport']:

        for protein_rec in records['IPGReport']['ProteinList']:

            if 'CDSList' in protein_rec:

                for cds in protein_rec['CDSList']:

                    cds_dict = {'acc': cds.attributes['accver'],
                                'start': cds.attributes['start'],
                                'stop': cds.attributes['stop'],
                                'strand': cds.attributes['strand'],
                                'p_score': 0}

                    # assign score according t0 the above defined criteria
                    for key in priority:

                        if cds_dict['acc'].startswith(key):
                            cds_dict['p_score'] = priority[key]

                        if cds_dict['p_score'] > 6.5:
                            best_genome = cds_dict  # Complete refseq

                            return best_genome

                    # special case: NZ_ records that map to "complete" WGS
                    # NZ_ records with a shorter ( <7 ) segment of trailing digits
                    # e.g. NZ_CP030158
                    if cds_dict['acc'].startswith('NZ_'):

                        pr, sf = split_accession(cds_dict['acc'])

                        # complete WGS genomes seem to always have < 7 trailing
                        # digits in their accession, so we pick as "complete"
                        # anything with less than 7
                        if len(sf.split('.')[0]) < 7:
                            cds_dict['p_score'] = 6.5

                            print '     Atypical NZ_: ' + cds_dict['acc'] + \
                                  ' upgraded to RefSeq complete genome.'

                    # finally, score gb WGS ("AAAA-AZZZ" prefixes)
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

                # So, to get all CDS regions available, simply continue into the
                # the next protein record.
                # (EXAMPLE: WP_150851050.1)
                continue

    else:

        print str(ortholog_acc) + " has no  ProteinList feature in its " \
                                 "IPGReport"
        log_file[ortholog_acc] = 'from genome_record_retrieval: no IPGReport'

        return None


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
    # -> genbank accessions do NOT contain '_' characters, refseq accessions do
    # https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.
    # refseq_accession_numbers_and_mole/?report=objectonly
    # http://www.nslc.wustl.edu/elgin/genomics/bio4342/1archives/2006/AccReference.pdf
    # https://www.ncbi.nlm.nih.gov/genbank/wgs/

    log_file = {}

    refseq = 0

    if '_' in nuc_record_acc:
        refseq = 1
    # for complete genome records, search other chromosomes or plasmids.
    # This can be performed retrieving the Assembly database  accession and
    # using it in ESearch

    # Download record information from Nucleotide
    nuc_rec = try_read_and_efetch(database='nucleotide',
                                  identifier=nuc_record_acc,
                                  ret_type=None, sleepy=sleepy,
                                  ret_mode='xml')

    if not nuc_rec:
        log_file[nuc_record_acc] = "initial nucleotide record couldn't be " \
                                   "retrieved by efetch"

        return None

    sp_name = prettyeze_spname(nuc_rec[0]['GBSeq_organism'])

    # if nuc_rec is part of a complete genome record (p_score >= 6), we need
    # to see if there are any plasmids or other unique genome records associated
    # with the organism, so we first attempt to esearch with a BioSample ID,
    # if its in nuc_rec. If not, we search with the BioProject ID. If that's
    # not in nuc_rec, or the search yields zero results, we just use the
    # record obtained in genome_record_retrieval
    if nuc_record_score >= 6:

        Assembly_found = False
        BioSample_found = False
        BioProject_found = False

        for element in nuc_rec[0]["GBSeq_xrefs"]:

            if element["GBXref_dbname"] == "Assembly":

                Assembly_found = True

                Assembly = element['GBXref_id'] + "[Assembly]"

            elif element["GBXref_dbname"] == "BioSample":

                BioSample_found = True

                BioSample = element["GBXref_id"] + "[BioSample]"

            else:

                if 'BioProject' in element['GBXref_dbname']:
                    BioProject = element['GBXref_id'] + '[BioProject]'

                    BioProject_found = True

        # in case more genome records exist. e.g. plasmids not associated with
        # the complete circular chromosome

        if Assembly_found:

            Id_List = try_read_and_esearch(database='nucleotide',
                                           term_val=Assembly,
                                           ID_type='acc',
                                           sleepy=sleepy)['IdList']

        elif BioSample_found:

            Id_List = try_read_and_esearch(database='nucleotide',
                                           term_val=BioSample,
                                           ID_type='acc',
                                           sleepy=sleepy)['IdList']

        elif BioProject_found:

            bioP_id_list = try_read_and_esearch(database='BioSample',
                                                term_val=BioProject,
                                                ID_type='acc',
                                                sleepy=sleepy)['IdList']

            if bioP_id_list:

                for Id in bioP_id_list:

                    BioSample = Id + '[Biosample]'

                    Id_List = try_read_and_esearch(database='nucleotide',
                                                   term_val=BioSample,
                                                   ID_type='acc',
                                                   sleepy=sleepy)['IdList']

                    if nuc_record_acc in Id_List:

                        break

                    else:

                        print '***Attempt to get all genomes for ' + \
                              nuc_record_acc + ' proved unsuccessful; Ortholog will be' \
                                               ' discarded'

                        return sp_name, None

            else:

                print '***Attempt to get all genomes for ' + \
                      nuc_record_acc + ' proved unsuccessful; Ortholog will be' \
                                       ' discarded'

                return sp_name, None

        else:

            print '***Attempt to get all genomes for ' + \
                  nuc_record_acc + ' proved unsuccessful; Ortholog will be' \
                                   ' discarded'

            return sp_name, None

        # Check genome records returned in Id_list. If the original genome record
        # is a refseq record, grab only refseq records. If the original
        # genome record is a genbank record, grab only genbank records.

        contigs = []

        for nuc_record in Id_List:

            if refseq == True and '_' in nuc_record:

                contigs.append(nuc_record)

            elif '_' not in nuc_record:

                contigs.append(nuc_record)

        return sp_name, contigs

    # for WGS records:

    elif nuc_record_score <= 5:

        # if the record is neither a primary record nor a master record, then
        # it must be a contig record --> a master record will be created using
        # the `create_master` function

        if "GBSeq_alt-seq" not in nuc_rec[0]:
            print '             ' + nuc_record_acc + ' is not a master record'

            # generate a XXX0000000-type master record accession from contig acc
            master_record_acc = create_master(nuc_record_acc)

            print '             --> Master record accession created: ' \
                  + master_record_acc

            # in principle, records missing "GBSeq_alt-seq" should be WGS
            # contigs (with an associated master), but not all are.
            # so first we check whether the generated master record can be
            # obtained;
            # grab the nucleotide associated with the master record
            try:

                handle = Entrez.efetch(db="nucleotide", id=master_record_acc, \
                                       rettype=None, retmode='xml')

            except:

                print '             --> Master record ' + master_record_acc + \
                      ' did not pan out'

                log_file[nuc_record_acc] = "Master record created didn't pan " \
                                           "out"

                # if we cannot get to the "master" we created
                # get species name, replacing spaces by underscores

                sp_name = prettyeze_spname(nuc_rec[0]['GBSeq_definition'])

                # return only the accession, as None [to discard it]
                return sp_name, None

            # if the master record accssion worked out, grab it
            master_rec = Entrez.read(handle)

        # if it was already a master record
        # assign the previously downloaded nucleotide info
        else:

            master_rec = nuc_rec

        # Double check that this is indeed a master record
        # ("GBSeq_alt-seq" must be a key of the record)
        if ("GBSeq_alt-seq" in master_rec[0]):

            # ranges deals with the fact that some master records do not hold
            # a contiguous set of contigs
            # for instance: the WGS mater record NZ_AUIV00000000.1 encompasses:
            # WGS_SCAFLD  NZ_AUIV01000001-NZ_AUIV01000002

            # WGS_SCAFLD  NZ_AUIV01000005-NZ_AUIV01000015

            # WGS_SCAFLD  NZ_AUIV01000020-NZ_AUIV01000026

            # WGS_SCAFLD  NZ_AUIV01000029-NZ_AUIV01000057

            # each interval is kept as an element in ranges
            ranges = []

            # process the record depending on whether it is genbank or refseq
            if refseq:

                # get first and last WGS_SCAFLD contig records listed in master
                # For non-contiguous records (example: NZ_AUIV00000000.1):
                # Each interval is presented as an element of the ranges list.
                # Then, several elements will contain the WGS_SCAFLD tag. So,
                # we iterate over the list and add those intervals to ranges
                for item in master_rec[0]["GBSeq_alt-seq"]:

                    if item['GBAltSeqData_name'] == 'WGS_SCAFLD':

                        # if this range contains more than one item
                        if "GBAltSeqItem_last-accn" in \
                                item["GBAltSeqData_items"][0].keys():

                            first = item["GBAltSeqData_items"][0] \
                                ["GBAltSeqItem_first-accn"]

                            last = item["GBAltSeqData_items"][0] \
                                ["GBAltSeqItem_last-accn"]
                        else:

                            first = item["GBAltSeqData_items"][0] \
                                ["GBAltSeqItem_first-accn"]

                            last = first

                        ranges.append([first, last])

            # if genbank
            else:
                # get first and last WGS contig records listed in master rec
                # the same procedure as for RefSeq, but looking for 'WGS' tags

                for item in master_rec[0]["GBSeq_alt-seq"]:
                    if item['GBAltSeqData_name'] == 'WGS':

                        # if this range contains more than one item
                        if "GBAltSeqItem_last-accn" in \
                                item["GBAltSeqData_items"][0]:

                            first = str(item["GBAltSeqData_items"][0] \
                                            ["GBAltSeqItem_first-accn"])

                            last = str(item["GBAltSeqData_items"][0] \
                                           ["GBAltSeqItem_last-accn"])
                        else:
                            
                            first = str(item["GBAltSeqData_items"][0] \
                                            ["GBAltSeqItem_first-accn"])

                            last = first

                        ranges.append([first, last])

            # get the list of contigs associated to the WGS record
            contig_list = []

            # for each interval in range
            for interval in ranges:
                # get the first and last accessions

                first = interval[0]
                last = interval[1]

                # get the prefix and the numbers for start and end
                prefix, first_num = split_accession(first)

                prefix, last_num = split_accession(last)

                suffix_len = len(first_num)

                # for each possible number between first and last, generate the
                # corresponding accession
                for number in range(int(first_num), int(last_num) + 1):
                    # get number from first to last, transform number to string,
                    # look at size difference between number and the "normal"
                    # suffix (numeric) length for these records
                    # pad with zeroes on the left to get the completed string
                    # corresponding to the numeric part
                    num_str = str(number)

                    padding = '0' * (suffix_len - len(num_str))

                    num_str = padding + num_str

                    # create the composite accession (prefix+suffix) and add
                    # to contig list
                    contig_list.append(prefix + num_str)

            # get species name, replacing spaces with underscores
            sp_name = prettyeze_spname(master_rec[0]['GBSeq_definition'])

            return sp_name, contig_list

        else:

            print '             --> Master record does not contain adequate' \
                  ' fields.'

            log_file[nuc_record_acc] = "Master record didn't contain adequate" \
                                       " fields"

            # if the master record that we generated does not work well
            # get species name, replacing spaces by underscores
            sp_name = prettyeze_spname(nuc_rec[0]['GBSeq_definition'])

            # return name, accession, as None [to discard it]
            return sp_name, None


######################################################################

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
    # e.g. "NZ_NFJD01000014" in "NZ_NFJD01000014.1"
    acc_preperiod = acc[:acc.find('.')]

    # get the part after the period in the accession
    # e.g. "NZ_NFJD01000014" in "NZ_NFJD01000014.1"
    acc_postperiod = acc[acc.find('.'):]

    # iterate through accession prefix (e.g. "NZ_NFJD01000014")
    for character in acc_preperiod:

        # if the character is a numeric digit, add a zero to master accession
        if character.isdigit():
            master += '0'

        # otherwise (e.g. NZ_NFJD) add the same character
        else:
            master += character

    # put together the master accession by combining new prefix and suffix
    master = master + acc_postperiod

    return master


######################################################################

def split_accession(acc):
    """Splits a genbank/refseq record into prefix (e.g. NZ_MDWE or MDWE) and
       numeric component (e.g. 0100092), and returns them as a list

    """
    # go through the record until we get the first numeric digit

    # store the prefix and use its length to get the numeric region

    prefix = ''
    for c in acc:
        if c.isdigit():
            break
        else:
            prefix = prefix + c

    suffix = acc[len(prefix):len(acc)]

    return [prefix, suffix]


######################################################################
def prettyeze_spname(spname):
    """Takes a nucleotide record, accesses its 'GBSeq_definition' tag to
       get the species name and then removes characters that might not work
       well with cgb (more appropriately, with programs that cgb calls, such
       as clustalw and BLAST)
    """

    mod_name = \
        spname.split(',')[0].replace(' ', '_').replace('.', '_').replace('-', '_')

    return (mod_name)

##############################################################################
