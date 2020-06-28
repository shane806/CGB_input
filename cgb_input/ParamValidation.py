from Bio import Entrez
import json, urllib2
from NCBI_funcs import try_read_and_efetch


def input_check(inputfile_name):
    """This function takes 1 input: A file name of the JSON file as a string.
        It will be the first function called when the main cgb_input function
        is called.

        The goal is to validate each of the 16 user-provided parameters
        involved in generating the genome accession numbers of the orthologs
        that ultimately are passed into the next phase of the CGB pipeline.

        A couple key points about this function:

            - If an improper/incompatible value(s) is(are) given. The function
              will either set a default value if the user's input is not
              necessary, or will instruct the user to revise their input and
              the program will terminate

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
    IGPs_checked = {}

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
            print "\n Your JSON file isn't structured properly somewhere. This" \
                  " type of error can often be tough to pinpoint." \
                  "\n A potentially time saving solution to this issue" \
                  " is to copy/paste the contents of your JSON file into a web" \
                  " based JSON formatter and/or validator such as:" \
                  "\n https://jsonformatter.curiousconcept.com/ " \
                  "\n to quickly identify what syntax error is" \
                  " preventing a successful read in."

            return None

        elif except_type == IOError:

            print "\n Remember, the file name passed into the script must include " \
                  " the .json extension. For example, If your file's name is" \
                  " 'input_f','input_f.json' would be the correct syntax used " \
                  " (assuming input_f is a JSON file, a necessary condition of the input " \
                  " file in order to be compatible with this program, and located in " \
                  " your current working directory), including the quotation marks, " \
                  " for the file to be properly recognized and read in by python."

            print "\n If the file is not in your current directory, you " \
                  " must use a relative or an absolute path."

            print "\n e.g. ../currentworkingdirectory/myjsonfile.json or" \
                  " home/User/otherdirectory/myjsonfile.json"
            return None

    CGB_parameters = in_file['cgb_parameters']

    print '\nWorking with parameters as follows:'

    ### *** internet connection test *** ###
    try:

        urllib2.urlopen("https://www.google.com/", timeout=1)

    except:

        print "It is necessary to have a working internet connection to use" \
              " this program. Connect to an available active network and restart" \
              " the program."

        return None

    ### *** NCBI account email address check *** ###
    if 'entrez_email' in CGB_parameters:

        user_email = str(CGB_parameters['entrez_email'])

        if "@" in user_email and "." in user_email:

            user_email_checked = user_email

            print "\n -NCBI account email address:", user_email_checked

            IGPs_checked['user_email'] = user_email_checked

        else:

            print "The email address provided under 'entrez_email' within" \
                  " the 'cgb_parameters' object of your JSON file is invalid. Please" \
                  " update that value to the email address used to register " \
                  "for an NCBI account at www.ncbi.nlm.nih.gov"

            return None
    else:

        print "It's necessary for a key, 'entrez_email' to be in the " \
              "CGB_parameters dictionary of your JSON file, and for that email to" \
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

            print "Update the entrez_apikey parameter in your JSON file to a valid" \
                  " api key."

            return None

        apikey_checked = apikey

        print "\n -NCBI account api key:", apikey_checked

    else:

        print "Non NCBI account api key specified. Setting parameter to none."

        apikey_checked = None

    IGPs_checked['apikey'] = apikey_checked

    IGPs = in_file['input_parameters']

    ### *** BLAST_eval check *** ###

    # isinstance() checks the type and returns a boolean
    if 'BLAST_eval' in IGPs:

        if isinstance(IGPs['BLAST_eval'], float) and IGPs['BLAST_eval']:

            BLAST_eval = IGPs['BLAST_eval']

            # if the e-value is pos. and less than 1, validate it
            if 0 <= BLAST_eval <= 1:

                BLAST_eval_checked = BLAST_eval

                print "\n -BLAST limiting e-value set to: ", BLAST_eval_checked

            # if not, return None
            else:

                print '\nThe "BLAST_eval" parameter reflects "the statistical ' \
                      'significance threshold for reporting matches against database' \
                      ' sequences." For the purposes of this program, the ' \
                      '"BLAST_eval" parameter should be a positive number ' \
                      'less than 1.'

                print "\nSee https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=" \
                      "Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp#expect for more " \
                      "details."

                return None

        # if user wants default
        elif not IGPs['BLAST_eval']:

            BLAST_eval_checked = 1e-10  # default

            print "\n -BLAST_eval parameter in your JSON set to null (default.)"

            print "\n -Setting parameter to default: 1e-10"

        # if not the right type, return None
        else:
            print '\nThe "BLAST_eval" parameter reflects "the statistical ' \
                  'significance threshold for reporting matches against database' \
                  ' sequences." For the purposes of this program, the ' \
                  '"BLAST_eval" parameter should be a positive number ' \
                  'less than 1.'

            print "\nSee https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=" \
                  "Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp#expect for more " \
                  "details."

            return None

    else:

        BLAST_eval_checked = 1e-10  # default

        print " -No BLAST_eval parameter in your JSON file. Setting parameter" \
              " to default: 1e-10"

    IGPs_checked['evalue'] = BLAST_eval_checked

    ### *** BLAST_dbase check *** ###
    if 'BLAST_dbase' in IGPs:

        databases = ['nr', 'refseq_protein', 'landmark', 'swissprot', \
                     'pataa', 'pdb', 'env_nr', 'tsa_nr']

        if isinstance(IGPs['BLAST_dbase'], unicode) and IGPs['BLAST_dbase']:

            BLAST_dbase = str(IGPs['BLAST_dbase'])

            # database options

            # if BLAST_dbase given is one of the databases available, validate it
            if BLAST_dbase in databases:

                BLAST_dbase_checked = BLAST_dbase

                print "\n -Database to be queried in search of orthologs:" \
                    , BLAST_dbase_checked

            # if not, print the available options and return None
            else:

                print "The 'BLAST_dbase' parameter reflects the NCBI protein " \
                      "database which will be used to find orthologs. Here's a list" \
                      " of the available databases."

                for database in databases:
                    print "\n -" + database

                    print "\nUpdate the 'BLAST_dbase' parameter in your JSON file with " \
                          "one of the above databases, or null to set to default."
                return None

        # if user wants default
        elif IGPs['BLAST_dbase'] == None:

            BLAST_dbase_checked = 'nr'  # default

            print "\n -The default BLAST database, nr, will be used."

        # if parameter not the right type
        else:

            print "The 'BLAST_dbase' parameter reflects the NCBI protein " \
                  "database which will be used to find orthologs. Here's a list" \
                  " of the available databases."

            for database in databases:
                print " \n" + database

            print "\nUpdate the 'BLAST_dbase' parameter in your JSON file with " \
                  "one of the above databases."

            return None

    else:

        BLAST_dbase_checked = 'nr'  # default

        print "No database specified. Setting BLAST_dbase parameter to " \
              "default: nr"

    IGPs_checked['blast_db'] = BLAST_dbase_checked

    ### *** BLAST max_hits check ***
    if 'max_hits' in IGPs:

        nhits = IGPs['max_hits']

        if nhits and isinstance(nhits, int):

            # if positive and less than 500, validate it
            if nhits > 0 and nhits < 1000:

                nhits_checked = nhits
                print "\n -Max number of hits to be returned in BLAST search:" \
                    , nhits_checked
            else:
                print 'This parameter limits the max number of BLAST hits returned, so' \
                      ' a positive integer is the only compatible type for this parameter.' \
                      ' Update the value for max_hits to a whole number between 1 and 1000.'
                return None

        elif nhits == None:

            nhits_checked = 50

            print "\N -Max number of hits to be returned in BLAST search set to " \
                  "default: 50"

        else:

            print 'This parameter limits the max number of BLAST hits returned, so' \
                  ' a positive integer is the only compatible type for this parameter.' \
                  ' Update the value for max_hits to a whole number between 1 and 100.'

            return None

    else:

        nhits_checked = 50

        print "no max_hits parameter in JSON file. Setting parameter to default" \
            , nhits_checked

    IGPs_checked['nhits'] = nhits_checked

    ### *** selected_taxon check *** ###

    if 'selected_taxon' in IGPs:

        if IGPs['selected_taxon'] == None:

            selected_taxon_checked = None  # default

            print "\n -Null (default) value for selected_taxon parameter given" \
                  " in JSON."

            print "\n  --Default setting will be used: None"

        # Define the taxon options available

        #  if selected_taxon is one of the options, validate it. Or if 'null'
        # given as value in JSON, validate default
        else:

            taxon_levels = ['species', 'genus', 'family', 'order', 'class', 'phylum']

            if IGPs['selected_taxon'] in taxon_levels:

                selected_taxon = str(IGPs['selected_taxon'])

                selected_taxon_checked = selected_taxon

                print '\n -Sampling results by taxonomic level: ' \
                      + selected_taxon_checked

            # if not validated, print the available options and return None
            else:

                print "\nThe taxonomic level indicated in selected_taxon is " \
                      "incompatible with this program."

                print "\nPlease update your JSON with one of the " \
                      "following:"

                for taxon in taxon_levels:
                    print " -" + taxon

                    return None

    # if parameter not in IGPs, set to default
    else:

        print '\n -Orthologs will not be sampled according to a taxonomic ' \
              'level.'

        selected_taxon_checked = None  # default

    IGPs_checked['selected_taxon'] = selected_taxon_checked
    
    # if ID_filter was set to True but None given for maxID, set to default
    if 'maxID' in IGPs:

        if isinstance(IGPs['maxID'], (float)):

            maxID = float(IGPs['maxID'])

            # if maxID a pos. number between 0 and 1, validate it, otherwise,
            # return none

            if maxID > 0 and maxID < 1:

                maxID_checked = maxID

                print "\n -Orthologs used in analysis will be less than " \
                      + str(maxID_checked * 100) + "% identical to TFs"

                IGPs_checked['maxID'] = maxID_checked

            else:

                print 'The maxID parameter must be a float number between zero ' \
                      'and one...'

                print "\nUpdate your JSON with a value for maxID within these" \
                      " constraints."

                return None

        elif not IGPs['maxID']:

            print '\n -maxID parameter in your JSON file set to null.'

            print "\n -Setting parameter to default: 0.8"

            maxID_checked = 0.8  # default

            IGPs_checked['maxID'] = maxID_checked

        else:

            print 'The maxID parameter must be a float number between zero ' \
                  'and one...'

            print "\nUpdate your JSON with a value for maxID within these" \
                  " constraints."

            return None
    else:

        print "\n -No maxID parameter in your JSON file."

        print "\n -Setting parameter to default: None"

        maxID_checked = 0.8  # default

        IGPs_checked['maxID'] = maxID_checked

    ### *** up_region and dw_region checks *** ###

    if 'up_region' in IGPs:

        if isinstance(IGPs['up_region'], int):
            up_region = IGPs['up_region']

            # if upregion between 25 and 1000, validate it, otherwise, return None
            if up_region >= 25 and up_region <= 1000:
                up_region_checked = up_region
                print "\n -Number of bases to check in the upstream region: ", \
                    up_region_checked
            else:
                print "the 'up_region' parameter takes only positive integer " \
                      "values from 25 to 1000. Update the value of 'up_region' " \
                      "in your JSON file such that it falls within these " \
                      "constraints"
                return None

        # if up_region set to None, validate it as default
        elif IGPs['up_region'] == None:
            print "\n -up_region parameter set to null in your JSON."
            print "\n -Setting parameter to default: 250"
            up_region_checked = 250

        # if up_region not the right type, return None
        else:
            print "the 'up_region' parameter takes only positive integer values from" \
                  " 25 to 1000. Update the value of 'up_region' in your JSON file such" \
                  " that it falls within these constraints"
            return None
    else:
        print "\n -No up_region parameter in your JSON file."
        print "\n -Setting parameter to default: 250"
        up_region_checked = 250

    IGPs_checked['up_region'] = up_region_checked

    if 'dw_region' in IGPs:
        if isinstance(IGPs['dw_region'], int):
            dw_region = IGPs['dw_region']

            # if dw_region between 10 and 100, validate it, otherwise, return None
            if dw_region >= 0 and dw_region <= 100:
                
                dw_region_checked = dw_region
                
                print "\n -Number of bases to check in the downstream region: ", \
                    dw_region_checked

            else:
                print "The 'dw_region' parameter takes only integers greater than or equal to zero. " \
                      "Update the value of 'dw_region' " \
                      "in your JSON file such that it falls within these " \
                      "constraints."
                return None

        # if dw_region set to None, validate as default
        elif IGPs['dw_region'] == None:
            print "\n -dw_region parameter in your JSON file set to null."
            print "\n -Setting parameter to default: 25"
            dw_region_checked = 25

        # if not the right type, return None
        else:

            print "the 'dw_region' parameter takes only positive integer values from" \
                  " 25 to 1000. Update the value of 'dw_region' in your JSON file such" \
                  " that it falls within these constraints"
            return None
    else:
        print "\n -No dw_region parameter in your JSON file."
        print "\n -Setting parameter to default: 25"
        dw_region_checked = 25

    IGPs_checked['dw_region'] = dw_region_checked

    ### *** tax_ID checks *** ###
    if 'tax_ID' in IGPs:

        if isinstance(IGPs['tax_ID'], int):

            tax_ID = str(IGPs['tax_ID'])

            # Taxonomic identifiers are numerical values ranging from 1 digit to
            # 7 digits, if tax_ID satisfies that req., validate it; otherwise,
            # return None
            if len(tax_ID) < 1 or len(tax_ID) > 7:
                print "The 'tax_ID' parameter reflects an NCBI taxonomic " \
                      "identifier which consists of all digits from 1 to 7 " \
                      "digits in length. Update the 'tax_ID' parameter in " \
                      "your JSON file such that it falls within these " \
                      "constraints."
                return None

            # Call the Entrez taxonomy database with the tax_ID given, if no
            # results returned, i.e. an empty list, return None
            Entrez.email = user_email_checked

            Entrez.apikey = apikey_checked

            tax_ID_check = try_read_and_efetch(database='taxonomy',
                                               identifier=tax_ID, ret_mode='xml', sleepy=.5,
                                               ret_type='text')
            if len(tax_ID_check) == 0:
                print "\n -The tax_ID parameter, " + tax_ID + ", in your JSON" \
                                                              " file isn't a valid taxonomic identifier."
                print " \n -Update your JSON file with a valid taxonomic " \
                      "identifier, or null, to set the default of None"
                return None

            else:
                tax_ID_checked = tax_ID

                IGPs_checked['tax_ID'] = tax_ID_checked

                print "\n -BLAST results will be restricted to those within the group" \
                      " specified by taxonomic identifier: ", tax_ID_checked

        # if tax_ID set to None, validate the default
        elif IGPs['tax_ID'] == None:
            tax_ID_checked = None

            print "\n -The tax_ID parameter set to null in your JSON file."

            print "\n -BLAST results will not be restricted to any taxonomic group."

            IGPs_checked['tax_ID'] = tax_ID_checked

        # if not the right type, return None
        else:
            print "The 'tax_ID' parameter reflects an NCBI taxonomic identifier " \
                  " which consists of all digits from 1 to 7 digits in length. Update" \
                  " the 'tax_ID' parameter in your JSON file such that it falls" \
                  " within these constraints."
    else:

        print "\n -No tax_ID parameter in your JSON file."

        print "\n -Setting parameter to default: None"

        tax_ID_checked = None

        IGPs_checked['tax_ID'] = tax_ID_checked

    ### *** min_cover check *** ###
    if 'min_cover' in IGPs:
        if isinstance(IGPs['min_cover'], (float, int)):
            min_cover = IGPs['min_cover']

            # if min_cover between .25 and 1, validate it, otherwise, return None
            if min_cover >= 0.25 and min_cover <= 1:

                min_cover_checked = min_cover

                print "\n -Minimum amino acid sequence coverage set to: " \
                      + str(min_cover_checked * 100) + "%"

            else:

                print "\nThe 'min_cover' parameter reflects a ratio of the BLAST" \
                      " alignment length to the length of the query amino acid sequence" \
                      " expressed as a float value between 0 and 1. Update the" \
                      " 'min_cover' parameter in your JSON file such that it" \
                      " falls within these constraints."

                return None
        # if min_cover set to None, validate the default
        elif IGPs['min_cover'] == None:

            min_cover_checked = 0.75

            print "\n -The min_cover parameter in your JSON file set to null."

            print "\n -Setting parameter to default: 0.75"

        # if not the right type, return None
        else:

            print "\nThe 'min_cover' parameter reflects a ratio of the BLAST" \
                  " alignment length to the length of the query amino acid sequence" \
                  " expressed as a float value between 0 and 1. Update the" \
                  " 'min_cover' parameter in your JSON file such that it" \
                  " falls within these constraints."

            return None
    else:
        print "\n -No min_cover parameter in your JSON file."
        print "\n -Setting parameter to default: 0.75"
        min_cover_checked = 0.75
    IGPs_checked['min_cover'] = min_cover_checked

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
                print "\nThe 'sleepy' parameter reflects an amount of time, " \
                      "in seconds, that the program will wait to query the NCBI " \
                      "database again in the event an HTTP error is returned. " \
                      "Update the 'sleepy' parameter such that it falls within" \
                      " these constraints (positive integer or float " \
                      "<= 5)."

                return None

        # if sleepy set to None, validate the default
        elif IGPs['sleepy'] == None:
            sleepy_checked = 0.5  # default
            print "\n -The sleepy parameter in your JSON file set to null."
            print "\n -Setting parameter to default: 0.5 seconds"

        # if not the right type, return None
        else:

            print "\nThe 'sleepy' parameter reflects an amount of time, " \
                  "in seconds, that the program will wait to query the NCBI " \
                  "database again in the event an HTTP error is returned. " \
                  "Update the 'sleepy' parameter such that it falls within" \
                  " these constraints (positive integer or float " \
                  "<= 5)."

            return None

    else:

        sleepy_checked = 0.5  # default

        print "\n -No sleepy parameter in your JSON file."

        print "\n -Setting parameter to default: 0.5 seconds"

    IGPs_checked['sleepy'] = sleepy_checked

    ### *** TF_family checks *** ###
    if 'TF_family' in IGPs:

        if isinstance(IGPs['TF_family'], unicode):

            TF_family_checked = str(IGPs['TF_family'])


        else:

            print "\nThe 'TF_family' parameter reflects a name for the TF to be " \
                  "analyzed. Given the absence of a user provided parameter, "\
                  "The default,'None', will be used"
        
            TF_family_checked = None


    else:

        print "\nThe 'TF_family' parameter reflects a name for the TF to be " \
              "analyzed. Given the absence of a user provided parameter, "\
              "The default,'None', will be used"
        
        TF_family_checked = None

    print "\n -TF family: " + TF_family_checked

    IGPs_checked['TF_family'] = TF_family_checked

    ### *** outputfile_name checks *** ###
    if 'outputfile_name' in IGPs:

        if isinstance(IGPs['outputfile_name'], unicode):

            outputfile_name = str(IGPs['outputfile_name'])

            # if outputfile_name ends with .json, validate it; otherwise, return
            # None
            if outputfile_name.endswith('.json'):

                outputfile_name_checked = outputfile_name

                print "\n -Name of output JSON file: " + outputfile_name_checked

            else:

                outputfile_name_checked = outputfile_name + '.json'

                print "\n -Name of output JSON file: " \
                      + outputfile_name_checked

        else:

            print "\n -Name for output file given in your JSON file is invalid."

            print "\n - Setting output file name to default: " \
                  "\n -cgb_input_generated.json"

            outputfile_name_checked = "cgb_input_generated.json"

    # if not the right type, return None
    else:

        print "\n -No outputfile_name parameter in your JSON file."

        print "\n -Setting output file name to default:" \
              " \n -cgb_input_generated.json"

        outputfile_name_checked = "cgb_input_generated.json"

    IGPs_checked['outputfile_name'] = outputfile_name_checked

    # Finally, add the CGB_parameters dict, a list of reference protein
    # accession numbers, and the whole 'motifs' object in the input file
    IGPs_checked['cgb_parameters'] = CGB_parameters

    TF_accessions = []

    for acc in in_file['motifs']:
        
        TF_accessions.append(acc['protein_accession'])

    IGPs_checked['TF_accessions'] = TF_accessions

    motif_object = in_file['motifs']

    IGPs_checked['motif_object'] = motif_object

    return IGPs_checked
