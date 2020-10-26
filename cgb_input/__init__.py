from cgb_input import cgb_input_main


def create_file(input_file):
	"""
	Runs the input file creation script and returns the final JSON
	"""
    
	final = cgb_input_main(input_file)

	return final
