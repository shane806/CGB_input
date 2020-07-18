import time
from cgb_input import create_file

input_file = 'pseudomonadales_test_erill.json'

start_time = time.time()

cgb_input_file = create_file(input_file)

overall_time = time.time() - start_time

print overall_time