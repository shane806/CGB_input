import os

def make_log_directories(TF_family):
    curr_dir_list = os.listdir(os.getcwd())

    if 'log_files' not in curr_dir_list:

        main_log_directory = os.getcwd() + "/log_files/"

        os.mkdir(main_log_directory)

    else:

        main_log_directory = os.getcwd() + "/log_files/"

    TF_log_dir = str(TF_family) + "_logs"

    if TF_log_dir not in os.listdir(main_log_directory):

        os.mkdir(main_log_directory + TF_log_dir)

    log_file_directory = main_log_directory + TF_log_dir + '/'
    
    if 'output_files' not in curr_dir_list:
        
        output_dir = os.getcwd() + "/output_files/"
        
        os.mkdir(output_dir)
        
    else:
        
        output_dir = os.getcwd() + '/output_files/'
        
    return log_file_directory, output_dir


