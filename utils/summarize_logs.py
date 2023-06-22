'''Summarize all SUMMA logs in a folder. Assumes all .txt files in folder are SUMMA logs.
Summary file is placed inside the log folder. Specifying a summary file name is optional.
Usage: python summarize_logs.py [log_folder] [name_of_summary_file.txt] [log file extension]'''

# Modules
import os
import re
import sys
import statistics as sts

# ----------------------
# Set defaults
summaryFile = '_log_summary.txt' # default, placed at the top of the log folder
ext = '.txt'

# Handle input arguments
if len(sys.argv) == 1: # sys.argv only contains the script name
    sys.exit('Error: no input folder specified')

else: # at least 2 elements in sys.argv; len(sys.argv) cannot be zero or we wouldn't be in this script

    # The first input argument specifies the folder where the log files are
    folder = sys.argv[1] # sys.argv values are strings by default so this is fine

    # Check if there are more arguments
    if len(sys.argv) == 3:

        # Assume the second argument is the name for the log file
        summaryFile = sys.argv[2] # string

    # No extra argument so no summary file name is specified
    elif len(sys.argv) == 4:

        # Assume the second argument is the name for the log file and the third is the file extension
        summaryFile = sys.argv[2] # string
        ext = sys.argv[3] # string

# End of input arguments
# ----------------------

# -------------
# Sub functions

# Define a function to grab the last line in a file
# See: https://stackoverflow.com/questions/136168/get-last-n-lines-of-a-file-similar-to-tail
def tail(folder, file, lines=1, _buffer=4098):

    """Tail a file and get X lines from the end"""

    # open the file
    with open(folder + '/' + file,'r') as f:

        # place holder for the lines found
        lines_found = []

        # block counter will be multiplied by buffer to get the block size from the end
        block_counter = -1

        # loop until we find X lines
        while len(lines_found) < lines:
            try:
                f.seek(block_counter * _buffer, os.SEEK_END)
            except IOError:  # either file is too small, or too many lines requested
                f.seek(0)
                lines_found = f.readlines()
                break

            lines_found = f.readlines()

            # decrement the block counter to get the next X bytes
            block_counter -= 1

    return lines_found[-lines:]

# Define a function to handle the three types of log results: success, SUMMA error, non-SUMMA error
def determine_output(folder,file,nLines=2):

    # initialize outputs
    success, summa, other = [0,0,0] # error flags
    time = -1 # times are only extracted for successful simulations - a default negative time allows filtering outside this function

    # get the file contents
    log_txt = tail(folder,file,nLines)

    # loop over the lines in the log text to find what happened
    msg = 'no line read' # initialize output
    for line in log_txt:

        # determine if the log contains a SUMMA statement
        if 'successfully' in line:

            # determine time taken
            more_lines = tail(folder,file,6) # extracting bottom 6 lines includes time in [h]
            time = float(re.sub("[^\d\.]", "", more_lines[0])) # hours should be in top entry; store as float

            # process output flags
            success = 1
            msg = 'success after ' + '{:.2f}'.format(time) + ' h \n' # message string
            return success, summa, other, msg, time # we know what happened, stop function call

        elif 'FATAL ERROR' in line:
            summa = 1
            msg = line
            return success, summa, other, msg, time # we know what happened, stop function call

    # if we reach this, no SUMMA termination statement was found
    other = 1
    msg = 'check SLURM logs - simulation terminated early at: ' + line

    return success, summa, other, msg, time

# End of sub functions
# --------------------

# -------------------
# Start of processing

# Remove the summar file if it exists
try:
    os.remove(folder + '/' + summaryFile)
except OSError:
    pass

# Find the .txt files in the folder
files = []
for file in os.listdir(folder):
    #if file.endswith(".txt"):
    if file.endswith(ext):
        files.append(file)

# Sort the list
files.sort()

# Initialize case counters
total_success, total_summa, total_other = [0,0,0]

# Initialize time list
computation_time = []

# Open the summary file
with open(folder + '/' + summaryFile, 'w') as sf:

    # Add a header
    sf.write('Summarizing log files in ' + folder + '\n \n')
    sf.write('Log files' + '\n')

    # Loop over the log files
    for file in files:

        size = os.path.getsize(folder + '/' + file)
        if (size==0): continue

        # Find the result contained in each log file
        success, summa, other, msg, time = determine_output(folder,file)  # default of using last 2 lines should suffice to catch all success/error cases

        # Increment the counters
        total_success += success
        total_summa   += summa
        total_other   += other

        # Save the computation time (time == -1 if error encountered)
        computation_time.append(time)

        # Add the file name and error message to the summary file
        sf.write(file + '\t' + msg)

    # Calculate percentages
    total       = total_success + total_summa + total_other
    pct_success = total_success / total * 100
    pct_summa   = total_summa / total * 100
    pct_other   = total_other / total * 100

    # Calculate computation time stats
    computation_time = [time for time in computation_time if time >= 0] # remove negative times that are associated with simulations that return an error
    st_min    = min(computation_time)
    st_max    = max(computation_time)
    st_mean   = sts.mean(computation_time)
    st_median = sts.median(computation_time)

    # add a statistical summary
    sf.write('\nSuccess stats' + '\n')
    sf.write('Success' + '\t \t \t \t {:.2f}% \n'.format(pct_success))
    sf.write('SUMMA error' + '\t \t \t {:.2f}% \n'.format(pct_summa))
    sf.write('Early termination' + '\t {:.2f}% \n'.format(pct_other))
    sf.write('\nTime needed for successful computations' + '\n')
    sf.write('Min time ' + '\t \t \t {:.2f} h \n'.format(st_min))
    sf.write('Median time ' + '\t \t {:.2f} h \n'.format(st_median))
    sf.write('Mean time ' + '\t \t \t {:.2f} h \n'.format(st_mean))
    sf.write('Max time ' + '\t \t \t {:.2f} h \n'.format(st_max))


    #sf.write('Min time ' + '\t \t \t \t' + str(min(computation_time)) +'h \n')
    #sf.write('Median time ' + '\t \t \t' + str(sts.median(computation_time)) +'h \n')
    #sf.write('Mean time ' + '\t \t \t \t' + str(sts.mean(computation_time)) +'h \n')
    #sf.write('Max time ' + '\t \t \t \t' + str(max(computation_time)) +'h \n')
# done, summary file is closed
# ----------------------------
