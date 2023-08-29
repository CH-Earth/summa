#!/bin/bash
# run as bash write_seff.bash /path/to/log.txt /path/to/output.txt

# Set the path to the log file
LOG_FILE="$1"

# Set the path to the output file
OUTPUT_FILE="$2"

# Create the output file if it doesn't exist, and overwrite it with the new output
touch "$OUTPUT_FILE"
> "$OUTPUT_FILE"

# Loop through each line in the log file with "log_slurmId_"
 grep "log_slurmId_" "$LOG_FILE" | while read -r line ; do 

    # Extract the job ID from the line
    job_id=$(echo "$line" | grep -oP 'log_slurmId_\K\d+(?=_arrayId_)')

    # Call the seff command on the job ID and store the output in a variable
    seff_output=$(seff "$job_id")

    # Extract the CPU efficiency percentage from the seff output and convert it to a decimal format
    cpu_efficiency=$(echo "$seff_output" | grep -oP 'CPU Efficiency: \K\d+\.\d+')
    cpu_efficiency=$(echo "scale=4; $cpu_efficiency/100" | bc)

    # Extract the Array Job ID from the seff output
    array_id=$(echo "$seff_output" | grep -oP 'Array Job ID: \K\d+_\d+' | awk -F'_' '{print $2}')

    # Extract the Job Wall-clock time from the seff output, and convert to hours
    wall_clock=$(echo "$seff_output" | grep -oP 'Job Wall-clock time: \K.*')
    if [[ "$wall_clock" == *-* ]]; then
        days=$(echo "$wall_clock" | awk -F'-' '{print $1}')
        hours=$(echo "$wall_clock" | awk -F'-' '{print $2}' | awk -F':' '{print $1}')
        minutes=$(echo "$wall_clock" | awk -F'-' '{print $2}' | awk -F':' '{print $2}')
    else
        days=0
        hours=$(echo "$wall_clock" | awk -F':' '{print $1}')
        minutes=$(echo "$wall_clock" | awk -F':' '{print $2}')
    fi
    total_hours=$(echo "scale=2; ($days*24)+$hours+($minutes/60)" | bc)

    # Call the sacct command on the job ID and store the output in a variable
    # Note, the user has to be changed to the user who ran the job
    sacct_output=$(sacct -j "$job_id" -u avanb -o NodeList)

    # Extract the node number from the NodeList field
    # Note, this assumes it is running on graXXXX
    node_list=$(echo "$sacct_output" | tail -n 1 | sed 's/^[[:space:]]*//')
    node_number=$(echo "$node_list" | sed 's/gra//')

    # Write the CPU efficiency, Array Job ID, Job Wall-clock time, and the node number to the output file
    echo "CPU Efficiency: $cpu_efficiency, Array ID: $array_id, Job Wall-clock time: $total_hours, Node Number: $node_number" >> "$OUTPUT_FILE"
 done
