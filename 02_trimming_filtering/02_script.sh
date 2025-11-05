#!/usr/bin/bash
###########################################################################################
#
# Authors: Nerea Barrio and Gabriel Coll
# Date: 2025-05-18
# Version: 1.0.1
#
#
# Description: Script that trims the fastq.gz files of the RNA-seq using Fastp. 
# The output will be a fastq.gz file for each trimmed read and also a .html and .json reports for the results of each sample.
# The script will create a 02_results folder in the output directory to store all the results and, inside it, one folder for the results of each sample, named with the sample's name.
# It will also create a 02_logs folder in the output directory to store the logs of the script.
# By default it will detect and eliminate automatically the adapters for paired-end reads, do not evaluate the duplication rate, cut the bases of the sliding window if the mean quality of that window is lower than 20, and discard the reads shorter than 30 bases
# The user can then add more options to the trimming process and modify the default values of the parameters.
# The following options of fastp will be used: 
# -i: input file with the forward read (fastq.gz)
# -I: input file with the reverse read (fastq.gz)
# -o: output file with the trimmed forward read (fastq.gz)
# -O: output file with the trimmed reverse read (fastq.gz)
# --failed_out: output file with the failed reads (fastq.gz)
# --detect_adapter_for_pe: detect automatically adapters for paired-end reads and eliminate them
# --dont_eval_duplication: do not evaluate the duplication rate (as it is a RNA-seq, we normally expect many duplications that do not mean something went wrong, so we will save time and resources if we do not evaluate them)
# --trim_front1: trim the front of the forward read in N number of bases (default: 0)
# --trim_front2: trim the front of the reverse read in N number of bases (default: 0)
# --trim_tail1: trim the tail of the forward read in N number of bases (default: 0)
# --trim_tail2: trim the tail of the reverse read in N number of bases (default: 0)
# --trim_poly_x: trim poly-X sequences (where X is A, C, G or T) (useful for polyA in mRNA and polyG that appears in Illumina sequencing)
# --cut_front: cut the front of the read (useful as the extremes usually have a lower quality)
# --cut_tail: cut the tail of the read (useful as the extremes usually have a lower quality)
# --cut_mean_quality: cut the bases of the sliding window if the mean quality of that window is lower than the specified value (20 is the default value)
# --cut_window_size: size of the sliding window (4 is the default value)
# --length_required: minimum length of the read to be kept (30 is the default value), if the read is shorter than this value, it will be discarded
# -j: output file with the json report (by default it will be the name of the sample without the extension followed by " fastp_report.json")
# -h: output file with the html report (by default it will be the name of the sample without the extension followed by " fastp_report.html")
# --thread: number of threads to use (useful to save time)
# --report_title: title that will appear in the html report (by default, it will be the name of the sample without the extension followed by " fastp report")
#
#
# Requirements: having fastp installed.
# It can be insalled from Bioconda: "mamba install -c bioconda fastp"
# Link to fastp: https://github.com/OpenGene/fastp
#
# Input parameters and arguments:
# # - i: Input directory where the FASTQ files are located
# # - o: Output directory (optional, default: current directory)
# # - d: Use this option if you want to store the failed reads in a separate file and specify the path and name of that file (e.g., ./failed_reads.fastq.gz)
# # - r: Number of bases to cut from the 5' end of the read forward (optional, default 0)
# # - R: Number of bases to cut from the 5' end of the read reverse (optional, default 0)
# # - e: Number of bases to cut from the 3' end of the read forward (optional, default 0)
# # - E: Number of bases to cut from the 3' end of the read reverse (optional, default 0)
# # - x: Use this option if you want to trim poly-X sequences 
# # - f: Use this option if you want to trim from 5' end
# # - t: Use this option if you want to trim from 3' end
# # - q: Quality threshold for trimming (optional, default: 20)
# # - w: Size of the sliding window (optional, default: 4)
# # - l: Minimum length of the read to be kept (optional, default: 30)
# # - T: Number of threads to use (optional, default: 1)
# # - h: Help message (optional, displays the usage of the script)
# # - v: Version of the script (optional)
#
#
# Output format: fastq.gz files, .html report and .json report
#
#
# Usage: ./02_script.sh -i </path/to/input/directory> -o </path/to/output/directory> -d <name_of_file_to_store_failed_reads.fastq.gz> -r <number_of_bases_to_cut_in_read_forward_5'> -R <number_of_bases_to_cut_in_read_reverse_5'> -e <number_of_bases_to_cut_in_read_forward_3'> -E <number_of_bases_to_cut_in_read_reverse_3'> -x -f -t -q <quality_threshold> -w <size_of_sliding_window> -l <minimum_length_of_read> -T <number_of_threads>
#
# Usage example: ./02_script.sh -i /data/2025/grado_biotech/gabriel.coll/GN_01/scripts/00_raw_data/00_results -o /data/2025/grado_biotech/gabriel.coll/scripts/GN_01/02_trimming_filtering -d ./failed_reads.fastq.gz -r 10 -R 10 -e 10 -E 10 -x -f -t -q 20 -w 4 -l 30 -T 6
#
# Expected output: A fastq.gz file for each trimmed read, a .html and .json report for each sample and a log directory with the output and error log files of the general execution and of each sample
#
############################################################################################

# First, define the variables needed
OUTPUT_DIRECTORY=$(pwd)
FAILED_OUT=""
POLY_X=""
CUT_FRONT=""
CUT_TAIL=""
QUALITY=20
WINDOW_SIZE=4
MIN_LENGTH=30
THREADS=1
TRIM_FRONT1=0
TRIM_FRONT2=0
TRIM_TAIL1=0
TRIM_TAIL2=0

# Define the version of the script
VERSION="1.0.1"

# Define the help message that will display if necessary as a function
help_message() {
 echo -e "\nScript that trims the fastq.gz files of the RNA-seq using Fastp."
 echo -e "The output will be a fastq.gz file for each trimmed read and also a .html and .json reports for the results of each sample. The script will create a 02_results folder in the output directory to store all the results and, inside it, one folder for the results of each sample, named with the sample's name. It will also create a 02_logs folder in the output directory to store the logs of the script."
 echo -e "By default it will detect and eliminate automatically the adapters for paired-end reads, do not evaluate the duplication rate, cut the bases of the sliding window (with default size of 4 bases) if the mean quality of that window is lower than 20, and discard the reads shorter than 30 bases."
 echo -e "The user can then add more options to the trimming process and modify the default values of the parameters."
 echo -e "\nInput options and arguments:"
 echo -e "\t-i: Input directory where the FASTQ files are located"
 echo -e "\t-o: Output directory (optional, default: current directory)"
 echo -e "\t-d: Use this option if you want to store the failed reads in a separate file and specify the path and name of that file (optional, e.g., ./failed_reads.fastq.gz)"
 echo -e "\t-r: Number of bases to cut from the 5' end of the read forward (optional, default 0)"
 echo -e "\t-R: Number of bases to cut from the 5' end of the read reverse (optional, default 0)"
 echo -e "\t-e: Number of bases to cut from the 3' end of the read forward (optional, default 0)"
 echo -e "\t-E: Number of bases to cut from the 3' end of the read reverse (optional, default 0)"
 echo -e "\t-x: Use this option if you want to trim poly-X sequences (optional)"
 echo -e "\t-f: Use this option if you want to trim from 5' end (optional)"
 echo -e "\t-t: Use this option if you want to trim from 3' end (optional)"
 echo -e "\t-q: Quality threshold for trimming (optional, default: 20)"
 echo -e "\t-w: Size of the sliding window (optional, default: 4)"
 echo -e "\t-l: Minimum length of the read to be kept (optional, default: 30)"
 echo -e "\t-T: Number of threads to use (optional, default: 1)"
 echo -e "\t-h: Help message (optional, displays the usage of the script)"
 echo -e "\t-v: Version of the script (optional)\n"
 echo -e "\nUsage: ./02_script.sh -i </path/to/input/directory> -o </path/to/output/directory> -d <name_of_file_to_store_failed_reads.fastq.gz> -r <number_of_bases_to_cut_in_read_forward_5'> -R <number_of_bases_to_cut_in_read_reverse_5'> -e <number_of_bases_to_cut_in_read_forward_3'> -E <number_of_bases_to_cut_in_read_reverse_3'> -x -f -t -q <quality_threshold> -w <size_of_sliding_window> -l <minimum_length_of_read> -T <number_of_threads>\n"
 echo -e "\nUsage example: ./02_script.sh -i /data/2025/grado_biotech/gabriel.coll/GN_01/scripts/00_raw_data/00_results -o /data/2025/grado_biotech/gabriel.coll/scripts/GN_01/02_trimming_filtering -d ./failed_reads.fastq.gz -r 10 -R 10 -e 10 -E 10 -x -f -t -q 20 -w 4 -l 30 -T 6\n"
}


# Check if the user has provided the required arguments with getopts
# If they have, assign the arguments to the corresponding variables
# If not, print an error message and exit
{
while getopts "i:o:d:r:R:e:E:xftq:w:l:T:hv" opt; do
 case $opt in 
 i) INPUT_DIRECTORY=$OPTARG ;;
 o) OUTPUT_DIRECTORY=$OPTARG ;;
 d) FAILED_OUT="--failed_out ${OPTARG}" ;;
 r) TRIM_FRONT1="${OPTARG}" ;;
 R) TRIM_FRONT2="${OPTARG}" ;;
 e) TRIM_TAIL1="${OPTARG}" ;;
 E) TRIM_TAIL2="${OPTARG}" ;;
 x) POLY_X="--trim_poly_x" ;;
 f) CUT_FRONT="--cut_front" ;;
 t) CUT_TAIL="--cut_tail" ;;
 q) QUALITY=$OPTARG ;;
 w) WINDOW_SIZE=$OPTARG ;;
 l) MIN_LENGTH=$OPTARG ;;
 T) THREADS=$OPTARG ;;
 h) help_message ;;
 v) echo -e "\nVersion: $VERSION\n" ;;
 ?) echo -e "\nInvalid option or missing argument. Use -h for help.\n" >&2
    sleep 0.1s # Command to give bash enough time to close the pipe between the subshell with tee
    exit 1 ;;
 esac
done


# The only essential argument is the input directory, so we check if the user has provided it
if [[ -z "$INPUT_DIRECTORY" ]]; then
 sleep 0.01s # This time is used to print the stderror after printing the stdout, because if not, the error is printed faster than the output (for example, if using just the -v option)
 echo -e "\nYou have not provided the input directory. Use -h for help.\n" >&2
 sleep 0.1s # Again, time for bash to close the pipe with the subshell of tee before exiting
 exit 1
fi


# Check if the input_directory exists and is a directory
# If not, print an error message and exit
# If it exists, check if the user has the necessary permissions to enter and read the content
# Lastly, if it exists and the user can read it, check if the directory contains fastq.gz files
if [[ ! -d "$INPUT_DIRECTORY" ]]; then
 sleep 0.01s #Same reason, to print stdout before stderror and not the other way around. This sleeps will appear before every stderror message
 echo -e "\nInput directory does not exist or is not a directory.\n" >&2
 sleep 0.1s # This command again. It would be optional in case we had not created a { } block redirected to a subshell
 exit 1
elif [[ ! -r "$INPUT_DIRECTORY" ]] || [[ ! -x "$INPUT_DIRECTORY" ]]; then
 sleep 0.01s
 echo -e "\nYou do not have the necessary permissions to enter and read the directory with the fastq files." >&2
 echo -e "Please, ask to change the permissions of the directory and try again.\n" >&2
 sleep 0.1s
 exit 1
else
 FILES_FOUND=$(find "$INPUT_DIRECTORY" -name "*.fastq.gz" | wc -l) # Find the fastq.gz files in the input directory and count the number of lines given as output with wc -l. Command seen in some class examples. More info on: https://www.ionos.es/digitalguide/servidores/configuracion/comando-de-linux-wc/
 if [[ "$FILES_FOUND" -eq 0 ]]; then # If it was not able to find any fastq.gz file, print an error message and exit
  sleep 0.01s
  echo -e "\nThe input directory does not contain any fastq.gz file.\n" >&2
  sleep 0.1s
  exit 1
 fi
fi


# Check if the output directory exists and is a directory
# If not, create it
if [[ ! -d "$OUTPUT_DIRECTORY" ]]; then
  sleep 0.01s
  echo -e "\nOutput directory does not exist. Creating it..."
  mkdir -p "$OUTPUT_DIRECTORY" # Create recursively also the parent directories if they do not exist
  sleep 2 # Time for the user to read it
  # Check if it was successfully created
  if [[ $? -ne 0 ]]; then
   echo -e "\nError creating the output directory.\n" >&2
   exit 1
  else
   echo -e "\nOutput directory successfully created.\n"
   sleep 2
  fi
fi

# Check if the user has the necessary permissions to access and modify the output directory
if [ ! -r "$OUTPUT_DIRECTORY" ] || [ ! -w "$OUTPUT_DIRECTORY" ] || [ ! -x "$OUTPUT_DIRECTORY" ]; then
 sleep 0.01s
 echo -e "\nYou do not have the necessary permissions to access and modify the output directory." >&2
 echo -e "Please, change the permissions of the output directory and try again.\n" >&2
 sleep 0.1s
 exit 1
fi

# Once created, the results will be stored in a 02_results folder
if [[ ! -d "${OUTPUT_DIRECTORY}/02_results" ]]; then
 sleep 0.01s
 echo -e "\nCreating the 02_results folder..."
 mkdir -p "${OUTPUT_DIRECTORY}/02_results"
 sleep 2
 # Check if it was successfully created
 if [[ $? -ne 0 ]]; then
  echo -e "\nError creating the 02_results folder.\n" >&2
  exit 1
 else
  echo -e "\n02_results folder successfully created.\n"
 fi
fi

# Now the first general part of the script is over, so redirect the standard output and standard error to the corresponding log files
} 2> >(tee ./02_logs_general.err) > >(tee ./02_logs_general.out) # Command to redirect standard output and standard input to different files and print them on the screen at the same time.
# >( ) This structure is called process substitution, and it is used to create a temporary file from the subshell that can be used as input or output for a command. In other words, it gives the subshell a file descriptor and treats it as a file.
# This process receives the standard error (2>) and output (>) respectively, and gives them as input to the subshell, which redirects each one to a different file.


# Create the logs directory that will be used later on for redirections
# It will be used to redirect the standard output and standard error of the general part of the script (before using fastp)
# Now we move it to the output directory
mkdir -p "${OUTPUT_DIRECTORY}"/02_logs/02_logs_general
mv ./02_logs_general.err "${OUTPUT_DIRECTORY}"/02_logs/02_logs_general
mv ./02_logs_general.out "${OUTPUT_DIRECTORY}"/02_logs/02_logs_general


echo -e "\nAll the checks have been passed successfully. Now, we will start the trimming process.\n"

# Perform the trimming with fastp with each sample of the input directory
for file in "${INPUT_DIRECTORY}"/*_1.fastq.gz; do
 SAMPLE_NAME=$(basename "$file" _1.fastq.gz) # Extract the name of the sample without the extension
 FILE_FW=$file # Path to the forward read
 FILE_RV="${INPUT_DIRECTORY}/${SAMPLE_NAME}_2.fastq.gz" # Path to the reverse read
 # Create a log file for each sample, where the output and error will be redirected
 mkdir -p "${OUTPUT_DIRECTORY}/02_logs/02_logs_${SAMPLE_NAME}"
 LOG_SAMPLE="${OUTPUT_DIRECTORY}/02_logs/02_logs_${SAMPLE_NAME}"
 # Start the analysis
 {
 echo -e "\nProcessing sample: $SAMPLE_NAME"
 sleep 2 # To give some time the user to read the message
 echo "Trimming the reads of $SAMPLE_NAME with fastp..."
 sleep 2
 echo "Creating the reports for $SAMPLE_NAME..."
 sleep 2
 fastp \
  -i "$FILE_FW" \
  -I "$FILE_RV" \
  -o "${OUTPUT_DIRECTORY}/02_results/${SAMPLE_NAME}_1.fastq.gz" \
  -O "${OUTPUT_DIRECTORY}/02_results/${SAMPLE_NAME}_2.fastq.gz" \
  "$FAILED_OUT" \
  --detect_adapter_for_pe \
  --dont_eval_duplication \
  --trim_front1 "$TRIM_FRONT1"\
  --trim_front2 "$TRIM_FRONT2"\
  --trim_tail1 "$TRIM_TAIL1"\
  --trim_tail2 "$TRIM_TAIL2"\
  "$POLY_X" \
  "$CUT_FRONT" \
  "$CUT_TAIL" \
  --cut_mean_quality "$QUALITY" \
  --cut_window_size "$WINDOW_SIZE" \
  --length_required "$MIN_LENGTH" \
  -j "${OUTPUT_DIRECTORY}/02_results/${SAMPLE_NAME}_fastp_report.json" \
  -h "${OUTPUT_DIRECTORY}/02_results/${SAMPLE_NAME}_fastp_report.html" \
  --thread "$THREADS" \
  --report_title "${SAMPLE_NAME}_fastp_report" 
 # Check if the trimming process was carried out successfully
 if [[ $? -ne 0 ]]; then
   echo -e "\nError trimming $SAMPLE_NAME.\n" >&2
   exit 1
  else
   echo -e "\nThe trimming of $SAMPLE_NAME was finished successfully.\n"
   sleep 2 # To give the user time to read the message 
  fi
 } 2> >(tee "${LOG_SAMPLE}/02_logs_${SAMPLE_NAME}.err") > >(tee "${LOG_SAMPLE}/02_logs_${SAMPLE_NAME}.out") # Redirect the standard output and standard error to the corresponding log files
 done


# Tell the user the process has finished
sleep 3
echo -e "\nAll fastq.gz files have been trimmed.\n"


# End of the script
