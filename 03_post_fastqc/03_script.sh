#!/usr/bin/bash
###########################################################################################
#
# Authors: Nerea Barrio and Gabriel Coll
# Date: 2025-05-12
# Version: 1.0.1
# 
# Description: Script to perform a quality control of trimmed and processed reads from RNA-Seq sequence data, using FastQC and MultiQC as tools
# 
# Requirements:
## - Having FastQC installed (it can be installed from the zipped file available in https://github.com/s-andrews/FastQC?tab=readme-ov-file#) (a Java-running environment and Perl script are also needed).
## - Having MultiQC installed (it can be installed from Bioconda: "mamba install multiqc").
##### Link to FastQC: https://github.com/s-andrews/FastQC?tab=readme-ov-file#
##### Link to MultiQC: https://github.com/MultiQC/MultiQC
##### FastQC commands used in this script: 
###### - fastqc -o <output_directory> <input_file>: Performs the quality control of the reads and generates an html file with the results
##### MultiQC commands used in this script: 
###### - multiqc -o <output_directory> -n <report_name> <input_directory>: Generates a report combining all the results delivered by FastQC
# 
#
# Input parameters and arguments:
# # - i: Directory (path) with RNA-seq reads, in fastq format or fastq.gz
# # - o: Output directory (optional, default: current directory)
# # - h: Help message (optional, displays the usage of the script)
# # - v: Version of the script (optional)
#
# Output format: html file
# 
# Usage: ./03_script.sh -i <directory_with_reads> -o </path/to/output/directory>
#
# Usage example: ./03_script.sh -i /data/reads -o /home/user/output_directory
#
# Expected output: An html file with the FastQC report of each sample, an html file with the MultiQC report and a log directory with the output and error log files of the general execution and of each sample
#
#########################################################################

# Define the output directory as the current directory by default
OUTPUT_DIRECTORY=$(pwd)

# Define the file with the reads as an empty variable for later use with getopts
READS=""

# Define the version of the script
VERSION="1.0.1"

# Define the output selected variable as 0 by default
OUTPUT_SELECTED=0


# Check if the user has provided the required arguments with getopts
# If they have, assign the arguments to the corresponding variables
# If not, print an error message and exit
{
while getopts "i:o:hv" opt; do
 case $opt in
 i) READS=$OPTARG ;;
 o) OUTPUT_SELECTED=1  # This variable is used to check if the user has provided the output directory
    OUTPUT_DIRECTORY=$OPTARG ;; 
 v) echo -e "\nVersion: $VERSION\n" ;;
 h) echo -e "\nScript to perform a quality control of trimmed and processed reads from RNA-Seq sequence data, from NCBI, using FastQC and MultiQC as tools"
    echo -e "\nInput options and arguments:"
    echo -e "\t-i: Directory (path) with the RNA-seq reads, in fastq format"
    echo -e "\t-o: Output directory (optional, default: current directory)"
    echo -e "\t-h: Help message (optional, displays the usage of the script)"
    echo -e "\t-v: Version of the script (optional)\n"
    echo -e "\nUsage: $0 -i <file_with_reads> -o </path/to/output/directory>\n";;
 ?) echo -e "\nInvalid option or missing argument. Use -h for help.\n" >&2
    sleep 0.25s # Command to give bash enough time to close the pipe between the subshell with tee before exiting
    exit 1 ;;
  esac
done

# Check if the user has provided the fastq file with the reads
# If not, it checks if the user has provided an output directory
# If not, it just exits, because it means that they have provided the -h or -v option, or both, but not the file with the reads, so the user seems not to want to download anything
# In case they have provided the output directory but not the fastq file, it prints an error message and exits
if [[ "$READS" == "" ]]; then
 if [[ "$OUTPUT_SELECTED" -ne 1 ]]; then
  sleep 0.25s # This line is necessary to have something before "exit 0" so that bash has enough time to close the pipe and tee does not wait forever to receive an input
  exit 0
 else
  echo -e "\nYou have not provided the fastq file with the reads. Use -h for help.\n" >&2
  sleep 0.25s # Again, time for bash to close the pipe with the subshell of tee
  exit 1
  fi
fi

# Check if the fastq file exists, is not empty and with the correct extension (.fastq.gz))
# If not, print an error message and exit
# If it exists, check if the user has the necessary permissions to read the file
# Lastly, if it exists and the user can read it, check if the file is, indeed, a fastq file, without errors (including unusual or out-of-place characters and sequence-quality length mismatches)
for file in "$READS"/*.fastq.gz ; do
  symb_path=$(readlink -f "$file") # New command not learned in class to get the absolute path of the original file the link symbolic points to, in case that "file" is a symbolic link.
  if [[ ! -f "$symb_path" ]]; then 
   echo -e "\nRNA-seq read file does not exist or is not a file.\n" >&2
   sleep 0.25s # Again, this command. It would be optional in case we had not created a { } block redirected to a subshell
   exit 1
  elif [[ ! -s "$symb_path" ]]; then
   echo -e "\nThe file provided with the reads is empty.\n" >&2
   sleep 0.25s 
   exit 1
  elif [[ ! "$file" == *.fastq.gz ]]; then
    echo -e "\nThe extension of the file provided is not .fastq.gz.\n" >&2
    sleep 0.25s
    exit 1
  elif [[ ! -r "$file" ]]; then
    echo "File is not readable. Fixing..." >&2
    (chmod +r "$file" && echo "Permissions fixed") || (echo "Could not change permissions" >&2 && exit 1)
  else
   format_var=0
   line_num=0
   zcat $file | head -n 4 | while read line; do
   line_num=$((line_num + 1))
    if (( $line_num % 4 == 1 )); then
      seq_id=$line
      if [[ ! $line =~ ^@ ]]; then
        format_var=1
      fi
    elif (( $line_num % 4 == 2 )); then
      seq_length=${#line}
      if [[ ! $line =~ ^[ACGTNactgn]+$ ]]; then
        format_var=1
      fi
    elif (( $line_num % 4 == 3 )); then
      if [[ ! $line =~ ^\+ ]]; then
        format_var=1
      fi
    elif (( $line_num % 4 == 4 )); then
      quality_length=${#line}
      if ! LC_ALL=C grep -q '^[!-~]*$' <<< "$line"; then # Check if the quality line contains only printable ASCII characters. The variable LC_ALL=C is used to set the locale to C, which is the default locale and uses ASCII as the character encoding. The grep takes all ASCII characters from ! to ~ and checks if the quality line contains only these characters. <<< is a here-string, used to pass strings. The -q (quiet) option is used to prevent grep from showing the matching patterns.
        format_var=1
      elif [[ $seq_length -ne $quality_length ]]; then
        format_var=2
      fi
    fi
    case $format_var in
    1) echo -e "\nInvalid format: $line" >&2
     echo -e "\nPlease provide a valid fastq file.\n" >&2
     sleep 0.25s
     exit 1 ;;
    2) echo -e "Sequence and quality length do not match. ID:  $seq_id\n" >&2
     sleep 0.25s
     exit 1 ;;
    esac
   done 
  fi
done


# Check if the output directory exists and is a directory
# If not, create it
if [[ ! -d "$OUTPUT_DIRECTORY" ]]; then
  echo -e "\nOutput directory does not exist. Creating it..."
  sleep 2 # To give some time to the user to read the message
  mkdir -p "$OUTPUT_DIRECTORY" # Create recursively also the parent directories if they do not exist
  echo -e "Output directory successfully created.\n"
fi


# Check if the user has the necessary permissions to access and modify the output directory
if [ ! -r "$OUTPUT_DIRECTORY" ] || [ ! -w "$OUTPUT_DIRECTORY" ] || [ ! -x "$OUTPUT_DIRECTORY" ]; then
 echo -e "\nYou do not have the necessary permissions to access and modify the output directory." >&2
 echo -e "Please, change the permissions of the output directory and try again.\n" >&2
 sleep 0.25s
 exit 1
fi
} 2> >(tee ./03_logs.err) > >(tee  ./03_logs.out) # Command to redirect standard output and standard input to different files and print them on the screen at the same time.
# >( ) This structure is called process substitution, and it is used to create a temporary file from the subshell that can be used as input or output for a command. 
# In other words, it gives the subshell a file descriptor and treats it as a file.
# This process receives the standard output (>) and error (2>) respectively, and gives them as input to the subshell, which redirects each one to a different file.


# Create the logs directory that will be used later on for redirections
mkdir -p "${OUTPUT_DIRECTORY}"/03_logs/03_logs_general
# Move the logs to the output directory
mv ./03_logs.err "${OUTPUT_DIRECTORY}"/03_logs/03_logs_general
mv ./03_logs.out "${OUTPUT_DIRECTORY}"/03_logs/03_logs_general

# Now we start with the analysis itself

# Create separate directories for the FastQC and MultiQC analyses
FASTQC_OUTPUT="$OUTPUT_DIRECTORY/03_results/fastqc_results"
MULTIQC_OUTPUT="$OUTPUT_DIRECTORY/03_results/multiqc_results"

mkdir -p "$FASTQC_OUTPUT"
mkdir -p "$MULTIQC_OUTPUT"


# Now we will perform the FastQC analysis with each of the files located in the input directory
for file in "$READS"/*.fastq.gz; do
 SAMPLE_NAME=$(basename "$file" .fastq.gz) # Extract the name of the sample without the extension
 # Create a log file for each sample, where the output and error will be redirected
 FASTQC_LOGS="${OUTPUT_DIRECTORY}"/03_logs/03_fastqc_logs/03_logs_"${SAMPLE_NAME}"
 mkdir -p "${FASTQC_LOGS}"
 {
 # Start the analysis 
 echo -e "\nProcessing file: $SAMPLE_NAME"
 echo "Checking the quality of $SAMPLE_NAME with fastqc..."
 sleep 2 # To give some time to the user to read the message
 echo -e "Generating the html file...\n"
 sleep 2
 # Then, apply fastqc
 # Use FastQC to perform the quality control of the raw reads, generating an html file with the results
 fastqc -o "$FASTQC_OUTPUT" "$file"
 # Check if the html report was created successfully
 if [[ $? -ne 0 ]]; then
  echo -e "\nError creating html report of $SAMPLE_NAME.\n" >&2
  exit 1
 else
  echo -e "\nFastqc html report of $SAMPLE_NAME created successfully.\n"
 fi
} 2> >(tee "${FASTQC_LOGS}"/03_logs_"${SAMPLE_NAME}".err) > >(tee "${FASTQC_LOGS}"/03_logs_"${SAMPLE_NAME}".out) # Redirect the standard output and standard error to the corresponding log files
done

# Now perform the MultiQC analysis
# First, create the log file where the standout and stderror of the multiqc analysis will be stored
MULTIQC_LOG="${OUTPUT_DIRECTORY}"/03_logs/03_multiqc_logs
mkdir -p "${MULTIQC_LOG}"
{
# Tell the user what is going to happen
echo -e "\nAll the samples have been processed through FastQC."
sleep 2
echo -e "\nNow, we will generate a MultiQC report with all the results.\n"
sleep 2
echo -e "Generating...\n"
sleep 2
# Use MultiQC to generate a report with the results of the quality control
multiqc -o "$MULTIQC_OUTPUT" -n multiqc_report.html "$FASTQC_OUTPUT"
# Check if the multiqc report was created successfully
if [[ $? -ne 0 ]]; then
 echo -e "\nError creating the MultiQC report.\n" >&2
 exit 1
else
 echo -e "\nMultiQC report created successfully in "${MULTIQC_OUTPUT}"."
 echo -e "The process is over."
fi
} 2> >(tee "${MULTIQC_LOG}"/03_multiqc_logs.err) > >(tee "${MULTIQC_LOG}"/03_multiqc_logs.out) # Redirect the standard output and standard error to the corresponding log files


# End of the script
