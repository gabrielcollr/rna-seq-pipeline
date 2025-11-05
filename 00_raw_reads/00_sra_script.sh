#!/usr/bin/bash
###########################################################################################
#
# Authors: Nerea Barrio and Gabriel Coll
# Date: 2025-05-05
# Version: 1.0.1
# 
# Description: Script to download RNA-Seq sequence data from NCBI using the sra-tools utility with a file containing a list of the accession numbers of the SRA files (e.g., SRR1234567)
# 
# Requirements:
## - Having sra-tools installed (they can be installed from Bioconda: "mamba install -c bioconda sra-tools").
##### Link to sra-tools: https://github.com/ncbi/sra-tools 
##### Sra-tools commands used in this script: 
###### - prefetch <accession_number> -O <output_directory>: Downloads SRA files from NCBI in .sra format in the specified output directory
###### - fasterq-dump <sra_file> --outdir <output_directory>: Converts SRA files to FASTQ format and saves them in the specified output directory
# 
#
# Input parameters and arguments:
# # - a: File with the accessions number of the SRA files (e.g., SRR1234567) in each line
# # - o: Output directory (optional, default: current directory)
# # - h: Help message (optional, displays the usage of the script)
# # - v: Version of the script (optional)
#
# Output format: FASTQ file
# 
# Usage: ./00_sra_script.sh -a <file_with_accession_numbers> -o </path/to/output/directory>
#
# Usage example: ./00_sra_script.sh -a accession_file.txt -o /home/user/output_directory
#
# Expected output: A fastq file for each accession number in the input file, and a log directory with the output and error log files of the general execution and of each sample
#
#########################################################################

# Define the output directory as the current directory by default
OUTPUT_DIRECTORY=$(pwd)

# Define the accession number as an empty variable for later use with getopts
ACCESSION_NUMBERS=""

# Define the version of the script
VERSION="1.0.1"

# Define the output selected variable as 0 by default
OUTPUT_SELECTED=0


# Check if the user has provided the required arguments with getopts
# If they have, assign the arguments to the corresponding variables
# If not, print an error message and exit
{
while getopts "a:o:hv" opt; do
 case $opt in
 a) ACCESSION_NUMBERS=$OPTARG ;;
 o) OUTPUT_SELECTED=1  # This variable is used to check if the user has provided the output directory
    OUTPUT_DIRECTORY=$OPTARG ;; 
 v) echo $VERSION ;;
 h) echo -e "Script to download RNA-Seq sequence data from NCBI using the sra-tools utility with the accession number of the SRA files (e.g., SRR1234567)"
    echo -e "\nInput options and arguments:"
    echo -e "\t-a: File with the accessions number of the SRA files (e.g., SRR1234567) in each line"
    echo -e "\t-o: Output directory (optional, default: current directory)"
    echo -e "\t-h: Help message (optional, displays the usage of the script)"
    echo -e "\t-v: Version of the script (optional)\n"
    echo -e "\nUsage: $0 -a <file_with_accession_numbers.txt> -o </path/to/output/directory>";;
 ?) echo "Invalid option or missing argument. Use -h for help." >>&2
    exit 1 ;;
  esac
done


# Check if the user has provided the file with the accession numbers
# If not, it checks if the user has provided an output directory
# If not, it just exits, because it means that they have provided the -h or -v option, or both, but not the file with the accession numbers, so the user seems not to want to download anything
# In case they have provided the output directory but not the file with the accession numbers, it prints an error message and exits
if [[ "$ACCESSION_NUMBERS" == "" ]]; then
 if [[ "$OUTPUT_SELECTED" -ne 1 ]]; then
  exit 0
 else
  echo "You have not provided the file with the accession numbers. Use -h for help." >>&2
  exit 1
  fi
fi


# Check if the accession number file exists and is not empty
# If not, print an error message and exit
# If it exists, check if the user has the necessary permissions to read the file
# Lastly, if it exists and the user can read it, check if the accession numbers are valid SRA accession numbers
if [[ ! -f "$ACCESSION_NUMBERS" ]]; then 
 echo "Accession number file does not exist or is not a file." >>&2
 exit 1
elif [[ ! -s "$ACCESSION_NUMBERS" ]]; then
 echo "The file provided with the accession numbers is empty." >>&2
 exit 1
elif [[ ! -r "$ACCESSION_NUMBERS" ]]; then
 echo "You do not have the necessary permissions to read the file with the accession numbers." >>&2
 echo "Please, change the permissions of the file and try again." >>&2
 exit 1
else
 while read line; do
  if [[ ! $line =~ ^SRR[0-9]+$ ]]; then
   echo "Invalid accession number: $line" >>&2
   echo "Please provide valid SRA accession numbers (e.g., SRR1234567)." >>&2
   exit 1
  fi
 done < "$ACCESSION_NUMBERS"
fi


# Check if the output directory exists and is a directory
# If not, create it
if [[ ! -d "$OUTPUT_DIRECTORY" ]]; then
  echo "Output directory does not exist. Creating it..."
  mkdir -p "$OUTPUT_DIRECTORY" # Create recursively also the parent directories if they do not exist
  echo "Output directory successfully created."
fi


# Check if the user has the necessary permissions to access and modify the output directory
if [ ! -r "$OUTPUT_DIRECTORY" ] || [ ! -w "$OUTPUT_DIRECTORY" ] || [ ! -x "$OUTPUT_DIRECTORY" ]; then
 echo "You do not have the necessary permissions to access and modify the output directory." >>&2
 echo "Please, change the permissions of the output directory and try again." >>&2
 exit 1
fi


# Create a "logs/00_sra_logs_general" directory inside the output directory
# It will be used to redirect the standard output and standard error of the script
mkdir -p "${OUTPUT_DIRECTORY}"/00_sra_logs/00_sra_logs_general


# Now the first general part of the script is over, so redirect the standard output and standard error to the corresponding log files
} > >(tee -a "${OUTPUT_DIRECTORY}"/00_sra_logs/00_sra_logs_general/00_sra_logs.out) 2> >(tee -a "${OUTPUT_DIRECTORY}"/00_sra_logs/00_sra_logs_general/00_sra_logs.err) # This command was not exactly seen in class. 
# It is used to redirect the standard output and standard error to different files, and also to the terminal.
# >( ) This structure is called process substitution, and it is used to create a temporary file from the subshell that can be used as input or output for a command. 
# This process receives the standard output (>) and error (2>) respectively, and gives them as input to the subshell, which redirects each one to a different file. 


# Download the RNA-Seq sequences data from NCBI using sra-tools, going through the list of accession numbers line by line
{
while read line; do
 # Create a variable with the name of the accession number of the sample
 SAMPLE="$line"
 # Create a log file for each accession number, where the oitput and error will be redirected
 mkdir -p "${OUTPUT_DIRECTORY}"/00_sra_logs/00_sra_logs_"${SAMPLE}"
 # Start the download
 echo "Downloading $line RNA-Seq sequence data from NCBI..."
 prefetch "$line" -O "$OUTPUT_DIRECTORY" # Download the SRA file to the output directory in .sra format
 fasterq-dump "$OUTPUT_DIRECTORY"/"$line".sra --outdir "$OUTPUT_DIRECTORY" # Convert the SRA file to FASTQ format and save it in the output directory
 echo "Removing temporary files created by sra-tools ($OUTPUT_DIRECTORY/$line.sra)..." 
 rm "$OUTPUT_DIRECTORY"/"$line".sra # Remove the temporary .sra file created by sra-tools
 echo "Temporary file $OUTPUT_DIRECTORY/$line.sra successfully removed."
 echo "$line sequence data successfully downloaded and saved in $OUTPUT_DIRECTORY as fastq files."
done < "$ACCESSION_NUMBERS"

# Tell the user the process has finished
echo "End of the script. All sequence data has been downloaded and saved in $OUTPUT_DIRECTORY"
} > >(tee -a "${OUTPUT_DIRECTORY}"/00_sra_logs/00_sra_logs_"${SAMPLE}"/00_sra_logs_"${SAMPLE}".out) 2> >(tee -a "${OUTPUT_DIRECTORY}"/00_sra_logs/00_sra_logs_"${SAMPLE}"/00_sra_logs_"${SAMPLE}".err) # We use again this sintax


# End of the script

# Warning: This script has never been tested
