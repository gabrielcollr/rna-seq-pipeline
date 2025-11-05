#!/usr/bin/bash
###########################################################################################
#
# Authors: Nerea Barrio and Gabriel Coll
# Date: 2025-05-11
# Version: 1.0.1
#
# Description: Script that creates symbolic links to the compressed FASTQ files downloaded with the sra-toolkit by our classmates
#
# Requirements: None
#
# Input parameters and arguments:
# # - i: Input directory where the FASTQ files are located
# # - o: Output directory (optional, default: current directory)
# # - h: Help message (optional, displays the usage of the script)
# # - v: Version of the script (optional)
#
# Output format: fastq.gz files
#
# Usage: ./00_symlink_script.sh -i </path/to/input/directory> -o </path/to/output/directory>
#
# Usage example: ./00_symlink_script.sh -i /data/2025/grado_biotech/marta.taboada/SRV_01/00-raw_data -o /data/2025/grado_biotech/gabriel.coll/GN_01/00_raw_data
#
# Expected output: A symbolic link to each fastq.gz file in the output directory, and a log directory with the output and error log files of the general execution and of each sample
#
############################################################################################


# Define the output directory as the current directory by default
OUTPUT_DIRECTORY=$(pwd)


# Define the output selected variable as 0 by default
OUTPUT_SELECTED=0


# Define the input directory as an empty variable for later use
INPUT_DIRECTORY=""

# Define the version of the script
VERSION="1.0.1"

# Define the help message that will display if necessary as a function
help_message() {
  echo "Script that creates symbolic links to the compressed FASTQ files downloaded with the sra-toolkit by our classmates"
  echo -e "\nInput options and arguments:"
  echo -e "\t-i: Input directory where the FASTQ files are located"
  echo -e "\t-o: Output directory (optional, default: current directory)"
  echo -e "\t-h: Help message (optional, displays the usage of the script)"
  echo -e "\t-v: Version of the script (optional)\n"
  echo -e "\nUsage: ./00_symlink_script.sh -i </path/to/input/directory> -o </path/to/output/directory>\n"
}


# Use getopts to check the options provided by the user and create the corresponding variables
{
while getopts "i:o:hv" opt; do
 case $opt in
  i) INPUT_DIRECTORY=$OPTARG ;;
  o) OUTPUT_SELECTED=1  # This variable is used to check if the user has provided the output directory
     OUTPUT_DIRECTORY=${OPTARG} ;;
  v) echo -e "\nVersion: $VERSION\n" ;;
  h) echo ""
     help_message ;;
  ?) echo -e "\nInvalid option or missing argument. Use -h for help.\n" >&2
     sleep 0.25s # Command to give bash enough time to close the pipe between the subshell with tee
     exit 1 ;;
 esac
done


# Check if the user has provided the directory with the fastq.gz files
# If not, it checks if the user has provided an output directory
# If not, it just exits, because it means that they have provided the -h or -v option, or both, but not the directory with the fastq.gz files, so the user seems not to want to create the symbolic links
# In case they have provided the output directory but not the directory with the fastq.gz files, it prints an error message and exits
if [[ "$INPUT_DIRECTORY" == "" ]]; then
 if [[ "$OUTPUT_SELECTED" -eq 0 ]]; then
  sleep 0.25s # This line is necessary to have something before "exit 0" so that bash has enough time to close the pipe and tee does not wait forever to receive an input
  exit 0
 else
  echo -e "\nYou have not provided the directory with the fastq.gz files. Use -h for help.\n" >&2
  sleep 0.25s # Again, time for bash to close the pipe with the subshell of tee
  exit 1
  fi
fi


# Check if the input_directory exists and is a directory
# If not, print an error message and exit
# If it exists, check if the user has the necessary permissions to enter and read the content
# Lastly, if it exists and the user can read it, check if the directory contains fastq.gz files
if [[ ! -d "$INPUT_DIRECTORY" ]]; then
 echo -e "\nInput directory does not exist or is not a directory.\n" >&2
 sleep 0.25s # Again this command. It would be optional in case we had not created a { } block redirected to a subshell
 exit 1
elif [[ ! -r "$INPUT_DIRECTORY" ]] || [[ ! -x "$INPUT_DIRECTORY" ]]; then
 echo -e "\nYou do not have the necessary permissions to enter and read the directory with the fastq files." >&2
 echo "Please, ask to change the permissions of the directory and try again." >&2
 sleep 0.25s
 exit 1
else
 FILES_FOUND=$(find "$INPUT_DIRECTORY" -name "*.fastq.gz" | wc -l) # Find the fastq.gz files in the input directory and count the number of lines given as output with wc -l. Command seen in some class examples. More info on: https://www.ionos.es/digitalguide/servidores/configuracion/comando-de-linux-wc/
 if [[ "$FILES_FOUND" -eq 0 ]]; then # If it was not able to find any fastq.gz file, print an error message and exit
  echo -e "\nThe input directory does not contain any fastq.gz file.\n" >&2
  sleep 0.25s
  exit 1
 fi
fi


# Check if the output directory exists and is a directory
# If not, create it
if [[ ! -d "$OUTPUT_DIRECTORY" ]]; then
  echo -e "\nOutput directory does not exist. Creating it..."
  mkdir -p "$OUTPUT_DIRECTORY" # Create recursively also the parent directories if they do not exist
  sleep 1
  # Check if it was successfully created
  if [[ $? -ne 0 ]]; then
   echo -e "\nError creating the output directory.\n" >&2
   exit 1
  else
   echo -e "\nOutput directory successfully created.\n"
  fi
fi

# Once created, the results will be stored in a 00_results folder
if [[ ! -d "${OUTPUT_DIRECTORY}/00_results" ]]; then
 echo -e "\nCreating the 00_results folder..."
 mkdir -p "${OUTPUT_DIRECTORY}/00_results"
 sleep 1
 # Check if it was successfully created
 if [[ $? -ne 0 ]]; then
  echo -e "\nError creating the 00_results folder.\n" >&2
  exit 1
 else
  echo -e "\n00_results folder successfully created.\n"
 fi
fi


# Check if the user has the necessary permissions to access and modify the output directory
if [ ! -r "$OUTPUT_DIRECTORY" ] || [ ! -w "$OUTPUT_DIRECTORY" ] || [ ! -x "$OUTPUT_DIRECTORY" ]; then
 echo -e "\nYou do not have the necessary permissions to access and modify the output directory." >&2
 echo "Please, change the permissions of the output directory and try again.\n" >&2
 sleep 0.25s
 exit 1
fi
# Now the first general part of the script is over, so redirect the standard output and standard error to the corresponding log files
} 2> >(tee ./00_logs.err) > >(tee ./00_logs.out) # Command to redirect standard output and standard input to different files and print them on the screen at the same time.
# >( ) This structure is called process substitution, and it is used to create a temporary file from the subshell that can be used as input or output for a command. In other words, it gives the subshell a file descriptor and treats it as a file.
# This process receives the standard output (>) and error (2>) respectively, and gives them as input to the subshell, which redirects each one to a different file.


# Create the logs directory that will be used for redirections
# Create a "logs/00_logs_general" directory inside the output directory
# It will be used to move the standard output and standard error of the script
mkdir -p "${OUTPUT_DIRECTORY}"/00_logs/00_logs_general
mv ./00_logs.err "${OUTPUT_DIRECTORY}"/00_logs/00_logs_general
mv ./00_logs.out "${OUTPUT_DIRECTORY}"/00_logs/00_logs_general

# Create symbolic links to the fastq.gz files in the output directory
# If the symbolic link already exists, it will be overwritten
# First, we extract the name of each fastq.gz file and create a log file for each one
for file in "$INPUT_DIRECTORY"/*.fastq.gz; do
 SAMPLE_NAME=$(basename "$file" .fastq.gz) # Extract the name of the sample without the extension
 NAME=$(basename "$file") # Extract the name of the file with the extension to create the symbolic link
 mkdir -p "${OUTPUT_DIRECTORY}/00_logs/00_logs_${SAMPLE_NAME}" # Create a log file for each sample, where the output and error will be redirected
 {
 echo -e "\nProcessing file: $NAME"
 echo "Creating symbolic link to $NAME..."
 sleep 2 # To give some time to the user to read the message
 # Then, we create the symbolic link
 ln -sf "$file" "$OUTPUT_DIRECTORY"/00_results/"$NAME"
 # Check if the symbolic link was created successfully
 if [[ $? -ne 0 ]]; then
  echo -e "\nError creating symbolic link to $NAME.\n" >&2
  exit 1
 else
  echo -e "\nSymbolic link to $NAME successfully created.\n"
  sleep 2 # To give the user time to read the message 
 fi
} 2> >(tee "${OUTPUT_DIRECTORY}"/00_logs/00_logs_"${SAMPLE_NAME}"/00_logs_"${SAMPLE_NAME}".err) > >(tee "${OUTPUT_DIRECTORY}"/00_logs/00_logs_"${SAMPLE_NAME}"/00_logs_"${SAMPLE_NAME}".out) # Redirect the standard output and standard error to the corresponding log files
done


# Tell the user the process has finished
sleep 3
echo -e "\nAll symbolic links have been created successfully in $OUTPUT_DIRECTORY.\n"


# End of the script
