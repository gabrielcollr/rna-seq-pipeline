#!/usr/bin/bash
###########################################################################################
#
# Authors: Nerea Barrio and Gabriel Coll
# Date: 2025-05-31
# Version: 1.0.1
# 
#
# Description: Script to perform an alignment of trimmed and processed reads from RNA-Seq sequence data, using the tool HISAT2, and a later quality control of the results with MultiQC as tool 
# 
# Requirements:
## - Having downloaded de reference genome and annotation files (GTF)
## - Having HISAT2 installed (it can be installed from Bioconda: "mamba install -c bioconda hisat2").
##### Link to Hisat2 for more information: https://github.com/DaehwanKimLab/hisat2
##### Hisat2 commands used in this script: 
###### - hisat2-build <reference_genome> <genome_index_prefix> # Creates the genome index
###### - hisat2 -x <genome_index> -1 <reads_1.fq> -2 <read2_2.fq> -S output.sam --summary-file <summary_text_file.txt> # Performs the alignment
## - Having MultiQC installed (it can be installed from Bioconda: "mamba install -c bioconda multiqc").
##### Link to MultiQC: https://github.com/MultiQC/MultiQC
##### MultiQC commands used in this script:
###### - multiqc -o <output_directory> -n <report_name> <input_directory>: Generates a report combining all the summary results delivered by Hisat2
## - Having samtools installed (it can be installed from Bioconda: "mamba install -c bioconda samtools").
##### Samtools is used to convert and sort the sam files to bam files before passing them to multiqc.
##### Link to samtools: https://github.com/samtools/samtools
##### Samtools commands used in this script:
###### - samtools view -b <input.sam> -o <output.bam>  # Converts the sam files to bam files
###### - samtools sort <input.bam> -o <output.sorted.bam> # Sorts the bam files by genome position
###### - samtools index <input.sorted.bam> # Creates the index for the sorted bam files
###### - samtools flagstat <input.sorted.bam> > <output.txt> # Generates a text file with the statistics of the alignment of the sorted bam files
###### - samtools idxstats <input.sorted.bam> > <output.txt> # Generates a text file with the statistics of the sorted bam files by chromosomes
#
#
# Input parameters and arguments:
# # - i: Directory (path) with trimmed RNA-seq reads, in fastq.gz format
# # - o: Output directory (optional, default: current directory)
# # - g: Reference genome file (FASTA format)
# # - a: Annotation file (GTF format)
# # - h: Help message (optional, displays the usage of the script)
# # - v: Version of the script (optional)
#
#
# Output format: sam files with the alignment results, summary text file of the alignment and html file with the MultiQC report
# 
#
# Usage: ./04_script.sh -i <directory_with_reads> -o </path/to/output/directory> -g <reference_genome> -a <annotation_file>
#
# Usage example: ./04_script.sh -i /data/reads -o /home/user/output_directory -g /data/reference_genome.fa -a /data/annotation.gtf
#
# Expected output: A sam file for each sample, a summary text file for each sample, a .html file with the MultiQC report and a log directory with the output and error log files of the general execution and of each sample
#
#########################################################################

# Define the output directory as the current directory by default 
OUTPUT_DIRECTORY=$(pwd)
# Define the rest of variables that will be used as empty
INPUT_DIRECTORY=""
REFERENCE_GENOME=""
ANNOTATION_FILE=""
# Define a variable that will check if the user has provided at least one option
OPTIONS_PROVIDED=0

# Define the version of the script
VERSION="1.0.1"

# Define the help message that will display if necessary as a function
help_message() {
 echo "Script to perform an alignment of trimmed and processed reads from RNA-Seq sequence data, using the tool HISAT2, and a later quality control of the results with MultiQC as tool"
 echo -e "\nInput options and arguments:"
 echo -e "\t-i: Directory (path) with trimmed RNA-seq reads, in fastq.gz format"
 echo -e "\t-o: Output directory (optional, default: current directory)"
 echo -e "\t-g: Reference genome file (FASTA format)"
 echo -e "\t-a: Annotation file (GTF format)"
 echo -e "\t-h: Help message (optional, displays the usage of the script)"
 echo -e "\t-v: Version of the script (optional)\n"
 echo -e "\nUsage: ./04_script.sh -i <directory_with_reads> -o </path/to/output/directory> -g <reference_genome> -a <annotation_file>\n"
 echo -e "\nUsage example: ./04_script.sh -i /data/reads -o /home/user/output_directory -g /data/reference_genome.fa -a /data/annotation.gtf\n"
}


# Use getopts to check the options provided by the user and create the corresponding variables
{
while getopts "i:o:g:a:hv" opt; do
 case $opt in
  i) INPUT_DIRECTORY=$OPTARG ;;
  o) OUTPUT_DIRECTORY=$OPTARG ;;
  g) REFERENCE_GENOME=$OPTARG ;;
  a) ANNOTATION_FILE=$OPTARG ;;
  v) echo -e "\nVersion: $VERSION\n" ;;
  h) echo ""
     help_message ;;
  ?) echo -e "\nInvalid option or missing argument. Use -h for help.\n" >&2
     sleep 0.1s # Command to give bash enough time to close the pipe between the subshell with tee
     exit 1 ;;
 esac
 OPTIONS_PROVIDED=1 # Check something was provided as option
done


# Now we start with the different check points

# Check if the user has provided any option or not
if [[ OPTIONS_PROVIDED -eq 0 ]]; then
 sleep 0.01s # This time is used to print the stderr after the stdout (for example if the -v option was given). By default the error is faster, so we need to slow it down a little bit
 echo -e "\nYou have not provided any option. Use -h for help.\n" >&2
 sleep 0.1s # Again, time for bash to close the pipe with the subshell of tee
 exit 1
fi


# Check if the user has provided the arguments needed
# To do so, if any of the variables INPUT_DIRECTORY, REFERENCE_GENOME or ANNOTATION_FILE is not empty, it checks if the rest of them are also not empty or something is missing
# If any of them is empty, it prints an error message and exits
# If all of them are not empty, it will continue and start the analysis
# If all of them are empty, it will just exit, because it means that they have provided the -h or -v option, or both, but not the arguments needed, so the user seems not to want to start the alignment
if [[ -n "$INPUT_DIRECTORY" || -n "$REFERENCE_GENOME" || -n "$ANNOTATION_FILE" ]]; then
 if [[ -z "$INPUT_DIRECTORY" ]]; then
  sleep 0.01s # Same as before, it will be used before any error message
  echo -e "\nYou have not provided the input directory with the fastq.gz files." >&2
  echo -e "Use -h for help and try again.\n" >&2
  sleep 0.1s # Again, time for bash to close the pipe with the subshell of tee
  exit 1
 elif [[ -z "$REFERENCE_GENOME" ]]; then
  sleep 0.01s # Same reason as before
  echo -e "\nYou have not provided the file with the reference genome." >&2
  echo -e "Use -h for help and try again.\n" >&2
  sleep 0.1s # This time will be added after any error message during this check points
  exit 1
 elif [[ -z "$ANNOTATION_FILE" ]]; then
  sleep 0.01s
  echo -e "\nYou have not provided the annotation file." >&2
  echo -e "Use -h for help and try again.\n" >&2
  sleep 0.1s
  exit 1
 fi
else
 sleep 0.1s # This line is necessary in order to have something before "exit 0" so that bash has enough time to close the pipe and tee does not wait forever to receive an input
 exit 0
fi


# Check if the input_directory exists and is a directory
# If not, print an error message and exit
# If it exists, check if the user has the necessary permissions to enter and read the content
# Lastly, if it exists and the user can read it, check if the directory contains fastq.gz files
if [[ ! -d "$INPUT_DIRECTORY" ]]; then
 sleep 0.01s
 echo -e "\nInput directory does not exist or is not a directory.\n" >&2
 sleep 0.1s # Again this command. It would be optional in case we had not created a { } block redirected to a subshell
 exit 1
elif [[ ! -r "$INPUT_DIRECTORY" ]] || [[ ! -x "$INPUT_DIRECTORY" ]]; then
 sleep 0.01s
 echo -e "\nYou do not have the necessary permissions to enter and read the directory with the fastq files." >&2
 echo "Please, ask to change the permissions of the directory and try again." >&2
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


# Check if the user has the necessary permissions to read the reference genome and annotation files
# First, check if the user can enter and read the directories where the files are located and then the files
# If he/she does, check if the files exist and are in the correct format
# If not, print an error message and exit
DIRECTORY_REF_GENOME=$(dirname "$REFERENCE_GENOME")
DIRECTORY_ANNOTATION=$(dirname "$ANNOTATION_FILE")
# Check if the directory with the reference genome exists and is a directory
if [[ ! -d "$DIRECTORY_REF_GENOME" ]]; then
 sleep 0.01s
 echo -e "\nThe directory with the reference genome file does not exist or is not a directory.\n" >&2
 echo -e "Please, provide a valid directory and try again.\n" >&2
 sleep 0.1s
 exit 1
# If it exists, check permissions to enter and read it
elif [[ ! -r "$DIRECTORY_REF_GENOME" || ! -x "$DIRECTORY_REF_GENOME" ]]; then
 sleep 0.01s
 echo -e "\nYou do not have the necessary permissions to enter and read the directory with the reference genome file." >&2
 echo -e "Please, ask to change the permissions of the directory and try again." >&2
 sleep 0.1s
 exit 1
# Check if the reference genome file exists
elif [[ ! -e "$REFERENCE_GENOME" ]]; then
 sleep 0.01s
 echo -e "\nThe reference genome file does not exist in the path provided.\n" >&2
 echo -e "Please, provide a valid path and try again.\n" >&2
 sleep 0.1s
 exit 1
# Check if the reference genome file is in the correct format
elif [[ ! "$REFERENCE_GENOME" == *.fa && ! "$REFERENCE_GENOME" == *.fasta ]]; then
 sleep 0.01s
 echo -e "\nThe reference genome file is not in the correct format." >&2
 echo -e "Please, provide a file in FASTA format and try again.\n" >&2
 sleep 0.1s
 exit 1
# Check permission to read the reference genome file
elif [[ ! -r "$REFERENCE_GENOME" ]]; then
 sleep 0.01s
 echo -e "\nYou do not have the necessary permissions to read the reference genome file." >&2
 echo -e "Please, ask to change the permissions of the file and try again." >&2
 sleep 0.1s
 exit 1
fi


# Now it repeats the exact same procedure with the annotation file
# Check if the directory with the annotation file exists and is a directory
if [[ ! -d "$DIRECTORY_ANNOTATION" ]]; then
 sleep 0.01s
 echo -e "\nThe directory with the annotation file does not exist or is not a directory.\n" >&2
 echo -e "Please, provide a valid directory and try again.\n" >&2
 sleep 0.1s
 exit 1
# If it exists, check permissions to enter and read it
elif [[ ! -r "$DIRECTORY_ANNOTATION" || ! -x "$DIRECTORY_ANNOTATION" ]]; then
 sleep 0.01s
 echo -e "\nYou do not have the necessary permissions to enter and read the directory with the annotation file." >&2
 echo -e "Please, ask to change the permissions of the directory and try again." >&2
 sleep 0.1s
 exit 1
# Check if the annotation file exists
elif [[ ! -e "$ANNOTATION_FILE" ]]; then
 sleep 0.01s
 echo -e "\nThe annotation file does not exist in the path provided.\n" >&2
 echo -e "Please, provide a valid path and try again.\n" >&2
 sleep 0.1s
 exit 1
# Check if the annotation file is in the correct format
elif [[ ! "$ANNOTATION_FILE" == *.gtf && ! "$ANNOTATION_FILE" == *.gff && ! "$ANNOTATION_FILE" == *.gff3 ]]; then
 sleep 0.01s
 echo -e "\nThe annotation file is not in the correct format." >&2
 echo -e "Please, provide a file in GTF format and try again.\n" >&2
 sleep 0.1s
 exit 1
# Check permission to read the reference genome file
elif [[ ! -r "$ANNOTATION_FILE" ]]; then
 sleep 0.01s
 echo -e "\nYou do not have the necessary permissions to read the annotation file." >&2
 echo -e "Please, ask to change the permissions of the file and try again." >&2
 sleep 0.1s
 exit 1
fi

# Check if the output directory exists and is a directory
# If not, create it
if [[ ! -d "$OUTPUT_DIRECTORY" ]]; then
  sleep 0.01s
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


# Check if the user has the necessary permissions to access and modify the output directory
if [ ! -r "$OUTPUT_DIRECTORY" ] || [ ! -w "$OUTPUT_DIRECTORY" ] || [ ! -x "$OUTPUT_DIRECTORY" ]; then
 sleep 0.01s
 echo -e "\nYou do not have the necessary permissions to access and modify the output directory." >&2
 echo -e "Please, change the permissions of the output directory and try again.\n" >&2
 sleep 0.1s
 exit 1
fi


# Once created, the results will be stored in a 04_results folder
if [[ ! -d "${OUTPUT_DIRECTORY}/04_results" ]]; then
 echo -e "\nCreating the 04_results folder..."
 mkdir -p "${OUTPUT_DIRECTORY}/04_results"
 sleep 1
 # Check if it was successfully created
 if [[ $? -ne 0 ]]; then
  echo -e "\nError creating the 04_results folder.\n" >&2
  exit 1
 else
  echo -e "\n04_results folder successfully created.\n"
 fi
fi

# The output of the alignment and multiqc analysis will be stored in different folders
# First we name them and then we create them
ALIGNMENT_OUTPUT="$OUTPUT_DIRECTORY/04_results/04_alignment_results"
MULTIQC_OUTPUT="$OUTPUT_DIRECTORY/04_results/04_multiqc_results"
echo -e "\nCreating the 04_alignment_results and 04_multiqc_results folders..."
sleep 1
mkdir -p "$ALIGNMENT_OUTPUT"
mkdir -p "$MULTIQC_OUTPUT"

# Now it starts with the alignment process
echo -e "\nAll the checks have been passed successfully. Now, we will start the alignment process.\n"
sleep 2 # Time for the user to read it

# Create the genome index that Hisat2 will use for the alignment
# To do so, we will use the hisat2-build command
# The genome index will be stored in a directory named genome_index inside the 04_results folder
echo -e "\nCreating the genome index...\n"
mkdir -p "${ALIGNMENT_OUTPUT}"/genome_index
# Check if it was already created
GENOME_INDEX=$(find "${ALIGNMENT_OUTPUT}"/genome_index -name "*.ht2" | wc -l) # Find the .ht2 files in the genome_index directory and count the number of lines given as output with wc
# If it was not, create it
if [[ "$GENOME_INDEX" -eq 0 ]]; then
 hisat2-build "$REFERENCE_GENOME" "${ALIGNMENT_OUTPUT}"/genome_index/genome_index
 # Check if it was successfully created
 if [[ $? -ne 0 ]]; then
  echo -e "\nError creating the genome index.\n" >&2
  exit 1
 else
  echo -e "\nGenome index successfully created in "${ALIGNMENT_OUTPUT}"/genome_index.\n"
  sleep 2 # Time for the user to read it
 fi
else
 echo -e "\nThe genome index already exists from previous analysis in the output directory.\n"
 echo -e "We are not going to create it again.\n"
 echo -e "Instead, we will use the existing one.\n"
 sleep 2 # Time for the user to read it 
fi

# The first general part of the script is over, so redirect the standard output and standard error to the corresponding log files
 } 2> >(tee ./04_logs.err) > >(tee ./04_logs.out) # Command to redirect standard output and standard input to different files and print them on the screen at the same time.
# >( ) This structure is called process substitution, and it is used to create a temporary file from the subshell that can be used as input or output for a command. 
# In other words, it gives the subshell a file descriptor and treats it as a file.
# This process receives the standard output (>) and error (2>) respectively, and gives them as input to the subshell, which redirects each one to a different file.


# Create the logs directory that will be used for redirections
# Create a "logs/04_logs_general" directory inside the output directory
# It will be used to move the standard output and standard error of the script
mkdir -p "${OUTPUT_DIRECTORY}"/04_logs/04_logs_general
mv ./04_logs.err "${OUTPUT_DIRECTORY}"/04_logs/04_logs_general
mv ./04_logs.out "${OUTPUT_DIRECTORY}"/04_logs/04_logs_general

# Now we will perform the alignment with each of the files located in the input directory
echo -e "\nNow, we will start the alignment process with each of the samples.\n"
sleep 2 # Time for the user to read it
for file in "$INPUT_DIRECTORY"/*_1.fastq.gz; do
 SAMPLE_NAME=$(basename "$file" _1.fastq.gz) # Extract the name of the sample without the extension
 FILE_FW=$file # Path to the forward read
 FILE_RV="${INPUT_DIRECTORY}/${SAMPLE_NAME}_2.fastq.gz" # Path to the reverse read
 # Create a log file for each sample, where the output and error will be redirected
 mkdir -p "${OUTPUT_DIRECTORY}"/04_logs/04_alignment_logs/04_logs_"${SAMPLE_NAME}"
 LOG_SAMPLE="${OUTPUT_DIRECTORY}"/04_logs/04_alignment_logs/04_logs_"${SAMPLE_NAME}"
 # Start the analysis
 {
 echo -e "\nProcessing sample: $SAMPLE_NAME"
 sleep 1 # To give some time the user to read the message
 echo "Aligning the reads of $SAMPLE_NAME with Hisat2..."
 sleep 1
 echo "Creating the summary file for $SAMPLE_NAME..."
 sleep 1
 hisat2 \
  -x "${ALIGNMENT_OUTPUT}"/genome_index/genome_index \
  -1 "$FILE_FW" \
  -2 "$FILE_RV" \
  -S "${ALIGNMENT_OUTPUT}"/"${SAMPLE_NAME}".sam \
  --summary-file "${ALIGNMENT_OUTPUT}"/"${SAMPLE_NAME}"_summary.txt
 # Check if the alignment process was carried out successfully
 if [[ $? -ne 0 ]]; then
  echo -e "\nError aligning $SAMPLE_NAME.\n" >&2
  exit 1
 else
  echo -e "\nThe alignment of $SAMPLE_NAME finished successfully.\n"
  sleep 1 # To give the user time to read the message
 fi
 # Now we will convert the sam file to bam file and sort it
 echo -e "\nConverting the sam file to bam file...\n"
 samtools view -b "${ALIGNMENT_OUTPUT}"/"${SAMPLE_NAME}".sam -o "${ALIGNMENT_OUTPUT}"/"${SAMPLE_NAME}".bam
  # Check if the conversion process was carried out successfully
 if [[ $? -ne 0 ]]; then
  echo -e "\nError converting $SAMPLE_NAME.\n" >&2
  exit 1
 else
  echo -e "\nThe conversion of $SAMPLE_NAME finished successfully.\n"
  sleep 1 # To give the user time to read the message
 fi
 # Sort the bam file by genome position
 echo -e "\nSorting the bam file...\n"
 samtools sort "${ALIGNMENT_OUTPUT}"/"${SAMPLE_NAME}".bam -o "${ALIGNMENT_OUTPUT}"/"${SAMPLE_NAME}"_sorted.bam
  # Check if the sorting process was carried out successfully
 if [[ $? -ne 0 ]]; then
  echo -e "\nError sorting $SAMPLE_NAME.\n" >&2
  exit 1
 else
  echo -e "\nSorting of $SAMPLE_NAME finished successfully.\n"
  sleep 1 # To give the user time to read the message
 fi
 # Create the index for the sorted bam file
 echo -e "\nCreating the index for the sorted bam file...\n"
 samtools index "${ALIGNMENT_OUTPUT}"/"${SAMPLE_NAME}"_sorted.bam
  # Check if the index was created successfully
 if [[ $? -ne 0 ]]; then
  echo -e "\nError creating the index in sample $SAMPLE_NAME.\n" >&2
  exit 1
 else
  echo -e "\nIndex of bam file of $SAMPLE_NAME created successfully.\n"
  sleep 1 # To give the user time to read the message
 fi
 # Create the flagstat file for the sorted bam file (summary of the bam file)
 echo -e "\nCreating the flagstat file for the sorted bam file...\n"
 samtools flagstat "${ALIGNMENT_OUTPUT}"/"${SAMPLE_NAME}"_sorted.bam > "${ALIGNMENT_OUTPUT}"/"${SAMPLE_NAME}"_flagstat.txt
 # Check if the flagstat file was created successfully
 if [[ $? -ne 0 ]]; then
  echo -e "\nError creating the flagstat file in sample $SAMPLE_NAME.\n" >&2
  exit 1
 else
  echo -e "\nFlagstat file of $SAMPLE_NAME created successfully.\n"
  sleep 1 # To give the user time to read the message
 fi
 # Create the idxstats file for the sorted bam file (summary of the bam file by chromosomes))
 echo -e "\nCreating the idxstats file for the sorted bam file...\n"
 samtools idxstats "${ALIGNMENT_OUTPUT}"/"${SAMPLE_NAME}"_sorted.bam > "${ALIGNMENT_OUTPUT}"/"${SAMPLE_NAME}"_idxstats.txt
 # Check if the idxstats file was created successfully
 if [[ $? -ne 0 ]]; then
  echo -e "\nError creating the idxstats file in sample $SAMPLE_NAME.\n" >&2
  exit 1
 else
  echo -e "\nIdxstats file of $SAMPLE_NAME created successfully.\n"
  sleep 1 # To give the user time to read the message
 fi
 # Delete the sam file to save space
 echo -e "\nDeleting the sam file to save space...\n"
 rm "${ALIGNMENT_OUTPUT}"/"${SAMPLE_NAME}".sam
 # Check if the sam file was deleted successfully
 if [[ $? -ne 0 ]]; then
  echo -e "\nError deleting the sam file in sample $SAMPLE_NAME.\n" >&2
  exit 1
 else
  echo -e "\nSam file of $SAMPLE_NAME deleted successfully.\n"
  sleep 1 # To give the user time to read the message
 fi
 # Delete also the bam file not sorted to save space
 echo -e "\nDeleting the bam file not sorted to save space...\n"
 rm "${ALIGNMENT_OUTPUT}"/"${SAMPLE_NAME}".bam
 # Check if the bam file was deleted successfully
 if [[ $? -ne 0 ]]; then
  echo -e "\nError deleting the bam file in sample $SAMPLE_NAME.\n" >&2
  exit 1
 else
  echo -e "\nBam file of $SAMPLE_NAME deleted successfully.\n"
  sleep 1 # To give the user time to read the message
 fi
 # Change the name of the _sorted.bam file to .bam to make it easier to use in the next script
 echo -e "\nChanging the name of the _sorted.bam file to .bam...\n"
 mv "${ALIGNMENT_OUTPUT}"/"${SAMPLE_NAME}"_sorted.bam "${ALIGNMENT_OUTPUT}"/"${SAMPLE_NAME}".bam
 # Check if the name was changed successfully
 if [[ $? -ne 0 ]]; then
  echo -e "\nError changing the name of the bam file in sample $SAMPLE_NAME.\n" >&2
  exit 1
 else
  echo -e "\nName of bam file of $SAMPLE_NAME changed successfully.\n"
  sleep 1 # To give the user time to read the message
 fi
 # Also, change the name of the _sorted.bam.bai file to .bam.bai to make it easier to use in the next script
 echo -e "\nChanging the name of the _sorted.bam.bai file to .bam.bai...\n"
 mv "${ALIGNMENT_OUTPUT}"/"${SAMPLE_NAME}"_sorted.bam.bai "${ALIGNMENT_OUTPUT}"/"${SAMPLE_NAME}".bam.bai
 # Check if the name was changed successfully
 if [[ $? -ne 0 ]]; then
  echo -e "\nError changing the name of the bam.bai file in sample $SAMPLE_NAME.\n" >&2
  exit 1
 else
  echo -e "\nName of bam.bai file of $SAMPLE_NAME changed successfully.\n"
  sleep 1 # To give the user time to read the message
 fi
 } 2> >(tee "${LOG_SAMPLE}"/04_logs_"${SAMPLE_NAME}".err) > >(tee "${LOG_SAMPLE}"/04_logs_"${SAMPLE_NAME}".out) # Redirect the standard output and standard error to the corresponding log files
done




# Now perform the MultiQC analysis
# First, create the log file where the standout and stderror of the multiqc analysis will be stored
MULTIQC_LOG="${OUTPUT_DIRECTORY}"/04_logs/04_multiqc_logs
mkdir -p "${MULTIQC_LOG}"
{
# Tell the user what is going to happen
echo -e "\nAll the samples have been successfullt aligned through Hisat2."
sleep 2
echo -e "\nNow, we will generate a MultiQC report with all the results.\n"
sleep 2
echo -e "Generating...\n"
sleep 2
# Use MultiQC to generate a report with the results of the quality control
multiqc -o "$MULTIQC_OUTPUT" -n multiqc_report.html "$ALIGNMENT_OUTPUT"
# Check if the multiqc report was created successfully
if [[ $? -ne 0 ]]; then
 echo -e "\nError creating the MultiQC report.\n" >&2
 exit 1
else
 echo -e "\nMultiQC report created successfully in "${MULTIQC_OUTPUT}"."
 echo -e "The process is over."
fi
} 2> >(tee "${MULTIQC_LOG}"/04_multiqc_logs.err) > >(tee "${MULTIQC_LOG}"/04_multiqc_logs.out) # Redirect the standard output and standard error to the corresponding log files


echo -e "\nEverything is finished.\n"


# End of the script
