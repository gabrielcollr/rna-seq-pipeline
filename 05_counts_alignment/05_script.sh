#!/usr/bin/bash
###########################################################################################
#
# Authors: Nerea Barrio and Gabriel Coll
# Date: 2025-05-28
# Version: 1.0.1
# 
# Description: Script to perform a quantification of reads per gene from RNA-Seq sequencing data, which has been previously processed (trimmed and aligned), using featureCounts as tool and the annotation file of the genome of interest.
# 
# Requirements:
## - Having featureCounts installed (it can be installed from Bioconda via the subreads package: "mamba install -c bioconda subread").

##### Link to featureCounts: https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html
##### FeatureCounts commands used in this script: 
###### - featureCounts -o <output_file> -a <reference_genome_annotation_file> <input_file_1> <input_file_2>...: Performs the read count and generates a tsv file with the results of all samples
###### Other parameters are displayed exactly the same way in the command line of this script (explained below)
# 
#
# Input parameters and arguments:
# # - i: Directory (path) with processed RNA-seq data, in BAM or SAM format
# # - o: Output directory (optional, default: current directory)
# # - a: Reference genome annotation file (GTF or GFF3 format)
# # - T: Followed by the number of threads to use (optional, default: 1)
# # - O: Assigns reads to all their overlapping meta-features
# # - p: Paired-end reads (optional, default: single-end reads)
# # - s: Strandness of the data, followed by 0 (unstranded reads), 1 (stranded reads) or 2 (reversely stranded reads) (optional, default: 0)
# # - h: Help message (optional, displays the usage of the script)
# # - v: Version of the script (optional)
#
# Output format: TSV file (read count file)
# 
# Usage: ./05_script.sh -i <directory_with_processed_RNAseq_data> -o </path/to/output/directory> -a <path/to/annotation/file>
#
# Usage example: ./05_script.sh -i /data/reads.bam -o /home/user/output_directory -a /home/user/annotation_file.gtf -O -s 1 -T 6
#
# Expected output: A .tsv file with the read-count file containing all samples, and a log directory with the output and error log files of the general execution and of the featureCounts analysis
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

# Define the variables regarding optional parameters as their default values:
THREADS=1
OVERLAP=0
STRAND=0
PAIRED=0

# Check if the user has provided the required arguments with getopts
# If they have, assign the arguments to the corresponding variables
# If not, print an error message and exit
{
while getopts "i:o:a:OpT:s:hv" opt; do
 case $opt in
 i) READS=$OPTARG ;;
 o) OUTPUT_SELECTED=1  # This variable is used to check if the user has provided the output directory
    OUTPUT_DIRECTORY=$OPTARG ;;
 a) REFGEN_0=$OPTARG
    REFGEN=$(readlink -f "$REFGEN_0") ;;
 T) THREADS=$OPTARG ;;
 O) OVERLAP=1 ;; # This variable is used to verify later on that the user has provided the -O option
 p) PAIRED=1 ;; # This variable is used to verify later on that the user has provided the -p option
 s) STRAND=$OPTARG ;;
 v) echo -e "\nVersion: $VERSION\n" ;;
 h) echo -e "\nScript to perform a quantification of reads per gene from RNA-Seq sequencing data, which has been previously processed (trimmed and aligned), using featureCounts as tool and the annotation file of the genome of interest"
    echo -e "\nInput options and arguments:"
    echo -e "\t-i: Directory (path) with the RNA-seq processed data, in BAM or SAM format"
    echo -e "\t-o: Output directory (optional, default: current directory)"
    echo -e "\t-a: Reference genome annotation file (GTF or GFF3 format)"
    echo -e "\t-T: Followed by the number of threads to use (optional, default: 1)"
    echo -e "\t-O: Assigns reads to all their overlapping meta-features"
    echo -e "\t-p: Paired-end reads (optional, default: single-end reads)"
    echo -e "\t-s: Strandness of the data, followed by 0 (unstranded reads), 1 (stranded reads) or 2 (reversely stranded reads) (optional, default: 0)"
    echo -e "\t-h: Help message (optional, displays the usage of the script)"
    echo -e "\t-v: Version of the script (optional)\n"
    echo -e "\nUsage: $0 Usage: ./05_script.sh -i <directory_with_processed_RNAseq_data> -o </path/to/output/directory> -a <path/to/annotation/file>\n";;
 ?) echo -e "\nInvalid option or missing argument. Use -h for help.\n" >&2
    sleep 0.25s # Command to give bash enough time to close the pipe between the subshell with tee before exiting
    exit 1 ;;
  esac
done

# Check if the user has provided the bam file with the reads and the reference genome annotation file
# If not, it checks if the user has provided an output directory
# If not, it just exits, because it means that they have provided the -h or -v option, or both, but none of the needed files, so the user seems not to want to download anything
# In case they have provided the output directory but not the necessary files, it prints an error message and exits
if [[ "$READS" == "" ]]; then
 if [[ "$OUTPUT_SELECTED" -ne 1 ]]; then
  if [[ "$REFGEN" == "" ]]; then
   sleep 0.25s # This line is necessary to have something before "exit 0" so that bash has enough time to close the pipe and tee does not wait forever to receive an input
   exit 0
  else
    echo -e "\nYou have not provided the bam/sam file(s) with the processed reads. Use -h for help.\n" >&2
    sleep 0.25s # Again, time for bash to close the pipe with the subshell of tee
    exit 1
  fi
 else
  echo -e "\nYou have not provided the reference genome annotation file. Use -h for help.\n" >&2
  sleep 0.25s # Again, time for bash to close the pipe with the subshell of tee
  exit 1
  fi
fi

if [[ "$REFGEN" == "" ]]; then
 if [[ "$OUTPUT_SELECTED" -ne 1 ]]; then
  if [[ "$READS" == "" ]]; then
   sleep 0.25s # This line is necessary to have something before "exit 0" so that bash has enough time to close the pipe and tee does not wait forever to receive an input
   exit 0
  else
    echo -e "\nYou have not provided the reference genome annotation file. Use -h for help.\n" >&2
    sleep 0.25s # Again, time for bash to close the pipe with the subshell of tee
    exit 1
  fi
 else
  echo -e "\nYou have not provided the bam/sam file(s) with the processed reads. Use -h for help.\n" >&2
  sleep 0.25s # Again, time for bash to close the pipe with the subshell of tee
  exit 1
  fi
fi

# Check if the bam file exists, is not empty and with the correct extension (.bam or .sam))
# If not, print an error message and exit
# If it exists, check if the user has the necessary permissions to read the file
for file in "$READS"/*.[bs]am ; do
  symb_path=$(readlink -f "$file") # New command not learned in class to get the absolute path of the original file the link symbolic points to, in case that "file" is a symbolic link.
  if [[ ! -f "$symb_path" ]]; then 
   echo -e "\nRNA-seq read file does not exist or is not a file.\n" >&2
   sleep 0.25s # Again this command. It would be optional in case we had not created a { } block redirected to a subshell
   exit 1
  elif [[ ! -s "$symb_path" ]]; then
   echo -e "\nThe file provided with the reads is empty.\n" >&2
   sleep 0.25s 
   exit 1
  elif [[ ! "$file" == *.[bs]am ]]; then
    echo -e "\nThe extension of the file provided is neither .bam nor .sam.\n" >&2
    sleep 0.25s
    exit 1
  elif [[ ! -r "$file" ]]; then
    echo "File is not readable. Fixing..." >&2
    (chmod +r "$file" && echo "Permissions fixed") || (echo "Could not change permissions" >&2 && exit 1)
  fi
done

# Check if the reference genome annotation file exists, is not empty and with the correct extension (.gtf or .gff3))
# If not, print an error message and exit
# If it exists, check if the user has the necessary permissions to read the file
symb_path=$(readlink -f "$REFGEN") # New command not learned in class to get the absolute path of the original file the link symbolic points to, in case that "file" is a symbolic link.
if [[ ! -f "$symb_path" ]]; then 
 echo -e "\nReference genome annotation file does not exist or is not a file.\n" >&2
 sleep 0.25s # Again this command. It would be optional in case we had not created a { } block redirected to a subshell
 exit 1
elif [[ ! -s "$symb_path" ]]; then
 echo -e "\nThe file provided with the reference genome is empty.\n" >&2
 sleep 0.25s 
 exit 1
elif [[ ! "$REFGEN" == *.gtf && ! "$REFGEN" == *.gff3 ]]; then
  echo -e "\nThe extension of the file provided is neither .gtf nor .gff3.\n" >&2
  sleep 0.25s
  exit 1
elif [[ ! -r "$REFGEN" ]]; then
  echo "File is not readable. Fixing..." >&2
  (chmod +r "$REFGEN" && echo "Permissions fixed") || (echo "Could not change permissions" >&2 && exit 1)
fi

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
} 2> >(tee ./05_logs.err) > >(tee  ./05_logs.out) # Command to redirect standard output and standard input to different files and print them on the screen at the same time.
# >( ) This structure is called process substitution, and it is used to create a temporary file from the subshell that can be used as input or output for a command. 
# In other words, it gives the subshell a file descriptor and treats it as a file.
# This process receives the standard output (>) and error (2>) respectively, and gives them as input to the subshell, which redirects each one to a different file.


# Create the logs directory that will be used later on for redirections
mkdir -p "${OUTPUT_DIRECTORY}"/05_logs/05_logs_general
mkdir -p "${OUTPUT_DIRECTORY}"/05_logs/05_logs_featurecounts
# Move the logs to the output directory
mv ./05_logs.err "${OUTPUT_DIRECTORY}"/05_logs/05_logs_general
mv ./05_logs.out "${OUTPUT_DIRECTORY}"/05_logs/05_logs_general

# Now we start with the analysis itself

# Create a directory for the results
FEATURECOUNTS_OUTPUT="$OUTPUT_DIRECTORY/05_results/"

mkdir -p "$FEATURECOUNTS_OUTPUT"


# Now we will perform the FeatureCounts analysis with each of the files located in the input directory

{
declare -a files_array=("$READS"/*.[bs]am) # Create an array with all the files in the input directory

# Start with the analysis 
printf "\nProcessing files:\n%s\n" "${files_array[@]}"
sleep 2 # To give some time to the user to read the message
echo -e "Generating the tsv file...\n"
sleep 2
# Then, apply featureCounts
# Use FeatureCounts to count the aligned reads per gene, generating a tsv file with the results
 if [[ $OVERLAP -eq 1 ]]; then
  if [[ $PAIRED -eq 1 ]]; then
    featureCounts -O -p -T $THREADS -s $STRAND -o "$FEATURECOUNTS_OUTPUT"/results.tsv -a $REFGEN "${files_array[@]}"
  else
    featureCounts -O -T $THREADS -s $STRAND -o "$FEATURECOUNTS_OUTPUT"/results.tsv -a $REFGEN "${files_array[@]}"
  fi
elif [[ $PAIRED -eq 1 ]]; then
  featureCounts -p -T $THREADS -s $STRAND -o "$FEATURECOUNTS_OUTPUT"/results.tsv -a $REFGEN "${files_array[@]}"
else
  featureCounts -T $THREADS -s $STRAND -o "$FEATURECOUNTS_OUTPUT"/results.tsv -a $REFGEN "${files_array[@]}"
 fi

 # Check if the tsv file was created successfully
 if [[ $? -ne 0 ]]; then
  echo -e "\nError creating tsv file.\n" >&2
  exit 1
 else
  echo -e "\nFeatureCounts tsv file of the samples created successfully.\n"
 fi
} 2> >(tee "${OUTPUT_DIRECTORY}"/05_logs/05_logs_featurecounts/05_logs_featurecounts.err) > >(tee "${OUTPUT_DIRECTORY}"/05_logs/05_logs_featurecounts/05_logs_featurecounts.out) # Redirect the standard output and standard error to the corresponding log files


# End of the script
