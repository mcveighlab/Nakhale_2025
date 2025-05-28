#!/bin/bash

# Specify resource requirements for SLURM scheduler
#SBATCH --job-name=BI_mouse2_Full_run
#SBATCH --mail-user=<enter_email>
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --cpus-per-task=5
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-18
#SBATCH --output=./Logs/BI_mouse2_Full_run%a.log
#SBATCH --partition=k2-medpri
#SBATCH --chdir=/<path_to_working_directory>

# Load modules
module load apps/fastqc/0.11.8
module load apps/trimgalore/0.6.10
module load apps/bowtie2/2.5.2/gcc-14.1.0
module load miRDeep2/2.0.1.3/conda-2024.06
module load fastx_toolkit/0.0.14

# activate anaconda environment containing miRDeep2 apps
source activate mirdeep2

# Specify the path to the metadata file containing the sample names
# NOTE myarray.txt is a two column .txt file, left column titled ArrayTaskID, right column titled SampleName
# ArrayTaskID column contains integers from 1 to the total number of samples to be processed
# SampleName column contains the fastq file names, with .fq extension removed.
config=<path_to_working_directory>/myarray.txt

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
SampleName=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
ShortName=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
echo "This is array task ${SLURM_ARRAY_TASK_ID}, the sample name is ${SampleName}, ${ShortName}"

# Use bowtie2 to map readfiles to the host genome, outputting consensus host-mapped reads as a fastq

echo "----- mapping ${SampleName}, ${ShortName} against host genome --------------------------------"
bowtie2 -N 1 -p 5 --un ${ShortName}_HOST_UNMAPPED.fq --al ${ShortName}_HOST_MAPPED.fq -x ../Genomes/GCF_000001635.27_GRCm39_genomic_CLEAN -U ./Read_Files/${SampleName}_trimmed.fq -S ${ShortName}_vs_HOST.sam

# Use bowtie2 to map host-UNmapped readfiles to the worm genome, outputting consensus host-unmapped worm-mapped reads as a fastq

echo "----- bowtie2 mapping host unmapped fastq files against WORM genome, outputting host-unmapped, wormmapped reads ----"
bowtie2 -p 5 --un ${ShortName}_HOST_UNMAPPED_WORM_UNMAPPED.fq --al ${ShortName}_HOST_UNMAPPED_WORM_MAPPED.fq -x ../Genomes/dimmitis_WSI_2.2 -U ${ShortName}_HOST_UNMAPPED.fq -S ${ShortName}_HOST_UNMAPPED_vs_WORM.sam

echo "---------- fastx_toolkit converting fastq to fasta host reads ---------"
fastq_to_fasta -i ${ShortName}_HOST_MAPPED.fq -o ${ShortName}_HOST_MAPPED.fa
fastq_to_fasta -i ${ShortName}_HOST_UNMAPPED_WORM_MAPPED.fq -o ${ShortName}_HOST_UNMAPPED_WORM_MAPPED.fa

# Collapse duplicate reads to produce a fasta formatted for input to a BLASTn search
echo "---------- collapsing HOST fasta reads ----------"
perl ../collapse_reads.pl ${ShortName}_HOST_MAPPED.fa Mmu > ${ShortName}_HOST_MAPPED_COLLAPSED.fa

echo "---------- collapsing WORM fasta reads ----------"
perl ../collapse_reads.pl ${ShortName}_HOST_UNMAPPED_WORM_MAPPED.fa Dim > ${ShortName}_HOST_UNMAPPED_WORM_MAPPED_COLLAPSED.fa

echo "---------- mirdeep2 mapper generating arf file vs host ---------"
# use mirdeep2 mapper.pl to map host-mapped against host genome, to generate host.arf
mapper.pl ${ShortName}_HOST_MAPPED.fq -e -h -m -p ../Genomes/GCF_000001635.27_GRCm39_genomic_CLEAN -t ${ShortName}_HOST.arf -v -s ${ShortName}discard_host.fa

echo "---------- mirdeep2 mapper generating arf file vs worm ---------"
# use mirdeep2 mapper.pl to map host-unmapped worm-mapped against worm genome, to generate worm.arf
mapper.pl ${ShortName}_HOST_UNMAPPED_WORM_MAPPED.fq -e -h -m -p ../Genomes/dimmitis_WSI_2.2 -t ${ShortName}_WORM.arf -v -s ${ShortName}discard_worm.fa

echo "------- mirdeep2.pl finding host miRNAs ${ShortName} ----------"
miRDeep2.pl ${ShortName}_HOST_MAPPED_COLLAPSED.fa ../Genomes/GCF_000001635.27_GRCm39_genomic_CLEAN.fa ${ShortName}_HOST.arf ../reference_miRNAs/Mouse_miRNAs_mature_29Oct24.fa ../reference_miRNAs/Mouse_miRNAs_mature_29Oct24.fa ../reference_miRNAs/Mouse_miRNAs_hairpin_29Oct24.fa -c -d

echo "---------- mirdeep2.pl finding WORM miRNAs ${ShortName} ----------"
miRDeep2.pl ${ShortName}_HOST_UNMAPPED_WORM_MAPPED_COLLAPSED.fa ../Genomes/dimmitis_WSI_2.2.fa ${ShortName}_WORM.arf ../reference_miRNAs/Nematode_miRNAs_mature_29Oct24.fa ../reference_miRNAs/Nematode_miRNAs_mature_29Oct24.fa ../reference_miRNAs/Nematode_miRNAs_hairpin_29Oct24.fa -c -d
