

#Analyse public datasets ENCODE or specific papers (GEO) to analyse epigenetic data using the T2T genome

####### Encode #############
#download H1

#get metadata from ENCODE (download >selected files from the filters) here only fastq files
#move to the metadata and download the information there

mkdir FastqFiles
mkdir BigWigFiles
mkdir BamFiles
mkdir LogFiles
mkdir Metadata

cd Metadata
#move ${TARGET}_Metadata.tsv ${TARGET}_Files.txt into the Metadata Folder
TARGET='HCT116'

#pulls out the download links based on the metadata file, not needed if pre filtered but its easier to filter the Metadata if encode has a lot of different Histone ChIPs that are not needed
awk -F'\t' 'NR==FNR {ids[$1]; next} {for (id in ids) if (id != "" && $0 ~ id) print $0}' ${TARGET}_Metadata.tsv ${TARGET}_Files.txt > filtered_${TARGET}_links.txt

#download the fastqFiles
xargs -n 1 curl --output-dir ../FastqFiles -O -L < "filtered_${TARGET}_links.txt"
#xargs -n 1 curl -O -L < filtered_H1_Control_links.txt


#for ChIP-seq renaming, inspect metadata to manually change columns if needed (or add in Control into column 23
awk -F'\t' '
NR==FNR && FNR > 1 {
    new_name = $1 "_" $11 "_" $23 "_" "SE" $37 "_" $35 ".fastq.gz"; 
    filename= $1 ".fastq.gz";
    print "mv " filename " " new_name;
    }' ${TARGET}_Metadata.tsv > rename_files_${TARGET}_ChIP.sh

mv rename_files_${TARGET}_ChIP.sh ../FastqFiles/
cd ../FastqFiles/
bash rename_files_${TARGET}_ChIP.sh
#Manually cat 2 replicates

#for control 
#not needed for most cell lines
awk -F'\t' '
NR==FNR {
# Construct the name first
    raw_name = $1 "_" $8 "_" $35;
    
    # Clean the whole thing: replace any whitespace with a single underscore
    gsub(/[[:space:]]+/, "_", raw_name);
    
    # Remove any trailing underscores that might have come from empty columns
    gsub(/_+$/, "", raw_name);
    
    mapping[$1] = raw_name; 
    next
}
{
    # For every link, extract the ENC ID (e.g., ENCFF000AYG)
    split($0, url_parts, "/");
    filename = url_parts[length(url_parts)];
    split(filename, name_parts, ".");
    file_id = name_parts[1];

    # If the ID exists in our mapping, print a rename command
    if (file_id in mapping) {
        new_name = mapping[file_id] ".fastq.gz";
        print "mv " filename " " new_name;
    }
}' H1_Metadata_Controls.txt filtered_H1_Control_links.txt > rename_files_H1_Control.sh


###### SRR download #########

SRR_Dataset=HNSCC

prefetch=/home/ec2-user/Programs/sratoolkit.3.0.0-ubuntu64/bin/prefetch
$prefetch --option-file SRR_Acc_List_HNSCC.txt

fasterq=/home/ec2-user/Programs/sratoolkit.3.0.0-ubuntu64/bin/fasterq-dump

xargs -P 4 -I {} sh -c "
  $fasterq --stdout {} | gzip > {}.fastq.gz && \
  rm -rf /home/ec2-user/Programs/sratoolkit.3.0.0-ubuntu64/SRAdeposit/sra/{}.sra
" < SRR_Acc_List_HNSCC.txt

#edit names using SraRunTable information

awk -F',' '
NR==FNR && FNR > 1 {
    # Remove quotes and carriage returns from the fields
    gsub(/"/, "", $1); gsub(/"/, "", $8); 
    gsub(/"/, "", $30);
    gsub(/\r/, "", $0);

    new_name = $1 "_" $8 "_" $30 ".fastq.gz"; 
    filename = $1 ".fastq.gz";
    
    if ($1 != "") {
        print "mv " filename " " new_name;
    }
}' "SraRunTable_${SRR_Dataset}.csv" > "rename_files_${SRR_Dataset}.sh"



# Standard Illumina Adapter

ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
CUTADAPT=~/miniconda3/envs/cutadapt/bin/cutadapt

#test case
zcat ENCFF023NCH_HCT116_H2AFZ_SE36_1.fastq.gz | head -n 400000 | gzip > test.fastq.gz

conda activate BioChannelKristjan
{
for q in *.fastq.gz
do
    # Define base name for cleaner logs and outputs
    BASE=${q%.fastq.gz}

    echo "Processing $BASE..."

    # 1. Trim with cutadapt
    # -m 15: discard reads shorter than 15bp after trimming
    # -q 20: quality trim ends at Phred 20
    $CUTADAPT -a $ADAPTER -j 8 -m 20 -q 30 $q 2> ${BASE}.TrimSummary.log | \

    # 2. Align with Bowtie2
    # -U - tells Bowtie2 to read from the pipe
    bowtie2 -t -p 30 -L 18 -N 0 --very-sensitive --phred33 \
        -x /home/ec2-user/Genomes/HumanT2Tv2_bowtie2/chm13v2.0/chm13v2.0 \
        -U - 2> ${BASE}.AlignmentSummary.log | \
    
    # 3. Convert and Sort
    # Using sambamba directly for speed; -l 0 for uncompressed stream to sort
    samtools view -@ 10 -u -h - | \
    sambamba sort -t 10 /dev/stdin -o ${BASE}_sorted.bam && \
    
    # 4. Mark Duplicates
    # This takes the sorted BAM we just created
    sambamba markdup -t 10 ${BASE}_sorted.bam ${BASE}_final.bam 2> ${BASE}.DuplicateSummary.log

    #Clean up the intermediate sorted BAM
    rm ${BASE}_sorted.bam
    rm ${BASE}_sorted.bam.bai
    mv ${BASE}*.log ../LogFiles
    mv ${BASE}_final.bam* ../BamFiles
done
}> Encode_${TARGET}_ChIP_pipeline.log 2>&1 &
disown

tail -f Encode_${TARGET}_ChIP_pipeline.log

conda deactivate

#T2T Bigwigs for inspection
for b in *.bam
do
BASE=${b%_final.bam}
bamCoverage -b ${BASE}_final.bam \
    -o ${BASE}.bw \
    --binSize 200 \
    --normalizeUsing CPM \
    --effectiveGenomeSize 3117292070 \
    --extendReads 200 \
    -p 60
done

#FeatureCount
/home/ec2-user/Programs/subread-2.0.3-Linux-x86_64/bin/featureCounts -T 40 -O \
    -M --fraction \
    -a databaseS01_cenSat_Annotation.SAF -F SAF \
    -o ../FeatureCountFiles/H1_centromere_counts.txt \
    *_final.bam


#
