# NGS1 Project Pipline:

## 1. Getting a sample data from NCBI, we got ours from these links:
  *It's a paired end read*

   [Read1](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR830/SRR830985/SRR830985_1.fastq.gz), [Read2](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR830/SRR830985/SRR830985_2.fastq.gz)

   And the refernce would be the **[human transcriptome](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.transcripts.fa.gz)**

**This was done using the following commands:**

   ###### 1. Making a new work diroctory:
```bash
mkdir ~/workdir/fqData2 && cd ~/workdir/fqData2
```

   ###### 2. Downloading the files & unzipping:
```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR830/SRR830985/SRR830985_1.fastq.gz
gunzip SRR830985_1.fastq
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR830/SRR830985/SRR830985_1.fastq.gz
gunzip SRR830985_2.fastq

```
   ###### 3. Checking that the files contains enough reads & subsetting:
*We want 5 million reads, and each read corrosponds to 4 lines in fastq files. Therefore, we need to extract 20 million lines to have 5 million reads
Also, the two files should have the same length*
```bash
wc -l SRR830985_1.fastq
wc -l SRR830985_2.fastq
```


   ###### 4- Subsetting 5 million reads & checking:
```bash
sed -n '1,20000000p' SRR830985_1.fastq > SRR830985_subset1.fastq
sed -n '1,20000000p' SRR830985_2.fastq > SRR830985_subset2.fastq
wc -l SRR830985_subset2.fastq
less SRR830985_subset2.fastq
```
*We want 5 million reads, and each read corrosponds to 4 lines in fastq files. Therefore, we need to extract 20 million lines to have 5 million reads
Also, the two files should have the same length*
```bash
wc -l SRR830985_1.fastq
wc -l SRR830985_2.fastq
```

## 2. Alignment:

 ###### 1. Downloading the reference (human trnscriptome) and renaming it:
 ```bash
cd ~/workdir/sample_data
 wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.transcripts.fa.gz
gunzip gencode.v33.transcripts.fa.gz
mv gencode.v33.transcripts.fa.gz trancriptome.fa
 ```
 ###### 2. Hisat alignment:
 
```bash
mkdir -p ~/workdir/hisat_align/hisatIndex && cd ~/workdir/hisat_align/hisatIndex
ln -s ~/workdir/sample_data/transcriptome.fa .
hisat2-build -p 1 transcriptome.fa transcriptome
cd ~/workdir/hisat_align
R1="$HOME/workdir/sample_data/SRR830985_subset1.fastq"
R2="$HOME/workdir/sample_data/SRR830985_subset2.fastq"
hisat2 -p 1 -x hisatIndex/transcriptome --dta --rna-strandness RF -1 $R1 -2 $R2 -S SRR.sam

```

## 3. Extracting the secondry alignments (and the primary for comparison):

 ###### 1. Secondary alignment:
 ```bash
samtools view -H -f 0x100 SRR.sam > SRR_secondry.sam && samtools view -f 0x100 SRR.sam >> SRR_secondry.sam 
 ```
 ###### 1. Primary alignment:
 ```bash
samtools view -H -F 0x100 SRR.sam > SRR_primary.sam && samtools view -F 0x100 SRR.sam >> SRR_primary.sam 
 ```
## 3. Getting the average GC content:
Before we calculate the GC content, we suspected that some reads might be repeated in the secondary alignment:
**By definition, a secondary alignment means the read was mapped to multiple places, so there might be unbalanced repetitions that could throw of our calculations.**

 ###### 1. Checking for repetitions in the secondary alignment file:
 ```bash
cat SRR_primary.sam | awk '{print $1}' | sort 	| uniq -c 	| sort -k2nr 	| awk '{printf("%s\t%s\n",$2,$1)}END{print}' > primary_reads_count.txt
cat SRR_secondry.sam | awk '{print $1}' | sort 	| uniq -c 	| sort -k2nr 	| awk '{printf("%s\t%s\n",$2,$1)}END{print}' > secondry_reads_count.txt
 ```
 
Then using:
 ```bash
less primary_reads_count.txt
```
**We notice that all reads are repeated twice. which makes sense for paired-end reads.
 Similarly:**
 ```bash
less secondry_reads_count.txt
```
**We notice that all reads exist only once!
To double check for the whole file, we ran the following commands:**
```bash
File=$"primary_reads_count.txt"
File2=$"secondry_reads_count.txt"

cat $File | awk '$2!="2"' | wc -l # Result is 1
cat $File2 | awk '$2!="1"' | wc -l # Result is 2196313
```
**-> So, now we have to figure which reads to filter in the repeatitions.**
We were concerned that there might be overlapping in the repetitions, one idea to resolve this was using an asmebler, which we did work on with stringtie, but then we get a better idea, that we should check first for overlaping ... So that's what we did.

###### 1. Checking for overlapping:
```bash
sed -n '2p' ~/workdir/DataSRR830985_subset1.fastq | wc -c # Getting the original count = 101bp per read (102 character)
samtools view ~/workdir/hisat_align/SRR_secondry.sam | awk '{print $1, length($10)}' | awk '{ seen[$1] += $2 } END { for (i in seen) print i, seen[i] }' > secondry_lengthes.txt
head secondry_lengthes.txt # Found lengthes over the expected total which is 101
```
**We noticed that the read length was constant, so the aligner always aligns the entire read, as opposed to our initial understanding that it aligns parts of the read ... We discovered this becasuse the results of the previous command were all multiples of 101.**
-> Our idea now was to extract unique reads using their Queryname in the first column of the same file, but that would mean we could lose the paired reads because they share the same Queryname just like the repeatitions ... **So the question becomes: how to differntiate repeatitions from mate reads?**

###### 2. Getting non-repeated sequences for GC content calculations:
First we extracted the unique first mate reads and stored them in a file, then we extracted the unique second mates and **added** them to the file.
```bash
samtools view -F 0x0040 SRR_secondry.sam | awk '{print $1, $10}' | sort | uniq  > secondry_uniq.sam # Unique first reads in a mate
samtools view -F 0x0080 SRR_secondry.sam | awk '{print $1, $10}' | sort | uniq  >> secondry_uniq.sam # Add unique second read in a mate to the file
```
The output file consists of two columns as follows:

Queryname  |  Actual sequence
------------------------------
980027648  | GTCATGGCAAT
           |
-> The above is just an example for illustration.           
###### 2. To count the G/C characters:
First we extract the pure sequence and discard the Queryname:
```bash
#Loop over lines in the file:
while IFS= read -r line; do string=${line#* } ; echo $string; done < secondry_uniq.sam > pure_2nd_seqs.txt

read C <<<$(tr -cd 'C' < pure_2nd_seqs.txt | wc -c)
read G <<<$(tr -cd 'G' < pure_2nd_seqs.txt | wc -c)
GC=$(($G+$C))
```
###### 3. To get avg GC content:
```bash
read lines filename <<< $(wc -l secondry_uniq.sam)
total=$(($lines*101))
GC=$(($GC*100))
echo "GC Content is $(($GC/$total))%" # GC Content is 50%
```
###### 4. Repeat for primary, but with no need for eliminating repetitions:
```bash
samtools view SRR_primary.sam | awk '{print $10}' > pure_1st_seqs.txt

read C <<<$(tr -cd 'C' < pure_1st_seqs.txt | wc -c)
read G <<<$(tr -cd 'G' < pure_1st_seqs.txt | wc -c)
GC=$(($G+$C))

#To get avg GC content:
read lines filename <<< $(wc -l SRR_primary.sam)
total=$(($lines*101))
GC=$(($GC*100))
echo "GC Content is $(($GC/$total))%" # GC Content is 47%
```
## 4. Getting the average mapping quality score:

 ##### 1. For Secondary:
 ```bash
 samtools view SRR_secondry.sam | awk '{sum+=$5} END { print "Mean MAPQ =",sum/NR}'
 #Output: Mean MAPQ = 0.985319
 ```
 ##### 1. For Primary:
 ``` bash
 samtools view SRR_primary.sam | awk '{sum+=$5} END { print "Mean MAPQ =",sum/NR}'
 #Output: Mean MAPQ = 26.4591
 ```
 

## 4. Getting the average sequncing quality score:

 ##### 1. For Secondary:
 ```bash
samtools view SRR_secondry.sam | perl -lane '@quals = split(undef, $F[10]); foreach $qual (@quals) {print ord($qual)-33}' | awk '{sum+=$0} END { print "Mean BASEQUAL =",sum/NR}
#Output: Mean BASEQUAL = 36.4889
 ```
 ##### 1. For Primary:
 ``` bash
samtools view SRR_primary.sam | perl -lane '@quals = split(undef, $F[10]); foreach $qual (@quals) {print ord($qual)-33}' | awk '{sum+=$0} END { print "Mean BASEQUAL =",sum/NR}'
#Output: Mean BASEQUAL = 36.1437
 ```
 
  
