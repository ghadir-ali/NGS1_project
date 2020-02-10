We discussed we Dr. Tamer that if we did extra work on the project we could get some bonus marks to substitute for any shortage in the overall grade of the course due to any other graded process.
-----------------------------------------------------------------

#### 1. Different machines:
**We had trouble working inside the VM due to space limitaions and so on, so we worked on ubuntu 18.04 as well:**
Of course, we had to install the software needed and handle dependicies ourselves, we used virtual enviroments in both cases.
At first we working in parallel and comparing the results which differed slightly, for example the overall aligment scores from the hisat:

| VM  | Ubuntu 18.04|
|-----|-------------|
| 93% |      87%    |

We had many errors and incompatabilities while using the Ubuntu, but we were able to resolve them and continue working on it. However, we couldn't continue to compare the results because the VM failed to go any further due to space limitaions dispeite increasing its working space, of course there was no comparison of the spapce dedicated to the VM and the space dedicated to a Ubuntu OS installed on the hard disk as dual boot.

#### 2. Assembly:
We first thought of the assembly as a solution to the repetitions problem in the secondary alignments, but we couldn't get a refernce, so we did a ref-free assembly of the secondary alignment file and another for the primary.
**> This step was used for comparison now**
We already had the bam files as an attempt to reduce the project space (please refer to the end of the main [pipeline](https://github.com/ghadir-ali/NGS1_project/blob/master/Project_pipeline.md).
S, we started by sorting the bam then ref-free assembly:
```bash
samtools sort SRR_primary.bam  -o SRR_sorted_primary.bam 
samtools sort SRR_secondry.bam  -o SRR_sorted_secondry.bam 
stringtie SRR_sorted_primary.bam --rf -l ref_free -o primary.gtf
stringtie SRR_sorted_secondry.bam --rf -l ref_free -o secondry.gtf
```
Counting transcripts in each gtf to check the success of the process, and to get a clue of what the results might look like compared to each other:
```bash
cat primary.gtf | grep -v "^@" | awk '$3=="transcript"' | wc -l # Output: 7274
cat secondry.gtf | grep -v "^@" | awk '$3=="transcript"' | wc -l # Output: 7107
```
**Output:**
| primary.gtf | secondry.gtf |
|-------------|--------------|
|    7274     |     7107     |

**Interpretation:**
It makes sense that the secondry would have less transcripts because not all reads get a secondry mapping position, this also indicates that assembler did in-fact account for repititions, because otherwise the secondry would have had a higher count.

Then we compared the resulting gtf files against each other:
```bash
conda deactivate # Because all our work so far has been in a virtenv called ngs1, now we need to make a different one.

conda create -n ngs-gtf python=3.6 anaconda
source activate ngs-gtf
conda install -c conda-forge pypy3.5
wget https://bootstrap.pypa.io/get-pip.py
pypy3 get-pip.py


```

The resulting files are uploaded on the repo in the assembly results folder, and below is the table from the stats file:
#-----------------| Sensitivity | Precision  |
        Base level:    20.9     |    24.0    |
        Exon level:    10.7     |    14.3    |
      Intron level:     3.7     |    35.0    |
Intron chain level:     2.9     |    25.4    |
  Transcript level:    16.8     |    16.6    |
       Locus level:    17.0     |    16.6    |

     Matching intron chains:      73
       Matching transcripts:    1179
              Matching loci:    1177

          Missed exons:    8080/9928	( 81.4%)
           Novel exons:    5583/7408	( 75.4%)
        Missed introns:    2796/2911	( 96.0%)
         Novel introns:     186/311	( 59.8%)
           Missed loci:    5280/6942	( 76.1%)
            Novel loci:    5391/7070	( 76.3%)

 Total union super-loci across all input datasets: 7050 
7107 out of 7107 consensus transcripts written in gffcmp.annotated.gtf (0 discarded as redundant)


#### 3. Differential Expression:
We wanted to compare the same pipeline on the whole paired end reads for the heart failure patient (HFP) and the normal case so we got both of the reads of each case:
  ##### a. Getting the data, renaming and checking:

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR830/SRR830985/SRR830985_1.fastq.gz # hfp_read1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR830/SRR830985/SRR830985_2.fastq.gz # hfp_read2
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR117/001/SRR1175541/SRR1175541_1.fastq.gz # normal_read1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR117/001/SRR1175541/SRR1175541_2.fastq.gz #normal_read2
gunzip SRR830985_1.fastq.gz && mv SRR830985_1.fastq hfp_1.fastq.gz
gunzip SRR830985_2.fastq.gz && mv SRR830985_2.fastq hfp_2.fastq.gz
gunzip SRR1175541_1.fastq.gz && mv SRR1175541_1.fastq.gz normal_1.fastq
gunzip SRR1175541_2.fastq.gz && mv SRR1175541_2.fastq.gz normal_2.fastq
ls -tral
```
  ##### b. Alignment:
  First, we tried the following code (using the same transcriptome reference and its idx):
  ```bash
  norm_R1=$"normal_1.fastq"
  norm_R2=$"normal_2.fastq"
  hisat2 -p 1 -x hisatIndex/transcriptome --dta --rna-strandness RF -1 $norm_R1 -2 $norm_R2 -S norm.sam
  ```
  * First issue was that we got an error that hisat2 isn't installed, because we were working inside the ngs-gtf environment. We resolved it"
  ```bash
  conda deactivate
  conda activate ngs1
  ```
  Then we tried the same code again, but this tim got this output:
  *Error: Read SRR1175541.1 1_11_1966/1 has more quality values than read characters.
terminate called after throwing an instance of 'int'
Aborted*

We searched on biostars and found that most likely this error occurs when the file is corrupted, so we tried the other two reads with the same command:
```bash
hfp_R1=$"hfp_1.fastq"
hfp_R2=$"hfp_2.fastq"
hisat2 -p 1 -x hisatIndex/transcriptome --dta --rna-strandness RF -1 $hfp_R1 -2 $hfp_R2 -S hfp.sam
```
It worked with following output:
*
13176801 reads; of these:
  13176801 (100.00%) were paired; of these:
    2165394 (16.43%) aligned concordantly 0 times
    5430495 (41.21%) aligned concordantly exactly 1 time
    5580912 (42.35%) aligned concordantly >1 times
    ----
    2165394 pairs aligned concordantly 0 times; of these:
      81534 (3.77%) aligned discordantly 1 time
    ----
    2083860 pairs aligned 0 times concordantly or discordantly; of these:
      4167720 mates make up the pairs; of these:
        3457853 (82.97%) aligned 0 times
        397110 (9.53%) aligned exactly 1 time
        312757 (7.50%) aligned >1 times
__86.88% overall alignment rate__

*
Which confirmed the possibility that the normal files were corrupted, so redownloaded them and repeated this process, and got the following output:
