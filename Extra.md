#### 1. We wanted to compare the same pipeline on the whole paired end reads for the heart failure patient (HFP) and the normal case so we got both of the reads of each case:

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
#### 1. We had trouble working inside the VM due to space limitaions and so on, so we worked on ubuntu 18.04 as well:
Of course, we had to install the software needed and handle dependicies ourselves, we used virtual enviroments in both cases.
At first we working in parallel and comparing the results which differed slightly, for example the overall aligment scores from the hisat:

VM | Ubuntu 18.04
-------------------
93% |  87% 

