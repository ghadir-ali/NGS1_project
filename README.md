# NGS1 Project Pipline:

## 1- Getting a sample data from NCBI, we got ours from these links:
*It's a paired end read*

...* [Read1.](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR830/SRR830985/SRR830985_1.fastq.gz)

...* [Read2.](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR830/SRR830985/SRR830985_2.fastq.gz)

And the refernce would be the **[human genome]**(ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz)

**This was done using the following commands:**

### 1- Making a new work diroctory:
```bash
mkdir ~/workdir/fqData2 && cd ~/workdir/fqData2
```

### 2- Downloading the files & unzipping:
```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR830/SRR830985/SRR830985_1.fastq.gz
gunzip SRR830985_1.fastq
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR830/SRR830985/SRR830985_1.fastq.gz
gunzip SRR830985_2.fastq

```

### 3- Checking that the files contains enough reads & subsetting:
*We want 5 million reads, and each read corrosponds to 4 lines in fastq files. Therefore, we need to extract 20 million lines to have 5 million reads
Also, the two files should have the same length*
```bash
wc -l SRR830985_1.fastq
wc -l SRR830985_2.fastq
```


## 2- Subsetting 5 million reads & checking:
```bash
sed -n '1,20000000p' SRR830985_1.fastq > SRR830985_subset1.fastq
sed -n '1,20000000p' SRR830985_2.fastq > SRR830985_subset2.fastq
wc -l SRR830985_subset2.fastq
less SRR830985_subset2.fastq
```
