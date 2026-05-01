# nanopore_QC

### Requirements
[Nextflow](https://www.nextflow.io/)

[Conda](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html)

### Setup
Clone and cd into repo
```
git clone https://github.com/graemefox/nanopore_QC.git
cd nanopore_QC
```

### Inputs

#### Input Data
The pipeline takes either BAM or FastQ files.
Pass either a single BAM/FAstQ or a list of BAMs/FastQs like 
(note the single quotes)

```
--bam '/path/to/my.bam'
# or
--bam '/path/to/*bam'
# or
--fastq '/path/to/*fq.gz'
# or
--fastq '/path/to/my.fq.gz'
```

#### Reference Genomes
You also pass a directory containing reference genomes to screen against (in .fa or .fa.gz).

For eg. download Human hg38, Mouse mm10 and A. fumigatus for use:
```
mkdir refs
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -P refs/
wget https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz -P refs/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/GCF_000002655.1_ASM265v1/GCF_000002655.1_ASM265v1_genomic.fna.gz -P refs/
mv refs/GCF_000002655.1_ASM265v1_genomic.fna.gz refs/Aspergillus_fumigatus_GCF_000002655.1_ASM265v1_genomic.fna.gz
```

Note - sequence IDs must be unique across all reference genomes used. 
Eg. Human and Mouse cannot both have a sequence called "chr1".
If necessary prepend sequences with the species (it is necessary with these examples)
```
zcat refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | sed 's/^>/>hg38_/' \
    | pigz > refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.tmp \
    && mv refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.tmp \
    refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

zcat refs/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz | sed 's/^>/>mm10_/' \
    | pigz > refs/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz.tmp \
    && mv refs/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz.tmp \
    refs/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz

zcat refs/Aspergillus_fumigatus_GCF_000002655.1_ASM265v1_genomic.fna.gz | sed 's/^>/>Af293_/' \
    | pigz > refs/Aspergillus_fumigatus_GCF_000002655.1_ASM265v1_genomic.fna.gz.tmp \
    && mv refs/Aspergillus_fumigatus_GCF_000002655.1_ASM265v1_genomic.fna.gz.tmp \
    refs/Aspergillus_fumigatus_GCF_000002655.1_ASM265v1_genomic.fna.gz
```

Filter to only main chromosome sequences, where appropriate
```
zcat refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
  | awk 'BEGIN{keep=0} /^>/{keep=($0 ~ /^>hg38_([0-9]+|X|Y|MT) /)} keep' \
  | pigz > refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.tmp \
  && mv refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.tmp refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

zcat refs/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz \
  | awk 'BEGIN{keep=0} /^>/{keep=($0 ~ /^>mm10_([0-9]+|X|Y|MT) /)} keep' \
  | pigz > refs/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz.tmp \
  && mv refs/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz.tmp refs/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
```

#### (optional) Download Test Data
```
git clone https://github.com/LooseLab/ROBIN_test_set_A.git
```
#### Choosing a profile for EPI2ME subworkflows
By default the EPI2ME workflows will try to use docker which won't work on the HPC.

Specify singularity:
```
--EPI2ME_profile singularity
```

#### Demo command
```
nextflow run main.nf --bam 'path/to/*.bam' --ref_dir path/to/refs/ --outdir results --threads 8
# or
nextflow run main.nf --fastq 'path/to/*.fq.gz' --ref_dir path/to/refs/ --outdir results --threads 8
```

To run with the demo data downloaded above
```
nextflow run main.nf --bam 'ROBIN_test_set_A/test_data_set/*.bam' --ref_dir refs/ --outdir results --threads 8
```

#### Downsample large data
For large datasets you may want to downsample the data for wf-alignment
(NanoPlot takes the entire dataset)

It works with FastQ and BAM inputs.

Adding "--downsample 0.1" runs wf-alignmemt with a random 10% of the data

```
nextflow run main.nf --downsample 0.1 --fastq 'path/to/*.fq.gz' --ref_dir path/to/refs/ --outdir results --threads 8
```

