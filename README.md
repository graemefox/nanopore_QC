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
Pass either a single BAM/FAstQ or a list of BAMs/FastQc like 
```
/path/to/*bam
#or
/path/to/*fq.gz
```

#### Reference Genomes
You also pass a directory containing reference genomes to screen against (in .fa or .fa.gz).

For eg. download hg38 and mm10 for use
```
mkdir refs
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -P refs/
wget https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz -P refs/
```

Note - sequence IDs must be unique across all reference genomes used. 
Eg. Human and Mouse cannot both have a sequence called "chr1".
If necessary prepend sequences with the species
```
zcat refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | sed 's/^>/>hg38_/' \
    | pigz > refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.tmp \
    && mv refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz.tmp \
    refs/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

zcat refs/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz | sed 's/^>/>mm10_/' \
    | pigz > refs/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz.tmp \
    && mv refs/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz.tmp \
    refs/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
```

Filter to only main chromosome sequences
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

#### Demo command
To run with the demo data
```
nextflow run main.nf --bams 'ROBIN_test_set_A/test_data_set/*.bam' --ref_dir refs/ --outdir results --threads 8
```
