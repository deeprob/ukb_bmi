# Variant Quality Control and Preparation
Whole exome sequencing variant call files (vcfs) downloaded from UKBiobank were prepared as follows:

1. Split multi-allelic records and left normalize variants - ./src/1_split_leftnorm.py
2. Annotate variants using VEP - ./src/5_vep_annot.sh
3. Add CADD annotations using vcfanno
4. Add loftee annotations
5. Add dbNSFP annotations using ANNOVAR
6. Filter for LoF and deleterious missense variant with intra-cohort freq<0.01 affecting protein coding genes


# Virtual environments

## Vep with loftee

### Conda env
```bash
foo@bar:~$ conda create -n loftee -c conda-forge -c anaconda -c bioconda perl ensembl-vep samtools perl-bioperl=1.7.2  perl-list-moreutils perl-dbd-sqlite perl-bio-bigfile -y
```

### Loftee installation
```bash
foo@bar:~$ git clone --single-branch --branch grch38 https://github.com/konradjk/loftee.git
```
**Had to change lines 263 to 267 in LoF.pm based on this commit: https://github.com/konradjk/loftee/commit/0f84baf09830d9bcf7e5fbc41e1bf04b8aaec2f6**


