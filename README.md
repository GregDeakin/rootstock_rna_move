# rootstock_rna_move

Test whether RNA has moved across the graft union.

Plan:
Identify snps in Gala
Identify snps in M27 and M116 
Identify snps unique between Gala:M27 and Gala:M116
Align RNAQseq reads to apple genome
Count SNPs in RNAseq 

The allele-seq pipeline may provide a convenient method for doing the majority of this work.

## Create maternal and paternal chromosomes 

###Assemble to reference
```shell
./bowtie.sh ~/projects/apple_rootstock/rootstock_genetics/gala/conc/phix_filtered.1  ~/projects/apple_rootstock//rootstock_genetics/gala/conc/phix_filtered.2 ~/projects/apple_rootstock/rootstock_genetics/ref/v1/Md ~/projects/apple_rootstock/rootstock_genetics/gala/analysis_v1/ gala_v1 250 500
```
###SAM to BAM and sort
```shell
./samtools.sh ~/projects/apple_rootstock/rootstock_genetics/gala/analysis/gala_v1.sam ~/projects/apple_rootstock/rootstock_genetics/gala/analysis/ gala_v1.bam gala_v1.sorted
```
###sort regions file hack
```shell
cat regions |xargs -I l touch l
ls -v chr* >regions
rm chr*
```
###index bam files
```shell
cat bam_files|xargs -I file samtools index file
```
####pileup
```shell
~/projects/apple_rootstock/scripts/pileup2.sh ~/projects/apple_rootstock/rootstock_genetics/ref/v1/Malus_x_domestica.v1.0-primary.pseudo.fa bam_files ~/projects/apple_rootstock/rna-seq/RNA_trans piledup.bcf regions

ls -v piledup.bcf.* > files
bcftools concat -O v -f files >rna.pileup.vcf
```
