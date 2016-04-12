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

##extra step in here...
bcftools view -h flt.all.vcf>all_new.vcf
awk '{gsub("[^ACTGN,]","N",$4); gsub("[^ACTGN,]","N",$5)}1' OFS="\t" flt.all.vcf |grep ^chr >>all_new.vcf

```

###phasing with Beagle
The parents of Gala are not necessary therefore a beagle ped file of:
gala 1 0 0 
is appropriate
A gala vcf was produced seperatly and at a differnt time from the rootstocks.
```shell
awk '{gsub("[^ACTGN,]","N",$4); gsub("[^ACTGN,]","N",$5)}1' OFS="\t" ../gala_all_piledup.vcf |grep ^chr >>gala_new.vcf
java -jar /home/deakig/projects/apple_rootstock/scripts/beagle.r1399.jar gtgl=gala_new.vcf usephase=true ped=../../beagle/gala.ped out=galabgl chrom=1 nthreads=8

java -jar /home/deakig/projects/apple_rootstock/scripts/beagle.r1399.jar gtgl=../bcf/all_new.vcf ped=pedigree.ped out=newbgl nthreads=16

```

### Make diploid genomes with vcf2diploid
```shell
java -jar ~/projects/apple_rootstock/vcf2diploid/vcf2diploid.jar -id m116 -chr ~/projects/apple_rootstock/rootstock_genetics/ref/v1/Malus_x_domestica.v1.0-primary.pseudo.fa -vcf newbgl.vcf -outDir ~/projects/apple_rootstock/allele/m116/newnew

java -jar ~/projects/apple_rootstock/vcf2diploid/vcf2diploid.jar -id m27 -chr ~/projects/apple_rootstock/rootstock_genetics/ref/v1/Malus_x_domestica.v1.0-primary.pseudo.fa -vcf newbgl.vcf -outDir ~/projects/apple_rootstock/allele/m27/newnew

java -jar ~/projects/apple_rootstock/vcf2diploid/vcf2diploid.jar -id gala -chr ~/projects/apple_rootstock/rootstock_genetics/ref/v1/Malus_x_domestica.v1.0-primary.pseudo.fa -vcf galabgl.vcf -outDir ~/projects/apple_rootstock/allele/gala

```

#### Indexing 

Create mat and pat folders, link to chromosomes and build bowtie index
```shell
mkdir mat pat
ln -s *maternal.fa mat/.
ln -s *paternal.fa pat/.
```
Then

```shell
bowtie2-build chr10_m116_paternal.fa,chr13_m116_paternal.fa,chr16_m116_paternal.fa,chr2_m116_paternal.fa,chr5_m116_paternal.fa,chr8_m116_paternal.fa,chr11_m116_paternal.fa,chr14_m116_paternal.fa,chr17_m116_paternal.fa,chr3_m116_paternal.fa,chr6_m116_paternal.fa,chr9_m116_paternal.fa,chr12_m116_paternal.fa,chr15_m116_paternal.fa,chr1_m116_paternal.fa,chr4_m116_paternal.fa,chr7_m116_paternal.fa m116_paternal_index
bowtie2-build chr10_m116_maternal.fa,chr11_m116_maternal.fa,chr12_m116_maternal.fa,chr13_m116_maternal.fa,chr14_m116_maternal.fa,chr15_m116_maternal.fa,chr16_m116_maternal.fa,chr17_m116_maternal.fa,chr1_m116_maternal.fa,chr2_m116_maternal.fa,chr3_m116_maternal.fa,chr4_m116_maternal.fa,chr5_m116_maternal.fa,chr6_m116_maternal.fa,chr7_m116_maternal.fa,chr8_m116_maternal.fa,chr9_m116_maternal.fa m116_maternal_index
bowtie2-build chr10_m27_maternal.fa,chr13_m27_maternal.fa,chr16_m27_maternal.fa,chr2_m27_maternal.fa,chr5_m27_maternal.fa,chr8_m27_maternal.fa,chr11_m27_maternal.fa,chr14_m27_maternal.fa,chr17_m27_maternal.fa,chr3_m27_maternal.fa,chr6_m27_maternal.fa,chr9_m27_maternal.fa,chr12_m27_maternal.fa,chr15_m27_maternal.fa,chr1_m27_maternal.fa,chr4_m27_maternal.fa,chr7_m27_maternal.fa m27_maternal_index
bowtie2-build chr10_m27_paternal.fa,chr13_m27_paternal.fa,chr16_m27_paternal.fa,chr2_m27_paternal.fa,chr5_m27_paternal.fa,chr8_m27_paternal.fa,chr11_m27_paternal.fa,chr14_m27_paternal.fa,chr17_m27_paternal.fa,chr3_m27_paternal.fa,chr6_m27_paternal.fa,chr9_m27_paternal.fa,chr12_m27_paternal.fa,chr15_m27_paternal.fa,chr1_m27_paternal.fa,chr4_m27_paternal.fa,chr7_m27_paternal.fa m27_paternal_index
bowtie2-build chr10_gala_maternal.fa,chr13_gala_maternal.fa,chr16_gala_maternal.fa,chr2_gala_maternal.fa,chr5_gala_maternal.fa,chr8_gala_maternal.fa,chr11_gala_maternal.fa,chr14_gala_maternal.fa,chr17_gala_maternal.fa,chr3_gala_maternal.fa,chr6_gala_maternal.fa,chr9_gala_maternal.fa,chr12_gala_maternal.fa,chr15_gala_maternal.fa,chr1_gala_maternal.fa,chr4_gala_maternal.fa,chr7_gala_maternal.fa gala_maternal_index
bowtie2-build chr10_gala_paternal.fa,chr13_gala_paternal.fa,chr16_gala_paternal.fa,chr2_gala_paternal.fa,chr5_gala_paternal.fa,chr8_gala_paternal.fa,chr11_gala_paternal.fa,chr14_gala_paternal.fa,chr17_gala_paternal.fa,chr3_gala_paternal.fa,chr6_gala_paternal.fa,chr9_gala_paternal.fa,chr12_gala_paternal.fa,chr15_gala_paternal.fa,chr1_gala_paternal.fa,chr4_gala_paternal.fa,chr7_gala_paternal.fa gala_paternal_index
```
