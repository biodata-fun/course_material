# IGSR - Variant Calling Course

We will use for this course the genomic data generated for [_Oryza sativa_ Japonica](https://plants.ensembl.org/Oryza_sativa/Info/Index) (rice) 

## Installation of dependencies
To install all the necessary software required in this course we will mainly use [Conda](https://docs.conda.io/en/latest/), which is a package and environment management system that is very convenient for installing all of the dependencies we will use.

1. First, create a new conda environment named `masters_jan2020`

        conda create --name masters_jan2020
2. Check that the new environment exists by doing

        conda env list
3. Activate envirnoment

        conda activate masters_jan2020
4. Install the necessary software

        conda install freebayes
        conda install -c bioconda bcftools
        conda install -c bioconda igv

## Variant calling
Variant calling is the process to identify variants from sequence data ([see](https://www.ebi.ac.uk/training-beta/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/#:~:text=What%20is%20variant%20calling%3F,creating%20BAM%20or%20CRAM%20files.)).
It starts with the sequencing data in the `FASTQ` format, these sequence reads need to aligned to the reference genome following the procedure explained in the first part of this course, that generated an alignment in the `BAM` file format. Once we have an alignment or multiple alignment files we can use a variant discovery tool to identify the germline variants. There are multiple variant calling tools available, the ones we have more experience with in our group are [SAMTools mpileup](http://samtools.github.io/bcftools/bcftools.html#mpileup), the [GATK suite](https://gatk.broadinstitute.org/hc/en-us) and [FreeBayes](https://github.com/ekg/freebayes). In this course we are going to use FreeBayes, as it is sensitive, precise and relatively simple to use. 

### Freebayes
FreeBayes is a haplotype-based variant detector, that uses a joint genotyping method capable of reporting variants on a single sample or on a cohort of multiple samples. It's going to be capable of detecting SNPs (single nucleotide polymorphisms), indels (short insertions and deletions) and MNPs (multi-nucleotide polymorphisms)

#### **Reference Genome**
FreeBayes needs the reference sequence in the `FASTA` format. In this section of the course we are going to use the same chromosome 10 sequence extracted from the *Oryza_sativa* (rice) genome that we used for the alignment part of the course.

#### **Using FreeBayes**
To run Freebayes you need to specify the ploidy of the genome being analysed, the FASTA reference sequence used for the alignment and the analhysis-ready BAM generated in the first section of the course. Once you have these prepared enter the following command in your terminal:

        freebayes -f /path/to/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --ploidy 2 /path/to/SAMEA2569438.chr10.reheaded.sorted.mark_duplicates.bam |bgzip -c > SAMEA2569438.chr10.vcf.gz &

This command pipes the output of FreeBayes to `bgzip`, which is a special compression/decompression program that is part of SAMTools. It is better to compress the VCF to make the file size smaller and also to use some of the `BCFTools` commands that are discussed later in this section.

#### **Understanding the FreeBayes output VCF**
After running FreeBayes, you will see a vcf file named `SAMEA2569438.chr10.vcf` containing the identified variants. The complete `VCF` specification with an explanation for each of the pieces of information in the file can be found [here](https://samtools.github.io/hts-specs/VCFv4.3.pdf). The most relevant sections for us are the header secion lines (start with `##`) and then the lines containing the variants. These will contain the following fields: 
| Col  | Field       | Brief description     |
| -----| ----------- | --------------------- | 
| 1    | CHROM       | Chromosome where the genetic variant was found   |
| 2    | POS         | Position in the chromosome where the genetic variant was found |
| 3    | ID          | SNP id |
| 4    | REF         | Reference allele |
| 5    | ALT         | Alternate allele |
| 6    | QUAL        | Variant quality |
| 7    | FILTER      | Filter string (`PASS` if it passsed all filters) |
| 8    | INFO        | Semicolon-separated series of variant additional information fields |
| 9    | GENOTYPE    | Genotype information (if present) | 

#### **Exploring the VCF file useing BCFTools** 
[BCFTools](http://samtools.github.io/bcftools/bcftools.html) is a suite of tools written in C that are quite efficient to manipulate files in the `VCF` format. Here we are going to see some of most useful commmands to manipulate the `VCF` we have just generated:

* Print the header section

        bcftools view -H SAMEA2569438.chr10.vcf.gz

You get:

        ##fileformat=VCFv4.1
        ##FILTER=<ID=PASS,Description="All filters passed">
        ##fileDate=20201117
        ##source=freeBayes v0.9.21
        ##reference=/hps/nobackup/production/reseq-info/ernesto/MASTER_COURSE/REFERENCE/Oryza_sativa.IRGSP-1.0.dna.toplevel.chr10.fa
        ##phasing=none
        ##commandline="freebayes -f /hps/nobackup/production/reseq-info/ernesto/MASTER_COURSE/REFERENCE/Oryza_sativa.IRGSP-1.0.dna.toplevel.chr10.fa /hps/nobackup/production/reseq-info/ernesto/MASTER_COURSE/ALIGNMENT/analysis_chr10/NEW_ALN_ONLY_CHR10_REF/POSTPROCESSING/SAMEA2569438.chr10.sorted.mark_duplicates.bam --ploidy 2"
        ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
        ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">
        ##INFO=<ID=DPB,Number=1,Type=Float,Description="Total read depth per bp at the locus; bases in reads overlapping / bases in haplotype">
        ##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
        ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
        ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">
        ##INFO=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count, with partial observations recorded fractionally">
        ##INFO=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observations, with partial observations recorded fractionally">
        ##INFO=<ID=PRO,Number=1,Type=Float,Description="Reference allele observation count, with partial observations recorded fractionally">
        ##INFO=<ID=PAO,Number=A,Type=Float,Description="Alternate allele observations, with partial observations recorded fractionally">
        ##INFO=<ID=QR,Number=1,Type=Integer,Description="Reference allele quality sum in phred">
        ##INFO=<ID=QA,Number=A,Type=Integer,Description="Alternate allele quality sum in phred">
        ##INFO=<ID=PQR,Number=1,Type=Float,Description="Reference allele quality sum in phred for partial observations">
        ##INFO=<ID=PQA,Number=A,Type=Float,Description="Alternate allele quality sum in phred for partial observations">
        ##INFO=<ID=SRF,Number=1,Type=Integer,Description="Number of reference observations on the forward strand">
        ##INFO=<ID=SRR,Number=1,Type=Integer,Description="Number of reference observations on the reverse strand">
        ##INFO=<ID=SAF,Number=A,Type=Integer,Description="Number of alternate observations on the forward strand">
        ##INFO=<ID=SAR,Number=A,Type=Integer,Description="Number of alternate observations on the reverse strand">
        ##INFO=<ID=SRP,Number=1,Type=Float,Description="Strand balance probability for the reference allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SRF and SRR given E(SRF/SRR) ~ 0.5, derived using Hoeffding's inequality">
        ##INFO=<ID=SAP,Number=A,Type=Float,Description="Strand balance probability for the alternate allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SAF and SAR given E(SAF/SAR) ~ 0.5, derived using Hoeffding's inequality">
        ##INFO=<ID=AB,Number=A,Type=Float,Description="Allele balance at heterozygous sites: a number between 0 and 1 representing the ratio of reads showing the reference allele to all reads, considering only reads from individuals called as heterozygous">
        ##INFO=<ID=ABP,Number=A,Type=Float,Description="Allele balance probability at heterozygous sites: Phred-scaled upper-bounds estimate of the probability of observing the deviation between ABR and ABA given E(ABR/ABA) ~ 0.5, derived using Hoeffding's inequality">
        ##INFO=<ID=RUN,Number=A,Type=Integer,Description="Run length: the number of consecutive repeats of the alternate allele in the reference genome">
        ##INFO=<ID=RPP,Number=A,Type=Float,Description="Read Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5, derived using Hoeffding's inequality">
        ##INFO=<ID=RPPR,Number=1,Type=Float,Description="Read Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5, derived using Hoeffding's inequality">
        ##INFO=<ID=RPL,Number=A,Type=Float,Description="Reads Placed Left: number of reads supporting the alternate balanced to the left (5') of the alternate allele">
        ##INFO=<ID=RPR,Number=A,Type=Float,Description="Reads Placed Right: number of reads supporting the alternate balanced to the right (3') of the alternate allele">
        ##INFO=<ID=EPP,Number=A,Type=Float,Description="End Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5, derived using Hoeffding's inequality">
        ##INFO=<ID=EPPR,Number=1,Type=Float,Description="End Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5, derived using Hoeffding's inequality">
        ##INFO=<ID=DPRA,Number=A,Type=Float,Description="Alternate allele depth ratio.  Ratio between depth in samples with each called alternate allele and those without.">
        ##INFO=<ID=ODDS,Number=1,Type=Float,Description="The log odds ratio of the best genotype combination to the second-best.">
        ##INFO=<ID=GTI,Number=1,Type=Integer,Description="Number of genotyping iterations required to reach convergence or bailout.">
        ##INFO=<ID=TYPE,Number=A,Type=String,Description="The type of allele, either snp, mnp, ins, del, or complex.">
        ##INFO=<ID=CIGAR,Number=A,Type=String,Description="The extended CIGAR representation of each alternate allele, with the exception that '=' is replaced by 'M' to ease VCF parsing.  Note that INDEL alleles do not have the first matched base (which is provided by default, per the spec) referred to by the CIGAR.">
        ##INFO=<ID=NUMALT,Number=1,Type=Integer,Description="Number of unique non-reference alleles in called genotypes at this position.">
        ##INFO=<ID=MEANALT,Number=A,Type=Float,Description="Mean number of unique non-reference allele observations per sample with the corresponding alternate alleles.">
        ##INFO=<ID=LEN,Number=A,Type=Integer,Description="allele length">
        ##INFO=<ID=MQM,Number=A,Type=Float,Description="Mean mapping quality of observed alternate alleles">
        ##INFO=<ID=MQMR,Number=1,Type=Float,Description="Mean mapping quality of observed reference alleles">
        ##INFO=<ID=PAIRED,Number=A,Type=Float,Description="Proportion of observed alternate alleles which are supported by properly paired read fragments">
        ##INFO=<ID=PAIREDR,Number=1,Type=Float,Description="Proportion of observed reference alleles which are supported by properly paired read fragments">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype">
        ##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
        ##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count">
        ##FORMAT=<ID=QR,Number=1,Type=Integer,Description="Sum of quality of the reference observations">
        ##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">
        ##FORMAT=<ID=QA,Number=A,Type=Integer,Description="Sum of quality of the alternate observations">
        ##bcftools_viewVersion=1.9+htslib-1.9
        ##bcftools_viewCommand=view -h SAMEA2569438.chr10.vcf; Date=Tue Nov 17 12:51:20 2020
        #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  unknown

* Print some SNPs:

        bcftools view -H -v snps SAMEA2569438.chr10.vcf.gz |less
You get:

        10      5893    .       C       T       55.2291 .       AB=0;ABP=0;AC=2;AF=1;AN=2;AO=2;CIGAR=1X;DP=2;DPB=2;DPRA=0;EPP=3.0103;EPPR=0;GTI=0;LEN=1;MEANALT=1;MQM=45.5;MQMR=0;NS=1;NUMALT=1;ODDS=7.37776;PAIRED=0;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=75;QR=0;RO=0;RPL=0;RPP=7.35324;RPPR=0;RPR=2;RUN=1;SAF=1;SAP=3.0103;SAR=1;SRF=0;SRP=0;SRR=0;TYPE=snp  GT:DP:RO:QR:AO:QA:GL    1/1:2:0:0:2:75:-6.72676,-0.60206,0
        10      5944    .       A       C       54.2577 .       AB=0;ABP=0;AC=2;AF=1;AN=2;AO=2;CIGAR=1X;DP=2;DPB=2;DPRA=0;EPP=7.35324;EPPR=0;GTI=0;LEN=1;MEANALT=1;MQM=47;MQMR=0;NS=1;NUMALT=1;ODDS=7.37776;PAIRED=0;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=76;QR=0;RO=0;RPL=1;RPP=3.0103;RPPR=0;RPR=1;RUN=1;SAF=1;SAP=3.0103;SAR=1;SRF=0;SRP=0;SRR=0;TYPE=snp    GT:DP:RO:QR:AO:QA:GL    1/1:2:0:0:2:76:-6.62961,-0.60206,0
        10      9454    .       G       T       37.074  .       AB=0.6;ABP=3.44459;AC=1;AF=0.5;AN=2;AO=3;CIGAR=1X;DP=5;DPB=5;DPRA=0;EPP=3.73412;EPPR=3.0103;GTI=0;LEN=1;MEANALT=1;MQM=51;MQMR=44;NS=1;NUMALT=1;ODDS=3.40948;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=102;QR=72;RO=2;RPL=0;RPP=9.52472;RPPR=7.35324;RPR=3;RUN=1;SAF=2;SAP=3.73412;SAR=1;SRF=1;SRP=3.0103;SRR=1;TYPE=snp     GT:DP:RO:QR:AO:QA:GL    0/1:5:2:72:3:102:-7.95337,0,-5.18999
        ...
* Print some INDELs

         bcftools view -H -v indels SAMEA2569438.chr10.vcf.gz |less
You get:

        10      16538   .       AGG     AG      0.014514        .       AB=0.666667;ABP=3.73412;AC=1;AF=0.5;AN=2;AO=2;CIGAR=1M1D1M;DP=3;DPB=2.33333;DPRA=0;EPP=7.35324;EPPR=5.18177;GTI=0;LEN=1;MEANALT=1;MQM=19;MQMR=32;NS=1;NUMALT=1;ODDS=5.77664;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=22;QR=40;RO=1;RPL=1;RPP=3.0103;RPPR=5.18177;RPR=1;RUN=1;SAF=1;SAP=3.0103;SAR=1;SRF=0;SRP=5.18177;SRR=1;TYPE=del       GT:DP:RO:QR:AO:QA:GL    0/1:3:1:40:2:22:-0.0253243,0,-2.23306
        10      52050   .       TAT     TT      7.82678 .       AB=0.285714;ABP=5.80219;AC=1;AF=0.5;AN=2;AO=2;CIGAR=1M1D1M;DP=7;DPB=6.33333;DPRA=0;EPP=7.35324;EPPR=6.91895;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=60;NS=1;NUMALT=1;ODDS=1.62193;PAIRED=1;PAIREDR=0.8;PAO=0;PQA=0;PQR=0;PRO=0;QA=78;QR=190;RO=5;RPL=2;RPP=7.35324;RPPR=3.44459;RPR=0;RUN=1;SAF=0;SAP=7.35324;SAR=2;SRF=2;SRP=3.44459;SRR=3;TYPE=del  GT:DP:RO:QR:AO:QA:GL    0/1:7:5:190:2:78:-5.29557,0,-15.3579
        10      53159   .       CAAAAAAAAAT     CAAAAAAAAT      57.2726 .       AB=0;ABP=0;AC=2;AF=1;AN=2;AO=2;CIGAR=1M1D9M;DP=2;DPB=1.81818;DPRA=0;EPP=3.0103;EPPR=0;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=0;NS=1;NUMALT=1;ODDS=7.37776;PAIRED=1;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=73;QR=0;RO=0;RPL=2;RPP=7.35324;RPPR=0;RPR=0;RUN=1;SAF=1;SAP=3.0103;SAR=1;SRF=0;SRP=0;SRR=0;TYPE=del  GT:DP:RO:QR:AO:QA:GL    1/1:2:0:0:2:73:-6.9311,-0.60206,0
        10      55082   .       CTGT    CTGTGT  44      .       AB=0;ABP=0;AC=2;AF=1;AN=2;AO=2;CIGAR=1M2I3M;DP=2;DPB=3;DPRA=0;EPP=7.35324;EPPR=0;GTI=0;LEN=2;MEANALT=1;MQM=60;MQMR=0;NS=1;NUMALT=1;ODDS=7.37776;PAIRED=1;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=59;QR=0;RO=0;RPL=1;RPP=3.0103;RPPR=0;RPR=1;RUN=1;SAF=1;SAP=3.0103;SAR=1;SRF=0;SRP=0;SRR=0;TYPE=ins        GT:DP:RO:QR:AO:QA:GL    1/1:2:0:0:2:59:-5.60384,-0.60206,0
        ...

* Print variants for a specific region

To fetch the variants located in a specific genomic region you need first to index the VCF, for this use `bcftools index`:

        bcftools index SAMEA2569438.chr10.vcf.gz

And then you can use `bcftools view` with the `-r` option to query a specific region:

        bcftools view -r 10:10000-20000 SAMEA2569438.chr10.vcf.gz |less

* Print some basic stats for the VCF file

We can use the `stats` command to generate a basic report on the number of variants in a VCF file:

        bcftools stats SAMEA2569438.chr10.vcf.gz |grep ^SN

We pipe the output of the `stats` command to the UNIX `grep` command to print only the lines starting with `SN`:

        SN      0       number of samples:      1
        SN      0       number of records:      105805
        SN      0       number of no-ALTs:      0
        SN      0       number of SNPs: 88827
        SN      0       number of MNPs: 7910
        SN      0       number of indels:       8340
        SN      0       number of others:       1092
        SN      0       number of multiallelic sites:   458
        SN      0       number of multiallelic SNP sites:       56

* Selecting the multiallelic SNPs

Use the following command to select the multiallelic SNPs:

        bcftools view -m3 -v snps SAMEA2569438.chr10.vcf.gz

And you get:

        10      134998  .       GCCA    CCCG,GCCG       185.768 .       AB=0.375,0.5;ABP=4.09604,3.0103;AC=1,1;AF=0.5,0.5;AN=2;AO=3,4;CIGAR=1X2M1X,3M1X;DP=8;DPB=8.75;DPRA=0,0;EPP=3.73412,3.0103;EPPR=5.18177;GTI=0;LEN=4,1;MEANALT=2,2;MQM=52.6667,60;MQMR=27;NS=1;NUMALT=2;ODDS=17.466;PAIRED=1,1;PAIREDR=1;PAO=1,0;PQA=39,0;PQR=0;PRO=0;QA=120,153;QR=40;RO=1;RPL=2,2;RPP=3.73412,3.0103;RPPR=5.18177;RPR=1,2;RUN=1,1;SAF=0,2;SAP=9.52472,3.0103;SAR=3,2;SRF=0;SRP=5.18177;SRR=1;TYPE=complex,snp   GT:DP:RO:QR:AO:QA:GL    1/2:8:1:40:3,4:120,153:-23.2698,-10.5583,-11.4412,-11.0232,0,-11.9036
        10      135116  .       G       A,C     106.967 .       AB=0.6,0.4;ABP=3.44459,3.44459;AC=1,1;AF=0.5,0.5;AN=2;AO=3,2;CIGAR=1X,1X;DP=5;DPB=5;DPRA=0,0;EPP=9.52472,7.35324;EPPR=0;GTI=0;LEN=1,1;MEANALT=2,2;MQM=60,55;MQMR=0;NS=1;NUMALT=2;ODDS=4.02222;PAIRED=0.666667,1;PAIREDR=0;PAO=0,0;PQA=0,0;PQR=0;PRO=0;QA=119,80;QR=0;RO=0;RPL=2,0;RPP=3.73412,7.35324;RPPR=0;RPR=1,2;RUN=1,1;SAF=1,2;SAP=3.73412,7.35324;SAR=2,0;SRF=0;SRP=0;SRR=0;TYPE=snp,snp GT:DP:RO:QR:AO:QA:GL    1/2:5:0:0:3,2:119,80:-16.7553,-6.96125,-6.05816,-10.1914,0,-9.58935
        10      297744  .       ATT     AG,AGT  167.268 .       AB=0.75,0.25;ABP=7.35324,7.35324;AC=1,1;AF=0.5,0.5;AN=2;AO=6,2;CIGAR=1M1D1X,1M1X1M;DP=8;DPB=6;DPRA=0,0;EPP=4.45795,3.0103;EPPR=0;GTI=0;LEN=2,1;MEANALT=2,2;MQM=60,60;MQMR=0;NS=1;NUMALT=2;ODDS=1.25872;PAIRED=0.833333,1;PAIREDR=0;PAO=0,0;PQA=0,0;PQR=0;PRO=0;QA=200,75;QR=0;RO=0;RPL=2,2;RPP=4.45795,7.35324;RPPR=0;RPR=4,0;RUN=1,1;SAF=4,1;SAP=4.45795,3.0103;SAR=2,1;SRF=0;SRP=0;SRR=0;TYPE=complex,snp     GT:DP:RO:QR:AO:QA:GL    1/2:8:0:0:6,2:200,75:-22.6737,-6.51804,-4.71186,-16.52,0,-15.918
        10      342811  .       TA      TCG,TC  180.859 .       AB=0.75,0.25;ABP=7.35324,7.35324;AC=1,1;AF=0.5,0.5;AN=2;AO=6,2;CIGAR=1M1I1X,1M1X;DP=8;DPB=11;DPRA=0,0;EPP=3.0103,3.0103;EPPR=0;GTI=0;LEN=2,1;MEANALT=2,2;MQM=60,60;MQMR=0;NS=1;NUMALT=2;ODDS=1.03749;PAIRED=1,1;PAIREDR=0;PAO=0,0;PQA=0,0;PQR=0;PRO=0;QA=213,74;QR=0;RO=0;RPL=3,2;RPP=3.0103,7.35324;RPPR=0;RPR=3,0;RUN=1,1;SAF=4,1;SAP=4.45795,3.0103;SAR=2,1;SRF=0;SRP=0;SRR=0;TYPE=complex,snp       GT:DP:RO:QR:AO:QA:GL    1/2:8:0:0:6,2:213,74:-23.7598,-6.42196,-4.61578,-17.7038,0,-17.1017

#### **Filtering the artifactual variants**

The process for identifiying variants is not perfect, and FreeBayes and in general all tools used to identify variants will report variants that are not real. These artifactual variants must be idenfitied and flagged so that users or tools using them do not take them into account, or treat them with caution in any subsequent analysis.

There are several filtering tools and strategies available for variant filtering, with varying degrees of complexity and sophistication. However, in this course we will use a very simple, yet effective approach, which consists of using the quality value assigned by FreeBayes as a proxy to estimate the likelihood of a variant being real. The lower the quality value, the less likely it is that a variant is real.

In this course, we will use `bcftools filter` with a hard cut-off value of `<=1` to flag the variants that have a low quality. For this enter the following in your terminal:

        bcftools filter -sQUALFILTER -e'QUAL<1' SAMEA2569438.chr10.vcf.gz -o SAMEA2569438.chr10.filt.vcf.gz -Oz

Where the string passed using the `-s` option will set the label used for the filtered lines in the 7th column of the VCF and `-Oz` is used in `bcftools` for generating the output VCF in a compressed format.

Now, use 'bcftools view' to check that the 7th column has 2 new labels: `QUALFILTER` and `PASS`.

        bcftools view -H SAMEA2569438.chr10.filt.vcf.gz |less

`-H` is used to skip the header section and only print the data lines:

        10      5893    .       C       T       55.2291 PASS    AB=0;ABP=0;AC=2;AF=1;AN=2;AO=2;CIGAR=1X;DP=2;DPB=2;DPRA=0;EPP=3.0103;EPPR=0;GTI=0;LEN=1;MEANALT=1;MQM=45.5;MQMR=0;NS=1;NUMALT=1;ODDS=7.37776;PAIRED=0;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=75;QR=0;RO=0;RPL=0;RPP=7.35324;RPPR=0;RPR=2;RUN=1;SAF=1;SAP=3.0103;SAR=1;SRF=0;SRP=0;SRR=0;TYPE=snp  GT:DP:RO:QR:AO:QA:GL    1/1:2:0:0:2:75:-6.72676,-0.60206,0
        10      9569    .       C       T       8.11386 PASS    AB=0.5;ABP=3.0103;AC=1;AF=0.5;AN=2;AO=2;CIGAR=1X;DP=4;DPB=4;DPRA=0;EPP=7.35324;EPPR=3.0103;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=43.5;NS=1;NUMALT=1;ODDS=1.69704;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=65;QR=81;RO=2;RPL=2;RPP=7.35324;RPPR=3.0103;RPR=0;RUN=1;SAF=2;SAP=7.35324;SAR=0;SRF=2;SRP=7.35324;SRR=0;TYPE=snp     GT:DP:RO:QR:AO:QA:GL    0/1:4:2:81:2:65:-4.96917,0,-6.07964
        10      9682    .       C       T       0.00381086      QUALFILTER      AB=0.222222;ABP=9.04217;AC=1;AF=0.5;AN=2;AO=2;CIGAR=1X;DP=9;DPB=9;DPRA=0;EPP=7.35324;EPPR=10.7656;GTI=0;LEN=1;MEANALT=1;MQM=25;MQMR=22.5714;NS=1;NUMALT=1;ODDS=7.03816;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=67;QR=262;RO=7;RPL=0;RPP=7.35324;RPPR=18.2106;RPR=2;RUN=1;SAF=2;SAP=7.35324;SAR=0;SRF=6;SRP=10.7656;SRR=1;TYPE=snp GT:DP:RO:QR:AO:QA:GL    0/1:9:7:262:2:67:-1.90251,0,-11.6176
        10      10021   .       G       A       51.1653 PASS    AB=0.75;ABP=5.18177;AC=1;AF=0.5;AN=2;AO=3;CIGAR=1X;DP=4;DPB=4;DPRA=0;EPP=3.73412;EPPR=5.18177;GTI=0;LEN=1;MEANALT=1;MQM=30;MQMR=19;NS=1;NUMALT=1;ODDS=6.50714;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=104;QR=33;RO=1;RPL=2;RPP=3.73412;RPPR=5.18177;RPR=1;RUN=1;SAF=1;SAP=3.73412;SAR=2;SRF=1;SRP=5.18177;SRR=0;TYPE=snp  GT:DP:RO:QR:AO:QA:GL    0/1:4:1:33:3:104:-6.16983,0,-0.679135
        ......

We can also print only the variants that have been filtered by doing:

        bcftools view -f QUALFILTER SAMEA2569438.chr10.filt.vcf.gz |less

And you get:

        10      9682    .       C       T       0.00381086      QUALFILTER      AB=0.222222;ABP=9.04217;AC=1;AF=0.5;AN=2;AO=2;CIGAR=1X;DP=9;DPB=9;DPRA=0;EPP=7.35324;EPPR=10.7656;GTI=0;LEN=1;MEANALT=1;MQM=25;MQMR=22.5714;NS=1;NUMALT=1;ODDS=7.03816;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=67;QR=262;RO=7;RPL=0;RPP=7.35324;RPPR=18.2106;RPR=2;RUN=1;SAF=2;SAP=7.35324;SAR=0;SRF=6;SRP=10.7656;SRR=1;TYPE=snp GT:DP:RO:QR:AO:QA:GL    0/1:9:7:262:2:67:-1.90251,0,-11.6176
        10      12889   .       T       A       0.884523        QUALFILTER      AB=0.222222;ABP=9.04217;AC=1;AF=0.5;AN=2;AO=2;CIGAR=1X;DP=9;DPB=9;DPRA=0;EPP=7.35324;EPPR=3.32051;GTI=0;LEN=1;MEANALT=1;MQM=48.5;MQMR=41.1429;NS=1;NUMALT=1;ODDS=1.4877;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=75;QR=263;RO=7;RPL=0;RPP=7.35324;RPPR=10.7656;RPR=2;RUN=1;SAF=0;SAP=7.35324;SAR=2;SRF=3;SRP=3.32051;SRR=4;TYPE=snp        GT:DP:RO:QR:AO:QA:GL    0/1:9:7:263:2:75:-4.31305,0,-20.0551
        10      16538   .       AGG     AG      0.014514        QUALFILTER      AB=0.666667;ABP=3.73412;AC=1;AF=0.5;AN=2;AO=2;CIGAR=1M1D1M;DP=3;DPB=2.33333;DPRA=0;EPP=7.35324;EPPR=5.18177;GTI=0;LEN=1;MEANALT=1;MQM=19;MQMR=32;NS=1;NUMALT=1;ODDS=5.77664;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=22;QR=40;RO=1;RPL=1;RPP=3.0103;RPPR=5.18177;RPR=1;RUN=1;SAF=1;SAP=3.0103;SAR=1;SRF=0;SRP=5.18177;SRR=1;TYPE=del       GT:DP:RO:QR:AO:QA:GL    0/1:3:1:40:2:22:-0.0253243,0,-2.23306
        ....

* How many variants have been filtered?

We can use the `stats` command together with the `-f` option to generate a report taken into account only the filtered variants:

        bcftools stats -f QUALFILTER SAMEA2569438.chr10.filt.vcf.gz |grep ^SN

And you get:

        SN      0       number of samples:      1
        SN      0       number of records:      5498
        SN      0       number of no-ALTs:      0
        SN      0       number of SNPs: 4851
        SN      0       number of MNPs: 416
        SN      0       number of indels:       211
        SN      0       number of others:       34
        SN      0       number of multiallelic sites:   25
        SN      0       number of multiallelic SNP sites:       8



#### **Exploring the identified variants using IGV**

The Integrative Genomics Viewer [IVG](http://software.broadinstitute.org/software/igv/) is a useful interactive tool that can be used to explore visually your genomic data. We are going to use it here to display the variants we have identified. In this example we will explore the variants in a specific region in chromosome 10.

First, open the `igv` viewer by going to your terminal and typing:

        igv

You will need to load in `IGV` the FASTA file containing the chromosome 10 sequence for rice, as this sequence is not included by default in `IGV`:

![load_genome_igv](https://www.ebi.ac.uk/~ernesto/IGSR/masters_IAMZ_jan2020/load_genome.png)

Look for you file and open it.

Now, load the `GTF` file containing the rice gene annotations for chromosome 10:

![load_annotation_igv](https://www.ebi.ac.uk/~ernesto/IGSR/masters_IAMZ_jan2020/load_annotation.png)

![annotation_view_igv](https://www.ebi.ac.uk/~ernesto/IGSR/masters_IAMZ_jan2020/annotation_view.png)

You will see the new track with genes annotated in chromosome 10

Now, load the `VCF` file containing the variants:
![load_variants_igv](https://www.ebi.ac.uk/~ernesto/IGSR/masters_IAMZ_jan2020/load_alignments.png)

![load_variants1_igv](https://www.ebi.ac.uk/~ernesto/IGSR/masters_IAMZ_jan2020/load_variants1.png)

Now, you can click on a particular variant (red vertical bar) to display information such as:

* Position
* Reference and alternate alleles
* Type of variant (SNPs or INDEL)
* Variant quality
* Filtering status
* Variant attributes
* etc ...

![load_variants2_igv](https://www.ebi.ac.uk/~ernesto/IGSR/masters_IAMZ_jan2020/load_variants2.png)

You can also click on the blue vertical bar to display genotype information and attributes:

![load_variants3_igv](https://www.ebi.ac.uk/~ernesto/IGSR/masters_IAMZ_jan2020/load_variants3.png)

Let's examine in detail a SNP and an INDEL variant. For this, enter the following genomic coordinate in your navigate box:

        10:16,532-16,568

Click on the variant on the left side of the screen:

![indel_example1_igv](https://www.ebi.ac.uk/~ernesto/IGSR/masters_IAMZ_jan2020/indel_example1.png)

You can see that this is a 1bp deletion (AGG->AG) that have an alternate allele_count=1. Which means it is an heterozygous variant, this can be confirmed by clicking on the genotype information bar:

![indel_example2_igv](https://www.ebi.ac.uk/~ernesto/IGSR/masters_IAMZ_jan2020/indel_example2.png)

Now, click on the variant on the right side of the screen:

![snp_example1_igv](https://www.ebi.ac.uk/~ernesto/IGSR/masters_IAMZ_jan2020/snp_example1.png)

You can see that this is a single nucleotide substitution (T->A) that have an alternate allele_count=2. Which means it is a homozygous variant, this can be confirmed by clicking on the genotype information bar:

![snp_example2_igv](https://www.ebi.ac.uk/~ernesto/IGSR/masters_IAMZ_jan2020/snp_example2.png)




