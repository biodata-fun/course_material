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

## Variant calling
Variant calling is the process to identify variants from sequence data ([see](https://www.ebi.ac.uk/training-beta/online/courses/human-genetic-variation-introduction/variant-identification-and-analysis/#:~:text=What%20is%20variant%20calling%3F,creating%20BAM%20or%20CRAM%20files.)).
It starts with the sequencing data in the `FASTQ` format, these sequence reads need to aligned to the reference genome following the procedure explained in the first part of this course, that generated an alignment in the `BAM` file format. Once we have an alignment or multiple alignment files we can use a variant discovery tool to identify the germline variants. There are multiple variant calling tools available, the ones we have more experience with in our group are [SAMTools mpileup](http://samtools.github.io/bcftools/bcftools.html#mpileup), the [GATK suite](https://gatk.broadinstitute.org/hc/en-us) and [FreeBayes](https://github.com/ekg/freebayes). In this course we are going to use FreeBayes, as it is sensitive, precise and relatively simple to use. 

### Freebayes
FreeBayes is a haplotype-based variant detector, that uses a joint genotyping method capable of reporting variants on a single sample or on a cohort of multiple samples. It's going to be capable of detecting SNPs (single nucleotide polymorphisms), indels (short insertions and deletions) and MNPs (multi-nucleotide polymorphisms)

#### **Reference Genome**
FreeBayes needs the reference sequence in the `FASTA` format. In this section of the course we are going to use the same chromosome 10 sequence extracted from the *Oryza_sativa* (rice) genome that we used for the alignment part of the course.

#### **Using FreeBayes**
To run Freebayes you need to specify the ploidy of the genome being analysed, the FASTA reference sequence used for the alignment and the analhysis-ready BAM generated in the first section of the course. Once you have these prepared enter the following command in your terminal:

        freebayes -f /path/to/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --ploidy 2 /path/to/SAMEA2569438.chr10.sorted.mark_duplicates.bam |bgzip -c > SAMEA2569438.chr10.vcf.gz &

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

#### **Filtering the artifactual variants**

The process for identifiying variants is not perfect, and FreeBayes and in general all tools used to identify variants will report variants that are not real. These artifactual variants must be idenfitied and flagged so that users or tools using them do not take them into account, or treat them with caution in any subsequent analysis.

There are several filtering tools and strategies available for variant filtering, with varying degrees of complexity and sophistication. However, in this course we will use a very simple, yet effective approach, which consists of using the quality value assigned by FreeBayes as a proxy to estimate the likelihood of a variant being real. The lower the quality value, the less likely it is that a variant is real.

In this course, we will use `bcftools` XXXX with a hard cut-off value of X to flag the variants that have a low quality. 
#### **Exploring the identified variants using IGV**

