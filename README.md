[![Build Status](https://travis-ci.org/ken0-1n/EBFilter.svg?branch=devel)](https://travis-ci.org/ken0-1n/EBFilter) 

# EBFilter (Empirical Bayesian Mutation Filtering)

## Introduction

EBFilter is a software for filtering false poistive somatic mutations in cancer genome sequencing data analysis.
EBFilter accepts a list of candidates somatic mutation alreadly narrowed down to some extents, and performs the main step of [EBCall](https://github.com/friend1ws/EBCall), 

1. estimate the parameters of the beta-binomial sequencing error model using multiple non-matched contorol sequencing data at the position of interest
2. get the predictive mismatch ratio derived from the above estimation and compare it with the observed mismatch ratio of tumor samples
3. if the mismatch ratio of the tumor sample is significantly deviated from the predicted mismatch ratio, then we idenitfy it as highly-likely somatic mutation

Therefore, you can use EBFilter along with your own mutation calling program 
or some popular mutation callers (e.g., MuTect, VarScan 2),
which we believe will reduce large parts of false positives with a slight computational cost.

## Paper

We would like you to kindly cite the following paper when you use this software; 

"[An empirical Bayesian framework for mutation detection from cancer genome sequencing data](http://nar.oxfordjournals.org/content/41/7/e89.long)", Shiraishi et al.,  
Nucleic Acids Research, 2013.


# Motivation

One major source of false positive somatic mutation in cancer genome sequencing data analysis is,
"broken symetry of mismatches between tumor and normal sequencing data that randomly occurs at the error-prone sites"
(see e.g., Figure 1, Shiraishi et al., NAR, 2013).
Therefore, one of the most important techniques for reducing false positive mutations is to check many control sequencing data (including non-matched normal samples), to see whether the candidates are artifacts produced at the error-prone sites, and many somatic mutation calling pipeline adopt this strategy.

EBCall uses 

1. estimate the parameters of the sequencing error model using multiple non-matched contorol sequencing data at the position of interest
2. get the predictive mismatch ratio obtained from the above estimation and compare it with the observed mismatch ratio of tumor samples
3. if the mismatch ratio of the tumor sample is significantly deviated from the predicted mismatch ratio, then we idenitfy it as highly-likely somatic mutation


EBCall is actually integrating many other steps (Fisher's exact test), and could not perfrom purely the above beta-binomial filteringt step.
Therefore, we decided to implement a software which just perform beta-binomial step for already collected candidates of somatic mutations.

Therefore, you can use this software after performing popular mutation callers (e.g., mutects, VarScan2,and so on),
or your own inhouse mutation caling program, which we believe will reduce large parts of false positives.



## Dependency

### Software
[samtools](http://www.htslib.org/)

### Python
EBFilter >= 0.2.2
Python (>= 3.7), `pysam`, `scipy`, `numpy`, `vcfpy` packages

EBFilter <= 0.2.1
Python (>= 2.7), `pysam`, `scipy`, `numpy`, `pyVCF` packages


## Install

```
git clone https://github.com/friend1ws/EBFilter.git
cd EBFilter
python setup.py build
python setup.py install
```

## Preparation
- add path to samtools.
- **target somatic mutation candidats**: the somatic mutation candidates (should be .vcf format).
- **target tumor sample**: the indexed bam file of the target tumor sample.
- **list of normal reference samples**: the list of paths to the indexed bam files for non-paired normal reference samples. Please name the text file as you like (e.g., myNormalRef.txt), and list the paths of .bam files as follows:  

		/home/yshira/ngs/data/sequence/normalreference1.bam
		/home/yshira/ngs/data/sequence/normalreference2.bam
		...
		/home/yshira/ngs/data/sequence/normalreference10.bam

## Commands
    EBFilter [-h] [--version] [-f {vcf,anno}] [-t thread_num]
                [-q mapping_qual_thres] [-Q base_qual_thres] [-ff filter_flags]
                [--loption] [--region REGION] [--debug]
                target.vcf target.bam controlBam_list.txt output.vcf
- **-f**: input mutation data format (indexed vcf format or annovar format) [default: vcf]
- **-q**: The threshold of mapping quality. The short reads wholse mapping quality are smaller than this value are skipped [default: 20].
- **-Q**: The threshold of base quality. The bases whose base quality are smaller than this value are skipped [default: 15].
- **-t**: The number of threads [default: 1].
- **--ff**: skip reads with mask bits set. [default: UNMAP,SECONDARY,QCFAIL,DUP] 
- **--loption**: If this option is on, -l option in samtools mpileup is used. In the default settings, EBFilter calculate the bases for each option repeatedly using samtools mpileup with -r option. This is suitable for investigating small number of mutation, However, using --loption is highly recomended for large number of mutations and will be effective when combining --region option below.
- **--region**: speficify genomic region for investigation.
- **--debug**: do not delete intermediate file (for mainly debugging).

When finished, the output.vcf includes the score calculated by validation by beta-binomial sequencing error model (**EB** tag).



## Test run

We provide a set of test data files in the testdata directory.
Type the following command after installing EBFilter:
	
	EBFilter testdata/input.vcf.gz testdata/tumor.bam testdata/list_normal_sample.txt output.vcf
	
Or for the annovar format:

	EBFilter -f anno testdata/input.anno testdata/tumor.bam testdata/list_normal_sample.txt output.anno
	
Then, compare the result with the golden-standard output in the testdata directory

	testdata/output.golden.vcf
	testdata/output.golden.anno
	



