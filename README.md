# EBFilter (Empirical Bayesian Mutation Filtering)

## Introduction

EBFilter is a software for filtering false poistive somatic mutations in cancer genome sequencing data analysis,
using the framework 


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

### Python
Python (>= 2.7), `pysam`, `scipy`, `numpy`, `pyVCF` packages

### Software
tabix, bgzip, blat

## Install

```
git clone https://github.com/friend1ws/EBFilter.git
cd EBFilter
python setup.py build
python setup.py install
```

## Preparation

All the input should be indexed bam files.  
- **target somatic mutation candidats**: the somatic mutation candidates (should be .vcf format).
- **target tumor sample**: the .bam file of the target tumor sample.
- **list of normal reference samples**: the list of paths to .bam files for non-paired normal reference samples. Please name the text file as you like (e.g., myNormalRef.txt), and list the paths of .bam files as follows:  

	/home/yshira/ngs/data/sequence/normalreference1.bam
	/home/yshira/ngs/data/sequence/normalreference2.bam
	...
	/home/yshira/ngs/data/sequence/normalreference10.bam


## Commands

    EBFilter target.vcf target.bam list_normal_sample.txt output.vcf
  
- **-q**: The threshold of mapping quality. The short reads wholse mapping quality are smaller than this value are skipped [defaut: 20].
- **-Q**: The threshold of base quality. The bases whose base quality are smaller than this value are skipped [default: 15].
- **-t**: The number of threads [default: 1].

When finished, the output.vcf includes the score calculated by validation by beta-binomial sequencing error model (**EB** tag).



## Test run

We provide a set of test data files in the EBCall-master/testdata directory and the result in the EBCall-master/testresult directory.   
Edit EBCall/testdata/list_normal_sample.txt to adjust the paths to the EBCall-master directory.

	/home/your_username/EBCall-master/testdata/normalreference1.bam
	/home/your_username/EBCall-master/testdata/normalreference2.bam
	...
	/home/your_username/EBCall-master/testdata/normalreference10.bam

Type the following command after setup EBCall and compiling C++ programs. 

	sh ebCall_v2.sh testdata/tumor.bam testdata/normal.bam testout testdata/list_normal_sample.txt

Result is stored under the testout directory.


