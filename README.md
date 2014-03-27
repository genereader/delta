# DELTA
a Distal Enhancer Locating Tool based on AdaBoost and shape features of chromatin modifications

## Introduction
In the post-genomic era, accurate functional annotation of genomic sequences, especially DNA regulatory elements, has become an urgent need to understand of the complex mechanisms of human genome. Recent large-scale chromatin states mapping efforts revealed characteristic chromatin modification signatures for various types of functional DNA elements, promoting the emergence of a series of supervised and unsupervised-based methods for DNA elements prediction. However, these methods suffered from two major issues: (1) incomplete feature extraction of chromatin signatures and (2) inconsistency of histone modification importance between cell types. Here, we address these issues by introducing three important shape parameters of chromatin modifications from probability theory. We developed a novel method DELTA (a Distal Enhancer Locating Tool based on AdaBoost and shape features of chromatin modifications). We not only show that DELTA outperformed previous methods in different datasets, but also find that histone modification importance for enhancer prediction in different cell types are highly correlated. In summary, our results give insight into the consistency of variable importance of chromatin modifications across cell types and provide an accurate tool for enhancer prediction based on chromatin signatures.

## Install
Please check the file 'INSTALL' in the distribution.

## Usage
	Usage: macs14 <-t tfile> [-n name] [-g genomesize] [options]
	Example: macs14 -t ChIP.bam -c Control.bam -f BAM -g h -n test -w --call-subpeaks

	--version
										show program's version number and exit
	-h, --help
										show this help message and exit
	-c CHIP_BEDS, --chip_bed=CHIP_BEDS
										ChIP-seq bed file of histone modifications
	-R, --read
										Read existing training and predicting data instead of generate from ChIP-seq (default: False)
	-E ENHANCER, --enhancer=ENHANCER
										BED file containing the enhancer loci
	-P PROMOTER, --promoter=PROMOTER
										BED file containing the promoter loci
	-g GENOME, --genome=GENOME
										Genome assembly should be one of the followings: dm3, mm9, hg17, hg18, hg19
	-b BIN_SIZE, --bin_size=BIN_SIZE
										Length of dividing bins (default: 100)
	-s STEP_SIZE, --step=STEP_SIZE
										Step size of sliding window, should be integer times of bin size (default: 2000)
	-w WIN_SIZE, --window_size=WIN_SIZE
										Length of sliding window, should be integer times of bin size (default: 4000)
	--iteration_number=ITER_NUM
										Number of iteration for AdaBoost (default: 100)
	--pvalue_threshold=P_THRES
										P-value threshold for enhancer prediction (default: 0.5)
	-o OUTPUT, --output=OUTPUT
										Output file name (default output file is "predicted_enhancer.bed")
## Parameters

## Notes
