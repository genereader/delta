# DELTA
a Distal Enhancer Locating Tool based on AdaBoost and shape features of chromatin modifications

## Introduction
Accurate identification of DNA regulatory elements has become an urgent need in the post-genomic era. Recently, a number of machine learning-based approaches have emerged to predict DNA elements based on chromatin signatures. Two major constraints to the development of this area include (1) incomplete feature extraction of chromatin signatures and (2) lack of multicell models for element prediction across cell types. Here, we introduce a set of shape parameters of chromatin modifications distributed throughout DNA elements. We show these parameters show high consistency in variable importance across cell types. By training on a combined dataset of p300- binding sites in six cell types, we construct a multicell model for systemic identification of functional enhancers based on chromatin modifications. We show our model can predict enhancers with high validation rate and low misclassification rate not only in six training cell types, but also in six independent test cell types. Although the aim of our multicell model is to make general and cell-type- independent predictions across multiple cell types, we found it still outperformed the state-of-art single-cell enhancer prediction approaches that are specifically trained and optimized on individual cell type. Our study presents a general framework for identifying DNA elements from multicell epigenetic data.

## Install
Please check the file 'INSTALL' in the distribution.

## Usage
	Usage: delta.py [-c chip_files] [-P promoter_loci] [-E enhancer_loci] [options]

	Example: delta.py -c H3K4me1.bed,H3K4me3.bed,H3K27ac.bed -E p300.bed -P tss.bed -g hg19

	--version
										Show program's version number and exit
	-h, --help
										Show this help message and exit
	-c CHIP_BEDS, --chip_bed=CHIP_BEDS
										ChIP-seq bed file of histone modifications
	-E ENHANCER, --enhancer=ENHANCER
										BED file containing the enhancer loci
	-P PROMOTER, --promoter=PROMOTER
										BED file containing the promoter loci
	-R, --read
										Read existing training and predicting data instead of 
										generate from ChIP-seq (default: False)
	-g GENOME, --genome=GENOME
										Genome assembly should be one of the followings: dm3, 
										mm9, hg17, hg18, hg19
	-b BIN_SIZE, --bin_size=BIN_SIZE
										Length of dividing bins (default: 100)
	-s STEP_SIZE, --step_size=STEP_SIZE
										Step size of sliding window, should be integer times 
										of bin size (default: 2000)
	-w WIN_SIZE, --window_size=WIN_SIZE
										Length of sliding window, should be integer times of 
										bin size (default: 4000)
	--iteration_number=ITER_NUM
										Number of iteration for AdaBoost (default: 100)
	--pvalue_threshold=P_THRES
										P-value threshold for enhancer prediction (default: 
										0.5)
	-o OUTPUT, --output=OUTPUT
										Output file name (default output file is 
										"predicted_enhancer.bed")
## Parameters
-c / --chip_bed

ChIP-seq files contain chromatin modifications mapping data. User should provide ChIP-seq files separated by comma, e.g. H3K4me1.bed,H3K4me3.bed,H3K27ac.bed.

The BED format is defined in "http://genome.ucsc.edu/FAQ/FAQformat#format1".

-R / --read

The "-R" option lets user read existing training and predicting data instead generate them from ChIP-seq files, which would be a time consuming process. WARNING: Use with care!!!, wrong training and predicting data could be load.

-s / --step

Step size of sliding window, should be integer times of bin size. A good step size should at least half length of the window size, because a small step size will produce redundant predictions as well as increase computing time.

--pvalue_threshold

P-value threshold for enhancer prediction. User could adjust number of predictions by tuning this parameter. 

## Output files

1.predicted_enhancer.bed is a BED format file containing the predicted enhancers. User should be aware that if the step size is smaller than window size, the predicted enhancers may be redundant. _uniq_ command should be used in this situation to remove repetitive predictions.

2.adaboost.R is a R script generated by delta.py for executing AdaBoost algorithm.

3.tmp_dir is a directory contains temporary files created by delta.py. It should not be removed until the entire training and prediction is done.

## License

Source code of DELTA is freely available for academic use. For commercial license please contact Dr. Chenggang Zhang (zhangcg@bmi.ac.cn).
