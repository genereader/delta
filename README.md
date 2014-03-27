# DELTA
a Distal Enhancer Locating Tool based on AdaBoost and shape features of chromatin modifications

## Options
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -c CHIP_BEDS, --chip_bed=CHIP_BEDS
                        ChIP-seq bed file of histone modifications
  -R, --read            Read existing training and predicting data instead of
                        generate from ChIP-seq (default: False)
  -E ENHANCER, --enhancer=ENHANCER
                        BED file containing the enhancer loci
  -P PROMOTER, --promoter=PROMOTER
                        BED file containing the promoter loci
  -g GENOME, --genome=GENOME
                        Genome assembly should be one of the followings: dm3,
                        mm9, hg17, hg18, hg19
  -b BIN_SIZE, --bin_size=BIN_SIZE
                        Length of dividing bins (default: 100)
  -s STEP_SIZE, --step=STEP_SIZE
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
