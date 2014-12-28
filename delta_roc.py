#!/usr/bin/env python

import os
import time
import subprocess
import numpy as np
from optparse import OptionParser
from pybedtools import *
from lib.features import *
from lib.process import *

def set_optparser():
	'''Options setter
	'''
	usage = '''usage: %prog [options]
	'''
	optparser = OptionParser(version='%prog 1.0.1',usage=usage,add_help_option=False)
	optparser.add_option("-h","--help",action="help",
		help="show this help message and exit")
	optparser.add_option('-c','--chip_bed',dest='chip_beds',type='string',
		help='ChIP-seq bed file of histone modifications', default='NA')
	optparser.add_option('-R','--read', action="store_true", dest='read',
		help='Read existing training and predicting data instead of generate from ChIP-seq (default: False)',default=False)
	optparser.add_option('-E','--enhancer',dest='enhancer',type='string',
		help='BED file containing the enhancer loci',default='NA')
	optparser.add_option('-P','--promoter',dest='promoter',type='string',
		help='BED file containing the promoter loci',default='NA')
	optparser.add_option('-B','--background_size',dest='back_size',type='int',
		help='Number of random genomic regions distal to known TSS (default: 10000)',default=10000)
	optparser.add_option('-g','--genome',dest='genome',type='string',
		help='Genome assembly should be one of the followings: dm3, mm9, hg17, hg18, hg19',default='hg19')
	optparser.add_option('-b','--bin_size',dest='bin_size',type='int',
		help='Length of dividing bins (default: 100)',default=100)
	optparser.add_option('-s','--step',dest='step_size',type='int',
		help='Step size of sliding window, should be integer times of bin size (default: 2000)',default=2000)
	optparser.add_option('-w','--window_size',dest='win_size',type='int',
		help='Length of sliding window, should be integer times of bin size (default: 4000)',default=4000)
	optparser.add_option('--iteration_number',dest='iter_num',type='int',
		help='Number of iteration for AdaBoost (default: 100)',default=100)
	optparser.add_option('--pvalue_threshold',dest='p_thres',type='float',
		help='P-value threshold for enhancer prediction (default: 0.5)',default=0.5)
	optparser.add_option('-o','--output',dest='output',type='string',
		help='Output file name (default output file is "predicted_enhancer.bed")',default='predicted_enhancer.bed')
	return optparser

def main():
	# Parsing options
	options,args = set_optparser().parse_args()
	tmp_dir = 'tmp_dir'
	
	# Path setting
	real_dir = os.path.dirname(os.path.realpath(__file__))
	
	# Generating training and predicting data
	if options.read == False:
		chipNames = options.chip_beds.split(',')
		# Check options
		if options.win_size%options.bin_size != 0:
			sys.exit('Window size should be integer times of bin size')
		if options.step_size%options.bin_size != 0:
			sys.exit('Step size of sliding window should be integer times of bin size')
		win = options.win_size/options.bin_size
		step = options.step_size/options.bin_size

		# Creat temporary directory
		if not os.path.exists(tmp_dir):
			os.makedirs(tmp_dir)
			print '@ ' + time.ctime(),
			print 'Folder "./' + tmp_dir + '" created.'
		
		# Calculating total count of reads
		print '@ ' + time.ctime(),
		print 'Calculating total count of reads in each BED file.'
		dictChip2Len = {}
		for chipName in chipNames:
			chipBed = BedTool(chipName)
			dictChip2Len[chipName] = len(open(chipName).readlines())

		# Train data
		if options.enhancer != 'NA':
			trainData = []
			enhancerBed = loci2bed(options.enhancer, 'enhancer', options.bin_size, options.win_size, tmp_dir)
			promoterBed = loci2bed(options.promoter, 'promoter', options.bin_size, options.win_size, tmp_dir)
			sampleSize = options.back_size
			shuffledBedName = shuffle_window(options.promoter, options.enhancer, sampleSize, options.win_size, options.genome, real_dir, tmp_dir)
			shuffledBed = loci2bed(shuffledBedName, 'shuffled', options.bin_size, options.win_size, tmp_dir)

			for chipName in chipNames:
				print '@ ' + time.ctime(),
				print 'Calculating coverage of ' + chipName + ' at training targets and background region.'
				# Loading ChIP-seq BED file
				chipBed = BedTool(chipName)
				# Line count of ChIP-seq for normalization
				lineCount = dictChip2Len[chipName]
				# Enhancer
				# Count calculation
				BedTool(enhancerBed).window(chipBed,w=0,c=True,output=os.path.join(tmp_dir,'enhancer_count.bed'))
				enhancerProfileMat = profile_target(os.path.join(tmp_dir,'enhancer_count.bed'), os.path.join(tmp_dir,'enhancer_temp_bin.bed'), lineCount, options.bin_size, options.win_size)
				# Calculate target position for each modification and add to position matrix
				#enhancerTargetPosList = [win/2]*len(enhancerProfileMat)
				# Smooth profiles
				#enhancerProfileMat = smooth_mat(enhancerProfileMat)
				# Calculate parameters
				#enhancerFeatureKurt, enhancerFeatureSkew, enhancerFeatureBimo, enhancerFeatureItst = shape(enhancerProfileMat, enhancerTargetPosList, win)

				# Promoter
				# Count calculation
				BedTool(promoterBed).window(chipBed,w=0,c=True,output=os.path.join(tmp_dir,'promoter_count.bed'))
				promoterProfileMat = profile_target(os.path.join(tmp_dir,'promoter_count.bed'), os.path.join(tmp_dir,'promoter_temp_bin.bed'), lineCount, options.bin_size, options.win_size)
				# Calculate target position for each modification and add to position matrix
				#promoterTargetPosList = [win/2]*len(promoterProfileMat)
				# Smooth profiles
				#promoterProfileMat = smooth_mat(promoterProfileMat)
				# Calculate parameters
				#promoterFeatureKurt, promoterFeatureSkew, promoterFeatureBimo, promoterFeatureItst = shape(promoterProfileMat, promoterTargetPosList, win)

				# Background window
				# Coverage calculation
				BedTool(shuffledBed).window(chipBed,w=0,c=True,output=os.path.join(tmp_dir,'shuffled_count.bed'))
				shuffledProfileMat = profile_target(os.path.join(tmp_dir,'shuffled_count.bed'), os.path.join(tmp_dir,'shuffled_temp_bin.bed'), lineCount, options.bin_size, options.win_size)
				# Calculate target position for each modification and add to position matrix
				#shuffledTargetPosList = [win/2]*len(shuffledProfileMat)
				# Smooth profiles
				#shuffledProfileMat = smooth_mat(shuffledProfileMat)
				# Calculate parameters
				#shuffledFeatureKurt, shuffledFeatureSkew, shuffledFeatureBimo, shuffledFeatureItst = shape(shuffledProfileMat, shuffledTargetPosList, win)

				if len(trainData) == 0:
					for i in range(0, len(enhancerProfileMat)):
						trainData.append(np.insert(enhancerProfileMat[i],0,1))
					for i in range(0, len(promoterProfileMat)):
						trainData.append(np.insert(promoterProfileMat[i],0,0))
					for i in range(0, len(shuffledProfileMat)):
						trainData.append(np.insert(shuffledProfileMat[i],0,0))
				else:
					for i in range(0, len(enhancerProfileMat)):
						trainData[i] = np.concatenate((trainData[i],enhancerProfileMat[i]))
					for i in range(0, len(promoterProfileMat)):
						trainData[len(enhancerProfileMat)+i] = np.concatenate((trainData[len(enhancerProfileMat)+i],promoterProfileMat[i]))
					for i in range(0, len(shuffledProfileMat)):
						trainData[len(enhancerProfileMat)+len(promoterProfileMat)+i] = np.concatenate((trainData[len(enhancerProfileMat)+len(promoterProfileMat)+i],shuffledProfileMat[i]))
			write_plain_format(trainData, 'trainData.txt')

if __name__ == '__main__':
	main()
