
import os
import numpy as np
from pybedtools import *


def loci2bed(lociName, targetName, bin_size, win_size, tmp_dir):
	infile = open(lociName, 'r')
	outfileName = os.path.join(tmp_dir,targetName+'_temp_bin.bed')
	outfile = open(outfileName, 'w')
	for ln in infile:
		Lln = ln.strip().split()
		chrom = Lln[0]
		start = int(Lln[1])
		end = int(Lln[2])
		strand = Lln[-1]
		if start > win_size/2:
			if strand == '-':
				for i in range(start+win_size/2, start-win_size/2, -bin_size):
					print >> outfile, chrom + '\t' + str(i-bin_size) + '\t' + str(i)
			else:
				for i in range(start-win_size/2, start+win_size/2, bin_size):
					print >> outfile, chrom + '\t' + str(i) + '\t' + str(i + bin_size)
	infile.close()
	outfile.close()
	binBed = BedTool(outfileName)
	return binBed

def smooth(x):
	window_len = 5
	x = np.array(x)
	s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
	w=np.hanning(window_len)
	y=np.convolve(w/w.sum(),s,mode='valid')
	return y

def smooth_mat(profileMat):
	profileMatNew = []
	for x in profileMat:
		y = smooth(x)
		margin = int((len(y)-len(x))/2)
		profileMatNew.append(y[margin:(len(y)-margin)])
	return profileMatNew

def normalization(chipFpkmMat, inputFpkmMat):
	chipFpkmNorMat = []
	i = 0
	for cc in chipFpkmMat:
		pp = inputFpkmMat[i]
		i += 1
		j = 0
		nn = []
		for c in cc:
			p = pp[j]
			j += 1
			nn.append(c-p)
		chipFpkmNorMat.append(nn)
	return chipFpkmNorMat

def bin_genome(genome, bin_size, win_size, tmp_dir):
	''' Binning genome
	'''
	# Get all chromosome names from BEDtools.
	chromsAll = eval('genome_registry.' + genome + '.keys()')
	# Filter random/unmapped chromosomes
	chromsCore = filter(lambda x: '_' not in x,chromsAll)
	# Get chromosome sizes
	chrom2size = eval('genome_registry.' + genome)
	# Genome bin file name
	binName = os.path.join(tmp_dir, genome + '_bins.bed')
	if not os.path.exists(binName):
		fout = open(binName,'w')
		for chrom in chromsCore:
			i = 0
			while i < (chrom2size[chrom][1]):
				print >> fout, chrom + '\t' + str(i) + '\t' + str(i + bin_size)
				i = i + bin_size
		fout.close()
	else:
		print binName + ' exists, ignore.'
	return binName

def shuffle_window(prom_file, enha_file, sampleSize, winSize, genome, real_dir, tmp_dir):
	outfileName = os.path.join(tmp_dir,'shuffled_window.bed')
	outfile = open(outfileName,'w')
	for i in range(0,sampleSize*2):
		ln = 'chrX\t0\t1\t'
		print >> outfile, ln
	outfile.close()
	aWinBed = BedTool(outfileName)
	shuffledWinBed = aWinBed.shuffle(genome=genome)
	tssBed = BedTool(prom_file)
	shuffledWinTss = shuffledWinBed.window(tssBed, w=winSize)
	enhaBed = BedTool(enha_file)
	shuffledWinEnha = shuffledWinBed.window(enhaBed, w=winSize)
	toFilterList = []
	for feature in shuffledWinTss:
		if feature[0]+'\t'+feature[1]+'\t'+feature[2] not in toFilterList:
			toFilterList.append(feature[0]+'\t'+feature[1]+'\t'+feature[2])
	for feature in shuffledWinEnha:
		if feature[0]+'\t'+feature[1]+'\t'+feature[2] not in toFilterList:
			toFilterList.append(feature[0]+'\t'+feature[1]+'\t'+feature[2])	
	outfileName = os.path.join(tmp_dir,'shuffled_window_filtered.bed')
	outfile = open(outfileName,'w')
	i = 0
	for feature in shuffledWinBed:
		if feature[0]+'\t'+feature[1]+'\t'+feature[2] not in toFilterList:
			print >> outfile, feature[0]+'\t'+feature[1]+'\t'+feature[2]
			i += 1
			if i >= sampleSize:
				break
	outfile.close()
	return outfileName

def profile_target(countBed, binBed, totReadCount, binSize, winSize):
	infile = open(countBed,'r')
	dictPosfpkm = {}
	for ln in infile:
		Lln = ln.strip().split()
		fpkm = float(Lln[3])*1000000000/(totReadCount*binSize)
		dictPosfpkm[(Lln[0],Lln[1],Lln[2])] = fpkm
	infile.close()
	windowList = []
	profileMat = []
	infile = open(binBed, 'r')
	i = 1
	for ln in infile:
		Lln = ln.strip().split()
		if i % (winSize/binSize) == 0:
			windowList.append(dictPosfpkm[(Lln[0],Lln[1],Lln[2])])
			profileMat.append(windowList)
			windowList = []
		else:
			windowList.append(dictPosfpkm[(Lln[0],Lln[1],Lln[2])])
		i += 1
	return profileMat

def profile_sliding_window(countBedName, totReadCount, win, step, binSize, tmp_dir):
	# Start/end positions of sliding windows
	plist = []
	# Entire profile matrix
	wlist = []
	# Open BED file
	countBed = open(countBedName, 'r')
	lines = countBed.readlines()
	def fpkm(x):
		return float(x.split()[3])*1000000000/(totReadCount*binSize)
	for i in range(0, len(lines), step):
		try:
			w = map(fpkm, lines[i:(i+win)])
			p = lines[i].split()[0:2] + lines[i+win-1].split()[2:3]
			wlist.append(w)
			plist.append(p)
		except IndexError:
			pass
	countBed.close()
	return wlist, plist

def target_position(profileMat, win_size):
	targetPosList = []
	for profile in profileMat:
		targetPosList.append(np.array(profile).argmax())
	return targetPosList

def write_plain_format(trainData, outfileName):
	outfile = open(outfileName, 'w')
	for i in range(0, len(trainData)):
		for j in range(0, len(trainData[i])-1):
			print >> outfile, str(trainData[i][j]) + '\t',
		print >> outfile, str(trainData[i][len(trainData[i])-1])
	outfile.close()
