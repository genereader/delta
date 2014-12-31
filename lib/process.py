
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
		start = start - start%bin_size
		strand = Lln[-1]
		if start > win_size/2:
			if strand == '-':
				for i in range(start+win_size/2+bin_size, start-win_size/2, -bin_size):
					print >> outfile, chrom + '\t' + str(i-bin_size) + '\t' + str(i)
			else:
				for i in range(start-win_size/2, start+win_size/2+bin_size, bin_size):
					print >> outfile, chrom + '\t' + str(i) + '\t' + str(i + bin_size)
	infile.close()
	outfile.close()
	binBed = BedTool(outfileName)
	return binBed

def smooth(x):
	window_len = 5
	#x = np.array(x)
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
		if i % (winSize/binSize+1) == 0:
			windowList.append(dictPosfpkm[(Lln[0],Lln[1],Lln[2])])
			profileMat.append(windowList)
			windowList = []
		else:
			windowList.append(dictPosfpkm[(Lln[0],Lln[1],Lln[2])])
		i += 1
	return profileMat

def shape(x, win):
	pos = win/2
	if len(x) == win + 1:
		x = np.asarray(x)
		y = smooth(x)
		margin = int((len(y)-len(x))/2)
		x = y[margin:(len(y)-margin)]
		h = x[pos]
		if x.mean() < 5 or sum(x[0:pos]!=0) <= 5 or sum(x[(pos+1):]!=0) <= 5:
			kurt = 0
			skew = 0
			bimo = 0
		else:
			# Kurtosis
			kurt = sum((np.array(range(0,len(x)))-len(x)/2)**4*x)*sum(x) / sum((np.array(range(0,len(x)))-len(x)/2)**2*x)**2
			# Skewness
			skew = np.sqrt(sum(x)) * sum((np.array(range(0,len(x)))-len(x)/2)**3*x) / sum((np.array(range(0,len(x)))-len(x)/2)**2*x)**1.5
			# Bimodality
			bimo = (skew**2 + 1) / kurt
	else:
		kurt = 0
		skew = 0
		bimo = 0
		h = 0
	return kurt, skew, bimo, h

def shape_mat(profileMat, win):
	# initialize parameters
	kurts = []
	skews = []
	bimos = []
	maxes = []
	pos = win/2
	for x in profileMat:
		x = np.asarray(x)
		y = smooth(x)
		margin = int((len(y)-len(x))/2)
		x = y[margin:(len(y)-margin)]
		h = x[pos]
		if x.mean() < 5 or sum(x[0:pos]!=0) <= 5 or sum(x[(pos+1):]!=0) <= 5:
			kurt = 0
			skew = 0
			bimo = 0
		else:
			# Kurtosis
			kurt = sum((np.array(range(0,len(x)))-len(x)/2)**4*x)*sum(x) / sum((np.array(range(0,len(x)))-len(x)/2)**2*x)**2
			# Skewness
			skew = np.sqrt(sum(x)) * sum((np.array(range(0,len(x)))-len(x)/2)**3*x) / sum((np.array(range(0,len(x)))-len(x)/2)**2*x)**1.5
			# Bimodality
			bimo = (skew**2 + 1) / kurt	
		maxes.append(h)
		kurts.append(kurt)
		skews.append(skew)
		bimos.append(bimo)
	return kurts, skews, bimos, maxes

def profile_sliding_window(countBedName, totReadCount, win, binSize, featureFileDir):
	# Open feature file
	featureFile = open(featureFileDir, 'w')
	# Open BED file
	fpkmVector = []
	countBed = open(countBedName, 'r')
	for ln in countBed:
		fpkm = float(ln.split()[3])*1000000000/(totReadCount*binSize)
		fpkmVector.append(fpkm)
	countBed.close()
	fpkmVector = np.asarray(fpkmVector)
	for i in range(0, len(fpkmVector)):
		kurt, skew, bimo, h = shape(fpkmVector[(i-win/2):(i+win/2+1)],win)
		featureFile.write(str(h)+'\t'+str(kurt)+'\t'+str(skew)+'\t'+str(bimo)+'\n')
	featureFile.close()

def profile_sliding_window_binvector(countBedName, totReadCount, win, binSize, featureFileDir):
	# Open feature file
	featureFile = open(featureFileDir, 'w')
	# Open BED file
	fpkmVector = []
	countBed = open(countBedName, 'r')
	for ln in countBed:
		fpkm = float(ln.split()[3])*1000000000/(totReadCount*binSize)
		fpkmVector.append(fpkm)
	countBed.close()
	for i in range(0, len(fpkmVector)):
		if len(fpkmVector[(i-win/2):(i+win/2+1)]) == win + 1:
			for fpkm in fpkmVector[(i-win/2):(i+win/2+1)]:
				featureFile.write(str(fpkm) + '\t')
		else:
			for fpkm in [0]*21:
				featureFile.write(str(fpkm) + '\t')
		featureFile.write('\n')
	featureFile.close()

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
