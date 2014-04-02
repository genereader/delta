#!/usr/bin/env python

import numpy as np

def shape(profileMat, targetPos, win):
	# initialize parameters
	kurts = []
	skews = []
	bimos = []
	maxes = []
	margin = 10
	halfwin = 21
	nonzero = 5
	i = 0

	for x in profileMat:
		pos = targetPos[i]
		i += 1
		h = max(x)

		if pos < margin or pos > win-margin-1 or sum(x[pos-margin:pos]!=0) <= nonzero or sum(x[pos+1:pos+margin+1]!=0) <= nonzero:
			# If 10 or less bins locating at either side of the target or 5 or less nonzero bins, set all shape parameters to be 0
			kurt = 0
			skew = 0
			bimo = 0
		else:
			x = x[pos-margin:pos+margin+1]
			# Calculation of Kurtosis
			kurt = sum((np.array(range(0,halfwin))-margin)**4*x)*sum(x) / sum((np.array(range(0,halfwin))-margin)**2*x)**2
			# Calculation of Skewness
			skew = np.sqrt(sum(x)) * sum((np.array(range(0,halfwin))-margin)**3*x) / sum((np.array(range(0,halfwin))-margin)**2*x)**1.5
			# Calculation of Bimodality
			bimo = (skew**2 + 1) / kurt
		
		maxes.append(h)
		kurts.append(kurt)
		skews.append(skew)
		bimos.append(bimo)

	return kurts, skews, bimos, maxes
