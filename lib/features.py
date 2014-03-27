#!/usr/bin/env python

import numpy as np

def shape(profileMat, targetPos, win):
	kurts = []
	skews = []
	bimos = []
	maxes = []
	def fwhm(x): return x > h/2
	i = 0
	for x in profileMat:
		pos = targetPos[i]
		i += 1
		h = max(x)
		maxes.append(h)
		if pos < 10 or pos > win-11 or sum(x[pos-10:pos]!=0) <= 5 or sum(x[pos+1:pos+11]!=0) <= 5:
			kurts.append(0)
			skews.append(0)
			bimos.append(0)
		else:
			x = x[pos-10:pos+11]
			kurt = sum((np.array(range(0,21))-10)**4*x)*sum(x) / sum((np.array(range(0,21))-10)**2*x)**2
			skew = np.sqrt(sum(x)) * sum((np.array(range(0,21))-10)**3*x) / sum((np.array(range(0,21))-10)**2*x)**1.5
			bimo = (skew**2 + 1) / kurt
			kurts.append(kurt)
			skews.append(skew)
			bimos.append(bimo)
	return kurts, skews, bimos, maxes
