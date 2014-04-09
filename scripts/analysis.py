#! /usr/bin/python

import pandas as pd
import numpy as np

in_dag = '/home/mbarga/Dropbox/#workbench/git/bayesnet/data/1000dag'
serial = '/home/mbarga/Dropbox/#workbench/git/bayesnet/data/output.txt'
parallel = '/home/mbarga/Dropbox/#workbench/git/bayesnet/data/output.txt'

#sensitivity = TP/P = TP/(TP+FN)
def sensitivity(tp, fn):
	denom = (tp + fn)
	if denom > 0:
		print "returning ", tp,
		print " / ", tp,
		print " + ", fn
		return float(tp) / denom
	else:
		return -1

#specificity = TN/N = TN/(FP+TN)
def specificity(tn, fp):
	denom = (fp + tn)
	if denom > 0:
		print "returning ", tn,
		print " / ", fp,
		print " + ", tn
		return float(tn) / denom
	else:
		return -1

# IMPORT ORIGINAL DAG
def get_true_network(output = False):
	f = open(in_dag, 'r')
	header = f.readline().split()
	p = int(header[2])
	gold = [[False for x in xrange(p)] for x in xrange(p)]

	child = 0
	for line in f:
		cleaned = line.split()
		size = len(cleaned)
		if size > 1:
			for col in xrange(size-1):
				parent = int(cleaned[col+1])
				gold[child][parent] = True
		child = child + 1
	f.close()
	gold = np.array(gold)

	print "Edges in true network: ", gold[np.where(gold > 0)].size
	if output is True:
		print gold
	return gold, p

# IMPORT OUTPUT G ADJACENCY MATRIX (serial)
def get_serial_network(output = False):
	G = np.loadtxt(serial, delimiter=',', dtype=np.bool)

	print "Edges in serial netowrk: ", G[np.where(G > 0)].size
	if output is True:
		print G
	return G

def main():
	true_network, p = get_true_network(output=False) # boolean
	serial_network = get_serial_network(output=False) # boolean
	#TODO assert that the sizes are equal?

	tp = 0
	tn = 0
	fp = 0
	fn = 0

	for i in xrange(p):
		for j in xrange(p):
			if true_network[i][j] == 1:
				if serial_network[i][j] == 1:
					tp = tp + 1
				else:
					fn = fn + 1
			else:
				if serial_network[i][j] == 0:
					tn = tn + 1
				else:
					fp = fp + 1

	print tp, fn, tn, fp
	sp = specificity(tn, fp)
	print "sp = ", sp
	sn = sensitivity(tp, fn)
	print "sn = ", sn

if __name__ == "__main__":
    main()
