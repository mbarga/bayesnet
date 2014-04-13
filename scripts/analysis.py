#! /usr/bin/python

import numpy as np

truenet =  '1000dag'
serial =   '1000serial'
parallel = '1000omp'

directory = '/home/mbarga/Dropbox/#workbench/git/bayesnet/data/'
in_dag = directory + truenet
serial = directory + serial
parallel = directory + parallel

# sensitivity = TP/P = TP/(TP+FN)
def sensitivity(tp, fn):
    denom = (tp + fn)
    if denom > 0:
        #print "returning ", tp,
        #print " / ", tp,
        #print " + ", fn
        return float(tp) / denom
    else:
        return -1

# specificity = TN/N = TN/(FP+TN)
def specificity(tn, fp):
    denom = (fp + tn)
    if denom > 0:
        #print "returning ", tn,
        #print " / ", fp,
        #print " + ", tn
        return float(tn) / denom
    else:
        return -1

# get values needed to compute specificity and sensitivity
def contingency(p, true_network, test_network):
    tp = 0
    tn = 0
    fp = 0
    fn = 0

    for i in xrange(p):
        for j in xrange(p):
            if true_network[i][j] == 1:
                if test_network[i][j] == 1:
                    tp = tp + 1
                else:
                    fn = fn + 1
            else:
                if test_network[i][j] == 0:
                    tn = tn + 1
                else:
                    fp = fp + 1

    return tp, tn, fp, fn

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
def get_est_network(f, output = False):
    G = np.loadtxt(f, delimiter=',', dtype=np.bool)

    print "Edges in estimated netowrk: ", G[np.where(G > 0)].size
    if output is True:
        print G
    return G

def main():
    true_network, p = get_true_network(output=False) # boolean
    serial_network = get_est_network(serial, output=False) # boolean
    parallel_network = get_est_network(parallel, output=False) # boolean
    #TODO assert that the sizes are equal?

    print
    print "SERIAL NETWORK"
    stp, stn, sfp, sfn = contingency(p, true_network, serial_network)
    sp = specificity(stn, sfp)
    sn = sensitivity(stp, sfn)
    print
    print "TP          FN          TN          FP          SP          SN"
    print '%12s%12s%12s%12s%12s%12s' % (str(stp).ljust(12), str(sfn).ljust(12), str(stn).ljust(12), str(sfp).ljust(12), str.format('{0:.4f}', sp).ljust(12), str.format('{0:.4f}', sn).ljust(12))
    print

    print "PARALLEL NETWORK"
    ptp, ptn, pfp, pfn = contingency(p, true_network, parallel_network)
    sp = specificity(ptn, pfp)
    sn = sensitivity(ptp, pfn)
    print "TP          FN          TN          FP          SP          SN"
    print '%12s%12s%12s%12s%12s%12s' % (str(ptp).ljust(12), str(pfn).ljust(12), str(ptn).ljust(12), str(pfp).ljust(12), str.format('{0:.4f}', sp).ljust(12), str.format('{0:.4f}', sn).ljust(12))
    print

if __name__ == "__main__":
    main()
