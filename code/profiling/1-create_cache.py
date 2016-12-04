#!/usr/bin/env python
import os
import sys
from optparse import OptionParser
import fnmatch


show_progress = True
parser = OptionParser("usage: %prog <dataset-name> ")
options, args = parser.parse_args()

if len(args) != 1:
    parser.error('Incorrect number of arguments')

dataset_name = sys.argv[1]

ddir="../input/"
flag=False
for i, fname in enumerate(os.listdir(os.path.join(ddir, dataset_name))):
    if fnmatch.fnmatch(fname, '*.properties'):
        prop_file = os.path.join(ddir, dataset_name, fname)
        cache_dir = os.path.join(ddir, dataset_name, "cache")
        cmd = "python -m cpa.profiling.cache -r %s %s %s" % \
              (prop_file, cache_dir, '""')
        bsub="bsub -q week -J %s -o %s.out " % \
              ("cache-" + dataset_name, "cache-" + dataset_name)

        print bsub + cmd
        flag=True


if not flag:
    parser.error('Dataset %s does not exist' % dataset_name)
