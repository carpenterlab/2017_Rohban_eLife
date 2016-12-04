#!/usr/bin/env python
import os
from cpa.profiling.cache import *
import sys
from optparse import OptionParser

show_progress = True
parser = OptionParser("usage: %prog <dataset-name> ")
options, args = parser.parse_args()

if len(args) != 1:
    parser.error('Incorrect number of arguments')

dataset_name = sys.argv[1]

if dataset_name == 'TA-MHP-B1':
    norm_methods = dict(
        {
            #            'ctrl_norm'  : "Image_Metadata_ASSAY_WELL_ROLE = 'CTRL'",
            'untreated_norm'  : "Image_Metadata_ASSAY_WELL_ROLE = 'Untreated'"
        })
elif dataset_name == 'TA-OE-B1':
    norm_methods = dict(
        {
            #        'lacz_norm'  : "Image_Metadata_GeneSymbol = 'LacZ'",
            #        'all_ctrl_norm'  : "Image_Metadata_ASSAY_WELL_ROLE = 'CTRL'"
            'untreated_norm'  : "Image_Metadata_ASSAY_WELL_ROLE = 'Untreated'"
        })
else:
    # Set a default
    print "Unknown dataset %s. Setting default values for norm_methods"%dataset_name
    norm_methods = dict({'scrambled_norm'  : "Image_Metadata_DetectRgt = 'Scrambled Ctrl'"})

datadir = '../input/' + dataset_name + '/'
import glob
try:
    properties_file = glob.glob(os.path.join(datadir, "*.properties"))[0]
except Exception, e: 
    print "Properties file not found. Exiting."
    print  e
    sys.exit(1)

################################################################

rename_normdir = True

cpa.properties.LoadFile(properties_file)
cache = Cache(os.path.join(datadir, 'cache'))

from cpa.profiling.normalization import *
norm_d = dict({#'std': StdNormalization
               'robust_std': RobustStdNormalization,
               'robust_linear': RobustLinearNormalization
               # 'dummy': DummyNormalization
           })
#norm_d = dict({'robust_linear': RobustLinearNormalization})

for norm_method in norm_methods.keys():
    for k, f in norm_d.iteritems():
        print "Computing normalization using %s:%s" % (norm_method, k)

        r = f(cache)
        r._create_cache(norm_methods[norm_method])
        
        if rename_normdir:
            try:
                print 'Renaming normalization dir...'
                os.rename(r.dir, r.dir + '_' +  norm_method)        
            except OSError, e:
                print e.args

