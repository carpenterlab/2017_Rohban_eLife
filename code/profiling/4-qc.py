#!/usr/bin/env python
import os
import cpa
import sys
from optparse import OptionParser

show_progress = True
parser = OptionParser("usage: %prog <dataset-name> ")
options, args = parser.parse_args()

if len(args) != 1:
    parser.error('Incorrect number of arguments')

dataset_name = sys.argv[1]


datadir = '../input/' + dataset_name + '/'
import glob
try:
    properties_file = glob.glob(os.path.join(datadir, "*.properties"))[0]
except Exception, e: 
    print "Properties file not found. Exiting."
    print  e
    sys.exit(1)

################################################################

cpa.properties.LoadFile(properties_file)


features = cpa.db.execute("""select Image_Metadata_Plate, Image_Metadata_Well, sum(Image_Metadata_QCFlag_isSaturated) as QCFlag_isSaturated, sum(Image_Metadata_QCFlag_isBlurry) as QCFlag_isBlurry from %s group by  Image_Metadata_Plate, Image_Metadata_Well""" % cpa.properties.image_table)

filename = datadir + 'qc.csv'

with open(filename, 'w') as f:
  print>>f, "Plate,Well,QCFlag_isSaturated,QCFlag_isBlurry"
  for i in features:
    print>>f, ",".join(str(j) for j in i ) 
  
  