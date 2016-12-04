#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Specify a dataset."
    exit
fi

lsfdir=/broad/hptmp/shsingh/.compute_`date +%s%N` \
dataset_name=$1

datadir=../input/${dataset_name}/
cache_dir=cache
properties_file=`find $datadir -name "*.properties"`
properties_file=`basename $properties_file`

curdir=`pwd`
OS=`uname`
if [ $OS == "Darwin" ]
then
    PYTHONBIN=/Users/shsingh/work/software/CPhomebrew/Cellar/cellprofiler-dev/1/cpdev/bin/python
    SED=gsed
    parallel_option=multiprocessing
else
    PYTHONBIN=python
    SED=sed
    #parallel_option=lsf-directory=$lsfdir
    parallel_option=multiprocessing
fi

# ---------- Do normalization  ----------

function run_profile_gen {
  echo 
  echo 
  set -x
  echo \#Computing profiles ${normalization_tag}:${normtype}:${method}
  csvname=well-summary-${profile_method}-${method}-${normalization_tag}-${normtype}.csv
  rm -rf $lsfdir
  cd ${datadir}${cache_dir}
  rm -rf ${normalization_tag}
  ln -s ${normalization_tag}_${normtype} ${normalization_tag}
  cd -
  cd $datadir
#  $PYTHONBIN -m cpa.profiling.${profile_method} -c -g -o $csvname  --${parallel_option} --njobs=300 --normalization=$normalization --method=${method} $properties_file $cache_dir Well
  $PYTHONBIN -m cpa.profiling.${profile_method} -c -g -o $csvname  --${parallel_option} --normalization=$normalization --method=${method} $properties_file $cache_dir Well
  rm -rf $lsfdir
  cd -
  echo 
}

# profile_method=profile_mean
# normalization=RobustLinearNormalization
# normalization_tag=robust_linear
# normtype=ctrl_norm
# method=mean
# run_profile_gen

# profile_method=profile_mean
# normalization=RobustLinearNormalization
# normalization_tag=robust_linear
# normtype=ctrl_norm
# method=median
# run_profile_gen

# profile_method=profile_mean
# normalization=DummyNormalization
# normalization_tag=dummy
# normtype=dummy_norm
# method=median
# run_profile_gen

# profile_method=profile_mean
# normalization=RobustLinearNormalization
# normalization_tag=robust_linear
# normtype=lacz_norm
# method=median
# run_profile_gen

profile_method=profile_mean
normalization=RobustStdNormalization
normalization_tag=robust_std
normtype=untreated_norm
method=median
run_profile_gen

profile_method=profile_mean
normalization=RobustLinearNormalization
normalization_tag=robust_linear
normtype=untreated_norm
method=median
run_profile_gen

# profile_method=profile_mean
# normalization=RobustLinearNormalization
# normalization_tag=robust_linear
# normtype=all_ctrl_norm
# method=median
# run_profile_gen

# profile_method=profile_mean
# normalization=StdNormalization
# normalization_tag=std
# normtype=all_ctrl_norm
# method=mean
# run_profile_gen





