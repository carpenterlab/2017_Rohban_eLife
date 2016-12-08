python -m cpa.profiling.cache -r ../../input/TA_OE_B1/TargetAccelerator.properties ../../input/TA_OE_B1/cache ""

python -m cpa.profiling.normalization --method=RobustStdNormalization \
	../../input/TA_OE_B1/TargetAccelerator.properties \
	../../input/TA_OE_B1/cache \
	"Image_Metadata_ASSAY_WELL_ROLE = 'Untreated'"

python -m cpa.profiling.profile_mean --method=median+mad \
 	--normalization RobustStdNormalization \
        -o ../../input/profiles/median.plus.mad.robstd.empty.well.profiles.csv \
        -c -g \
        ../../input/TA_OE_B1/TargetAccelerator.properties ../../input/TA_OE_B1/cache Well

python -m cpa.profiling.profile_mean --method=cellcount \
	--normalization DummyNormalization \
	-o ../../input/profiles/counts.well.csv \
	-c -g \
	../../input/TA_OE_B1/TargetAccelerator.properties ../../input/TA_OE_B1/cache Well
