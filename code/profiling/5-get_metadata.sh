#!/bin/bash
dataset_name="TA-POP-B1" 
#geneset_name="all_genes"
geneset_name="pop_genes"

#------------------------------------------------------------
ctrls_name="all_controls"
db_name="2011_07_13_TargetAccelerator"

datadir=../data/${dataset_name}/
csvfile=${datadir}cell_cnt-${ctrls_name}-${geneset_name}.csv


if [ $dataset_name = "TA-POP-B1" ]
then
  DBNAME='U2OS_reimage_Merged_2012_08_15_Analysis_Per_Image'
  CMD1='SELECT Image_Metadata_Well as Well,Image_Metadata_Plate as Plate, Image_Metadata_ASSAY_WELL_ROLE as Role, Image_Metadata_CompoundName_GeneSymbol as Gene, Image_Metadata_Source as RNAi, sum(Image_Count_Cells) as CellCnt, Image_Metadata_Performance_U2OS_AffyExp_ge_6 as Affy,Image_Metadata_Performance_ATARiS_hp_cscore as ataris,Image_Metadata_Performance_Conc_uM_IE_ge_07 as conc_ie, Image_Metadata_Performance_shRNA_CoreGeneSignature as cgs_shrna,Image_Metadata_Performance_CoreGeneSignature_A549 as cgs_a549,Image_Metadata_Performance_CoreGeneSignature_PC3  as cgs_pc3,Image_Metadata_Performance_CoreGeneSignature_VCAP  as cgs_vcap FROM '

  if [ $geneset_name = "all_genes" ]
  then
#    CMD2=' WHERE Image_Metadata_Type="shRNA" AND Image_Metadata_CompoundName_GeneSymbol<>"" AND Image_Metadata_TimePoint_Hours=144 group by Image_Metadata_Well, Image_Metadata_Plate order by sum(Image_Count_Cells)'
    CMD2=' WHERE Image_Metadata_Type="shRNA" AND Image_Metadata_CompoundName_GeneSymbol<>"" AND Image_Metadata_TimePoint_Hours=144 group by Image_Metadata_Well, Image_Metadata_Plate'
  elif [ $geneset_name = "pop_genes" ]
  then
	CMD2=' WHERE Image_Metadata_Type="shRNA" AND Image_Metadata_ProofOfPrincipleAnalysis="POPanalysis" AND Image_Metadata_CompoundName_GeneSymbol<>"" AND Image_Metadata_TimePoint_Hours=144 group by Image_Metadata_Well, Image_Metadata_Plate'
	
  else 
	  echo "No such geneset_name:" ${geneset_name}
	  exit
  fi
fi

CMD=${CMD1}${DBNAME}${CMD2}

sqlfile=.query_`date +%s%N`.sql
echo $CMD > $sqlfile
mysql -A -u cpuser -pcPus3r -h imgdb02 ${db_name} < $sqlfile | tr '\t' ',' > $csvfile
#mysql -A -u cpuser -pcPus3r -h imgdb02 ${db_name} < $sqlfile 
#rm $sqlfile
