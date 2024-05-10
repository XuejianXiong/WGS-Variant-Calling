# WGS-Variant-Calling

This repotary is about WGS variant calling and annotation using GATK.


Download pre-packaged data source from Funcotator

gatk FuncotatorDataSourceDownloader --germline --validate-integrity \
						--extract-after-download --hg38


cd scripts/

chmod +x variant_calling.sh
./variant_calling.sh

chmod +x variant_filtering_annotation.sh
./variant_filtering_annotation.sh


Useful variant Data Sources:
1) GNOMAD
2) dbSNP
3) ClinVar
4) COSMIC


