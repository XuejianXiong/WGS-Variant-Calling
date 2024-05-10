# WGS-Variant-Calling

This repotary is about WGS variant calling and annotation using GATK.


1) Download pre-packaged data source from Funcotator

	cd tools/

	gatk FuncotatorDataSourceDownloader --germline --validate-integrity --extract-after-download --hg38


2) Run scripts

	cd scripts/

	chmod +x variant_calling.sh ./variant_calling.sh

	chmod +x variant_filtering_annotation.sh ./variant_filtering_annotation.sh


3) Useful variant Data Sources:

	GNOMAD

	dbSNP

	ClinVar

	COSMIC


