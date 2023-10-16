

wget --directory-prefix=/data6/deepro/ukb_bmi/0_data_preparation_and_download/enrichment/data/go http://release.geneontology.org/2023-07-27/ontology/subsets/goslim_agr.obo # http://release.geneontology.org/2021-02-01/ontology/go-basic.obo

wget --directory-prefix=/data6/deepro/ukb_bmi/0_data_preparation_and_download/enrichment/data/go http://release.geneontology.org/2023-07-27/annotations/goa_human.gaf.gz # http://release.geneontology.org/2021-02-01/annotations/goa_human.gaf.gz

gunzip /data6/deepro/ukb_bmi/0_data_preparation_and_download/enrichment/data/go/goa_human.gaf.gz

# wget  --directory-prefix=/data6/deepro/ukb_bmi/0_data_preparation_and_download/enrichment/data/hpo https://raw.githubusercontent.com/tanghaibao/goatools/main/notebooks/data/hpo/hp.obo

# wget  --directory-prefix=/data6/deepro/ukb_bmi/0_data_preparation_and_download/enrichment/data/hpo https://raw.githubusercontent.com/tanghaibao/goatools/main/notebooks/data/hpo/hpo.annotation.tab