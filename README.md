I. TCGA_HNSC input file
tcga_hnsc_tpm_os.matrix.txt
##Original RSEM data was downloaded at:
http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/HNSC/20160128/gdac.broadinstitute.org_HNSC.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0.tar.gz


II. GSE65858 input file
GSE65858_exp.matrix.txt
##original GSE65858 input file was downloaded at:
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65858/matrix/GSE65858_series_matrix.txt.gz


III. to generate 15-gene signature score from TCGA-HNSC dataset:

tcga_hnsc_tpm_random_survival_forest.r

IV. to evaluate 15-gene signature score in GSE65858 dataset:

gse65858_microarray_random_survival_forest_15gene_risk.R

V. to generate NMF clustering result from TCGA-HNSC dataset:

tcga_hnsc_top_k_for_loops_nmf_nrun100.r

VI. to evaluate NMF clustering result in GSE65858 dataset:

limma_analysis_gse65858_candidate_gene_from_tcga_hnsc_top500_mad.R

VII. to generate ICI clustering result from TCGA-HNSC dataset:

TCGA_HNSC_cibersort.r

VIII. to generate ICI clustering result from GSE65858 dataset:

GSE65858_ciber_infiltration.r
