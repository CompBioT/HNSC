1. TCGA-HNSC Original RSEM data was downloaded at:
http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/HNSC/20160128/gdac.broadinstitute.org_HNSC.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0.tar.gz

2. GSE65858 original file was downloaded at:
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65858/matrix/GSE65858_series_matrix.txt.gz

3. our processed TCGA-HNSC and GSE65858 were in "input" folder.

4. Script to generate 15-gene signature score from TCGA-HNSC dataset:

tcga_hnsc_tpm_random_survival_forest.r

5. to evaluate 15-gene signature score in GSE65858 dataset:

gse65858_microarray_random_survival_forest_15gene_risk.R

6. to generate NMF clustering result from TCGA-HNSC dataset:

tcga_hnsc_top_k_for_loops_nmf_nrun100.r

7. to evaluate NMF clustering result in GSE65858 dataset:

limma_analysis_gse65858_candidate_gene_from_tcga_hnsc_top500_mad.R

8. to generate ICI clustering result from TCGA-HNSC dataset:

TCGA_HNSC_cibersort.r

9. to generate ICI clustering result from GSE65858 dataset:

GSE65858_ciber_infiltration.r
