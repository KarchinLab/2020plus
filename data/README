Gene lists
==========

The `gene_lists` folder contains text files 'oncogenes.txt' and 'tsgs.txt', which contain the training set
of previously labeled oncogenes and tumor suppressor genes, respectively.

Covariates
==========

'mutsigcv_gene_features.txt' contains features from MutsigCV about replication timing, expression from CCLE, and HiC.

Interaction network
===================

'biogrid_stats.txt' contains network features from BioGrid about the gene degree (number of connected genes) and gene betweeness centrality.

Mutation processing
===================

The following description is relevant to the processing of mutations found in:

Tokheim C, Papdopoulis N, Kinzler KW, Vogelstein B, Karchin R (2016) Evaluating the evaluation of cancer driver genes. Submitted `[bioRxiv preprint] <http://biorxiv.org/content/early/2016/06/23/060426>`_

Mutation cleaning
-----------------
Out mutation data set contains small somatic mutations in coding regions of genes across thousands of cancer samples. We observed that the two previous studies where we originally obtained our mutation data (TUSON (http://elledgelab.med.harvard.edu/wpcontent/uploads/2013/11/Mutation_Dataset.txt.zip) and Mutsig (http://www.tumorportal.org/load/data/per_ttype_mafs/PanCan.maf)) did not always filter out all duplicate samples, non-cancer samples, or non-systematic studies (i.e. not whole-exome nor whole-genome). Our study focused on whole-exome or -genome sequencing of cancer, and, as such, samples which could not be attributed to cancer or did not comprehensively report all mutations on the exome were excluded. We manually examined cancer samples that did not report any silent mutations. We backtracked the original study reporting the samples, and examined whether the study only reported non-silent mutations. If true, we removed those samples from our data set. In a similar manner, we manually excluded samples that were discovery-prevalence screens, which do not report the whole exome. Also, we excluded samples which were not firmly cancer. This included samples with "NS" ('not specified') tumor types, and pre-cancerous cysts (Pancreatic IPMNs, MCN, SCA, or SPN). Duplicate samples were detected by either having the same sample ID (done automatically) or manually excluding samples identified by having a high proportion of exactly the same mutations among two samples.

Note on sample IDs
------------------
Many sample IDs come from the original study and thus are chosen by the author. Unlike the standardized The Cancer Genome Atlas (TCGA) sample IDs, they may be arbitrary combinations of letters and numbers. Often they are either just a number, or a biolgical abbreviation plus a number (e.g. 587226 or IPMN11, respectively). To fully understand the sample ID naming convention, one would need to go back to the original study for each ID, although this was outside the scope of our filtering step. 

Resulting filter
----------------
Samples manually filtered by the above process are shown in "banned.samples.txt". The final mutation file excludes those samples in the "banned.samples.txt", hypermutators (>1000 mutations), mutations with mappability warning, and duplicate samples based on the exact same sample ID. Therefore filtered samples will not appear in the "Tumor_Sample" column in the final cleaned data set (available at http://karchinlab.org/data/Protocol/pancan-mutation-set-from-Tokheim-2016.txt.gz).
