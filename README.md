# PAM-clustering-automation-you-just-built
CLI-driven, Linux-friendly pipeline for PAM clustering on microbiome count tables. It imports TSV/CSV (biom header handled), optionally applies TSS scaling, computes JSD distances, and picks optimal k via CH index (silhouette tie-break). Runs PAM, plots PCoA, outputs Top-10 barplots (per-cluster + combined), and saves PNGs, TSVs, Excel reports
