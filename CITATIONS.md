# vaximile: Citations

## [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

> Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820.

## [nf-core](https://pubmed.ncbi.nlm.nih.gov/32055031/)

> Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x.

## Pipeline tools

### Quality control and reporting

- [fastp](https://pubmed.ncbi.nlm.nih.gov/30423086/) — Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018;34(17):i884-i890.
- [MultiQC](https://pubmed.ncbi.nlm.nih.gov/27312411/) — Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016;32(19):3047-8.
- [SAMtools](https://pubmed.ncbi.nlm.nih.gov/19505943/) — Li H, Handsaker B, Wysoker A, et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009;25(16):2078-9.
- [somalier](https://pubmed.ncbi.nlm.nih.gov/32664994/) — Pedersen BS, Bhetariya PJ, Brown J, et al. Somalier: rapid relatedness estimation for cancer and germline studies using efficient genotype sketches. Genome Med. 2020;12(1):62.

### Alignment and pre-processing

- [Minibwa](https://arxiv.org/abs/2606.15357) — Li H, Homer N. Fast genomic read alignment with minibwa. arXiv:2606.15357. 2026.
- [GATK](https://pubmed.ncbi.nlm.nih.gov/20644199/) — McKenna A, Hanna M, Banks E, et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010;20(9):1297-303.

### Somatic and germline variant calling

- [Mutect2](https://www.biorxiv.org/content/10.1101/861054v1) — Benjamin D, Sato T, Cibulskis K, et al. Calling Somatic SNVs and Indels with Mutect2. bioRxiv. 2019.
- [Strelka2](https://pubmed.ncbi.nlm.nih.gov/30013048/) — Kim S, Scheffler K, Halpern AL, et al. Strelka2: fast and accurate calling of germline and somatic variants. Nat Methods. 2018;15(8):591-594.
- [Manta](https://pubmed.ncbi.nlm.nih.gov/26647377/) — Chen X, Schulz-Trieglaff O, Shaw R, et al. Manta: rapid detection of structural variants and indels for germline and cancer sequencing applications. Bioinformatics. 2016;32(8):1220-2.
- [DeepVariant](https://pubmed.ncbi.nlm.nih.gov/30247488/) — Poplin R, Chang PC, Alexander D, et al. A universal SNP and small-indel variant caller using deep neural networks. Nat Biotechnol. 2018;36(10):983-987.
- [DeepSomatic](https://www.biorxiv.org/content/10.1101/2024.08.16.608331v1) — Park J, Cook DE, Chang PC, et al. DeepSomatic: Accurate somatic small variant discovery for multiple sequencing technologies. bioRxiv. 2024.

### Annotation and copy number

- [Ensembl VEP](https://pubmed.ncbi.nlm.nih.gov/27268795/) — McLaren W, Gil L, Hunt SE, et al. The Ensembl Variant Effect Predictor. Genome Biol. 2016;17(1):122.
- bam-readcount — Khanna A, Larson DE, Srivatsan SN, et al. Bam-readcount - rapid generation of basepair-resolution sequence metrics. J Open Source Softw. 2022;7(69):3722.
- [ASCAT](https://pubmed.ncbi.nlm.nih.gov/20837533/) — Van Loo P, Nordgard SH, Lingjærde OC, et al. Allele-specific copy number analysis of tumors. Proc Natl Acad Sci USA. 2010;107(39):16910-5.

### HLA typing

- [OptiType](https://pubmed.ncbi.nlm.nih.gov/25143287/) — Szolek A, Schubert B, Mohr C, et al. OptiType: precision HLA typing from next-generation sequencing data. Bioinformatics. 2014;30(23):3310-6.
- [HLA-HD](https://pubmed.ncbi.nlm.nih.gov/28419628/) — Kawaguchi S, Higasa K, Shimizu M, Yamada R, Matsuda F. HLA-HD: An accurate HLA typing algorithm for next-generation sequencing data. Hum Mutat. 2017;38(7):788-797.

### RNA quantification and fusion calling

- [STAR](https://pubmed.ncbi.nlm.nih.gov/23104886/) — Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013;29(1):15-21.
- [kallisto](https://pubmed.ncbi.nlm.nih.gov/27043002/) — Bray NL, Pimentel H, Melsted P, Pachter L. Near-optimal probabilistic RNA-seq quantification. Nat Biotechnol. 2016;34(5):525-7.
- [salmon](https://pubmed.ncbi.nlm.nih.gov/28263959/) — Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. Salmon provides fast and bias-aware quantification of transcript expression. Nat Methods. 2017;14(4):417-419.
- [Arriba](https://pubmed.ncbi.nlm.nih.gov/33500328/) — Uhrig S, Ellermann J, Walther T, et al. Accurate and efficient detection of gene fusions from RNA sequencing data. Genome Res. 2021;31(3):448-460.
- [STAR-Fusion](https://www.biorxiv.org/content/10.1101/120295v1) — Haas BJ, Dobin A, Stransky N, et al. STAR-Fusion: Fast and Accurate Fusion Transcript Detection from RNA-Seq. bioRxiv. 2017.

### Neoantigen prediction

- [pVACtools](https://pubmed.ncbi.nlm.nih.gov/31907209/) — Hundal J, Kiwala S, McMichael J, et al. pVACtools: A Computational Toolkit to Identify and Visualize Cancer Neoantigens. Cancer Immunol Res. 2020;8(3):409-420.

## Software packaging and containerisation

- [Anaconda](https://anaconda.com) — Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016.
- [Bioconda](https://pubmed.ncbi.nlm.nih.gov/29967506/) — Grüning B, Dale R, Sjödin A, et al. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018;15(7):475-476.
- [BioContainers](https://pubmed.ncbi.nlm.nih.gov/28379341/) — da Veiga Leprevost F, Grüning B, Alves Aflitos S, et al. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017;33(16):2580-2582.
- [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241) — Merkel D. Docker: lightweight Linux containers for consistent development and deployment. Linux Journal. 2014;2014(239):2.
- [Singularity](https://pubmed.ncbi.nlm.nih.gov/28494014/) — Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017;12(5):e0177459.
