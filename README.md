# cottonSaltGRN
Scripts used to analyze the 11 public transcriptome datasets in cotton related to salt and salt-alkali responses and tolerance

1. [transcriptome.md](transcriptome.md): Download public RNA-seq dataset, preprocessing, reference genome mapping, and transcript quantification
2. [batchEffectAnalysis.R](batchEffectAnalysis.R): Preliminerary analysis of TPM table, PCA and batch effect analysis
3. [WGCNA.R](WGCNA.R): Weighted Gene Coexpression Network Analysis
4. [GO.R](GO.R): GO enrichment analysis
5. [Pheatmap.R](Pheatmap.R): Heatmap and other analysis of GhGDHs expression patterns.

Other data used:
*  11 transcriptome datas were downloaded in NCBI under BioProject accession numbers [PRJNA531727](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA531727), [PRJNA601953](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA601953), [PRJNA623201](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA623201), [PRJNA253112](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA253112), [PRJNA722118](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA722118), [PRJNA248163](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA248163), [PRJNA482027](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA482027), [PRJNA485838](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA485838), [PRJNA532694](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA532694), [PRJNA490626](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA490626), and [PRJNA919499](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA919499).
* The raw sequence data generated have been deposited in China National Center for Bioinformation under GSA: [CRA016187](https://ngdc.cncb.ac.cn/gsa/search?searchTerm=CRA016187)
 
