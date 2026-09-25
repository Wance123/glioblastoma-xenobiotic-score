suppressPackageStartupMessages(library(MCPcounter))
root <- Sys.getenv('JOB49_ROOT', unset='.')
geo <- Sys.getenv('GEO_LOCAL_ROOT')
bulk <- file.path(Sys.getenv('JOB46_ROOT'),'reanalysis_v2','data')
out <- file.path(root,'results','virtual_validation')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
genes <- read.delim(file.path(root,'input','MCPcounter_genes.txt'),check.names=FALSE,stringsAsFactors=FALSE)
probes <- read.delim(file.path(root,'input','MCPcounter_probesets.txt'),header=FALSE,stringsAsFactors=FALSE)
# KYNU overlaps the tested Hallmark set; CSF1R overlaps the original five-gene
# marker panel. Exclude both from the monocytic reference.
genes <- genes[!(genes[['Cell population']]=='Monocytic lineage' & genes[['HUGO symbols']] %in% c('KYNU','CSF1R')),]
inputs <- list(
  TCGA=file.path(bulk,'TCGA_GBM_log_expression.tsv.gz'),
  CGGA693=file.path(bulk,'CGGA693_GBM_log_expression.tsv.gz'),
  CGGA325=file.path(bulk,'CGGA325_GBM_log_expression.tsv.gz'),
  GSE16011=file.path(geo,'GSE16011_genesymbol.txt.gz'),
  GSE149009=file.path(geo,'GSE149009_genesymbol.txt.gz'))
for (cohort in names(inputs)) {
  d <- read.delim(gzfile(inputs[[cohort]]),check.names=FALSE,row.names=1)
  d <- d[!duplicated(rownames(d)),,drop=FALSE]
  mat <- as.matrix(d)
  storage.mode(mat) <- 'double'
  if (cohort %in% c('GSE16011','GSE149009')) mat <- t(mat)
  scores <- MCPcounter.estimate(mat,featuresType='HUGO_symbols',probesets=probes,genes=genes)
  write.csv(data.frame(sample=colnames(scores),t(scores),check.names=FALSE),
            file.path(out,paste0(cohort,'_MCPcounter.csv')),row.names=FALSE)
  cat(cohort,nrow(mat),'genes',ncol(mat),'samples\n')
}
write.csv(genes,file.path(out,'MCPcounter_genes_without_KYNU_CSF1R.csv'),row.names=FALSE)
