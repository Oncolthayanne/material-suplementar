dir.create("C:/Users/usuario/Desktop/Modelo R", recursive = TRUE)
.libPaths("C:/Users/usuario/Desktop/Modelo R")


BiocManager::install("TCGAWorkflow", force = TRUE)
BiocManager::install("TCGAWorkflowData", force = TRUE)
BiocManager::install("TCGAbiolinks", force = TRUE)
install.packages("lme4")

BiocManager::install("sesameData")
BiocManager::install("sesame")
BiocManager::install("reactome.db")
BiocManager::install("ReactomePA")
install.packages("tidyverse")
BiocManager::install("limma")
BiocManager::install("EDASeq")
BiocManager::install("edgeR")

BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", force = TRUE)
BiocManager::install("TCGAWorkflow", force = TRUE)
BiocManager::install(c("sesameData", "sesame", "reactome.db", "ReactomePA", "limma", "EDASeq", "edgeR"), force = TRUE)


#### A parte acima pode variar dependendo do computador, trata-se da instalacao de programas #########


# Carregando pacotes necessários
library(TCGAWorkflowData)  # Pacote para acessar dados do TCGA
library(DT)                # Pacote para visualização de dados em tabelas interativas
library(TCGAbiolinks)      # Pacote para acesso a dados do TCGA via Bioconductor
library(SummarizedExperiment)  # Pacote para manipulação de dados em matrizes nomeadas
library(EDASeq)
library(edgeR)



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("RTCGAToolbox")
# Biblioteca necessária
library(TCGAbiolinks)

# Query e download de dados de metilação para TCGA-BRCA (câncer de mama)
query_met_brca <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450"
)
GDCdownload(query_met_brca)
met_brca_450k <- GDCprepare(
  query = query_met_brca,
  summarizedExperiment = TRUE
)


# Query e download de dados de expressão gênica para TCGA-BRCA (câncer de mama)
query_exp_brca <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)
GDCdownload(query_exp_brca)
exp_brca <- GDCprepare(
  query = query_exp_brca
)


# Query e download de dados de variação no número de cópias mascaradas para TCGA-BRCA (câncer de mama)
query_masked_cnv_brca <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Copy Number Variation",
  data.type = "Masked Copy Number Segment",
  sample.type = c("Primary Tumor")
)
GDCdownload(query_masked_cnv_brca)
masked_cnv_brca <- GDCprepare(query_masked_cnv_brca)



# Query e download de dados de expressão gênica para TCGA-BRCA (câncer de mama)
query_exp_brca <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)
GDCdownload(query_exp_brca)
exp_brca <- GDCprepare(query_exp_brca)

# Exemplo de uso dos dados
# Obtendo matriz de expressão gênica
data <- assay(exp_brca)
# Visualizando como tabela interativa
datatable(
  data = data[1:10,], 
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = TRUE
)

# Obtendo informações sobre genes
genes.info <- rowRanges(exp_brca)

# Obtendo informações sobre amostras
sample.info <- colData(exp_brca)
datatable(
  data = as.data.frame(sample.info), 
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = FALSE
)


# Obtendo dados clínicos indexados para amostras de TCGA-BRCA (câncer de mama)
brca_clin <- GDCquery_clinic(
  project = "TCGA-BRCA", 
  type = "Clinical"
)



# Query e download de dados clínicos diretamente dos arquivos XML para TCGA-BRCA (câncer de mama)
query <- GDCquery(
  project = "TCGA-BRCA",
  data.format = "bcr xml",
  data.category = "Clinical"
) 
GDCdownload(query)
clinical <- GDCprepare_clinic(
  query = query, 
  clinical.info = "patient"
)



# Query e download de dados clínicos diretamente dos arquivos XML para TCGA-BRCA (câncer de mama)
query <- GDCquery(
  project = "TCGA-BRCA",
  data.format = "bcr xml",
  data.category = "Clinical"
) 
GDCdownload(query)
clinical <- GDCprepare_clinic(
  query = query, 
  clinical.info = "patient"
)

# Visualizando dados clínicos como tabela interativa
datatable(
  data = clinical, 
  options = list(scrollX = TRUE, keys = TRUE), 
  rownames = FALSE
)



# Query e download de dados clínicos diretamente dos arquivos XML para TCGA-BRCA (câncer de mama)
query <- GDCquery(
  project = "TCGA-BRCA",
  data.format = "bcr xml",
  data.category = "Clinical"
) 
GDCdownload(query)

# Preparando dados clínicos relacionados a medicamentos
clinical_drug <- GDCprepare_clinic(
  query = query, 
  clinical.info = "drug"
)

# Visualizando dados clínicos relacionados a medicamentos como tabela interativa
datatable(
  data = clinical_drug,
  options = list(scrollX = TRUE, keys = TRUE),
  rownames = FALSE
)

# Preparando dados clínicos relacionados à radioterapia
clinical_radiation <- GDCprepare_clinic(
  query = query, 
  clinical.info = "radiation"
)

# Visualizando dados clínicos relacionados à radioterapia como tabela interativa
datatable(
  data = clinical_radiation,
  options = list(scrollX = TRUE, keys = TRUE),
  rownames = FALSE
)


# Query e download de dados clínicos diretamente dos arquivos XML para TCGA-BRCA (câncer de mama)
query <- GDCquery(
  project = "TCGA-BRCA",
  data.format = "bcr xml",
  data.category = "Clinical"
) 
GDCdownload(query)

# Visualizando dados clínicos relacionados à radioterapia
clinical_radiation <- GDCprepare_clinic(
  query = query, 
  clinical.info = "radiation"
)

clinical_radiation |>
  datatable(
    options = list(scrollX = TRUE, keys = TRUE), 
    rownames = FALSE
  )

# Preparando dados clínicos administrativos
clinical_admin <- GDCprepare_clinic(
  query = query, 
  clinical.info = "admin"
)

clinical_admin |>
  datatable(
    options = list(scrollX = TRUE, keys = TRUE), 
    rownames = FALSE
  )

# Query para obter dados de variação somática no DNA para TCGA-BRCA
query_mutation <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open",
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(query_mutation)
maf <- GDCprepare(query_mutation)



# Carregar a biblioteca DT para criar tabelas interativas
library(DT)

# Visualizando os dados de variação somática no DNA em tabela interativa
maf[1:10,] |>
  datatable(
    options = list(scrollX = TRUE, keys = TRUE), 
    rownames = FALSE
  )



# Consultando subtipos de câncer para TCGA-BRCA
brca_subtypes <- TCGAquery_subtype(
  tumor = "brca"
)

# Visualizando subtipos de câncer
datatable(
  brca_subtypes[1:10,], 
  options = list(scrollX = TRUE, keys = TRUE), 
  rownames = FALSE
)

library(RTCGAToolbox)  # Pacote para acessar dados do Firehose do TCGA

# Obter a última data de execução do TCGA Firehose
lastRunDate <- getFirehoseRunningDates()[1]  # Obtém a última data de execução disponível

# Obter dados de metilação de DNA, RNAseq2 e dados clínicos para TCGA-BRCA
brca_data <- getFirehoseData(
  dataset = "BRCA",
  runDate = lastRunDate, 
  gistic2Date = getFirehoseAnalyzeDates(1),
  Methylation = FALSE,  
  clinical = TRUE,
  RNASeq2GeneNorm  = FALSE, 
  Mutation = TRUE,
  fileSizeLimit = 10000
)


# Obter dados de mutação para TCGA-BRCA
brca_mut <- getData(brca_data, "Mutation")

# Obter dados clínicos para TCGA-BRCA
brca_clin <- getData(brca_data, "clinical")

# Download dos resultados GISTIC (Genomic Identification of Significant Targets in Cancer)
lastanalyzedate <- getFirehoseAnalyzeDates(1)
gistic_brca <- getFirehoseData(
  dataset = "BRCA",
  GISTIC = TRUE, 
  gistic2Date = lastanalyzedate
)

# Download dos resultados GISTIC para TCGA-BRCA
lastanalyzedate <- getFirehoseAnalyzeDates(1)
gistic_brca <- getFirehoseData(
  dataset = "BRCA",
  GISTIC = TRUE, 
  gistic2Date = lastanalyzedate
)

# Obter resultados GISTIC (todas as alterações genômicas por gene)
gistic_allbygene_brca <- getData(
  object = gistic_brca, 
  type = "GISTIC", 
  platform = "AllByGene"
)

# Obter resultados GISTIC (alterações genômicas por gene com limiar)
gistic_thresholdedbygene_brca <- getData(
  object = gistic_brca, 
  type = "GISTIC", 
  platform = "ThresholdedByGene"
)

# Visualizar as primeiras linhas dos resultados GISTIC (todas as alterações genômicas por gene) 

head(gistic_allbygene_brca) 

# Visualizar as primeiras linhas dos resultados GISTIC (alterações genômicas por gene com limiar) 

head(gistic_thresholdedbygene_brca)

####### # Visualizar as primeiras linhas dos resultados GISTIC (todas as alterações genômicas por gene)
head(gistic_allbygene)

####### Visualizar as primeiras linhas dos resultados GISTIC (alterações genômicas por gene com limiar)
head(gistic_thresholedbygene)



# Carregar pacote maftools para análise de mutações somáticas
library(maftools)

# Recuperar dados do pacote TCGAWorkflowData.
maf


# Preparar os dados clínicos para análise de câncer de mama (TCGA-BRCA)
brca_clin <- GDCquery_clinic(project = "TCGA-BRCA", type = "Clinical")


# Criar um nome de coluna padronizado para Tumor_Sample_Barcode colnames(brca_clin)[grep("submitter_id", colnames(brca_clin))] <- "Tumor_Sample_Barcode"


# Criar uma variável binária para status de sobrevivência (1 para morto, 0 para vivo) brca_clin$Overall_Survival_Status <- 1 # morto brca_clin$Overall_Survival_Status[which(brca_clin$vital_status != "Dead")] <- 0

# Definir o tempo (`time`) com base nos dias até a morte (`days_to_death`) ou nos dias até o último acompanhamento (`days_to_last_follow_up`) para pacientes vivos
brca_clin$time <- brca_clin$days_to_death
brca_clin$time[is.na(brca_clin$days_to_death)] <- brca_clin$days_to_last_follow_up[is.na(brca_clin$days_to_death)]


# Criar objeto para usar no maftools
maf <- read.maf(
  maf = maf, 
  clinicalData = brca_clin, 
  isTCGA = TRUE
)

# Gerar sumário do MAF (Mutational Annotation Format)
plotmafSummary(
  maf = maf,
  rmOutlier = TRUE,
  addStat = 'median',
  dashboard = TRUE
)

# Adicionar análise de mutações utilizando o maftools
# Criar um oncoplot para visualizar as 10 principais mutações
oncoplot(
  maf = maf,
  top = 10,
  legendFontSize = 8,
  clinicalFeatures = c("tissue_or_organ_of_origin")
)

# Análise de sobrevivência baseada apenas no gene da TP53
plot <- mafSurvival(
  maf = maf,
  genes = "TP53",
  time = 'time',
  Status = 'Overall_Survival_Status',
  isTCGA = TRUE
)

# Análise de sobrevivência baseada na combinação de mutações dos genes CTNNB1 e TP53
plot <- mafSurvival(
  maf = maf,
  genes = c("CTNNB1", "TP53"),
  time = 'time',
  Status = 'Overall_Survival_Status',
  isTCGA = TRUE
)

# Verificar se BiocManager está instalado, caso contrário, instalá-lo
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Usar BiocManager para instalar o pacote org.Hs.eg.db
BiocManager::install("org.Hs.eg.db")

# Carregar o pacote
library(org.Hs.eg.db)




library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)

library(org.Hs.eg.db)
# Análise de Enriquecimento de Caminhos para Genes Mutados
# Carregar o pacote ReactomePA para enriquecimento de caminho
library(ReactomePA)

# Extrair os nomes únicos dos genes
unique_genes <- unique(maf@data$Hugo_Symbol)


# Converter para um data frame
unique_genes_df <- data.frame(Hugo_Symbol = unique_genes)

# Usar bitr para converter os símbolos de genes para IDs ENTREZ
gene_list <- bitr(unique_genes_df$Hugo_Symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)



# Realizando a análise de enriquecimento
enriched_pathways <- enrichPathway(
  gene = gene_list$ENTREZID,
  organism = "human"
)

# Plotando os resultados
dotplot(enriched_pathways)


# Instale ou carregue o pacote clusterProfiler
install.packages("clusterProfiler")
library(clusterProfiler)

# Análise de enriquecimento para KEGG
kegg_enrichment <- enrichKEGG(gene = gene_list$ENTREZID)

# Análise de enriquecimento para processos GO
go_enrichment <- enrichGO(gene = gene_list$ENTREZID, OrgDb = org.Hs.eg.db, ont = "MF")

dotplot(kegg_enrichment, showCategory = 20)  # Mostra os top 20 termos enriquecidos

dotplot(go_enrichment, showCategory = 20)  # Mostra os top 20 termos enriquecidos



#================EXEMPLO COM MENOS AMOSTRAS==================

# Query e download de dados de expressão gênica para TCGA-BRCA (exemplo com 20 amostras)
query_exp_brca <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)
# Obter apenas as primeiras 20 amostras para tornar o exemplo mais rápido
query_exp_brca$results[[1]] <- query_exp_brca$results[[1]][1:20,]
GDCdownload(query_exp_brca)
exp_brca <- GDCprepare(
  query = query_exp_brca
)

# Pré-processamento dos dados de expressão
exp_brca_preprocessed <- TCGAanalyze_Preprocessing(
  object = exp_brca,
  cor.cut = 0.6,    
  datatype = "unstranded",
  filename = "BRCA_IlluminaHiSeq_RNASeqV2.png"
)

# Normalização dos dados de expressão
exp_normalized <- TCGAanalyze_Normalization(
  tabDF = exp_brca_preprocessed,
  geneInfo = TCGAbiolinks::geneInfoHT,
  method = "gcContent"
)

# Filtragem dos dados de expressão normalizados
exp_filtered <- TCGAanalyze_Filtering(
  tabDF = exp_normalized,
  method = "quantile",
  qnt.cut = 0.25
)

# Para realizar uma análise de expressão diferencial (DEA), precisamos de dois grupos.
# Vamos criar dois grupos fictícios de amostras para exemplificação.
# Dividir as amostras de forma aleatória em dois grupos para DEA
set.seed(123)  # Definir semente para reprodutibilidade
sample_indices <- sample(1:ncol(exp_filtered), ncol(exp_filtered) / 2)
group1 <- exp_filtered[, sample_indices]
group2 <- exp_filtered[, -sample_indices]

# Análise de expressão diferencial (DEA) entre os dois grupos fictícios
diff_expressed_genes <- TCGAanalyze_DEA(
  mat1 = group1,
  mat2 = group2,
  Cond1type = "Group1",
  Cond2type = "Group2",
  fdr.cut = 0.01,
  logFC.cut = 1,
  method = "glmLRT"
)

# Número de genes diferencialmente expressos (DEG)
nrow(diff_expressed_genes)

#-------------------  4.2 EA: enrichment analysis --------------------
# Realizando a análise de enriquecimento
ansEA <- TCGAanalyze_EAcomplete(
  TFname = "DEA genes Group1 Vs Group2", 
  RegulonList = diff_expressed_genes$gene_name
)

# Visualizando a análise de enriquecimento em um gráfico de barras
TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = NULL,
  GOBPTab = ansEA$ResBP,
  nRGTab = diff_expressed_genes$gene_name,
  nBar = 20
)

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = NULL,
  GOCCTab = ansEA$ResCC,
  nRGTab = diff_expressed_genes$gene_name,
  nBar = 20
)


TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = NULL,
  GOMFTab = ansEA$ResMF,
  nRGTab = diff_expressed_genes$gene_name,
  nBar = 20
)

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = NULL,
  PathTab = ansEA$ResPat,
  nRGTab = rownames(diff_expressed_genes),
  nBar = 20
)


#O ESTUDO DAS VIAS ESTA SUPRIMIDO PORQUE PRECISA VERIFICAR NAS IMAGENS (KEGG, GO, ETC) GERADAS PARA COLOCAR NO SCRIPT

#================PARA VALER, SÓ RODA EM UMA MÁQUINA POTENTE==================


# Query e download de dados de expressão gênica para TCGA-BRCA
query_exp_brca <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)

# Download dos dados
GDCdownload(query_exp_brca)

# Preparação dos dados
exp_brca <- GDCprepare(
  query = query_exp_brca
)

# Pré-processamento dos dados de expressão
exp_brca_preprocessed <- TCGAanalyze_Preprocessing(
  object = exp_brca,
  cor.cut = 0.6,    
  datatype = "unstranded",
  filename = "BRCA_IlluminaHiSeq_RNASeqV2.png"
)

# Normalização dos dados de expressão
exp_normalized <- TCGAanalyze_Normalization(
  tabDF = exp_brca_preprocessed,
  geneInfo = TCGAbiolinks::geneInfoHT,
  method = "gcContent"
)

# Filtragem dos dados de expressão normalizados
exp_filtered <- TCGAanalyze_Filtering(
  tabDF = exp_normalized,
  method = "quantile",
  qnt.cut = 0.25
)

# Para realizar uma análise de expressão diferencial (DEA), precisamos de dois grupos.
# Vamos criar dois grupos fictícios de amostras para exemplificação.
# Dividir as amostras de forma aleatória em dois grupos para DEA
set.seed(123)  # Definir semente para reprodutibilidade
sample_indices <- sample(1:ncol(exp_filtered), ncol(exp_filtered) / 2)
group1 <- exp_filtered[, sample_indices]
group2 <- exp_filtered[, -sample_indices]

# Análise de expressão diferencial (DEA) entre os dois grupos fictícios
diff_expressed_genes <- TCGAanalyze_DEA(
  mat1 = group1,
  mat2 = group2,
  Cond1type = "Group1",
  Cond2type = "Group2",
  fdr.cut = 0.01,
  logFC.cut = 1,
  method = "glmLRT"
)

# Número de genes diferencialmente expressos (DEG)
nrow(diff_expressed_genes)

#-------------------  4.2 EA: enrichment analysis --------------------
# Realizando a análise de enriquecimento
ansEA <- TCGAanalyze_EAcomplete(
  TFname = "DEA genes Group1 Vs Group2", 
  RegulonList = diff_expressed_genes$gene_name
)

# Visualizando a análise de enriquecimento em um gráfico de barras
TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = NULL,
  GOBPTab = ansEA$ResBP,
  nRGTab = diff_expressed_genes$gene_name,
  nBar = 20
)

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = NULL,
  GOCCTab = ansEA$ResCC,
  nRGTab = diff_expressed_genes$gene_name,
  nBar = 20
)

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = NULL,
  GOMFTab = ansEA$ResMF,
  nRGTab = diff_expressed_genes$gene_name,
  nBar = 20
)

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = NULL,
  PathTab = ansEA$ResPat,
  nRGTab = rownames(diff_expressed_genes),
  nBar = 20
)






