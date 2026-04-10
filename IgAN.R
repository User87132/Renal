## microbiome IgAN RYU LAB
#Package install----
#BiocManager::install("ALDEx2")
#BiocManager::install("EnhancedVolcano")

library(tidyverse)
library(ALDEx2)
library(EnhancedVolcano)

# 1. Braken 데이터 불러오기 및 준비
setwd("~/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/data")
IgAN_inhouse = read.csv("genus_braken_merged", sep = "\t")
IgAN_inhouse = IgAN_inhouse[,-2]

# colnames(IgAN_inhouse) 에서 끝 이름이 bracken_frac 인 열들 모두 제거하기
# 열 이름에 'bracken_frac'이 포함된 열 제거
IgAN_inhouse <- IgAN_inhouse[, !grepl("bracken_frac", colnames(IgAN_inhouse))]
rownames(IgAN_inhouse) <- IgAN_inhouse$name
IgAN_inhouse$name <- NULL
IgAN_inhouse = IgAN_inhouse[,-1]
colnames(IgAN_inhouse)

#1. DA Volcano plot----
colnames(IgAN_inhouse) <- gsub("_9606_filtered.bracken_num", "", colnames(IgAN_inhouse))
colnames(IgAN_inhouse)
IgAN_inhouse_0 <- IgAN_inhouse

IgAN_inhouse = IgAN_inhouse_0[,c("C1","C3","C4","C5","I1","I2","I3","I4","I5","I6","I7")]
IgAN_condition = c(rep("Control", 4), rep("IgAN", 7))

#mc.samples
aldex128 <- aldex(IgAN_inhouse, IgAN_condition, mc.samples=128, test="t", effect=TRUE, 
                    include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE, gamma=NULL) #mc.samples: Monte Carlo 반복 수
aldex256 <- aldex(IgAN_inhouse, IgAN_condition, mc.samples=256, test="t", effect=TRUE, 
                  include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE, gamma=NULL)
common <- rownames(aldex128)
cor_effect_128_256 <- cor(aldex128[common,"effect"], aldex256[common,"effect"], use="complete.obs")
cor_q_128_256 <- cor(aldex128[common,"we.eBH"], aldex256[common,"we.eBH"], use="complete.obs")

rm(list=c("common","aldex128","aldex256","cor_effect_128_256","cor_q_128_256"))


#IgAN_aldex
IgAN_aldex <- aldex(IgAN_inhouse, IgAN_condition, mc.samples=128, test="t", effect=TRUE, 
                  include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE, gamma=NULL) #mc.samples: Monte Carlo 반복 수
IgAN_aldex_output <- data.frame(IgAN_aldex)

#filtering
cutoff_data <- IgAN_aldex_output[
  IgAN_aldex_output$wi.eBH < 0.5 & abs(IgAN_aldex_output$effect) > 0.5, 
]

selected_labels <- rownames(cutoff_data)
IgAN_inhouse_filtered <- IgAN_inhouse[rownames(IgAN_inhouse) %in% selected_labels, ]
common_rows <- intersect(rownames(IgAN_aldex_output), rownames(IgAN_inhouse_filtered))
IgAN_aldex_filtered <- IgAN_aldex_output[common_rows, ]

IgAN_aldex_filt <- rownames(IgAN_aldex_filtered)[head(order(abs(IgAN_aldex_filtered$effect), decreasing=TRUE), 10)]


#중요 두 균 표시----
key_microb <- c("Escherichia","Salmonella")
IgAN_inhouse_key <- IgAN_inhouse[rownames(IgAN_inhouse) %in% key_microb, ]
common_keys <- intersect(rownames(IgAN_aldex_output), rownames(IgAN_inhouse_key))

IgAN_aldex_filt <- c(IgAN_aldex_filt,common_keys)

#visualization
EnhancedVolcano(
  IgAN_aldex_output,
  lab = rownames(IgAN_aldex_output),
  x = 'effect',                # 효과 크기
  y = 'wi.eBH',                # p-value
  xlim = c(-4.7,3),
  ylim = c(0, 3),
  xlab = "Effect size",
  ylab = "-log10(p-adj)",
  pCutoff = 0.5,              # p-value cutoff
  FCcutoff = 0.5,             # 효과 크기 cutoff
  title = "Differential Abundance Analysis",
  subtitle = "Control vs IgAN",
  caption = "Data: Aldex2",
  pointSize = 3.5,
  col=c('gray','#00028C','gray','#FF0B55'),
  labSize = 5.5,
  cutoffLineCol = "#405D72",      # 가로/세로 기준선 색
  cutoffLineWidth = 0.8,      # 선 두께(원하면)
  cutoffLineType = "dotted",
  selectLab = IgAN_aldex_filt,  # 모든 라벨 표시
  drawConnectors = TRUE,                # 연결선 추가
  widthConnectors = 0.5,                # 연결선 두께
  colConnectors = "black",               # 연결선 색상
  arrowheads = TRUE,                   # 화살표 비활성화
  lengthConnectors = unit(0.02, "npc"),  # 연결선 길이
  max.overlaps = Inf
)
#p-value only: pCutoff는 통과했는데 FCcutoff는 못 넘은 점
#log2FC only: FCcutoff는 통과했는데 pCutoff는 못 넘은 점

ggsave("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/DA_microbiome_IgAN.png", dpi=1000, dev='png', height=8, width=8, units="in", bg = "white")

#중요 두 균 표시
ggsave("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/DA_key_microbiome_IgAN.png", dpi=1000, dev='png', height=8, width=8, units="in", bg = "white")



#2. Microbiome Proportion----
setwd("~/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/data")
IgAN_inhouse_percentage = read.csv("genus_braken_merged", sep = "\t")
IgAN_inhouse_percentage = IgAN_inhouse_percentage[,-c(2,3)]
colnames(IgAN_inhouse_percentage)

# 열 이름에 'bracken_num'이 포함된 열 제거
IgAN_inhouse_percentage <- IgAN_inhouse_percentage[, !grepl("bracken_num", colnames(IgAN_inhouse_percentage))]

IgAN_inhouse_percentage = IgAN_inhouse_percentage[,-c(2:7)]
IgAN_inhouse_percentage = IgAN_inhouse_percentage[,-c(3,7)]
IgAN_inhouse_percentage = IgAN_inhouse_percentage[,-c(6:9)]
IgAN_inhouse_percentage = IgAN_inhouse_percentage[,-c(13:29)]

rownames(IgAN_inhouse_percentage) <- IgAN_inhouse_percentage$name
IgAN_inhouse_percentage$name <- NULL
colnames(IgAN_inhouse_percentage)


# IgAN_inhouse_percentage 행이 모두 0.01 미만일 경우 제거 
IgAN_inhouse_percentage_filtered <- IgAN_inhouse_percentage[rowSums(IgAN_inhouse_percentage >= 0.01) > 0, ]

library(tidyr)
library(dplyr)
library(tibble)
colnames(IgAN_inhouse_percentage_filtered) <- gsub("_9606_filtered.bracken_frac", "", colnames(IgAN_inhouse_percentage_filtered))
colnames(IgAN_inhouse_percentage_filtered)

IgAN_inhouse_percentage_filtered <- IgAN_inhouse_percentage_filtered1

plot_data <- IgAN_inhouse_percentage_filtered %>%
  as.data.frame() %>%
  rownames_to_column("Microbe") %>%
  pivot_longer(-Microbe, names_to = "Sample", values_to = "Proportion")


#remotes::install_github("sqjin/CellChat")
library(CellChat)
library(ggplot2)

ggplot(plot_data, aes(x = Sample, y = Proportion, fill = Microbe)) +
  geom_bar(stat = "identity", position = "fill", alpha = 0.8) +
  scale_y_continuous(labels = scales::percent) +  # y축을 퍼센트로 표시
  labs(
    y = "Proportion",
    fill = "Microbe"
  ) + 
  scale_fill_manual(values = scPalette(39)) + 
  theme_minimal(base_size = 14) + 
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1.6),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 17),
    legend.position = "none",
    legend.text = element_text(size = 15),
    panel.grid = element_blank()  # 격자 제거
  )

ggsave("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/Prop_microbiome_IgAN.png", dpi=1000, dev='png', height=6.5, width=11, units="in", bg = "white")


#3. Microbiome Boxplot----
library(ggplot2)
library(ggpubr)

rownames(IgAN_inhouse_percentage_filtered) # 29개

for (i in 1:3) {
  microbe_name <- rownames(IgAN_inhouse_percentage_filtered)[i]
  microbe_data <- IgAN_inhouse_percentage_filtered[microbe_name, ]
  plot_data_microbe <- data.frame(
    Abundance = as.numeric(microbe_data),
    Group = IgAN_condition
  )
  p <- ggplot(plot_data_microbe, aes(x = Group, y = Abundance, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black", width = 0.4) +
    geom_jitter(width = 0.2, alpha = 0.8, size = 3, aes(color = Group)) +
    stat_compare_means(
      method = "t.test",
      label = "p.format",
      comparisons = list(c("Control", "IgAN")),
      vjust = -0.003,
      size = 7
    ) +
    labs(
      title = paste(microbe_name),
      subtitle = "Relative Abundance",
      x = "Group"
    ) +
    scale_fill_manual(values = c("Control" = "#424242", "IgAN" = "#e4a426")) +
    scale_color_manual(values = c("Control" = "#424242", "IgAN" = "#e4a426")) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
      plot.subtitle = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 19, face = "bold"),
      axis.text.x = element_text(size = 19, face = "bold"),
      axis.text.y = element_text(size = 19, face = "bold"),
      legend.position = "none"
    );p
  
  ggsave(
    filename = paste0("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/Box_microbiome_IgAN_", microbe_name, ".png"),
    plot = p,
    dpi = 1000,
    device = "png",
    height = 6.7,
    width = 6,
    units = "in",
    bg = "white"
  )
}



# 1. 모든 P-value 한 번에 모으기
library(tidyverse)

p_results <- data.frame(
  taxon = rownames(IgAN_inhouse_percentage_filtered),
  p_value = NA
)

for(i in 1:nrow(IgAN_inhouse_percentage_filtered)) {
  microbe_data <- IgAN_inhouse_percentage_filtered[i, ]
  test <- t.test(as.numeric(microbe_data) ~ IgAN_condition)
  p_results$p_value[i] <- test$p.value
}

# 2. FDR 보정 (BH 방법)
p_results$q_value <- p.adjust(p_results$p_value, method = "BH")

# 3. 결과 확인 (P < 0.05 & q < 0.05)
p_results %>% 
  filter(p_value < 0.05) %>% 
  arrange(p_value) %>%
  print()










#4. DEG Volcano plot----
setwd("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/data")

#BiocManager::install("Glimma")
#BiocManager::install("GEOquery")
library(limma)
library(Glimma)
library(edgeR)
library(R.utils)
library(GEOquery)
library(gplots)
library(limma)
library(edgeR)

only_count = read.table(file = "Inhouse_counts.txt", header = T, row.names = "Geneid")
only_count = only_count[,-c(1:5)]
colnames(only_count)

only_count1 = only_count[,-c(8)]
only_count1 = only_count1[,-c(11)]
only_count1 = only_count1[,-c(36)]
only_count1 = only_count1[,-c(36)]
only_count1 = only_count1[,-c(36)]
only_count1 = only_count1[,-c(34,35)]
colnames(only_count1)


colnames(only_count1) = c("A1", "A2", "A3", "A4", "A5", "A6", "C1", "C3", "C4", "C5", 
                          "F1", "F2", "F3", "F4", "I1", "I2", "I3", "I4", "I5", "I6", "I7", "M10", 
                          "M11", "M12", "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9")

#protein-coding 이외 filtering - biotype
#protein-coding인데 ENGS로 표시된 경우 포함
#BiocManager::install("biomaRt")
library(biomaRt)

mart <- useMart("ensembl")
datasets <- listDatasets(mart)
mart <- useDataset("hsapiens_gene_ensembl",mart)
anno <- getBM(c("ensembl_gene_id","hgnc_symbol","gene_biotype"), 
              filters="ensembl_gene_id", values=rownames(only_count1), mart=mart)
unique(anno$gene_biotype)

nrow(anno)
nrow(only_count1[grepl("^ENSG", rownames(only_count1)),])

anno <- anno[anno$gene_biotype == "protein_coding",]
anno$gene_name <- paste0(anno$ensembl_gene_id, "protein_coding")

ensg <- grepl("^ENSG", rownames(only_count1))

keep <- !ensg | (ensg & rownames(only_count1) %in% anno$ensembl_gene_id)  #ENSG아닌 행 or ENSG행 중에서 protein_coding (둘다)

only_count1 <- only_count1[keep, , drop = FALSE] #[행,열,]
only_count <- only_count1
rm(only_count1)

#행의 합계가 10이상인 경우만 filtering
nrow(only_count)
only_count <- only_count[rowSums(only_count) >= 10,]
#40819 -> 27700


diseases_group <- c("AKI", "AKI", "AKI", "AKI", "AKI", "AKI", "Control", "Control", "Control", "Control",
                    "CKD", "CKD", "CKD", "AKI", "IgAN", "IgAN", "IgAN", "IgAN", "IgAN", "IgAN", "IgAN", "CKD", 
                    "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD")
ID = c("A1", "A2", "A3", "A4", "A5", "A6", "C1", "C3", "C4", "C5", 
       "F1", "F2", "F3", "F4", "I1", "I2", "I3", "I4", "I5", "I6", "I7", "M10", 
       "M11", "M12", "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9")

diseases_group <- factor(diseases_group, levels = c("Control", "IgAN", "AKI", "CKD"))

library(edgeR)
d0 <- DGEList(counts = only_count, group = diseases_group, samples = ID) #edgeR 형식 객체 #lib.size = 샘플의 전체 카운트 합
d0 <- calcNormFactors(d0) #총량+조성 차이 보정

#CPM: 시퀀싱 깊이(lib.size) 보정한 상대 발현 기준
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d0 <- d0[-drop,] 

plotMDS(d0, col = as.numeric(diseases_group))

mm <- model.matrix(~0+diseases_group)
y <- voom(d0, mm, plot = T)
fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(diseases_groupIgAN - diseases_groupControl , levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
top.table$adj.P.Val
length(which(top.table$adj.P.Val < 0.05))


# count값 80%이상 filtering
zero <- which(rowMeans(only_count == 0) >= 0.8)
d0$counts <- d0$counts[-zero,]

keep.exprs <- filterByExpr(d0, group=diseases_group)
d0$counts <- d0$counts[keep.exprs,]
d0$samples$lib.size <- colSums(d0$counts)

d <- d0
rm(d0)

#데이터 정규화
##회의안건 TMM, 라이브러리 크기로 조정##
d <- calcNormFactors(d, method = "TMM")

design <- model.matrix(~0 + diseases_group)
colnames(design) <- levels(diseases_group)  # Control, IgAN, AKI, CKD
design

contr.matrix <- makeContrasts(
  IgANvsControl = IgAN-Control,  #IgAN, Control 비교
  levels = design)

#카운트 데이터 이분산성 제거
##회의안건## voom 사용 (라이브러리 사이즈 압축보정)
par(mfrow=c(1,2))
v <- voom(d, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
res_efit <- topTable(efit,number=Inf,adjust="fdr")

#7. 이 과정에서 v 다시 생성함. AKI, CKD 제거함.

## volcanoplot
library(EnhancedVolcano)
max(res_efit$adj.P.Val)

# logFC 및 P.Value에 기반하여 색상 설정
keyvals <- ifelse(res_efit$adj.P.Val > 0.031 | abs(res_efit$logFC) < 1, '#A5A5A5', 
                  ifelse(res_efit$logFC < 0, '#1f77b4', '#ff7f0e'))

# 이름 설정
names(keyvals)[keyvals == '#A5A5A5'] <- 'Non-significant'
names(keyvals)[keyvals == '#1f77b4'] <- 'Downregulated'
names(keyvals)[keyvals == '#ff7f0e'] <- 'Upregulated'

EnhancedVolcano(res_efit,
                lab = rownames(res_efit),
                x = 'logFC',
                y = 'P.Value',
                title = NULL,
                subtitle = NULL,
                pCutoff = 1e-2, 
                FCcutoff = 1,
                pointSize = 1.5,
                labSize = 5.7,
                xlim = c(-4.7, 4.7),
                ylim = c(0, 7),
                axisLabSize = 21,
                caption = NULL,
                legendLabSize = 20,
                colCustom = keyvals
)

ggsave("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/DEG_volcano_IgAN.png", dpi=1000, dev='png', height=8, width=9, units="in", bg = "white")


#5. GO pathway----
library(clusterProfiler)
library(enrichplot)
require(DOSE)
library(tidyverse)
#BiocManager::install("org.Hs.eg.db")
library('org.Hs.eg.db', character.only = TRUE)
library(EnhancedVolcano)

FB_list <- res_efit$logFC
df <- res_efit %>% 
  rownames_to_column(var = "gene")
names(FB_list) <- df$gene
FB_list = sort(FB_list, decreasing = TRUE)

set.seed(716)
gse_FB <- gseGO(geneList=FB_list, 
                ont ="all", 
                keyType = "SYMBOL",     #change keytype here, only symbol works for this data
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE, 
                OrgDb = org.Hs.eg.db, 
                pAdjustMethod = "BH", eps = 0) #BH; Benjamini-Hochberg method, eps: 작은 p-value 더 낮게 추정.
IgAN_GO_Res = gse_FB@result
#pvalue: BH 다중검정 보정 이전 원래 p값. #p.adjust: pvalue에 대해 BH 방식으로 계산된 FDR 기반 보정 p값. #qvalue: 각 결과가 유의하다고 판단될 때 최소 FDR(거짓양성) 수준.

nrow(IgAN_GO_Res)
IgAN_GO_Res <- IgAN_GO_Res[IgAN_GO_Res$p.adjust<0.05,]

saveRDS(gse_FB, file = "/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/data/GO_gse_FB_IgAN.rds")

#이후 과정은 'IgAN_GSE175759.R' file에.


#6. IgAN, IgAN_GEO Common pathway bar plot----
#file: IgAN_GSE175759.R
IgAN_GO_Res_filt

#7.GO genes Heat map----
#IgAN 관련 NES Top3 Pathway, Top10 genes heat map
top3 <- IgAN_GO_Res_filt[1:3,]

#top3$core_enrichment #Leading edge genes = "pathway enrichment 신호를 주도한 핵심 유전자 집합
core_list <- strsplit(as.character(top3$core_enrichment), "/")

picked <- character(0)
picked_by_pathway <- vector("list", length(core_list))

#아래 필터링한 res_efit과 겹치는 유전자 개수 30개 이하이기 때문에 head 제외.
for (i in seq_along(core_list)) {
  genes_i <- unique(core_list[[i]])
  genes_i <- genes_i[!genes_i %in% picked]   # 앞에서 뽑힌 중복 제거
  picked_by_pathway[[i]] <- genes_i
  picked <- c(picked, picked_by_pathway[[i]])
}

picked_by_pathway
picked_all <- unlist(picked_by_pathway, use.names = FALSE)
length(picked_all) #86


#limma-voom(v, res_efit)
sig_genes <- rownames(res_efit)[res_efit$adj.P.Val < 0.05 & abs(res_efit$logFC) >= 1]
sig_genes <- sig_genes[
  order(abs(res_efit[sig_genes, "t"]), decreasing = TRUE)
]
sig_genes <- intersect(sig_genes, picked_all)
length(sig_genes) #18
sig_genes <- intersect(sig_genes, rownames(v$E))
length(sig_genes) #18

IgAN_sig_genes_data <- res_efit[rownames(res_efit) %in% sig_genes,]
saveRDS(IgAN_sig_genes_data, "/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/data/IgAN_sig_genes_data.rds")


#Heatmap----
IgAN_sig_genes_data <- readRDS("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/data/IgAN_sig_genes_data.rds")
sig_genes <- rownames(IgAN_sig_genes_data)

library(ComplexHeatmap)
library(circlize)

#v: AKI, CKD 포함된 상태에서 만든 voom 결과
#v 다시 생성
#앞 부분은 res_efit만 이용해서 괜찮음.
diseases_group <- c("AKI", "AKI", "AKI", "AKI", "AKI", "AKI", "Control", "Control", "Control", "Control",
                    "CKD", "CKD", "CKD", "AKI", "IgAN", "IgAN", "IgAN", "IgAN", "IgAN", "IgAN", "IgAN", "CKD", 
                    "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD")
diseases_group <- factor(diseases_group, levels = c("Control", "IgAN", "AKI", "CKD"))
keep <- diseases_group %in% c("Control", "IgAN")
diseases_group <- droplevels(diseases_group[keep])


#########두 번째에는 skip 가능
d2 <- d[, keep, keep.lib.sizes = FALSE]
d2 <- calcNormFactors(d2, method = "TMM")

design2 <- model.matrix(~0 + diseases_group)
colnames(design2) <- levels(diseases_group)
v <- voom(d2, design2, plot = TRUE)
saveRDS(v, "/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/data/voom_IgAN.rds")
#================================================================================================


v <- readRDS("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/data/voom_IgAN.rds")
mat <- v[sig_genes, , drop = FALSE]
group_vec <- diseases_group
names(group_vec) <- colnames(v$E)              # v$E 컬럼명과 매칭
group_vec <- group_vec[colnames(mat)]     # expr_mat 컬럼 순서로 재정렬
ord <- order(group_vec)
mat <- mat[, ord, drop = FALSE]
group_vec <- group_vec[ord]

png("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/Heatmap_IgAN.png", width = 8, height = 6, units = "in", res = 600)
expr_mat_scaled <- t(scale(t(v$E[sig_genes,])))
expr_mat_scaled <- expr_mat_scaled[, ord, drop = FALSE]
Heatmap(expr_mat_scaled,
        col = colorRamp2(c(-2,0,2), c('#3396D3',"white",'#F67400')),
        column_split = group_vec,
        name = "Z-score",
        column_title_gp = gpar(fontsize=14))
dev.off()



#8. IgAN genes(18) & IgAN_GEO genes(102) Common genes----
#file: IgAN_GSE175759.R
Com_genes <- intersect(IgAN_sig_genes_data, IgAN_GEO_top3_genes)
length(Com_genes) #3. #"FCGR2B" "PRKCB"  "IGHG4" 


#9. Top10 pathway 공통 유전자 correlation----
IgAN_GO_Res_filt
core_list <- strsplit(as.character(IgAN_GO_Res_filt$core_enrichment), "/")
picked <- character(0)
picked_by_pathway <- vector("list", length(core_list))


for (i in seq_along(core_list)) {
  genes_i <- unique(core_list[[i]])
  genes_i <- genes_i[!genes_i %in% picked]   # 앞에서 뽑힌 중복 제거
  picked_by_pathway[[i]] <- genes_i
  picked <- c(picked, picked_by_pathway[[i]])
}

picked_by_pathway
picked_all <- unlist(picked_by_pathway, use.names = FALSE)
length(picked_all) #278



diseases_group <- c("AKI", "AKI", "AKI", "AKI", "AKI", "AKI", "Control", "Control", "Control", "Control",
                    "CKD", "CKD", "CKD", "AKI", "IgAN", "IgAN", "IgAN", "IgAN", "IgAN", "IgAN", "IgAN", "CKD", 
                    "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD", "CKD")
ID = c("A1", "A2", "A3", "A4", "A5", "A6", "C1", "C3", "C4", "C5", 
       "F1", "F2", "F3", "F4", "I1", "I2", "I3", "I4", "I5", "I6", "I7", "M10", 
       "M11", "M12", "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9")

diseases_group <- factor(diseases_group, levels = c("Control", "IgAN", "AKI", "CKD"))
keep <- diseases_group %in% c("Control", "IgAN")
diseases_group <- droplevels(diseases_group[keep])

design <- model.matrix(~0 + diseases_group)
colnames(design) <- levels(diseases_group)
design

contr.matrix <- makeContrasts(
  IgANvsControl = IgAN-Control,  #IgAN, Control 비교
  levels = design)

v <- readRDS("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/data/voom_IgAN.rds")
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
res_efit <- topTable(efit,number=Inf,adjust="fdr")

sig_genes <- rownames(res_efit)[res_efit$adj.P.Val < 0.05 & abs(res_efit$logFC) >= 1]
sig_genes <- sig_genes[
  order(abs(res_efit[sig_genes, "t"]), decreasing = TRUE)
]
sig_genes <- intersect(sig_genes, picked_all)
length(sig_genes) #75
sig_genes <- intersect(sig_genes, rownames(v$E))
length(sig_genes) #75


Com_genes <- intersect(sig_genes, sig_genes_GEO)
length(Com_genes) #13

diff <- setdiff(sig_genes, Com_genes)
all_genes <- c(diff[1:17], Com_genes) 

mat <- v$E[Com_genes, , drop = FALSE]
group_vec <- diseases_group
names(group_vec) <- colnames(v$E)              # v$E 컬럼명과 매칭
group_vec <- group_vec[colnames(mat)]     # expr_mat 컬럼 순서로 재정렬
group_vec <- factor(group_vec, levels = c("Control", "IgAN"))
desired_order <- names(group_vec)
ord <- order(group_vec)
mat <- mat[, ord, drop = FALSE]
group_vec <- group_vec[ord]
group_vec <- group_vec[desired_order]

png("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/Heatmap_IgAN.png", width = 8, height = 6, units = "in", res = 600)
expr_mat_scaled <- t(scale(t(v$E[Com_genes,])))
expr_mat_scaled <- expr_mat_scaled[, desired_order, drop = FALSE]
Heatmap(expr_mat_scaled,
        col = colorRamp2(c(-2,0,2), c('#3396D3',"white",'#F67400')),
        column_split = group_vec,
        name = "Z-score",
        column_order = desired_order,
        column_title_gp = gpar(fontsize=14))
dev.off()



#Final heatmap with micorbiome abundance----
library(ComplexHeatmap)
library(circlize)
library(grid)

taxa1 <- "Escherichia"
taxa2 <- "Salmonella"

#microbiome_ratio <- t(microbiome_ratio)
stopifnot(all(colnames(mat) %in% colnames(microbiome_ratio)))
sum_vec <- as.numeric(microbiome_ratio[taxa1, colnames(mat)] + microbiome_ratio[taxa2, colnames(mat)])
names(sum_vec) <- colnames(mat)

# 4) ord 순서로 재정렬 (유전자 heatmap과 동일한 샘플 순서)
sum_vec <- sum_vec[ord]

# 5) (선택) 값이 0~작은 값 많으면 log1p 권장
sum_vec_plot <- log1p(sum_vec)   # 필요 없으면 sum_vec 그대로 써도 됨

# 6) 1줄짜리 매트릭스로 만들기
sum_mat <- matrix(sum_vec_plot, nrow = 1)
rownames(sum_mat) <- paste0(taxa1, "+", taxa2)
colnames(sum_mat) <- names(sum_vec_plot)

# 7) 색상 범위 잡기(데이터 기반)
rng <- range(sum_mat, na.rm = TRUE)
col_sum <- colorRamp2(c(rng[1], rng[2]), c("#FAF3F3", "#FF0000"))

# 8) 기존 유전자 Z-score heatmap
expr_mat_scaled <- t(scale(t(v$E[Com_genes, , drop=FALSE])))
expr_mat_scaled <- expr_mat_scaled[,desired_order, drop = FALSE]
expr_mat_scaled <- expr_mat_scaled[order(rownames(expr_mat_scaled)), , drop = FALSE]
group_vec2 <- factor(group_vec[desired_order], levels = c("Control","IgAN"))
ro <- seq_len(nrow(expr_mat_scaled))
names(ro) <- rownames(expr_mat_scaled)

ht_gene <- Heatmap(expr_mat_scaled,
                   col = colorRamp2(c(-2,0,2), c('#3396D3',"white",'#F67400')),
                   column_split = group_vec2,
                   name = "Z-score",
                   column_order = desired_order,
                   row_order = ro,
                   column_title_gp = gpar(fontsize=14))

ht_micro <- Heatmap(sum_mat,
                    name = "Sum",
                    col = col_sum,
                    show_row_names = TRUE,
                    show_column_names = FALSE,
                    cluster_rows = FALSE,
                    cluster_columns = FALSE,
                    column_order = desired_order,
                    height = unit(4, "mm"))

# 10) 세로로 붙이기
png("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/Heatmap_IgAN_with_microbiome.png", width = 8, height = 6, units = "in", res = 600)
draw(ht_gene %v% ht_micro)
dev.off()









rna_normalized  = v$E
rna_normalized <- as.data.frame(t(rna_normalized))

setwd("~/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/data")
IgAN_inhouse_percentage = read.csv("genus_braken_merged", sep = "\t")
IgAN_inhouse_percentage = IgAN_inhouse_percentage[,-c(2,3)]
IgAN_inhouse_percentage <- IgAN_inhouse_percentage[, !grepl("bracken_num", colnames(IgAN_inhouse_percentage))]
IgAN_inhouse_percentage = IgAN_inhouse_percentage[,-c(2:7)]
IgAN_inhouse_percentage = IgAN_inhouse_percentage[,-c(3,7)]
IgAN_inhouse_percentage = IgAN_inhouse_percentage[,-c(6:9)]
IgAN_inhouse_percentage = IgAN_inhouse_percentage[,-c(13:29)]
rownames(IgAN_inhouse_percentage) <- IgAN_inhouse_percentage$name
IgAN_inhouse_percentage$name <- NULL
IgAN_inhouse_percentage_filtered <- IgAN_inhouse_percentage[rowSums(IgAN_inhouse_percentage >= 0.01) > 0, ] #frac이 0.01보다 큰 샘플이 최소 1개
colnames(IgAN_inhouse_percentage_filtered) <- gsub("_9606_filtered.bracken_frac", "", colnames(IgAN_inhouse_percentage_filtered))
microbiome_ratio = IgAN_inhouse_percentage_filtered
microbiome_ratio <- as.data.frame(t(microbiome_ratio))

#Target gene filtering
rna_normalized <- rna_normalized[,colnames(rna_normalized) %in% Com_genes]

#target microbiome filtering
microbiome_ratio <- microbiome_ratio[,colnames(microbiome_ratio) %in% c('Escherichia', 'Salmonella')]

#RNA-seq와 microbiome 데이터를 합치기
combined_data <- cbind(rna_normalized, microbiome_ratio)
combined_data <- as.matrix(combined_data)

#cor, p-value 계산
#install.packages("psych", dependencies=TRUE)
library(psych)
ncol(rna_normalized)
ncol(combined_data)
cor_results <- corr.test(combined_data[,1:13], combined_data[,14:15], 
                         method="spearman")

# 4. 결과
print(cor_results$r)  # correlation
print(cor_results$p)  # p-value

cor_mat <- cor_results$r
cor_test <- cor_results$p

ID <- colnames(IgAN_inhouse_percentage_filtered)
condition <- c(rep("Control", 4), rep("IgAN", 7))
sample_data <- data.frame(SampleID = ID, Condition = condition)


library(ggplot2)
library(dplyr)

gene_of_interest <- "IRF4"
taxa_of_interest <- "Escherichia"
#Escherichia
#Salmonella

# 숫자 벡터 추출 (drop=TRUE로 벡터로)
gene_vec <- rna_normalized[, gene_of_interest, drop = TRUE]
taxa_vec <- microbiome_ratio[, taxa_of_interest, drop = TRUE]

# (중요) 샘플 순서/길이 체크
stopifnot(length(gene_vec) == nrow(sample_data))
stopifnot(length(taxa_vec) == nrow(sample_data))

plot_df <- data.frame(
  Gene = as.numeric(gene_vec),
  Taxa = as.numeric(taxa_vec),
  Condition = sample_data$Condition
) %>%
  filter(!is.na(Gene), !is.na(Taxa), Taxa > 0)

# corr.test 결과에서 해당 gene-taxa 값 꺼내기
correlation <- cor_mat[gene_of_interest, taxa_of_interest]
p_value <- cor_test[gene_of_interest, taxa_of_interest]

ggplot(plot_df, aes(x = Gene, y = Taxa, color = Condition, fill = Condition)) +
  geom_point(size = 2, alpha = 0.8, shape = 21, stroke = 1.0) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.1, alpha = 0.12, fullrange = TRUE) +
  annotate(
    "text",
    x = mean(range(plot_df$Gene, na.rm = TRUE)),           # 가운데(범위 기준)
    y = 0.8,
    label = paste0("R = ", signif(correlation, 3), "\nP = ", signif(p_value, 3)),
    size = 5, hjust = 0.5, vjust = 1, color = "black"
  ) +
  labs(
    x = paste0(gene_of_interest, " expression"),
    y = paste0(taxa_of_interest, " abundance"),
    color = "Condition"
  ) +
  scale_color_manual(values = c("Control" = "#3396D3", "IgAN" = "#F67400")) +
  scale_fill_manual(values = c("Control" = "#3396D3", "IgAN" = "#F67400")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  filename = paste0("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/Correlation_", gene_of_interest,"_",taxa_of_interest,".png"),
  plot = last_plot(),          # 방금 그린 ggplot
  width = 6, height = 4,       # inch
  dpi = 300
)


#10. Escherichia graph----
base <- "/home/eunhye/raw_data1/eunhye/IgAN_renal_disease/DEH_IgAN_CODE/Inhouse/extract"
for (i in 1:7) {
  f <- file.path(base, sprintf("I%d_9606_filtered.bracken", i))
  assign(paste0("I", i), read.delim(f, stringsAsFactors = FALSE), envir = .GlobalEnv)
  df_name <- paste0("I",i)
  df <- get(df_name)
  df <- df[grepl('Escherichia',df$name),]
  df$sample <- df_name
  assign(df_name, df)
}
IgAN_Escherichia <- do.call(rbind, mget(paste0("I", 1:7)))
colnames(IgAN_Escherichia)[1] <- 'Microbe'


library(ggplot2)
unique(IgAN_Escherichia$Microbe)
IgAN_Escherichia$Microbe <- factor(
  IgAN_Escherichia$Microbe,
  levels = c("Escherichia albertii", "Escherichia coli", "Escherichia fergusonii")
)

IgAN_Escherichia$sample <- factor(
  IgAN_Escherichia$sample,
  levels = rev(unique(IgAN_Escherichia$sample))
)


IgAN_Escherichia0 <- aggregate(
  fraction_total_reads ~ Microbe,
  data = IgAN_Escherichia,
  FUN = mean,
  na.rm = TRUE
)
colnames(IgAN_Escherichia0)[2] <- "fraction_total_reads_mean"

ggplot(IgAN_Escherichia0, aes(x = "", y = fraction_total_reads_mean, fill = Microbe)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c(
    "Escherichia albertii"   = "steelblue",
    "Escherichia coli"       = "#FF5656",
    "Escherichia fergusonii" = "#8E7DBE"
  )) +
  theme_void() +
  labs(fill = "Microbe")
ggsave("/home/eunhye/raw_data1/eunhye/IgAN_renal_disease/DEH_IgAN_CODE/figure/Escherichia_plot_IgAN.png", dpi=1000, dev='png', height=9, width=12, units="in", bg = "white")



