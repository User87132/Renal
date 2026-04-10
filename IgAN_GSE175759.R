#0. pakages----
library(tidyverse)
library(ALDEx2)
library(EnhancedVolcano)

setwd("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/data")

#1. DA Volcano plot----
#IgAN_GSE175759 data
IgAN_GSE175759 = read.csv("GSE175759_9606_extract_genus_braken_merged", sep = "\t")
IgAN_GSE175759 = IgAN_GSE175759[,-c(2,3)]

#remove 'bracken_frac'
IgAN_GSE175759 <- IgAN_GSE175759[, !grepl("bracken_frac", colnames(IgAN_GSE175759))]

rownames(IgAN_GSE175759) <- IgAN_GSE175759$name
IgAN_GSE175759$name <- NULL

colnames(IgAN_GSE175759)
condition <- c(rep("IgAN", 46), rep("Control", 22))

#ALDEx2
#고처리량 시퀀싱 count 데이터의 DEG 분석 패키지.
#유전자 수를 그대로 비교하지 않고 샘플 안에서 서로의 비율이 얼마나 달라졌는지를 반복 계산해, 변화는 작지만 우연히 유의해 보이는 신호를 제외하는 차등발현 분석 방법.
#DESeq2/edgeR 패키지(절대 count를 모형화)


#mc.samples
aldex128 <- aldex(IgAN_GSE175759, condition, mc.samples=128, test="t", effect=TRUE, 
                  include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE, gamma=NULL) #mc.samples: Monte Carlo 반복 수
aldex256 <- aldex(IgAN_GSE175759, condition, mc.samples=256, test="t", effect=TRUE, 
                  include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE, gamma=NULL)
common <- rownames(aldex128)
cor_effect_128_256 <- cor(aldex128[common,"effect"], aldex256[common,"effect"], use="complete.obs")
cor_effect_128_256
cor_q_128_256 <- cor(aldex128[common,"we.eBH"], aldex256[common,"we.eBH"], use="complete.obs")
cor_q_128_256
#1에 가까우면 128로 설정.

rm(list=c("common","aldex128","aldex256","cor_effect_128_256","cor_q_128_256"))
IgAN_aldex <- aldex(IgAN_GSE175759, condition, mc.samples=128, test="t", effect=TRUE, 
                    include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE, gamma=NULL)
aldex_output <- data.frame(IgAN_aldex)


#filtering
cutoff_data <- aldex_output[
  aldex_output$wi.eBH < 0.2 & abs(aldex_output$effect) > 0.5, 
]
selected_labels <- rownames(cutoff_data)
IgAN_GSE175759_filtered <- IgAN_GSE175759[rownames(IgAN_GSE175759) %in% selected_labels, ]
common_rows <- intersect(rownames(aldex_output), rownames(IgAN_GSE175759_filtered))
aldex_filtered <- aldex_output[common_rows, ]
aldex_filt <- rownames(aldex_filtered)[head(order(abs(aldex_filtered$effect), decreasing=TRUE), 10)]


#중요 두 균 표시----
key_microb <- c("Escherichia","Salmonella")
IgAN_GSE175759_key <- IgAN_GSE175759[rownames(IgAN_GSE175759) %in% key_microb, ]
common_keys <- intersect(rownames(aldex_output), rownames(IgAN_GSE175759_key))

aldex_filt <- c(aldex_filt,common_keys)


#visualization
EnhancedVolcano(
  aldex_output,
  lab = rownames(aldex_output),
  x = 'effect',                # 효과 크기
  y = 'wi.eBH',                # p-value
  xlim = c(-1,1),
  ylim = c(0, 1.5),
  xlab = "Effect size",
  ylab = "-log10(p-adj)",
  pCutoff = 0.5,              # p-value cutoff
  FCcutoff = 0.5,             # 효과 크기 cutoff
  title = "Differential Abundance Analysis",
  subtitle = "Control vs IgAN(GEO)",
  caption = "Data: Aldex2",
  pointSize = 3.5,
  col=c('gray','#00028C','gray','#FF0B55'),
  labSize = 5.5,
  cutoffLineCol = "#405D72",      # 가로/세로 기준선 색
  cutoffLineWidth = 0.8,      # 선 두께(원하면)
  cutoffLineType = "dotted",
  selectLab = aldex_filt,  # 모든 라벨 표시
  drawConnectors = TRUE,                # 연결선 추가
  widthConnectors = 0.5,                # 연결선 두께
  colConnectors = "black",               # 연결선 색상
  arrowheads = TRUE,                   # 화살표 비활성화
  lengthConnectors = unit(0.02, "npc"),  # 연결선 길이
  max.overlaps = Inf
)
#p-value only: pCutoff는 통과했는데 FCcutoff는 못 넘은 점
#log2FC only: FCcutoff는 통과했는데 pCutoff는 못 넘은 점

ggsave("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/DA_microbiome_IgAN_GSE175759.png", dpi=1000, dev='png', height=8, width=8, units="in", bg = "white")

#중요 두 균 표시
ggsave("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/DA_key_microbiome_IgAN_GSE175759.png", dpi=1000, dev='png', height=8, width=8, units="in", bg = "white")




#2. Microbiome Proportion----
setwd("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/data")
IgAN_GSE175759_percentage = read.csv("GSE175759_9606_extract_genus_braken_merged", sep = "\t")
IgAN_GSE175759_percentage = IgAN_GSE175759_percentage[,-c(2,3)]
IgAN_GSE175759_percentage <- IgAN_GSE175759_percentage[, !grepl("bracken_num", colnames(IgAN_GSE175759_percentage))]
rownames(IgAN_GSE175759_percentage) <- IgAN_GSE175759_percentage$name
IgAN_GSE175759_percentage$name <- NULL

IgAN_GSE175759_percentage_filtered <- IgAN_GSE175759_percentage[rowSums(IgAN_GSE175759_percentage >= 0.01) > 0, ]

library(tidyr)
library(dplyr)
colnames(IgAN_GSE175759_percentage_filtered) <- gsub(".genus_bracken_frac", "", colnames(IgAN_GSE175759_percentage_filtered))

IgAN_GSE175759_percentage_filtered <- IgAN_GSE175759_percentage_filtered1

plot_data <- IgAN_GSE175759_percentage_filtered %>%
  as.data.frame() %>%
  rownames_to_column("Microbe") %>%
  pivot_longer(-Microbe, names_to = "Sample", values_to = "Proportion")


IgAN_factor <- as.factor(c("SRR14678579","SRR14678580","SRR14678581","SRR14678582","SRR14678583","SRR14678584",
                           "SRR14678585","SRR14678586","SRR14678587","SRR14678588","SRR14678589","SRR14678590","SRR14678591","SRR14678592",
                           "SRR14678593","SRR14678594","SRR14678595","SRR14678596","SRR14678597","SRR14678598","SRR14678599","SRR14678600",
                           "SRR14678533","SRR14678534","SRR14678535","SRR14678536","SRR14678537","SRR14678538","SRR14678539","SRR14678540",
                           "SRR14678541","SRR14678542","SRR14678543","SRR14678544","SRR14678545","SRR14678546","SRR14678547","SRR14678548",
                           "SRR14678549","SRR14678550","SRR14678551","SRR14678552","SRR14678553","SRR14678554","SRR14678555","SRR14678556",
                           "SRR14678557","SRR14678558","SRR14678559","SRR14678560","SRR14678561","SRR14678562","SRR14678563","SRR14678564",
                           "SRR14678565","SRR14678566","SRR14678567","SRR14678568","SRR14678569","SRR14678570","SRR14678571","SRR14678572",
                           "SRR14678573","SRR14678574","SRR14678575","SRR14678576","SRR14678577","SRR14678578"))
plot_data$Sample <- factor(plot_data$Sample, levels = IgAN_factor)

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
    axis.text.y = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 17),
    legend.position = "right",
    legend.text = element_text(size = 15),
    panel.grid = element_blank()  # 격자 제거
  )

ggsave("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/Prop_microbiome_IgAN_GSE175759.png", dpi=1000, dev='png', height=6, width=14, units="in", bg = "white")


#3. Microbiome Boxplot----
library(ggplot2)
library(ggpubr)

rownames(IgAN_GSE175759_percentage_filtered)

for (i in 20:21) {
  microbe_name <- rownames(IgAN_GSE175759_percentage_filtered)[i]
  microbe_data <- IgAN_GSE175759_percentage_filtered[microbe_name, ]
  plot_data_microbe <- data.frame(
    Abundance = as.numeric(microbe_data),
    Group = condition
  )
  p <- ggplot(plot_data_microbe, aes(x = Group, y = Abundance, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black", width = 0.4) +  # 박스 폭 좁게
    geom_jitter(width = 0.2, alpha = 0.8, size = 2, aes(color = Group)) +
    stat_compare_means(
      method = "t.test",
      label = "p.format",
      comparisons = list(c("Control", "IgAN")),
      vjust = 0.09,
      size = 7
    ) +
    labs(
      title = paste(microbe_name),
      subtitle = "Relative Abundance",
      x = "Group"
    ) +
    scale_fill_manual(values = c("Control" = "#424242", "IgAN" = "#A252E8")) +
    scale_color_manual(values = c("Control" = "#424242", "IgAN" = "#A252E8")) +
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
    filename = paste0("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/Box_microbiome_IgAN_GSE175759_", microbe_name, ".png"),
    plot = p,
    dpi = 1000,
    device = "png",
    height = 6.7,
    width = 6,
    units = "in",
    bg = "white"
  )
}


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

IgAN_GSE175759 = read.table(file = "GSE175759_counts.txt", header = T, row.names = "Geneid")
only_count = IgAN_GSE175759[,-c(1:5)]
rm(IgAN_GSE175759)
colnames(only_count)
condition <- c(rep("I", 46), rep("C", 22))

#protein-coding 이외 filtering - biotype
#protein-coding인데 ENGS로 표시된 경우 포함
#BiocManager::install("biomaRt")
library(biomaRt)

mart <- useMart("ensembl")
datasets <- listDatasets(mart)
mart <- useDataset("hsapiens_gene_ensembl",mart)
anno <- getBM(c("ensembl_gene_id","hgnc_symbol","gene_biotype"), 
              filters="ensembl_gene_id", values=rownames(only_count), mart=mart)
unique(anno$gene_biotype)

nrow(anno)
nrow(only_count[grepl("^ENSG", rownames(only_count)),])

anno <- anno[anno$gene_biotype == "protein_coding",]
anno$gene_name <- paste0(anno$ensembl_gene_id, "protein_coding")

ensg <- grepl("^ENSG", rownames(only_count))
keep <- !ensg | (ensg & rownames(only_count) %in% anno$ensembl_gene_id)  #ENSG아닌 행 or ENSG행 중에서 protein_coding (둘다)
only_count <- only_count[keep, , drop = FALSE] #[행,열,]

diseases_group <- c(rep("IgAN", 46), rep("Control", 22))

IgAN_IDs <- paste0("IgAN", 1:46)       # IgAN 1번부터 46번
Control_IDs <- paste0("Control", 1:22)   # Control 1번부터 22번
ID <- c(IgAN_IDs, Control_IDs)
colnames(only_count) = ID
colnames(only_count)


###최소 행 합계가 10개인 유전자 << 이부분 진행하기 !!!!!!
nrow(only_count) #40819 -> 30278
only_count <- only_count[rowSums(only_count) >= 10,]

new_levels <- c("Control", "IgAN")
diseases_group <- factor(diseases_group, levels = new_levels)
levels(factor(diseases_group))

#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
library(edgeR)

d0 <- DGEList(counts = only_count, group = diseases_group, samples = ID)
d0 <- calcNormFactors(d0)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

rm(d0)

group_colors <- c("#424242", "#A353E8")  # 그룹별 색상 지정
col <- group_colors[as.numeric(diseases_group)]  # diseases_group에 따른 색상 매핑
png("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/MDS plot.png", width = 4.5, height = 4.5, units = 'in', res = 1000)
plotMDS(d, col = col, main = "MDS Plot Control and IgAN")
dev.off()


mm <- model.matrix(~0+diseases_group)
y <- voom(d, mm, plot = T)
fit <- lmFit(y, mm)


## IgAN and Control
contr <- makeContrasts(diseases_groupIgAN - diseases_groupControl , levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)

head(top.table, 20)
top.table$adj.P.Val
length(which(top.table$adj.P.Val < 0.05))


# count값 80%이상 filtering
zero <- which(rowMeans(only_count == 0) >= 0.8)
d$counts <- d$counts[-zero,]

keep.exprs <- filterByExpr(d, group=diseases_group)
d$counts <- d$counts[keep.exprs,]
d$samples$lib.size <- colSums(d$counts)



#데이터 정규화
##회의안건 TMM, 라이브러리 크기로 조정##
d <- calcNormFactors(d, method = "TMM")

design <- model.matrix(~0+diseases_group)
colnames(design)
colnames(design) <- c('Control','IgAN')

contr.matrix <- makeContrasts(
  IgANvsControl = IgAN-Control, 
  levels = design)


#카운트 데이터 이분산성 제거
par(mfrow=c(1,2))
v <- voom(d, design, plot=T)
#saveRDS(v, file = "/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/data/voom_IgAN_GEO.rds")
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
res_efit <- topTable(efit,number=Inf,adjust="BH")


#volcanoplot
library(EnhancedVolcano)
#logFC 및 P.Value에 기반하여 색상 설정
keyvals <- ifelse(res_efit$adj.P.Val > 0.031 | abs(res_efit$logFC) < 1, '#A5A5A5', 
                  ifelse(res_efit$logFC < 0, '#E8C326', '#A353E8'))
names(keyvals)[keyvals == '#A5A5A5'] <- 'Non-significant'
names(keyvals)[keyvals == '#E8C326'] <- 'Downregulated'
names(keyvals)[keyvals == '#A353E8'] <- 'Upregulated'

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
                xlim = c(-3, 3),
                ylim = c(0, 24.5),
                axisLabSize = 21,
                caption = NULL,
                legendLabSize = 20,
                colCustom = keyvals
)

ggsave("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/DEG_volcano_IgAN_GSE175759.png", dpi=1000, dev='png', height=7, width=8, units="in", bg = "white")



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
#NES: Normalized enrichment scores
saveRDS(gse_FB, file = "/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/data/GO_gse_FB_IgAN_GSE175759.rds")


#visualization
gse_FB_GSE <- readRDS("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/data/GO_gse_FB_IgAN_GSE175759.rds")
IgAN_GO_Res_GSE <- gse_FB_GSE@result
View(IgAN_GO_Res_GSE)

IgAN_GSE_top10 <- IgAN_GO_Res_GSE %>% slice_head(n=10)

IgAN_GSE_top10 <- IgAN_GSE_top10 %>% mutate(NES_sign = ifelse(NES<0,"neg","pos"))
ggplot(IgAN_GSE_top10, aes(NES, fct_reorder(Description, NES), fill = NES_sign)) + 
  geom_barh(stat='identity', color="#000000", width = 0.65) + 
  scale_fill_manual(values = c(pos = "#85409D", neg = "#FFC400", guide = "none"))+
  theme_minimal() + 
  ylab(NULL) + 
  theme(axis.text.x = element_text(size = 23, colour="black"),
        axis.text.y = element_text(size = 20, colour="black",lineheight = 0.7),
        axis.title.x = element_text(size = 23, colour="black"),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.5),
        legend.position="none")+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 53))

ggsave("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/Bar_IgAN_GSE175759.png", dpi=1000, dev='png', height=6.7, width=13, units="in", bg = "white")





#6. IgAN & IgAN _GSE175759 Common pathway----
#1) Common pathway - NES Top 10(figure X)
gse_FB <- readRDS("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/data/GO_gse_FB_IgAN.rds")
IgAN_GO_Res <- gse_FB@result
common_pathway <- merge(
  IgAN_GO_Res[, c("Description", "ID", "ONTOLOGY", "NES", "p.adjust")],
  IgAN_GO_Res_GSE[, c("Description", "ID", "ONTOLOGY", "NES", "p.adjust")],
  by = "Description",
  suffixes = c("_IgAN", "_IgAN_GSE")
)

common_pathway <- common_pathway %>%
  dplyr::arrange(dplyr::desc(abs(NES_IgAN)), dplyr::desc(abs(NES_IgAN_GSE)))

common_pathway <- common_pathway$Description
common_pathway <- common_pathway[1:10]


#(1) IgAN Bar plot
library(clusterProfiler)
library(ggplot2)
library(ggstance)

IgAN_GO_Res_filt <- IgAN_GO_Res[IgAN_GO_Res$Description %in% common_pathway,]
nrow(IgAN_GO_Res) #431
nrow(IgAN_GO_Res_filt)

IgAN_GO_Res_filt <- IgAN_GO_Res_filt %>% mutate(NES_sign = ifelse(NES<0,"neg","pos"))
ggplot(IgAN_GO_Res_filt, aes(NES, fct_reorder(Description, NES), fill = NES_sign)) + 
  geom_barh(stat='identity', color="#000000", width = 0.65) + 
  scale_fill_manual(values = c(pos = '#FF6600', neg = '#0D5EA6', guide = "none"))+
  theme_minimal() + 
  ylab(NULL) + 
  theme(axis.text.x = element_text(size = 18, colour="black"),
        axis.text.y = element_text(size = 16.3, colour="black",lineheight = 0.3),
        axis.title.x = element_text(size = 18, colour="black"),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.5),
        legend.position="none")

ggsave("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/Bar_com_IgAN.png", dpi=1000, dev='png', height=5, width=11, units="in", bg = "white")


#(2) IgAN_GSE175759 Bar plot
IgAN_GO_Res_GSE_filt <- IgAN_GO_Res_GSE[IgAN_GO_Res_GSE$Description %in% common_pathway,]
nrow(IgAN_GO_Res_GSE) #431
nrow(IgAN_GO_Res_GSE_filt)

IgAN_GO_Res_GSE_filt <- IgAN_GO_Res_GSE_filt %>% mutate(NES_sign = ifelse(NES<0,"neg","pos"))

ggplot(IgAN_GO_Res_GSE_filt, aes(NES, fct_reorder(Description, NES), fill = NES_sign)) + 
  geom_barh(stat='identity', color="#000000", width = 0.65) + 
  scale_fill_manual(values = c(pos = "#85409D", neg = "#FFC400", guide = "none"))+
  theme_minimal() + 
  ylab(NULL) + 
  theme(axis.text.x = element_text(size = 18, colour="black"),
        axis.text.y = element_text(size = 16.3, colour="black",lineheight = 0.3),
        axis.title.x = element_text(size = 18, colour="black"),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.5),
        legend.position="none")

ggsave("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/Bar_com_IgAN_GSE175759.png", dpi=1000, dev='png', height=5, width=11, units="in", bg = "white")


#2) Common pathway - IgAN, NES Top 10----
#IgAN pathway, NES filtering => FIGURE.3
common_pathway <- merge(
  IgAN_GO_Res[, c("Description", "ID", "ONTOLOGY", "NES", "p.adjust")],
  IgAN_GO_Res_GSE[, c("Description", "ID", "ONTOLOGY", "NES", "p.adjust")],
  by = "Description",
  suffixes = c("_IgAN", "_IgAN_GSE")
)
common_pathway[,1:2]

IgAN_core_GOs <- c(
  # Core Adaptive Immune (1-5)
  "GO:0002253",  # activation of immune response
  "GO:0002250",  # adaptive immune response  
  "GO:0002460",  # Ig superfamily recombination
  "GO:0050851",  # antigen receptor signaling
  "GO:0042113",  # B cell activation
  
  # B cell/Ig Production (6-10)
  "GO:0050853",  # B cell receptor signaling
  "GO:0006959",  # humoral immune response
  "GO:0071735",  # IgG complex
  "GO:0019814",  # immunoglobulin complex
  "GO:0001816",  # cytokine production
  
  # Immune Signaling (11-15)
  "GO:0002429",  # immune response-activating receptor
  "GO:0002757",  # immune response-activating signaling
  "GO:0002252",  # immune effector process
  "GO:0046649",  # lymphocyte activation
  "GO:0030098",  # lymphocyte differentiation
  
  # Leukocyte Effector (16-20)
  "GO:0002449",  # lymphocyte mediated immunity
  "GO:0002821",  # positive regulation adaptive immune
  "GO:0002824",  # Ig superfamily positive regulation
  "GO:0001819",  # positive regulation cytokine production
  "GO:0050778",  # positive regulation immune response
  
  # Additional Core (21-24)
  "GO:0045089",  # positive regulation innate immune
  "GO:0032615",  # IL-12 production
  "GO:0002521",  # leukocyte differentiation
  "GO:0001909"   # leukocyte mediated cytotoxicity
)

common_pathway <- common_pathway[common_pathway$ID_IgAN %in% IgAN_core_GOs,]
common_pathway <- common_pathway$Description
#IgAN 관련 공통 pathway 24개


#(1) IgAN Bar plot----
IgAN_GO_Res_filt <- IgAN_GO_Res[IgAN_GO_Res$Description %in% common_pathway,]
nrow(IgAN_GO_Res) #431
nrow(IgAN_GO_Res_filt)

IgAN_GO_Res_filt <- IgAN_GO_Res_filt %>%
  arrange(desc(abs(NES))) %>%
  slice_head(n = 10)

IgAN_GO_Res_filt <- IgAN_GO_Res_filt %>% mutate(NES_sign = ifelse(NES<0,"neg","pos"))
ggplot(IgAN_GO_Res_filt, aes(NES, fct_reorder(Description, NES), fill = NES_sign)) + 
  geom_barh(stat='identity', color="#000000", width = 0.65) + 
  scale_fill_manual(values = c(pos = '#FF6600', neg = '#0D5EA6', guide = "none"))+
  theme_minimal() + 
  ylab(NULL) + 
  theme(axis.text.x = element_text(size = 23, colour="black"),
        axis.text.y = element_text(size = 20, colour="black",lineheight = 0.8),
        axis.title.x = element_text(size = 23, colour="black"),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.5),
        legend.position="none")+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 50))

ggsave("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/Bar_IgANcom_IgAN.png", dpi=1000, dev='png', height=7, width=12, units="in", bg = "white")


#(2) IgAN_GSE175759 Bar plot----
IgAN_GO_Res_GSE_filt <- IgAN_GO_Res_GSE[IgAN_GO_Res_GSE$Description %in% common_pathway,]
nrow(IgAN_GO_Res_GSE) #431
nrow(IgAN_GO_Res_GSE_filt)

IgAN_GO_Res_GSE_filt <- IgAN_GO_Res_GSE_filt %>%
  arrange(desc(abs(NES))) %>%
  slice_head(n = 10)

IgAN_GO_Res_GSE_filt <- IgAN_GO_Res_GSE_filt %>% mutate(NES_sign = ifelse(NES<0,"neg","pos"))

ggplot(IgAN_GO_Res_GSE_filt, aes(NES, fct_reorder(Description, NES), fill = NES_sign)) + 
  geom_barh(stat='identity', color="#000000", width = 0.65) + 
  scale_fill_manual(values = c(pos = "#85409D", neg = "#FFC400", guide = "none"))+
  theme_minimal() + 
  ylab(NULL) + 
  theme(axis.text.x = element_text(size = 23, colour="black"),
        axis.text.y = element_text(size = 20, colour="black",lineheight = 0.8),
        axis.title.x = element_text(size = 23, colour="black"),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.5),
        legend.position="none")+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 50))

ggsave("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/Bar_IgANcom_IgAN_GSE175759.png", dpi=1000, dev='png', height=7, width=12, units="in", bg = "white")




#7.GO genes Heat map----
#Top10 pathway 공통 유전자 13개 correlation 결과
IgAN_GO_Res_GSE_filt
core_list <- strsplit(as.character(IgAN_GO_Res_GSE_filt$core_enrichment), "/")
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
length(picked_all) #338

diseases_group <- c(rep("IgAN", 46), rep("Control", 22))
diseases_group <- factor(diseases_group, levels = c("Control", "IgAN"))

design <- model.matrix(~0 + diseases_group)
colnames(design) <- levels(diseases_group)
design

contr.matrix <- makeContrasts(
  IgANvsControl = IgAN-Control,  #IgAN, Control 비교
  levels = design)

v <- readRDS("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/data/voom_IgAN_GEO.rds")
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
res_efit <- topTable(efit,number=Inf,adjust="fdr")

sig_genes_GEO <- rownames(res_efit)[res_efit$adj.P.Val < 0.05 & abs(res_efit$logFC) >= 1]
sig_genes_GEO <- sig_genes_GEO[
  order(abs(res_efit[sig_genes_GEO, "t"]), decreasing = TRUE)
]
sig_genes_GEO <- intersect(sig_genes_GEO, picked_all)
length(sig_genes_GEO) #117
sig_genes_GEO <- intersect(sig_genes_GEO, rownames(v$E))
length(sig_genes_GEO) #117


Com_genes <- intersect(sig_genes, sig_genes_GEO)
length(Com_genes) #13

diff <- setdiff(sig_genes_GEO, Com_genes)

all_genes <- c(diff[1:17], Com_genes) #Com_genes(13) + diff(17)

mat <- v$E[Com_genes, , drop = FALSE]
group_vec <- diseases_group
names(group_vec) <- colnames(v$E)              # v$E 컬럼명과 매칭
group_vec <- group_vec[colnames(mat)]     # expr_mat 컬럼 순서로 재정렬
group_vec <- factor(group_vec, levels = c("Control", "IgAN"))
desired_order <- names(group_vec)
ord <- order(group_vec)
mat <- mat[, ord, drop = FALSE]
group_vec <- group_vec[desired_order]

png("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/Heatmap_IgAN_GSE175759.png", width = 8, height = 6, units = "in", res = 600)
expr_mat_scaled <- t(scale(t(v$E[Com_genes,])))
expr_mat_scaled <- expr_mat_scaled[, desired_order, drop = FALSE]
Heatmap(expr_mat_scaled,
        col = colorRamp2(c(-2,0,2), c('#FFCD39',"white",'#A353E8')),
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

microbiome_ratio = IgAN_GSE175759_percentage_filtered
microbiome_ratio <- as.data.frame(microbiome_ratio)

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
expr_mat_scaled <- t(scale(t(v$E[Com_genes,])))
expr_mat_scaled <- expr_mat_scaled[, desired_order, drop = FALSE]
expr_mat_scaled <- expr_mat_scaled[order(rownames(expr_mat_scaled)), , drop = FALSE]
sum_mat <- sum_mat[, desired_order, drop = FALSE]
ro <- seq_len(nrow(expr_mat_scaled))
names(ro) <- rownames(expr_mat_scaled)

ht_gene <- Heatmap(expr_mat_scaled,
                   col = colorRamp2(c(-2,0,2), c('#FFCD39',"white",'#A353E8')),
                   column_split = group_vec,
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
png("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/Heatmap_IgAN_with_microbiome_GSE175759.png", width = 8, height = 5, units = "in", res = 600)
draw(ht_gene %v% ht_micro)
dev.off()


#9. Microbiome&Common genes correlation----
#Top10 pathway 공통 유전자 13개 correlation 결과
rna_normalized  = v$E
rna_normalized <- as.data.frame(t(rna_normalized))

setwd("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/data")
IgAN_GSE175759_percentage = read.csv("GSE175759_9606_extract_genus_braken_merged", sep = "\t")
IgAN_GSE175759_percentage = IgAN_GSE175759_percentage[,-c(2,3)]
IgAN_GSE175759_percentage <- IgAN_GSE175759_percentage[, !grepl("bracken_num", colnames(IgAN_GSE175759_percentage))]
rownames(IgAN_GSE175759_percentage) <- IgAN_GSE175759_percentage$name
IgAN_GSE175759_percentage$name <- NULL
IgAN_GSE175759_percentage_filtered <- IgAN_GSE175759_percentage[rowSums(IgAN_GSE175759_percentage >= 0.01) > 0, ]
colnames(IgAN_GSE175759_percentage_filtered) <- gsub(".genus_bracken_frac", "", colnames(IgAN_GSE175759_percentage_filtered))
microbiome_ratio = IgAN_GSE175759_percentage_filtered
microbiome_ratio <- as.data.frame(t(microbiome_ratio))


library(ggplot2)
library(dplyr)
IgAN_IDs <- paste0("IgAN", 1:46)       # IgAN 1번부터 46번
Control_IDs <- paste0("Control", 1:22)   # Control 1번부터 22번

# 원하는 ID 벡터 생성
ID <- c(IgAN_IDs, Control_IDs)
colnames(IgAN_GSE175759_percentage_filtered) <- ID
condition <- c(rep("IgAN", 46), rep("Control", 22))
sample_data <- data.frame(SampleID = ID, Condition = condition)

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
ncol(combined_data)
cor_results <- corr.test(combined_data[,1:13], combined_data[,14:15], 
                         method="spearman")

# 4. 결과
print(cor_results$r)  # correlation
print(cor_results$p)  # p-value

cor_mat <- cor_results$r
cor_test <- cor_results$p

gene_of_interest <- "TLR7"
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
    y = 1,
    label = paste0("R = ", signif(correlation, 3), "\nP = ", signif(p_value, 3)),
    size = 5, hjust = 0.5, vjust = 1, color = "black"
  ) +
  labs(
    x = paste0(gene_of_interest, " expression"),
    y = paste0(taxa_of_interest, " abundance"),
    color = "Condition"
  ) +
  scale_color_manual(values = c("Control" = "#E8C223", "IgAN" = "#A252E8")) +
  scale_fill_manual(values = c("Control" = "#E8C223", "IgAN" = "#A252E8")) +
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







#10. Escherichia species circle graph----
base <- "/home/eunhye/raw_data1/eunhye/IgAN_renal_disease/DEH_IgAN_CODE/GSE175759_bracken/species"
#SRR14678533~SRR14678578 :IgAN Sample

for (i in 33:78) {
  f <- file.path(base, sprintf("SRR146785%d.species_bracken", i))
  df <- read.delim(f, stringsAsFactors = FALSE)
  df_name <- paste0("SRR146785",i)
  df <- df[grepl('Escherichia',df$name),]
  if (nrow(df) == 0) next
  df$sample <- df_name
  colnames(df)[1] <- "Microbe"
  assign(df_name, df)
}
objs <- paste0("SRR146785", 33:78)
dfs <- mget(intersect(objs, ls()), envir = .GlobalEnv)
IgAN_Escherichia <- do.call(rbind, dfs)
unique(IgAN_Escherichia$sample)

library(ggplot2)
unique(IgAN_Escherichia$Microbe)

aa <- objs
aa <- as.data.frame(aa)
aa$ID <- paste0("I",1:46)
IgAN_Escherichia$IgAN_sample <- aa$ID[match(IgAN_Escherichia$sample, aa$aa)]
IgAN_Escherichia$IgAN_sample <- factor(
  IgAN_Escherichia$IgAN_sample,
  levels = rev(unique(IgAN_Escherichia$IgAN_sample))
)

IgAN_Escherichia1 <- aggregate(
  fraction_total_reads ~ Microbe,
  data = IgAN_Escherichia,
  FUN = mean,
  na.rm = TRUE
)
colnames(IgAN_Escherichia1)[2] <- "fraction_total_reads_mean"

library(ggplot2)

ggplot(IgAN_Escherichia1, aes(x = "", y = fraction_total_reads_mean, fill = Microbe)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c(
    "Escherichia albertii"   = "steelblue",
    "Escherichia coli"       = "#FF5656",
    "Escherichia fergusonii" = "#8E7DBE",
    "Escherichia marmotae" = "#FFCF50"
  )) +
  theme_void() +
  labs(fill = "Microbe")

ggsave("/home/eunhye/raw_data1/eunhye/IgAN_renal_disease/DEH_IgAN_CODE/figure/Escherichia_plot_IgAN_GSE.png", dpi=1000, dev='png', height=9, width=12, units="in", bg = "white")






