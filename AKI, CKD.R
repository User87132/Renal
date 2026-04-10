## microbiome AKI, CKD RYU LAB
#Package install----
#BiocManager::install("ALDEx2")
#BiocManager::install("EnhancedVolcano")

library(tidyverse)
library(ALDEx2)
library(EnhancedVolcano)

#1. AKI, CKD Data----
# Braken 데이터 불러오기 및 준비
setwd("~/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/data")
IgAN_inhouse = read.csv("genus_braken_merged", sep = "\t")
IgAN_inhouse = IgAN_inhouse[,-2]

# colnames(IgAN_inhouse) 에서 끝 이름이 bracken_frac 인 열들 모두 제거하기
# 열 이름에 'bracken_frac'이 포함된 열 제거
IgAN_inhouse <- IgAN_inhouse[, !grepl("bracken_frac", colnames(IgAN_inhouse))]
aa = IgAN_inhouse$name
IgAN_inhouse = IgAN_inhouse[,-1]
rownames(IgAN_inhouse) = aa

head(IgAN_inhouse)

#1. AKI, CKD DA volcano plot----
colnames(IgAN_inhouse) <- gsub("_9606_filtered.bracken_num", "", colnames(IgAN_inhouse))
colnames(IgAN_inhouse)
IgAN_inhouse_0 <- IgAN_inhouse
colnames(IgAN_inhouse_0)

AKI_inhouse = IgAN_inhouse_0[, c("C1", "C3", "C4", "C5", "A1","A2","A3","A4","A5","A6","F4")]
CKD_inhouse = IgAN_inhouse_0[, c("C1", "C3", "C4", "C5", "M1","M2","M3","M4","M5","M6","M7","M8","M9","M10","M11","M12","F1","F2","F3")]

AKI_condition = c(rep("Control", 4), rep("AKI", 7))
CKD_condition = c(rep("Control", 4), rep("CKD", 15))

#mc.samples
aldex128 <- aldex(CKD_inhouse, CKD_condition, mc.samples=128, test="t", effect=TRUE, 
                  include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE, gamma=NULL) #mc.samples: Monte Carlo 반복 수
aldex256 <- aldex(CKD_inhouse, CKD_condition, mc.samples=256, test="t", effect=TRUE, 
                  include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE, gamma=NULL)
common <- rownames(aldex128)
cor_effect_128_256 <- cor(aldex128[common,"effect"], aldex256[common,"effect"], use="complete.obs")
cor_q_128_256 <- cor(aldex128[common,"we.eBH"], aldex256[common,"we.eBH"], use="complete.obs")
cor_effect_128_256
cor_q_128_256
rm(list=c("common","aldex128","aldex256","cor_effect_128_256","cor_q_128_256"))


AKI_aldex <- aldex(AKI_inhouse, AKI_condition, mc.samples=128, test="t", effect=TRUE, 
                   include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE, gamma=NULL) #mc.samples: Monte Carlo 반복 수
CKD_aldex <- aldex(CKD_inhouse, CKD_condition, mc.samples=128, test="t", effect=TRUE, 
                   include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE, gamma=NULL) #mc.samples: Monte Carlo 반복 수

#AKI
AKI_aldex_output <- data.frame(AKI_aldex)

#필터링
cutoff_data <- AKI_aldex_output[
  AKI_aldex_output$wi.eBH < 0.5 & abs(AKI_aldex_output$effect) > 0.5, 
]

selected_labels <- rownames(cutoff_data)

#AKI_inhouse_filtered는 rownames(AKI_inhouse)가 selected_labels인 것들이다. 
AKI_inhouse_filtered <- AKI_inhouse[rownames(AKI_inhouse) %in% selected_labels, ]

## rownames(aldex_output) 에서 rownames(AKI_inhouse_filtered) 와 같은 행만 추출하기 
# 공통 행 이름 추출
common_rows <- intersect(rownames(AKI_aldex_output), rownames(AKI_inhouse_filtered))
AKI_aldex_filtered <- AKI_aldex_output[common_rows, ]
head(AKI_aldex_filtered) # 필터링된 데이터 확인

AKI_aldex_filt <- rownames(AKI_aldex_filtered)[head(order(abs(AKI_aldex_filtered$effect), decreasing=TRUE), 10)]


#중요 두 균 표시----
key_microb <- c("Escherichia","Salmonella")
AKI_inhouse_key <- AKI_inhouse[rownames(AKI_inhouse) %in% key_microb, ]
common_keys <- intersect(rownames(AKI_aldex_output), rownames(AKI_inhouse_key))

AKI_aldex_filt <- c(AKI_aldex_filt,common_keys)


### 라벨 표시
EnhancedVolcano(
  AKI_aldex_output,
  lab = rownames(AKI_aldex_output),
  x = 'effect',                # 효과 크기
  y = 'wi.eBH',                # p-value
  ylim = c(0, 4),
  xlab = "Effect size",
  ylab = "-log10(p-adj)",
  pCutoff = 0.5,              # p-value cutoff
  FCcutoff = 0.5,             # 효과 크기 cutoff
  title = "Differential Abundance Analysis",
  subtitle = "Control vs AKI",
  caption = "Data: Aldex2",
  pointSize = 3.5,
  col=c('gray','#00028C','gray','#FF0B55'),
  labSize = 5.5,
  cutoffLineCol = "#405D72",      # 가로/세로 기준선 색
  cutoffLineWidth = 0.8,      # 선 두께(원하면)
  cutoffLineType = "dotted",
  selectLab = AKI_aldex_filt,  # 모든 라벨 표시
  drawConnectors = TRUE,                # 연결선 추가
  widthConnectors = 0.5,                # 연결선 두께
  colConnectors = "black",               # 연결선 색상
  arrowheads = TRUE,                   # 화살표 비활성화
  lengthConnectors = unit(0.02, "npc"),  # 연결선 길이
  max.overlaps = Inf
)

ggsave("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/DA_microbiome_AKI.png", dpi=1000, dev='png', height=8, width=8, units="in", bg = "white")

#중요 두 균 표시
ggsave("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/DA_key_microbiome_AKI.png", dpi=1000, dev='png', height=8, width=8, units="in", bg = "white")

#CKD
CKD_aldex_output <- data.frame(CKD_aldex)

#필터링
cutoff_data <- CKD_aldex_output[
  CKD_aldex_output$wi.eBH < 0.5 & abs(CKD_aldex_output$effect) > 0.5, 
]

selected_labels <- rownames(cutoff_data)

#CKD_inhouse_filtered는 rownames(CKD_inhouse)가 selected_labels인 것들이다. 
CKD_inhouse_filtered <- CKD_inhouse[rownames(CKD_inhouse) %in% selected_labels, ]

## rownames(aldex_output) 에서 rownames(CKD_inhouse_filtered) 와 같은 행만 추출하기 
# 공통 행 이름 추출
common_rows <- intersect(rownames(CKD_aldex_output), rownames(CKD_inhouse_filtered))
CKD_aldex_filtered <- CKD_aldex_output[common_rows, ]
head(CKD_aldex_filtered) # 필터링된 데이터 확인

CKD_aldex_filt <- rownames(CKD_aldex_filtered)[head(order(abs(CKD_aldex_filtered$effect), decreasing=TRUE), 10)]

#중요 두 균 표시----
key_microb <- c("Escherichia","Salmonella")
CKD_inhouse_key <- CKD_inhouse[rownames(CKD_inhouse) %in% key_microb, ]
common_keys <- intersect(rownames(CKD_aldex_output), rownames(CKD_inhouse_key))

CKD_aldex_filt <- c(CKD_aldex_filt,common_keys)

### 라벨 표시
EnhancedVolcano(
  CKD_aldex_output,
  lab = rownames(CKD_aldex_output),
  x = 'effect',                # 효과 크기
  y = 'wi.eBH',                # p-value
  xlim = c(-2,2),
  ylim = c(0, 1),
  xlab = "Effect size",
  ylab = "-log10(p-adj)",
  pCutoff = 0.5,              # p-value cutoff
  FCcutoff = 0.5,             # 효과 크기 cutoff
  title = "Differential Abundance Analysis",
  subtitle = "Control vs CKD",
  caption = "Data: Aldex2",
  pointSize = 3.5,
  col=c('gray','#00028C','gray','#FF0B55'),
  labSize = 5.5,
  cutoffLineCol = "#405D72",      # 가로/세로 기준선 색
  cutoffLineWidth = 0.8,      # 선 두께(원하면)
  cutoffLineType = "dotted",
  selectLab = CKD_aldex_filt,  # 모든 라벨 표시
  drawConnectors = TRUE,                # 연결선 추가
  widthConnectors = 0.5,                # 연결선 두께
  colConnectors = "black",               # 연결선 색상
  arrowheads = TRUE,                   # 화살표 비활성화
  lengthConnectors = unit(0.02, "npc"),  # 연결선 길이
  max.overlaps = Inf
)

ggsave("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/DA_microbiome_CKD.png", dpi=1000, dev='png', height=8, width=8, units="in", bg = "white")

#중요 두 균 표시
ggsave("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/DA_key_microbiome_CKD.png", dpi=1000, dev='png', height=8, width=8, units="in", bg = "white")


#2. Microbiome Proportion----
#AKI
setwd("~/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/data")
AKI_inhouse_percentage = read.csv("genus_braken_merged", sep = "\t")
AKI_inhouse_percentage = AKI_inhouse_percentage[,-c(2,3)]
colnames(AKI_inhouse_percentage)

# 열 이름에 'bracken_num'이 포함된 열 제거
AKI_inhouse_percentage <- AKI_inhouse_percentage[, !grepl("bracken_num", colnames(AKI_inhouse_percentage))]

colnames(AKI_inhouse_percentage) <- gsub("_9606_filtered.bracken_frac", "", colnames(AKI_inhouse_percentage))

rownames(AKI_inhouse_percentage) <- AKI_inhouse_percentage$name
AKI_inhouse_percentage <- AKI_inhouse_percentage[, c("C1", "C3", "C4", "C5", "A1","A2","A3","A4","A5","A6","F4")]
AKI_inhouse_percentage$name <- NULL
head(AKI_inhouse_percentage)


# AKI_inhouse_percentage 행이 모두 0.01 미만일 경우 제거 
AKI_inhouse_percentage_filtered <- AKI_inhouse_percentage[rowSums(AKI_inhouse_percentage >= 0.01) > 0, ]

library(tidyr)
library(dplyr)
library(tibble)

plot_data <- AKI_inhouse_percentage_filtered %>%
  as.data.frame() %>%
  rownames_to_column("Microbe") %>%
  pivot_longer(-Microbe, names_to = "Sample", values_to = "Proportion")

unique(plot_data$Sample)
plot_data$Sample <- factor(plot_data$Sample, levels = c("C1", "C3", "C4", "C5", "A1", "A2", "A3", "A4", "A5", "A6", "F4"))

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
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1.6),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 17),
    legend.position = "right",
    legend.text = element_text(size = 15),
    panel.grid = element_blank()  # 격자 제거
  )


ggsave("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/Prop_microbiome_AKI.png", dpi=1000, dev='png', height=6.5, width=11, units="in", bg = "white")


#CKD
CKD_inhouse_percentage = read.csv("genus_braken_merged", sep = "\t")
CKD_inhouse_percentage = CKD_inhouse_percentage[,-c(2,3)]
colnames(CKD_inhouse_percentage)

# 열 이름에 'bracken_num'이 포함된 열 제거
CKD_inhouse_percentage <- CKD_inhouse_percentage[, !grepl("bracken_num", colnames(CKD_inhouse_percentage))]

colnames(CKD_inhouse_percentage) <- gsub("_9606_filtered.bracken_frac", "", colnames(CKD_inhouse_percentage))

rownames(CKD_inhouse_percentage) <- CKD_inhouse_percentage$name
CKD_inhouse_percentage <- CKD_inhouse_percentage[, c("C1", "C3", "C4", "C5", "M1","M2","M3","M4","M5","M6","M7","M8","M9","M10","M11","M12","F1","F2","F3")]
CKD_inhouse_percentage$name <- NULL
head(CKD_inhouse_percentage)


# CKD_inhouse_percentage 행이 모두 0.01 미만일 경우 제거 
CKD_inhouse_percentage_filtered <- CKD_inhouse_percentage[rowSums(CKD_inhouse_percentage >= 0.01) > 0, ]

library(tidyr)
library(dplyr)
library(tibble)

plot_data <- CKD_inhouse_percentage_filtered %>%
  as.data.frame() %>%
  rownames_to_column("Microbe") %>%
  pivot_longer(-Microbe, names_to = "Sample", values_to = "Proportion")
plot_data$Sample <- factor(plot_data$Sample, levels = c("C1", "C3", "C4", "C5", "M1","M2","M3","M4","M5","M6","M7","M8","M9","M10","M11","M12","F1","F2","F3"))

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
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1.6),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 17),
    legend.position = "right",
    legend.text = element_text(size = 15),
    panel.grid = element_blank()  # 격자 제거
  )

ggsave("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/Prop_microbiome_CKD.png", dpi=1000, dev='png', height=6.5, width=11, units="in", bg = "white")


#3. Microbiome Boxplot----
library(ggplot2)
library(ggpubr)

#AKI
rownames(AKI_inhouse_percentage_filtered) # 29개
AKI_condition <- factor(AKI_condition, levels = c("Control", "AKI"))
for (i in 1:29) {
  microbe_name <- rownames(AKI_inhouse_percentage_filtered)[i]
  microbe_data <- AKI_inhouse_percentage_filtered[microbe_name, ]
  plot_data_microbe <- data.frame(
    Abundance = as.numeric(microbe_data),
    Group = AKI_condition
  )
  p <- ggplot(plot_data_microbe, aes(x = Group, y = Abundance, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black", width = 0.4) +
    geom_jitter(width = 0.2, alpha = 0.8, size = 3, aes(color = Group)) +
    stat_compare_means(
      method = "t.test",
      label = "p.format",
      comparisons = list(c("Control", "AKI")),
      vjust = -0.003,
      size = 7
    ) +
    labs(
      title = paste(microbe_name),
      subtitle = "Relative Abundance",
      x = "Group"
    ) +
    scale_fill_manual(values = c("Control" = "#424242", "AKI" = "#BF4F6F")) +
    scale_color_manual(values = c("Control" = "#424242", "AKI" = "#BF4F6F")) +
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
    filename = paste0("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/Box_microbiome_AKI_", microbe_name, ".png"),
    plot = p,
    dpi = 1000,
    device = "png",
    height = 6.7,
    width = 6,
    units = "in",
    bg = "white"
  )
}


#CKD
rownames(CKD_inhouse_percentage_filtered) # 29개
CKD_condition <- factor(CKD_condition, levels = c("Control", "CKD"))
for (i in 1:3) {
  microbe_name <- rownames(CKD_inhouse_percentage_filtered)[i]
  microbe_data <- CKD_inhouse_percentage_filtered[microbe_name, ]
  plot_data_microbe <- data.frame(
    Abundance = as.numeric(microbe_data),
    Group = CKD_condition
  )
  p <- ggplot(plot_data_microbe, aes(x = Group, y = Abundance, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black", width = 0.4) +
    geom_jitter(width = 0.2, alpha = 0.8, size = 3, aes(color = Group)) +
    stat_compare_means(
      method = "t.test",
      label = "p.format",
      comparisons = list(c("Control", "CKD")),
      vjust = -0.003,
      size = 7
    ) +
    labs(
      title = paste(microbe_name),
      subtitle = "Relative Abundance",
      x = "Group"
    ) +
    scale_fill_manual(values = c("Control" = "#424242", "CKD" = "#1F77B4")) +
    scale_color_manual(values = c("Control" = "#424242", "CKD" = "#1F77B4")) +
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
    filename = paste0("/home/eunhye/raw_data1/eunhye/IgAN renal disease/DEH_IgAN_CODE/figure/Box_microbiome_CKD_", microbe_name, ".png"),
    plot = p,
    dpi = 1000,
    device = "png",
    height = 6.7,
    width = 6,
    units = "in",
    bg = "white"
  )
}


