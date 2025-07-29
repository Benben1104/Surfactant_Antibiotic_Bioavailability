library(ggplot2)
library(dplyr)
library(tidyr)
library(car)
library(RColorBrewer)

d1 <- read.csv("TET_mRNA_fold_change_1.csv")
head(d1)
d1$gene <- factor(d1$gene, level = c("marA", "acrA", "acrB", "ompF"))

d_marA <- subset(d1, gene == "marA")
residuals <- (aov(log2FC ~ group, data = d_marA))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_marA)
aov_marA <- aov(log2FC ~ group, data = d_marA)
summary(aov_marA)
TukeyHSD(aov_marA)

d_acrA <- subset(d1, gene == "acrA")
residuals <- (aov(log2FC ~ group, data = d_acrA))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_acrA)
aov_acrA <- aov(log2FC ~ group, data = d_acrA)
summary(aov_acrA)
TukeyHSD(aov_acrA)

d_acrB <- subset(d1, gene == "acrB")
residuals <- (aov(log2FC ~ group, data = d_acrB))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_acrB)
aov_acrB <- aov(log2FC ~ group, data = d_acrB)
summary(aov_acrB)
TukeyHSD(aov_acrB)

d_ompF <- subset(d1, gene == "ompF")
residuals <- (aov(log2FC ~ group, data = d_ompF))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_ompF)
aov_ompF <- aov(log2FC ~ group, data = d_ompF)
summary(aov_ompF)
TukeyHSD(aov_ompF)

ggplot(d1, aes(gene, log2FC, group = group, color = group))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  geom_point(stat = "summary", fun = "mean", alpha = 1.0, size = 6,
             position = position_dodge(width = 0.6))+
  geom_errorbar(stat = "summary",
                fun.min = function(x) mean(x)-sd(x),
                fun.max = function(x) mean(x)+sd(x),
                width = 0, size = 2.25,
                position = position_dodge(width = 0.6))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1)+
  labs(x = NULL,
       y = expression(Log[2]("FC")),
       title = expression("TET"))+
    scale_y_continuous(limits = c(-2.5, 5))+
    scale_fill_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
    scale_color_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
    theme(axis.text.x = element_text(size = 18, face = "bold.italic", angle = 0))+
    theme(axis.text.y = element_text(size = 18, face = "bold"))+
    theme(axis.title.y = element_text(size = 20, vjust = 2))+
    theme(legend.text = element_text(size = 18, face = "bold"))+
    theme(legend.title = element_blank())+
    theme(legend.position = "bottom")+
    theme(title = element_text(size = 18))

ggsave("TET_mRNA_fold_change_1.png", width = 675/90, height = 450/90, dpi = 600, units = "in")



d2 <- read.csv("TET_mRNA_fold_change_2.csv")
head(d2)
d2$gene <- factor(d2$gene, level = c("cdsA", "pgsA", "pgpA", "pssA", "psd", "pldA"))

d_cdsA <- subset(d2, gene == "cdsA")
residuals <- (aov(log2FC ~ group, data = d_cdsA))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_cdsA)
aov_cdsA <- aov(log2FC ~ group, data = d_cdsA)
summary(aov_cdsA)
TukeyHSD(aov_cdsA)

d_pgsA <- subset(d2, gene == "pgsA")
residuals <- (aov(log2FC ~ group, data = d_pgsA))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_pgsA)
aov_pgsA <- aov(log2FC ~ group, data = d_pgsA)
summary(aov_pgsA)
TukeyHSD(aov_pgsA)

d_pgpA <- subset(d2, gene == "pgpA")
residuals <- (aov(log2FC ~ group, data = d_pgpA))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_pgpA)
aov_pgpA <- aov(log2FC ~ group, data = d_pgpA)
summary(aov_pgpA)
TukeyHSD(aov_pgpA)

d_pssA <- subset(d2, gene == "pssA")
residuals <- (aov(log2FC ~ group, data = d_pssA))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_pssA)
aov_pssA <- aov(log2FC ~ group, data = d_pssA)
summary(aov_pssA)
TukeyHSD(aov_pssA)

d_psd <- subset(d2, gene == "psd")
residuals <- (aov(log2FC ~ group, data = d_psd))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_psd)
aov_psd <- aov(log2FC ~ group, data = d_psd)
summary(aov_psd)
TukeyHSD(aov_psd)

d_pldA <- subset(d2, gene == "pldA")
residuals <- (aov(log2FC ~ group, data = d_pldA))$residuals
shapiro.test(residuals)
leveneTest(log2FC ~ group, data = d_pldA)
aov_pldA <- aov(log2FC ~ group, data = d_pldA)
summary(aov_pldA)
TukeyHSD(aov_pldA)

ggplot(d2, aes(gene, log2FC, group = group, color = group))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))+
  geom_point(stat = "summary", fun = "mean", alpha = 1.0, size = 6,
             position = position_dodge(width = 0.6))+
  geom_errorbar(stat = "summary",
                fun.min = function(x) mean(x)-sd(x),
                fun.max = function(x) mean(x)+sd(x),
                width = 0, size = 2.25,
                position = position_dodge(width = 0.6))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1)+
  labs(x = NULL,
       y = expression(Log[2]("FC")),
       title = expression("TET"))+
  scale_y_continuous(limits = c(-5, 1))+
  scale_fill_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  scale_color_manual(values = c("Ctrl" = "#C0C0C0", "DTAC" = "#bd1e2f", "SDS" = "#1f4e9f", "TX-100" = "#7367BE"))+
  theme(axis.text.x = element_text(size = 18, face = "bold.italic", angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size = 20, vjust = 2))+
  theme(legend.text = element_text(size = 18, face = "bold"))+
  theme(legend.title = element_blank())+
  theme(legend.position = "bottom")+
  theme(title = element_text(size = 18))

ggsave("TET_mRNA_fold_change_2.png", width = 1000/90, height = 450/90, dpi = 600, units = "in")
