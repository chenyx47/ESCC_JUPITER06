# ---------------------------------------------------------------
## Scripts to Reproduce the Analysis Figures
## Author: Yanxing Chen
## Version: R 4.0.1
# ---------------------------------------------------------------
dir.create("/data3/cyx/JS001_ESCC/R_Project/Running-Hub")
setwd("/data3/cyx/JS001_ESCC/R_Project/Running-Hub")
{
  library(forestplot) #2.0.1
  library(survival) #3.2-13
  library(survminer) #0.4.8
  library(ggplot2) #3.3.5
  library(ggrepel) #0.8.2
  library(ComplexHeatmap) #2.10.0
  library(circlize) #0.4.12
  library(ggpubr) #0.4.0
  library(patchwork) #1.1.0
  library(ggsci) #2.9
  library(psych) #2.0.7
  library(ggcorrplot) #0.1.3
  library(gghalves) #0.1.3
  source("./Function.R")
}


ESCC_Jupiter06_Clinical_SupplementaryTable2<-
  readRDS("./SourceData/ESCC_Jupiter06_Clinical_SupplementaryTable2.RDS")

# ---------------------------------------------------------------
## Scripts to Reproduce the Association between TMBs and the benefits----
# ---------------------------------------------------------------

dir.create("./Figure/Figure2",recursive = TRUE)
ESCC_Jupiter06_Immunogenicity<-
  ESCC_Jupiter06_Clinical_SupplementaryTable2[
    c(1:21)
  ]

##Figure2a----

Result<-
  SubgroupHR_Plot(ESCC_Jupiter06_Immunogenicity,
                  which(names(ESCC_Jupiter06_Immunogenicity)=="Total-TMB"),
                  c("High","Low"))

pdf("./Figure/Figure2/tTMB_PFS_HR_Cut80(SupplementaryFigure).pdf",
    width = 10,height = 8)
print(Result$PFS)
dev.off()

pdf("./Figure/Figure2/tTMB_OS_HR_Cut80(Figure2A).pdf",
    width = 10,height = 8)
print(Result$OS)
dev.off()

##Figure2b----

Result<-
  SubgroupHR_Plot(ESCC_Jupiter06_Immunogenicity,
                  which(names(ESCC_Jupiter06_Immunogenicity)=="Clonal-TMB"),
                  c("High","Low"))

pdf("./Figure/Figure2/clonalTMB_PFS_HR_Cut70(SupplementaryFigure).pdf",
    width = 10,height = 8)
print(Result$PFS)
dev.off()

pdf("./Figure/Figure2/clonalTMB_OS_HR_Cut70(Figure2B).pdf",
    width = 10,height = 8)
print(Result$OS)
dev.off()

##Supplementary(APOBEC)----

Result<-
  SubgroupHR_Plot(ESCC_Jupiter06_Immunogenicity,
                  which(names(ESCC_Jupiter06_Immunogenicity)=="APOBECstatus"),
                  c("High","Low"))

pdf("./Figure/Figure2/APOBEC_PFS_HR_Cut70(SupplementaryFigure).pdf",
    width = 10,height = 8)
print(Result$PFS)
dev.off()

pdf("./Figure/Figure2/APOBEC_OS_HR_Cut70(SupplementaryFigure).pdf",
    width = 10,height = 8)
print(Result$OS)
dev.off()

##Figure2d----

Result<-
  SubgroupHR_Plot(ESCC_Jupiter06_Immunogenicity,
                  which(names(ESCC_Jupiter06_Immunogenicity)=="ccTMBstatus"),
                  c("High","Low"))

pdf("./Figure/Figure2/ccTMB_PFS_HR_Cut70(SupplementaryFigure).pdf",
    width = 10,height = 8)
print(Result$PFS)
dev.off()

pdf("./Figure/Figure2/ccTMB_OS_HR_Cut70(Figure2d).pdf",
    width = 10,height = 8)
print(Result$OS)
dev.off()

InputGene<-"ccTMBstatus"

{
  InputData_Select<-
    ESCC_Jupiter06_Immunogenicity[c("SUBJID","PFS","PFS.status.","OS","OS.status.",
                                    "Grouping",InputGene)]%>%unique()
  InputData_Select$Grouping<-
    factor(InputData_Select$Grouping,
           levels = c("Chemotherapy+Placebo","Chemotherapy+Toripalimab"))
  names(InputData_Select)[7]<-"Item"
  InputData_Select$Item[
    InputData_Select$Item%in%c("Low")
  ]<-"Control"
  InputData_Select$Item[
    InputData_Select$Item%in%c("High")
  ]<-"Pos"
  InputData_Select<-
    subset(InputData_Select,Item%in%c("Control","Pos"))
  InputData_Select$Item<-
    factor(InputData_Select$Item,
           levels = c("Control","Pos"))
  
  SurvFit<-
    survfit(Surv(`PFS`,`PFS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Pos"))
  PFS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  PFS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`PFS`,`PFS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Pos"))%>%
    summary()
  Plot1<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Pos"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+anti-PD-1"),
               pval = paste0("Median PFS:\n",
                             PFS_JS001_Median," Months vs ",PFS_Placebo_Median," Months"),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  
  SurvFit<-
    survfit(Surv(`PFS`,`PFS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Control"))
  PFS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  PFS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`PFS`,`PFS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Control"))%>%
    summary()
  Plot2<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Control"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+anti-PD-1"),
               pval = paste0("Median PFS:\n",
                             PFS_JS001_Median," Months vs ",PFS_Placebo_Median," Months"),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  SurvFit<-
    survfit(Surv(`OS`,`OS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Pos"))
  OS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  OS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`OS`,`OS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Pos"))%>%
    summary()
  Plot3<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Pos"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+anti-PD-1"),
               pval = paste0("Median OS:\n",
                             OS_JS001_Median," Months vs ",OS_Placebo_Median," Months"),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  SurvFit<-
    survfit(Surv(`OS`,`OS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Control"))
  OS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  OS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`OS`,`OS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Control"))%>%
    summary()
  Plot4<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Control"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+anti-PD-1"),
               pval = paste0("Median OS:\n",
                             OS_JS001_Median," Months vs ",OS_Placebo_Median," Months"),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  TotalPlot<-
    Plot1$plot+ylab("Progress-free Survival")+xlab("Time (Months)")+
    Plot2$plot+ylab("Progress-free Survival")+xlab("Time (Months)")+
    Plot3$plot+ylab("Overall Survival")+xlab("Time (Months)")+
    Plot4$plot+ylab("Overall Survival")+xlab("Time (Months)")
  print(TotalPlot)
}
pdf("./Figure/Figure2/ccTMB_Curve_Type2.pdf",
    width = 10,height = 8)
print(TotalPlot)
dev.off()

##Figure2e----

Immunogenicity_Compare<-
  Interation_Screening_SV_BV_ForKvalue(
    ESCC_Jupiter06_Immunogenicity,
    Set_Column=which(names(ESCC_Jupiter06_Immunogenicity)%in%c("Total-TMB",
                                                               "Clonal-TMB",
                                                               "ccTMBstatus",
                                                               "APOBECstatus")),
    Set_Comparison_A="High",
    Set_Comparison_B="Low"
  )
Immunogenicity_Compare$OS_Inter_Status<-
  ifelse(Immunogenicity_Compare$OS.Inter.k>0,
         "Pos","Neg")
Immunogenicity_Compare<-
  Immunogenicity_Compare[
    order(Immunogenicity_Compare$OS.Inter.k,decreasing = T),
  ]
Plot<-
  ggplot(Immunogenicity_Compare,
         aes(x=Items,y=OS.Inter.k,size=-log(OS.Inter.Pval),
             color=OS_Inter_Status))+
  geom_point()+
  geom_errorbar(data = Immunogenicity_Compare,aes(ymin = OS.Inter.k.Down, ymax=OS.Inter.k.Up), #误差条表示95%的置信区间
                width=0.1, #误差条末端短横线的宽度
                position=position_dodge(0.9), 
                color="black",
                alpha = 0.7,
                size=0.5)+
  coord_flip()+
  scale_color_manual(values=c("darkred","darkblue"))+
  scale_x_discrete(limits=Immunogenicity_Compare$Items)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size=12,color="black"),
        axis.title.x = element_text(size=15),
        legend.position = "top")+
  geom_hline(yintercept = 0,linetype="dashed")+
  ylab("Interation k (OS)");Plot
pdf("./Figure/Figure2/TMB_Comparison_OS(Figure2e).pdf",
    width = 3,height = 6)
print(Plot)
dev.off()

Immunogenicity_Compare<-
  Immunogenicity_Compare[
    order(Immunogenicity_Compare$PFS.Inter.k),
  ]
Immunogenicity_Compare$PFS_Inter_Status<-
  ifelse(Immunogenicity_Compare$PFS.Inter.k>0,
         "Pos","Neg")
Immunogenicity_Compare$PFS_Inter_Status<-
  gsub("Neg","Favorable",Immunogenicity_Compare$PFS_Inter_Status)
Plot<-
  ggplot(Immunogenicity_Compare,
         aes(x=Items,y=PFS.Inter.k,size=-log(PFS.Inter.Pval),
             color=PFS_Inter_Status))+
  geom_point()+
  geom_errorbar(data = Immunogenicity_Compare,aes(ymin = PFS.Inter.k.Down, ymax=PFS.Inter.k.Up), #误差条表示95%的置信区间
                width=0.1, #误差条末端短横线的宽度
                position=position_dodge(0.9), 
                color="black",
                alpha = 0.7,
                size=0.5)+
  coord_flip()+
  scale_color_manual(values=c("darkred","darkblue"))+
  scale_x_discrete(limits=rev(Immunogenicity_Compare$Items))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size=12,color="black"),
        axis.title.x = element_text(size=15),
        legend.position = "top")+
  geom_hline(yintercept = 0,linetype="dashed")+
  ylab("Interation k (PFS)");Plot
pdf("./Figure/Figure2/TMB_Comparison_PFS(SupplementaryFigure).pdf",
    width = 3,height = 6)
print(Plot)
dev.off()


##GC cohort for validation----

GC_IO_Clinical_TMB<-
  readRDS("SourceData/Mono_IO_Cohort/GC_IO_CombinatorialBiomarker.RDS")%>%
  .[c("WES_ID","ORR","TMB","ccTMB")]

Plot1<-
  ggplot(subset(GC_IO_Clinical_TMB,ORR!="UK"),
         aes(x=ORR,y=log10(TMB+1),fill=ORR))+
  geom_violin()+geom_boxplot(width=0.2)+
  geom_signif(comparisons = list(c("Non-Response",
                                   "Response")))+
  ylab("log10(TMB+1)")+
  scale_fill_manual(values=c("#CD9C40","#31854F"))+
  theme(axis.title.x = element_blank(),
        axis.text = element_text(color="black"),
        legend.position = "top");Plot1


Plot2<-
  ggplot(subset(GC_IO_Clinical_TMB,ORR!="UK"),
         aes(x=ORR,y=log10(ccTMB+1),fill=ORR))+
  geom_violin()+geom_boxplot(width=0.2)+
  geom_signif(comparisons = list(c("Non-Response",
                                   "Response")))+
  ylab("log10(ccTMB+1)")+
  scale_fill_manual(values=c("#CD9C40","#31854F"))+
  theme(axis.title.x = element_blank(),
        axis.text = element_text(color="black"),
        legend.position = "top");Plot2
Plot_Response<-
  Plot1+Plot2
ggsave(Plot_Response,
       filename = "./Figure/Figure2/MonoICI_GC_ccTMB.pdf",
       width = 6,height = 5)

# ---------------------------------------------------------------
## Scripts to Reproduce the Association between HLA genotype and the benefits----
# ---------------------------------------------------------------

dir.create("./Figure/Figure3",recursive = TRUE)

##Figure3a----

All_HLAI_Freq<-
  table(ESCC_Jupiter06_Immunogenicity$HLA_Remain)%>%
  as.data.frame()
Plot<-
  ggdonutchart(All_HLAI_Freq,
               "Freq",
               label = "Var1",                               
               fill = "Var1",                            
               color = "white",                                
               palette = c("#FFEBC9","#D79771","#B05B3B","#753422") 
  )
ggsave(Plot,
       filename = "./Figure/Figure3/HLAI_AllCount_Freq(Figure3a).pdf",
       width = 4,height = 4)

Germline_HLAII_Freq<-
  table(ESCC_Jupiter06_Immunogenicity$Beta_Remain)%>%
  as.data.frame()
Plot<-
  ggdonutchart(Germline_HLAII_Freq,
               "Freq",
               label = "Var1",                               
               fill = "Var1",                            
               color = "white",                                
               palette = c("#FBE1C0","#CC8F74",
                           "#A65D46","#6C3627") 
  )
ggsave(Plot,
       filename = "./Figure/Figure3/HLAclassII_Status_Freq(Figure3a).pdf",
       width = 4,height = 4)

###Supplementary Figure (Frequency of HLA heterogeneity Respectively)----

HLAA_Freq<-
  table(ESCC_Jupiter06_Immunogenicity$HLA_A_Status)%>%
  as.data.frame()%>%
  mutate(Var1=factor(Var1,levels = c("Hetero","Somatic-Homo","Germline-Homo","Unknown")))
Plot<-
  ggdonutchart(HLAA_Freq,
               "Freq",
               label = "Var1",
               fill = "Var1",
               color = "white",
               palette = c("#D8E3E7","#51C4D3","#126E82","lightgrey")
  )+theme(legend.title = element_blank(),
          axis.text.x = element_blank())
ggsave(Plot,
       filename = "./Figure/Figure3/HLAA_Alt_Freq.pdf",
       width = 4,height = 4)

HLAB_Freq<-
  table(ESCC_Jupiter06_Immunogenicity$HLA_B_Status)%>%
  as.data.frame()%>%
  mutate(Var1=factor(Var1,levels = c("Hetero","Somatic-Homo","Germline-Homo","Unknown")))
Plot<-
  ggdonutchart(HLAB_Freq,
               "Freq",
               label = "Var1",
               fill = "Var1",
               color = "white",
               palette = c("#D8E3E7","#51C4D3","#126E82","lightgrey")
  )+theme(legend.title = element_blank(),
          axis.text.x = element_blank())
ggsave(Plot,
       filename = "./Figure/Figure3/HLAB_Alt_Freq.pdf",
       width = 4,height = 4)

HLAC_Freq<-
  table(ESCC_Jupiter06_Immunogenicity$HLA_C_Status)%>%
  as.data.frame()%>%
  mutate(Var1=factor(Var1,levels = c("Hetero","Somatic-Homo","Germline-Homo","Unknown")))
Plot<-
  ggdonutchart(HLAC_Freq,
               "Freq",
               label = "Var1",
               fill = "Var1",
               color = "white",
               palette = c("#D8E3E7","#51C4D3","#126E82","lightgrey")
  )+theme(legend.title = element_blank(),
          axis.text.x = element_blank())
ggsave(Plot,
       filename = "./Figure/Figure3/HLAC_Alt_Freq.pdf",
       width = 4,height = 4)

##Figure3b----

ESCC_Jupiter06_Immunogenicity$TotalRemain<-
  ifelse(ESCC_Jupiter06_Immunogenicity$HLA_Remain>=5,
         ">=5","<5")
TotalRemain_CoxResult<-
  SubgroupHR_Plot(ESCC_Jupiter06_Immunogenicity,
                  Set_Column = which(names(ESCC_Jupiter06_Immunogenicity)=="TotalRemain"),
                  Set_Comparison_Group = c(">=5","<5"))

HLAII_CoxResult<-
  SubgroupHR_Plot(ESCC_Jupiter06_Immunogenicity,
                  Set_Column = which(names(ESCC_Jupiter06_Immunogenicity)=="Beta_Remain_Status"),
                  Set_Comparison_Group = c(">=5","<5"))

###PFS
HR_Plot_Combine_PFS<-
  rbind(TotalRemain_CoxResult$Data_PFS,
        HLAII_CoxResult$Data_PFS)

HR_Plot_Combine_PFS[c(1,4),
                    c(2:7)]<-NA
HR_Plot_Combine_PFS[c(1,4),1]<-c("Residual HLA-I alleles",
                                 "Residual HLA-IIbeta alleles")
HR_Plot_Combine_PFS[1,c(2,6,7)]<-
  c("HR(95% CI)","P-value","Freq")
PFS_Plot<-
  forestplot(as.matrix(HR_Plot_Combine_PFS[c(1,7,6)]),
             HR_Plot_Combine_PFS$HR_exact,HR_Plot_Combine_PFS$Low_CI,HR_Plot_Combine_PFS$Up_CI,
             xlog = TRUE,
             zero = 1, 
             xticks = c(0.20,1,5),
             colgap = unit(10,"mm"),
             graphwidth=unit(50,"mm"),
             lineheight = unit(1,"cm"),
             graph.pos = 3,
             col = fpColors(lines="#0177BF", box="#0177BF",zero = "BLACK"),
             boxsize = 0.4,
             ci.vertices = F,
             lty.ci = 1,
             lwd.ci = 1.5,
             lwd.zero=0.4,
             lty.zero=3,
             txt_gp=fpTxtGp(label = gpar(cex=0.8),
                            ticks = gpar(cex=0.8)), is.summary = c(rep(FALSE,19)),
             fn.ci_norm = fpDrawNormalCI,
             hrzl_lines=list("1" = gpar(lwd=1, col="black")));PFS_Plot
pdf("./Figure/Figure3/HLAhetero_HRforrest_PFS(Figure3b).pdf",
    width = 10,height = 10)
print(PFS_Plot)
dev.off()

###OS
HR_Plot_Combine_OS<-
  rbind(TotalRemain_CoxResult$Data_OS,
        HLAII_CoxResult$Data_OS)

HR_Plot_Combine_OS[c(1,4),
                   c(2:7)]<-NA
HR_Plot_Combine_OS[c(1,4),1]<-c("Residual HLA-I alleles",
                                "Residual HLA-IIbeta alleles")
HR_Plot_Combine_OS[1,c(2,6,7)]<-
  c("HR(95% CI)","P-value","Freq")
OS_Plot<-
  forestplot(as.matrix(HR_Plot_Combine_OS[c(1,7,6)]),
             HR_Plot_Combine_OS$HR_exact,HR_Plot_Combine_OS$Low_CI,HR_Plot_Combine_OS$Up_CI,
             xlog = TRUE,
             zero = 1, 
             xticks = c(0.2,1,5),
             colgap = unit(10,"mm"),
             graphwidth=unit(50,"mm"),
             lineheight = unit(1,"cm"),
             graph.pos = 3,
             col = fpColors(lines="#0177BF", box="#0177BF",zero = "BLACK"),
             boxsize = 0.4,
             ci.vertices = F,
             lty.ci = 1,
             lwd.ci = 1.5,
             lwd.zero=0.4,
             lty.zero=3,
             txt_gp=fpTxtGp(label = gpar(cex=0.8),
                            ticks = gpar(cex=0.8)), is.summary = c(rep(FALSE,19)),
             fn.ci_norm = fpDrawNormalCI,
             hrzl_lines=list("1" = gpar(lwd=1, col="black")));OS_Plot
pdf("./Figure/Figure3/HLAhetero_HRforrest_OS(Figure3b).pdf",
    width = 10,height = 10)
print(OS_Plot)
dev.off()

###Supplementary Figure (Germline+Somatic Respectively)----

ESCC_Jupiter06_Immunogenicity_HLA_A_HR<-
  SubgroupHR_Plot(ESCC_Jupiter06_Immunogenicity,
                  Set_Column = which(names(ESCC_Jupiter06_Immunogenicity)=="HLA_A_Status"),
                  Set_Comparison_Group = c("Hetero","Germline-Homo","Somatic-Homo"))
pdf("./Figure/ImmunogenicFeature/CI_Score_OS_HR_Binary(Figurec).pdf",
    width = 15,height = 5)
print(ESCC_Jupiter06_Immunogenicity_HLA_A_HR$OS)
dev.off()
pdf("./Figure/ImmunogenicFeature/CI_Score_PFS_HR_Binary(Figured).pdf",
    width = 15,height = 5)
print(ESCC_Jupiter06_Immunogenicity_HLA_A_HR$PFS)
dev.off()

ESCC_Jupiter06_Immunogenicity_HLA_B_HR<-
  SubgroupHR_Plot(ESCC_Jupiter06_Immunogenicity,
                  Set_Column = which(names(ESCC_Jupiter06_Immunogenicity)=="HLA_B_Status"),
                  Set_Comparison_Group = c("Hetero","Germline-Homo","Somatic-Homo"))
pdf("./Figure/ImmunogenicFeature/CI_Score_OS_HR_Binary(Figurec).pdf",
    width = 15,height = 5)
print(ESCC_Jupiter06_Immunogenicity_HLA_B_HR$OS)
dev.off() 
pdf("./Figure/ImmunogenicFeature/CI_Score_PFS_HR_Binary(Figured).pdf",
    width = 15,height = 5)
print(ESCC_Jupiter06_Immunogenicity_HLA_B_HR$PFS)
dev.off()

ESCC_Jupiter06_Immunogenicity_HLA_C_HR<-
  SubgroupHR_Plot(ESCC_Jupiter06_Immunogenicity,
                  Set_Column = which(names(ESCC_Jupiter06_Immunogenicity)=="HLA_C_Status"),
                  Set_Comparison_Group = c("Hetero","Germline-Homo","Somatic-Homo"))
pdf("./Figure/ImmunogenicFeature/CI_Score_OS_HR_Binary(Figurec).pdf",
    width = 15,height = 5)
print(ESCC_Jupiter06_Immunogenicity_HLA_C_HR$OS)
dev.off() 
pdf("./Figure/ImmunogenicFeature/CI_Score_PFS_HR_Binary(Figured).pdf",
    width = 15,height = 5)
print(ESCC_Jupiter06_Immunogenicity_HLA_C_HR$PFS)
dev.off()

ESCC_Jupiter06_HLAI_Binary<-
  apply(ESCC_Jupiter06_Immunogenicity[which(names(ESCC_Jupiter06_Immunogenicity)%in%c("HLA_A_Status",
                                                                                      "HLA_B_Status",
                                                                                      "HLA_C_Status"))],2,
        function(x){
          recode(x,
                 "Germline-Homo"="Homo",
                 "Somatic-Homo"="Homo",
                 "Hetero"="Hetero")
        })%>%
  cbind(ESCC_Jupiter06_Immunogenicity[c(1:7)],
        .)

HLAIA_Individual_ForrestPlot<-
  SubgroupHR_Plot(ESCC_Jupiter06_HLAI_Binary,
                  Set_Column = which(names(ESCC_Jupiter06_HLAI_Binary)%in%c("HLA_A_Status")),
                  Set_Comparison_Group = c("Hetero","Homo"))
HLAIB_Individual_ForrestPlot<-
  SubgroupHR_Plot(ESCC_Jupiter06_HLAI_Binary,
                  Set_Column = which(names(ESCC_Jupiter06_HLAI_Binary)%in%c("HLA_B_Status")),
                  Set_Comparison_Group = c("Hetero","Homo"))
HLAIC_Individual_ForrestPlot<-
  SubgroupHR_Plot(ESCC_Jupiter06_HLAI_Binary,
                  Set_Column = which(names(ESCC_Jupiter06_HLAI_Binary)%in%c("HLA_C_Status")),
                  Set_Comparison_Group = c("Hetero","Homo"))

HR_Plot_Combine<-
  rbind(ESCC_Jupiter06_Immunogenicity_HLA_A_HR$Data_OS,
        HLAIA_Individual_ForrestPlot$Data_OS[3,],
        ESCC_Jupiter06_Immunogenicity_HLA_B_HR$Data_OS,
        HLAIB_Individual_ForrestPlot$Data_OS[3,],
        ESCC_Jupiter06_Immunogenicity_HLA_C_HR$Data_OS,
        HLAIC_Individual_ForrestPlot$Data_OS[3,])

HR_Plot_Combine[c(1,6,11),
                c(2:7)]<-NA
HR_Plot_Combine[c(1,6,11),1]<-c("HLA-A",
                                "HLA-B",
                                "HLA-C")
HR_Plot_Combine[1,c(2,6,7)]<-
  c("HR(95% CI)","P-value","Freq")
OS_Plot<-
  forestplot(as.matrix(HR_Plot_Combine[c(1,7,6)]),
             HR_Plot_Combine$HR_exact,HR_Plot_Combine$Low_CI,HR_Plot_Combine$Up_CI,
             xlog = TRUE,
             zero = 1, 
             xticks = c(0.2,1,5),
             colgap = unit(10,"mm"),
             graphwidth=unit(50,"mm"),
             lineheight = unit(1,"cm"),
             graph.pos = 3,
             col = fpColors(lines="#0177BF", box="#0177BF",zero = "BLACK"),
             boxsize = 0.4,
             ci.vertices = F,
             lty.ci = 1,
             lwd.ci = 1.5,
             lwd.zero=0.4,
             lty.zero=3,
             txt_gp=fpTxtGp(label = gpar(cex=0.8),
                            ticks = gpar(cex=0.8)), is.summary = c(rep(FALSE,19)),
             fn.ci_norm = fpDrawNormalCI,
             hrzl_lines=list("1" = gpar(lwd=1, col="black")));OS_Plot
pdf("./Figure/Figure3/HLAI_HRforrest_OS(Supplementary).pdf",
    width = 10,height = 10)
print(OS_Plot)
dev.off()

HR_Plot_Combine<-
  rbind(ESCC_Jupiter06_Immunogenicity_HLA_A_HR$Data_PFS,
        HLAIA_Individual_ForrestPlot$Data_PFS[3,],
        ESCC_Jupiter06_Immunogenicity_HLA_B_HR$Data_PFS,
        HLAIB_Individual_ForrestPlot$Data_PFS[3,],
        ESCC_Jupiter06_Immunogenicity_HLA_C_HR$Data_PFS,
        HLAIC_Individual_ForrestPlot$Data_PFS[3,])

HR_Plot_Combine[c(1,6,11),
                c(2:7)]<-NA
HR_Plot_Combine[c(1,6,11),1]<-c("HLA-A",
                                "HLA-B",
                                "HLA-C")
HR_Plot_Combine[1,c(2,6,7)]<-
  c("HR(95% CI)","P-value","Freq")
PFS_Plot<-
  forestplot(as.matrix(HR_Plot_Combine[c(1,7,6)]),
             HR_Plot_Combine$HR_exact,HR_Plot_Combine$Low_CI,HR_Plot_Combine$Up_CI,
             xlog = TRUE,
             zero = 1, 
             xticks = c(0.2,1,5),
             colgap = unit(10,"mm"),
             graphwidth=unit(50,"mm"),
             lineheight = unit(1,"cm"),
             graph.pos = 3,
             col = fpColors(lines="#0177BF", box="#0177BF",zero = "BLACK"),
             boxsize = 0.4,
             ci.vertices = F,
             lty.ci = 1,
             lwd.ci = 1.5,
             lwd.zero=0.4,
             lty.zero=3,
             txt_gp=fpTxtGp(label = gpar(cex=0.8),
                            ticks = gpar(cex=0.8)), is.summary = c(rep(FALSE,19)),
             fn.ci_norm = fpDrawNormalCI,
             hrzl_lines=list("1" = gpar(lwd=1, col="black")));PFS_Plot
pdf("./Figure/Figure3/HLAI_HRforrest_PFS(Supplementary).pdf",
    width = 10,height = 10)
print(PFS_Plot)
dev.off()

##Figure3c----

Kvalue_Result_2<-
  Interation_Screening_SV_BV_ForKvalue(ESCC_Jupiter06_Immunogenicity,
                                       Set_Column=which(names(ESCC_Jupiter06_Immunogenicity)%in%c("HLA_A_Status",
                                                                                                  "HLA_B_Status",
                                                                                                  "HLA_C_Status")),
                                       Set_Comparison_A="Hetero",
                                       Set_Comparison_B=c("Somatic-Homo","Germline-Homo"))

Kvalue_Result<-
  rbind(Kvalue_Result_2)

Kvalue_Result_ForPlot<-
  cbind(reshape2::melt(Kvalue_Result[c(1,3,7)],
                       id="Items")%>%
          mutate(variable=recode(variable,"OS.Inter.k"="OS","PFS.Inter.k"="PFS")),
        reshape2::melt(Kvalue_Result[c(1,4,8)],
                       id="Items")%>%.[c(3)],
        reshape2::melt(Kvalue_Result[c(1,5,9)],
                       id="Items")%>%.[c(3)],
        reshape2::melt(Kvalue_Result[c(1,6,10)],
                       id="Items")%>%.[c(3)])
names(Kvalue_Result_ForPlot)[c(3:6)]<-
  c("Inter-K","Inter-P","Upper95","Lower95")

Kvalue_Result_ForPlot$Items<-
  c("HLA-A Status",
    "HLA-B Status",
    "HLA-C Status")
Kvalue_Result_ForPlot$Items<-
  gsub(" Status","",Kvalue_Result_ForPlot$Items)

Plot<-
  ggplot(Kvalue_Result_ForPlot, aes(x=Items, y=`Inter-K`,fill=variable)) +
  # coord_flip()+
  # geom_violin(trim=FALSE,color="white") + #绘制小提琴图
  geom_point(data = Kvalue_Result_ForPlot,aes(x=Items, y=`Inter-K`,color=variable,size=-(log10(`Inter-P`))),
             pch=19,position=position_dodge(0.9))+ #绘制均值为点图
  geom_errorbar(data = Kvalue_Result_ForPlot,aes(ymin = Lower95, ymax=Upper95,color=variable), #误差条表示95%的置信区间
                width=0.1, #误差条末端短横线的宽度
                position=position_dodge(0.9), 
                alpha = 0.7,
                size=0.5)+
  # scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  scale_x_discrete(limits=c("HLA-B",
                            "HLA-C",
                            "HLA-A"))+
  scale_color_manual(values=c("#798777","#7D5A50"))+
  theme(axis.text.x=element_text(colour="black",size=20), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(colour="black",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(colour="black",  #设置图例的子标题的字体属性
                                 size=16),
        legend.title=element_text(colour="black", #设置图例的总标题的字体属性
                                  size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        legend.position = "top")+  #不显示网格线
  ylab("Inter-K")+xlab("")+
  geom_hline(yintercept = 0,linetype="dashed");Plot
ggsave(Plot,
       filename = "./Figure/Figure3/HLA_AlleleHetero_Kvalue_Plot(Figure3c).pdf",
       width = 10,height = 4)

##Figure3d----

ESCC_JS001_HLAI_SuperType_Inter<-
  readRDS("./SourceData/HLAsupertype/ESCC_JS001_HLAI_SuperType_Inter_PFS.RDS")
Plot<-
  ggplot(ESCC_JS001_HLAI_SuperType_Inter,
         aes(x=Items,y=PFS.Inter.k,size=-log(PFS.Inter.Pval),
             color=PFS_Inter_Status))+
  geom_point()+
  geom_errorbar(data = ESCC_JS001_HLAI_SuperType_Inter,aes(ymin = PFS.Inter.k.Down, ymax=PFS.Inter.k.Up), #误差条表示95%的置信区间
                width=0.1, #误差条末端短横线的宽度
                position=position_dodge(0.9), 
                color="black",
                alpha = 0.7,
                size=0.5)+
  # coord_flip()+
  scale_color_manual(values=c("darkblue","darkred"))+
  scale_x_discrete(limits=ESCC_JS001_HLAI_SuperType_Inter$Items%>%rev())+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size=12,color="black"),
        axis.title.y = element_text(size=15),
        legend.position = "top")+
  geom_hline(yintercept = 0,linetype="dashed")+
  ylab("Interation k (PFS)");Plot
pdf("./Figure/Figure3/HLAI_Supertype_Compare_PFS(Figure3d).pdf",
    width = 6,height = 3)
print(Plot)
dev.off()

##Figure3e----

ESCC_JS001_HLAI_SuperType_Inter<-
  readRDS("./SourceData/HLAsupertype/ESCC_JS001_HLAI_SuperType_Inter_OS.RDS")
Plot<-
  ggplot(ESCC_JS001_HLAI_SuperType_Inter,
         aes(x=Items,y=OS.Inter.k,size=-log(OS.Inter.Pval),
             color=OS_Inter_Status))+
  geom_point()+
  geom_errorbar(data = ESCC_JS001_HLAI_SuperType_Inter,aes(ymin = OS.Inter.k.Down, ymax=OS.Inter.k.Up), #误差条表示95%的置信区间
                width=0.1, #误差条末端短横线的宽度
                position=position_dodge(0.9), 
                color="black",
                alpha = 0.7,
                size=0.5)+
  # coord_flip()+
  scale_color_manual(values=c("darkblue","darkred"))+
  scale_x_discrete(limits=ESCC_JS001_HLAI_SuperType_Inter$Items%>%rev())+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size=12,color="black"),
        axis.title.y = element_text(size=15),
        legend.position = "top")+
  geom_hline(yintercept = 0,linetype="dashed")+
  ylab("Interation k (OS)");Plot
pdf("./Figure/Figure3/HLAI_Supertype_Compare_OS(Figure3e).pdf",
    width = 6,height = 3)
print(Plot)
dev.off()


##Figure3f----

B62_Detail_Freq<-
  table(ESCC_Jupiter06_Immunogenicity$B62_Detailed)%>%
  as.data.frame()

Plot<-
  ggdonutchart(B62_Detail_Freq,
               "Freq",
               label = "Var1",
               fill = "Var1",
               color = "white",
               palette = c("#F37030","#BA272D","#248C9A")
  )+theme(axis.text = element_blank(),
          legend.title = element_blank())
ggsave(Plot,
       filename = "./Figure/Figure3/B62supertype_Frequency.pdf",
       width = 4,height = 4)

##Figure3g----

InputGene<-"B62_Detailed"
{
  InputData_Select<-
    ESCC_Jupiter06_Immunogenicity[c("SUBJID","PFS","PFS.status.","OS","OS.status.",
                                    "Grouping",InputGene)]%>%unique()
  
  InputData_Select$Grouping<-
    factor(InputData_Select$Grouping,
           levels = c("Chemotherapy+Placebo","Chemotherapy+Toripalimab"))
  names(InputData_Select)[7]<-"Item"
  # InputData_Select$Item[
  #   InputData_Select$Item%in%c("Low")
  # ]<-"Control"
  # InputData_Select$Item[
  #   InputData_Select$Item%in%c("High")
  # ]<-"Pos"
  # InputData_Select<-
  #   subset(InputData_Select,Item%in%c("Control","Pos"))
  # InputData_Select$Item<-
  #   factor(InputData_Select$Item,
  #          levels = c("Control","Pos"))
  
  SurvFit<-
    survfit(Surv(`OS`,`OS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="nonB62"))
  OS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  OS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`OS`,`OS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="nonB62"))%>%
    summary()
  Plot1<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="nonB62"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median OS:\n",
                             OS_JS001_Median," Months vs ",OS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  
  SurvFit<-
    survfit(Surv(`OS`,`OS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="B62-B1501"))
  OS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  OS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`OS`,`OS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="B62-B1501"))%>%
    summary()
  Plot2<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="B62-B1501"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median OS:\n",
                             OS_JS001_Median," Months vs ",OS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  SurvFit<-
    survfit(Surv(`OS`,`OS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="B62-nonB1501"))
  OS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  OS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`OS`,`OS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="B62-nonB1501"))%>%
    summary()
  Plot3<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="B62-nonB1501"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median OS:\n",
                             OS_JS001_Median," Months vs ",OS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  
  TotalPlot<-
    Plot1$plot+ylab("Overall Survival")+xlab("Time (Months)")+
    Plot2$plot+ylab("Overall Survival")+xlab("Time (Months)")+
    Plot3$plot+ylab("Overall Survival")+xlab("Time (Months)")
  print(TotalPlot)
}
pdf("./Figure/Figure3/B62_Detail_OScurve(Figure3g).pdf",
    width = 18,height = 5)
print(TotalPlot)
dev.off()

{
  InputData_Select<-
    ESCC_Jupiter06_Immunogenicity[c("SUBJID","PFS","PFS.status.","OS","OS.status.",
                                    "Grouping",InputGene)]%>%unique()
  
  InputData_Select$Grouping<-
    factor(InputData_Select$Grouping,
           levels = c("Chemotherapy+Placebo","Chemotherapy+Toripalimab"))
  names(InputData_Select)[7]<-"Item"
  # InputData_Select$Item[
  #   InputData_Select$Item%in%c("Low")
  # ]<-"Control"
  # InputData_Select$Item[
  #   InputData_Select$Item%in%c("High")
  # ]<-"Pos"
  # InputData_Select<-
  #   subset(InputData_Select,Item%in%c("Control","Pos"))
  # InputData_Select$Item<-
  #   factor(InputData_Select$Item,
  #          levels = c("Control","Pos"))
  
  SurvFit<-
    survfit(Surv(`PFS`,`PFS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="nonB62"))
  PFS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  PFS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`PFS`,`PFS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="nonB62"))%>%
    summary()
  Plot1<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="nonB62"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median PFS:\n",
                             PFS_JS001_Median," Months vs ",PFS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  
  SurvFit<-
    survfit(Surv(`PFS`,`PFS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="B62-B1501"))
  PFS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  PFS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`PFS`,`PFS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="B62-B1501"))%>%
    summary()
  Plot2<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="B62-B1501"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median PFS:\n",
                             PFS_JS001_Median," Months vs ",PFS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  SurvFit<-
    survfit(Surv(`PFS`,`PFS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="B62-nonB1501"))
  PFS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  PFS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`PFS`,`PFS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="B62-nonB1501"))%>%
    summary()
  Plot3<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="B62-nonB1501"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median PFS:\n",
                             PFS_JS001_Median," Months vs ",PFS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  
  TotalPlot<-
    Plot1$plot+ylab("Progression-free Survival")+xlab("Time (Months)")+
    Plot2$plot+ylab("Progression-free Survival")+xlab("Time (Months)")+
    Plot3$plot+ylab("Progression-free Survival")+xlab("Time (Months)")
  print(TotalPlot)
}
pdf("./Figure/Figure3/B62_Detail_PFScurve(SupplementaryFigure).pdf",
    width = 18,height = 5)
print(TotalPlot)
dev.off()

HLAIA_Individual_ForrestPlot<-
  SubgroupHR_Plot(InputData =  ESCC_Jupiter06_Immunogenicity,
                  Set_Column = which(names(ESCC_Jupiter06_Immunogenicity)=="B62_Detailed"),
                  Set_Comparison_Group = c("nonB62","B62-nonB1501","B62-B1501"))

pdf("./Figure/Figure3/B62_Detail_OS_Forest(Supplementary).pdf",
    width = 10,height = 10)
print(HLAIA_Individual_ForrestPlot$OS)
dev.off()

pdf("./Figure/Figure3/B62_Detail_PFS_Forest(Supplementary).pdf",
    width = 10,height = 10)
print(HLAIA_Individual_ForrestPlot$PFS)
dev.off()

###Race B62

HLAB_Supertype_Freq_Select_Freq<-
  readRDS("./SourceData/B62Freq/HLAB_Supertype_Freq_Race.RDS")
Plot<-
  ggplot(HLAB_Supertype_Freq_Select_Freq,
         aes(x=Race,y=value,fill=variable))+
  geom_bar(stat = "identity",
           position = "stack",
           width = 0.7)+
  scale_fill_manual(values = c("#F37030","#BA272D"))+
  theme_classic()+
  ylab("Proportion")+
  theme(axis.text = element_text(color="black",size = 12),
        legend.title = element_blank(),
        legend.position = "top",
        axis.text.x = element_text(angle = 30,
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size=12))+
  scale_x_discrete(limits=unique(HLAB_Supertype_Freq_Select_Freq$Race));Plot

pdf("./Figure/Figure3/HLAB_Supertype_RaceFreq(Supplementary).pdf",
    width = 8,height = 6)
print(Plot)
dev.off()

###B62-Immune from TCGA

# TCGA_HLAgenotype<-
#   openxlsx::read.xlsx("TCGA_HLA_benchmark_20200810.xlsx")
# HLAB_Supertype<-
#   read.delim("HLABsupertype_2008_reformat.txt")
# 
# TCGA_HLAgenotype_Split<-
#   limma::strsplit2(TCGA_HLAgenotype$HLA_ClassI,"\\|")
# TCGA_HLAgenotype_Split<-
#   cbind(TCGA_HLAgenotype,
#         TCGA_HLAgenotype_Split)
# TCGA_HLAgenotype_Split[-1]<-
#   apply(TCGA_HLAgenotype_Split[-1],2,function(x){gsub("\\:","",x)})
# TCGA_HLAgenotype_Split$Supertype<-
#   apply(TCGA_HLAgenotype_Split[c(20:28)],1,
#         function(x){y=ifelse(any(x%in%HLAB_Supertype$Allele[HLAB_Supertype$Supertype=="B62"]),
#                              ifelse(any(x=="B*1501"),"B62-B1501","B62-nonB1501"),
#                              "Non-B62")})
# 
# TCGA_ImmuneLandscape_Score<-
#   read.delim("/data2/Public/TCGAxena/TCGA_Immune/TCGA_ImmuneLandscape_Score.txt")
# TCGA_ImmuneLandscape_Score<-
#   TCGA_ImmuneLandscape_Score[c(1,2,5,13,29,55)]
# TCGA_ImmuneLandscape_Score$T.Cells.CD8<-
#   TCGA_ImmuneLandscape_Score$T.Cells.CD8*
#   TCGA_ImmuneLandscape_Score$Leukocyte.Fraction
# TCGA_HLAgenotype_Split<-
#   merge(TCGA_HLAgenotype_Split,
#         TCGA_ImmuneLandscape_Score,
#         by.x="Sample_Barcode",
#         by.y="TCGA.Participant.Barcode")
# saveRDS(TCGA_HLAgenotype_Split,
#         "TCGA_HLAgenotype_Split_Immune.RDS")

TCGA_HLAgenotype_Split<-
  readRDS("SourceData/B62Freq/TCGA_HLAgenotype_Split_Immune.RDS")

p1<-ggplot(subset(TCGA_HLAgenotype_Split,
                  Supertype%in%c("Non-B62","B62-nonB1501")),
           aes(x=1,y = IFN.gamma.Response,fill=Supertype))+
  geom_split_violin(trim = T,colour=NA)+
  geom_boxplot(color="white",width=0.25)+
  scale_fill_manual(values = c("#D48852","#3D5057"))+
  theme_classic()+
  # mytheme+
  ylab("IFN-gamma Response")+xlab("Type")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top",
        legend.title = element_blank());p1
p2<-ggplot(subset(TCGA_HLAgenotype_Split,
                  Supertype%in%c("B62-B1501","B62-nonB1501"))%>%
             mutate(Supertype=factor(Supertype,levels = c("B62-nonB1501","B62-B1501"))),
           aes(x=1,y = IFN.gamma.Response,fill=Supertype))+
  geom_split_violin(trim = T,colour=NA)+
  geom_boxplot(color="white",width=0.25)+
  scale_fill_manual(values = c("#D48852","#4E8F7F"))+
  theme_classic()+
  # mytheme+
  ylab("IFN-gamma Response")+xlab("Type")+
  theme(axis.text = element_text(size = 12,color = "black"),
        axis.title = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top",
        legend.title = element_blank());p2
Plot_IFN<-p1+p2
pdf("Figure/Figure3/B62_IFNg_TCGA.pdf",
    width = 6,height = 6)
print(Plot_IFN)
dev.off()

p1<-ggplot(subset(TCGA_HLAgenotype_Split,
                  Supertype%in%c("Non-B62","B62-nonB1501")),
           aes(x=1,y = T.Cells.CD8,fill=Supertype))+
  geom_split_violin(trim = T,colour=NA)+
  geom_boxplot(color="white",width=0.25)+
  scale_fill_manual(values = c("#D48852","#3D5057"))+
  theme_classic()+
  # mytheme+
  ylab("T.Cells.CD8")+xlab("Type")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top",
        legend.title = element_blank());p1
p2<-ggplot(subset(TCGA_HLAgenotype_Split,
                  Supertype%in%c("B62-B1501","B62-nonB1501"))%>%
             mutate(Supertype=factor(Supertype,levels = c("B62-nonB1501","B62-B1501"))),
           aes(x=1,y = T.Cells.CD8,fill=Supertype))+
  geom_split_violin(trim = T,colour=NA)+
  geom_boxplot(color="white",width=0.25)+
  scale_fill_manual(values = c("#D48852","#4E8F7F"))+
  theme_classic()+
  # mytheme+
  ylab("T.Cells.CD8")+xlab("Type")+
  theme(axis.text = element_text(size = 12,color = "black"),
        axis.title = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top",
        legend.title = element_blank());p2
Plot_CD8<-p1+p2
pdf("Figure/Figure3/B62_CD8_TCGA.pdf",
    width = 6,height = 6)
print(Plot_CD8)
dev.off()

p1<-ggplot(subset(TCGA_HLAgenotype_Split,
                  Supertype%in%c("Non-B62","B62-nonB1501")),
           aes(x=1,y = Leukocyte.Fraction,fill=Supertype))+
  geom_split_violin(trim = T,colour=NA)+
  geom_boxplot(color="white",width=0.25)+
  scale_fill_manual(values = c("#D48852","#3D5057"))+
  theme_classic()+
  # mytheme+
  ylab("Leukocyte.Fraction")+xlab("Type")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top",
        legend.title = element_blank());p1
p2<-ggplot(subset(TCGA_HLAgenotype_Split,
                  Supertype%in%c("B62-B1501","B62-nonB1501"))%>%
             mutate(Supertype=factor(Supertype,levels = c("B62-nonB1501","B62-B1501"))),
           aes(x=1,y = Leukocyte.Fraction,fill=Supertype))+
  geom_split_violin(trim = T,colour=NA)+
  geom_boxplot(color="white",width=0.25)+
  scale_fill_manual(values = c("#D48852","#4E8F7F"))+
  theme_classic()+
  # mytheme+
  ylab("Leukocyte.Fraction")+xlab("Type")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top",
        legend.title = element_blank());p2
Plot_Leu<-p1+p2
pdf("Figure/Figure3/B62_Leukocyte_TCGA.pdf",
    width = 6,height = 6)
print(Plot_Leu)
dev.off()

#  ---------------------------------------------------------------
#Immunogenicity Feature-based classification----
#  ---------------------------------------------------------------

dir.create("./Figure/ImmunogenicFeature/")
ESCC_Jupiter06_Immunogenicity$HLA_B_Status[
  ESCC_Jupiter06_Immunogenicity$HLA_B_Status=="Unknown"
]<-NA
ESCC_Jupiter06_Immunogenicity$ccTMB_Status_Binary<-
  ifelse(ESCC_Jupiter06_Immunogenicity$ccTMBstatus%in%c("High"),
         1,0)
ESCC_Jupiter06_Immunogenicity$HLA_B_Status_Binary<-
  ifelse(ESCC_Jupiter06_Immunogenicity$HLA_B_Status%in%c("Hetero"),
         1,0)
ESCC_Jupiter06_Immunogenicity$Beta_Remain_Status_Binary<-
  ifelse(ESCC_Jupiter06_Immunogenicity$Beta_Remain_Status==">=5",
         1,0)
ESCC_Jupiter06_Immunogenicity$B62_Detailed_Binary<-
  ifelse(ESCC_Jupiter06_Immunogenicity$B62_Detailed=="B62-nonB1501",
         1,0)
ESCC_Jupiter06_Immunogenicity$HLAdiversity_Binary<-
  ifelse(ESCC_Jupiter06_Immunogenicity$HLA_B_Status_Binary==1|
           ESCC_Jupiter06_Immunogenicity$Beta_Remain_Status_Binary==1,
         1,0)
ESCC_Jupiter06_Immunogenicity$RiskCombine<-
  apply(ESCC_Jupiter06_Immunogenicity[which(names(ESCC_Jupiter06_Immunogenicity)%in%c("HLAdiversity_Binary",
                                                                                      "B62_Detailed_Binary",
                                                                                      "ccTMB_Status_Binary"))],
        1,sum)

ESCC_Jupiter06_Immunogenicity$RiskCombine_Status<-
  ifelse(ESCC_Jupiter06_Immunogenicity$RiskCombine>1,
         "immune feature-favorable","immune feature-unfavorable")

##Immunogenicity Feature-based classification----

BarPlot_Data<-
  table(ESCC_Jupiter06_Immunogenicity$RiskCombine)%>%
  as.data.frame()%>%
  mutate(Status=ifelse(as.numeric(as.character(.[,1]))>1,
                       "Immune features-favorable",
                       "Immune features-unfavorable"))
Plot<-
  ggplot(BarPlot_Data,
         aes(x=Var1, y=Freq, fill = Status))+
  geom_bar(stat = "identity",width = 0.5)+
  scale_fill_npg()+
  theme_bw()+
  theme(axis.text = element_text(size=15,color="black"),
        axis.title = element_text(size = 15),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 15),
        legend.title = element_blank())+
  ylab("Patients")+
  xlab("Count of favorable immune features");Plot
pdf("./Figure/ImmunogenicFeature/BarPlot(Figurea).pdf",
    width = 5,height = 5)
print(Plot)
dev.off()

##Immunogenicity Feature-based classification Frequency----

IM_Frequency<-
  table(ESCC_Jupiter06_Immunogenicity$RiskCombine_Status)%>%
  as.data.frame()
Plot<-
  ggdonutchart(IM_Frequency,
               "Freq",
               label = "Var1",
               fill = "Var1",
               color = "white",
               palette = c("#E64B35","#4DBBD5")
  )
ggsave(Plot,
       filename = "./Figure/ImmunogenicFeature/Proportion_Circle(Figureb).pdf",
       width = 5,height = 5)

##Immunogenicity Feature-based classification and Benefits----

HRresult<-
  SubgroupHR_Plot(InputData=ESCC_Jupiter06_Immunogenicity,
                  Set_Column = which(names(ESCC_Jupiter06_Immunogenicity)=="RiskCombine_Status"),
                  Set_Comparison_Group = c("immune feature-favorable","immune feature-unfavorable"))
pdf("./Figure/ImmunogenicFeature/CI_Score_OS_HR_Binary(Figurec).pdf",
    width = 15,height = 5)
print(HRresult$OS)
dev.off()

pdf("./Figure/ImmunogenicFeature/CI_Score_PFS_HR_Binary(Figured).pdf",
    width = 15,height = 5)
print(HRresult$PFS)
dev.off()

InputGene<-"RiskCombine_Status"
{
  InputData_Select<-
    ESCC_Jupiter06_Immunogenicity[
      c("SUBJID","PFS","PFS.status.","OS","OS.status.",
        "Grouping",InputGene)]%>%unique()
  
  InputData_Select$Grouping<-
    factor(InputData_Select$Grouping,
           levels = c("Chemotherapy+Placebo","Chemotherapy+Toripalimab"))
  names(InputData_Select)[7]<-"Item"
  InputData_Select$Item[
    InputData_Select$Item%in%c("immune feature-favorable")
  ]<-"Control"
  InputData_Select$Item[
    InputData_Select$Item%in%c("immune feature-unfavorable")
  ]<-"Pos"
  InputData_Select<-
    subset(InputData_Select,Item%in%c("Control","Pos"))
  InputData_Select$Item<-
    factor(InputData_Select$Item,
           levels = c("Control","Pos"))
  
  SurvFit<-
    survfit(Surv(`PFS`,`PFS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Pos"))
  PFS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  PFS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`PFS`,`PFS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Pos"))%>%
    summary()
  Plot1<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Pos"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median PFS:\n",
                             PFS_JS001_Median," Months vs ",PFS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  
  SurvFit<-
    survfit(Surv(`PFS`,`PFS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Control"))
  PFS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  PFS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`PFS`,`PFS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Control"))%>%
    summary()
  Plot2<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Control"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median PFS:\n",
                             PFS_JS001_Median," Months vs ",PFS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  SurvFit<-
    survfit(Surv(`OS`,`OS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Pos"))
  OS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  OS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`OS`,`OS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Pos"))%>%
    summary()
  Plot3<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Pos"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median OS:\n",
                             OS_JS001_Median," Months vs ",OS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  SurvFit<-
    survfit(Surv(`OS`,`OS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Control"))
  OS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  OS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`OS`,`OS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Control"))%>%
    summary()
  Plot4<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Control"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median OS:\n",
                             OS_JS001_Median," Months vs ",OS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  TotalPlot<-
    Plot1$plot+ylab("Progress-free Survival")+xlab("Time (Months)")+
    Plot2$plot+ylab("Progress-free Survival")+xlab("Time (Months)")+
    Plot3$plot+ylab("Overall Survival")+xlab("Time (Months)")+
    Plot4$plot+ylab("Overall Survival")+xlab("Time (Months)")
  print(TotalPlot)
}
pdf("./Figure/ImmunogenicFeature/CI_Score_Curve_Binary.pdf",
    width = 10,height = 8)
print(TotalPlot)
dev.off()

##Interaction Comparison with ccTMB----

Kvalue_Result_ccTMB<-
  Interation_Screening_SV_BV_ForKvalue(ESCC_Jupiter06_Immunogenicity,
                                       Set_Column=which(names(ESCC_Jupiter06_Immunogenicity)%in%c("ccTMBstatus")),
                                       Set_Comparison_A="High",
                                       Set_Comparison_B=c("Low"))


Kvalue_Result_ImmunogenicFeature<-
  Interation_Screening_SV_BV_ForKvalue(ESCC_Jupiter06_Immunogenicity,
                                       Set_Column=which(names(ESCC_Jupiter06_Immunogenicity)%in%c("RiskCombine_Status")),
                                       Set_Comparison_A="immune feature-favorable",
                                       Set_Comparison_B=c("immune feature-unfavorable"))

Interaction_Pvalue_ForImmunogenicFeature<-
  rbind(Kvalue_Result_ccTMB,
        Kvalue_Result_ImmunogenicFeature)

Interaction_Pvalue_ForImmunogenicFeature<-
  Interaction_Pvalue_ForImmunogenicFeature[c("Items","OS.Inter.Pval","PFS.Inter.Pval")]%>%
  reshape2::melt()
Interaction_Pvalue_ForImmunogenicFeature$Items<-
  c("ccTMB","Immunogenic Features",
    "ccTMB","Immunogenic Features")

Plot<-
  ggplot(Interaction_Pvalue_ForImmunogenicFeature,
         aes(x=variable,
             y=-log10(value),
             fill=Items))+
  geom_bar(stat = "identity",position = "dodge")+
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.text = element_text(size=12,color="black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15,color="black"),
        legend.text = element_text(size=12))+
  scale_fill_jama()
ggsave("./Figure/ImmunogenicFeature/Interaction_Pvalue.pdf",
       plot = Plot,width = 4,height = 5)

##Immunogenicity Feature-based classification and 12-month PFS rate----

SurvFit1<-
  survfit(Surv(`PFS`,`PFS.status.`)~1,
          data = subset(ESCC_Jupiter06_Immunogenicity,
                        RiskCombine_Status=="immune feature-favorable"&Grouping=="Chemotherapy+Toripalimab"))%>%
  summary()
SurvFit2<-
  survfit(Surv(`PFS`,`PFS.status.`)~1,
          data = subset(ESCC_Jupiter06_Immunogenicity,
                        RiskCombine_Status=="immune feature-favorable"&Grouping=="Chemotherapy+Placebo"))%>%
  summary()
SurvFit3<-
  survfit(Surv(`PFS`,`PFS.status.`)~1,
          data = subset(ESCC_Jupiter06_Immunogenicity,
                        RiskCombine_Status=="immune feature-unfavorable"&Grouping=="Chemotherapy+Toripalimab"))%>%
  summary()
SurvFit4<-
  survfit(Surv(`PFS`,`PFS.status.`)~1,
          data = subset(ESCC_Jupiter06_Immunogenicity,
                        RiskCombine_Status=="immune feature-unfavorable"&Grouping=="Chemotherapy+Placebo"))%>%
  summary()

PFS_12monthRate<-
  data.frame(Rate=c(SurvFit1$surv[SurvFit1$time==SurvFit1$time[SurvFit1$time<12]%>%tail(n=1)],
                    SurvFit2$surv[SurvFit2$time==SurvFit2$time[SurvFit2$time<12]%>%tail(n=1)],
                    SurvFit3$surv[SurvFit3$time==SurvFit3$time[SurvFit3$time<12]%>%tail(n=1)],
                    SurvFit4$surv[SurvFit4$time==SurvFit4$time[SurvFit4$time<12]%>%tail(n=1)]),
             Upper=c(SurvFit1$upper[SurvFit1$time==SurvFit1$time[SurvFit1$time<12]%>%tail(n=1)],
                     SurvFit2$upper[SurvFit2$time==SurvFit2$time[SurvFit2$time<12]%>%tail(n=1)],
                     SurvFit3$upper[SurvFit3$time==SurvFit3$time[SurvFit3$time<12]%>%tail(n=1)],
                     SurvFit4$upper[SurvFit4$time==SurvFit4$time[SurvFit4$time<12]%>%tail(n=1)]),
             Lower=c(SurvFit1$lower[SurvFit1$time==SurvFit1$time[SurvFit1$time<12]%>%tail(n=1)],
                     SurvFit2$lower[SurvFit2$time==SurvFit2$time[SurvFit2$time<12]%>%tail(n=1)],
                     SurvFit3$lower[SurvFit3$time==SurvFit3$time[SurvFit3$time<12]%>%tail(n=1)],
                     SurvFit4$lower[SurvFit4$time==SurvFit4$time[SurvFit4$time<12]%>%tail(n=1)]),
             Item=c("CI-High & antiPD-1+Chemo","CI-High & Chemo",
                    "CI-Low & antiPD-1+Chemo","CI-Low & Chemo"))

Plot<-
  ggplot(PFS_12monthRate)+
  geom_bar(aes(x=Item,y=Rate,fill=Item),
           stat = "identity",position = "dodge",width = 0.7,color="black")+
  geom_errorbar(aes(x=Item,ymin=Lower,ymax=Upper),color="black",
                stat = "identity",position = "dodge",width=0.3)+
  geom_point(aes(x=Item,y=Rate),color="black",size=2)+
  scale_x_discrete(limits=rev(PFS_12monthRate$Item))+
  # coord_flip()+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.y = element_text(size=12,color="black"),
        axis.text.x = element_text(size=12,color="black",
                                   angle = 45,hjust=1),
        axis.title.y = element_text(size=15),
        panel.grid = element_blank())+
  guides(fill=FALSE)+
  scale_fill_manual(values = c("#D36810","#395863",
                               "#F3A667","#79A5B3"))+
  ylab("12-Month PFS Rate");Plot

pdf("./Figure/ImmunogenicFeature/12monthPFSrate_CIscore(FigureE).pdf",
    width = 4,height = 10)
print(Plot)

##GC-monoIO for validating the Immunogenicity-based classification----

GC_IO_Clinical_CombinatorialBiomarker<-
  readRDS("./SourceData/Mono_IO_Cohort/GC_IO_CombinatorialBiomarker.RDS")%>%
  .[c("WES_ID","ORR","ImmunogenicityCombine_Binary")]

GC_IO_Immunogenecity_Combine_ORR<-
  table(subset(GC_IO_Clinical_CombinatorialBiomarker,ORR!="UK")[c("ImmunogenicityCombine_Binary","ORR")])%>%
  prop.table(margin = 1)%>%as.data.frame()%>%subset(.,ORR=="Response")
GC_IO_Immunogenecity_Combine_ORR$ImmunogenicityCombine_Binary<-
  factor(GC_IO_Immunogenecity_Combine_ORR$ImmunogenicityCombine_Binary,
         levels=c("Unfavorable","Favorable"))
GC_IO_Immunogenecity_Combine_ORR_Table<-
  table(subset(GC_IO_Clinical_CombinatorialBiomarker,ORR!="UK")[c("ImmunogenicityCombine_Binary","ORR")])
GC_IO_Immunogenecity_Combine_ORR_TestResult<-
  chisq.test(table(subset(GC_IO_Clinical_CombinatorialBiomarker,ORR!="UK")[c("ImmunogenicityCombine_Binary","ORR")]))

Plot<-
  ggplot(GC_IO_Immunogenecity_Combine_ORR,
         aes(x=ImmunogenicityCombine_Binary,y=Freq,fill=ImmunogenicityCombine_Binary))+
  geom_bar(stat = "identity",position = "dodge",width = 0.7)+
  ylab("ORR")+
  scale_fill_jama()+
  theme_classic()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title = element_text(size=15),
        axis.text = element_text(size=12,color = "black"),
        axis.title.x = element_blank())+
  annotate(geom = "text",x="Unfavorable",y=0.2,
           label=paste0("Chi-squared Test\nP-value=",
                        round(GC_IO_Immunogenecity_Combine_ORR_TestResult$p.value,digits = 3)))+
  annotate(geom = "text",x="Unfavorable",y=GC_IO_Immunogenecity_Combine_ORR$Freq[2],
           label=paste0(round(GC_IO_Immunogenecity_Combine_ORR$Freq[2],digits = 3)%>%scales::percent(),
                        "\n(",
                        GC_IO_Immunogenecity_Combine_ORR_Table[2,2],"/",
                        GC_IO_Immunogenecity_Combine_ORR_Table[2,2]+GC_IO_Immunogenecity_Combine_ORR_Table[2,1],
                        ")"))+
  annotate(geom = "text",x="Favorable",y=GC_IO_Immunogenecity_Combine_ORR$Freq[1],
           label=paste0(round(GC_IO_Immunogenecity_Combine_ORR$Freq[1],digits = 3)%>%scales::percent(),
                        "\n(",
                        GC_IO_Immunogenecity_Combine_ORR_Table[1,2],"/",
                        GC_IO_Immunogenecity_Combine_ORR_Table[1,2]+GC_IO_Immunogenecity_Combine_ORR_Table[1,1],
                        ")"))
ggsave(plot=Plot,
       filename = "./Figure/ImmunogenicFeature/GCIO_ImmunogenicityFeatures_Validation.pdf",
       width = 4,height = 5)



# ---------------------------------------------------------------
## Scripts to Reproduce the Association between Oncogenic Events and the benefits----
# ---------------------------------------------------------------

dir.create("./Figure/Figure4")

##Figure4a----

ESCC_Jupiter06_Mutation<-
  ESCC_Jupiter06_Clinical_SupplementaryTable2[
    c(1:7,22:56)
  ]

ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub<-
  Interation_Screening_SV_BV(
    InputData = ESCC_Jupiter06_Mutation,
    Set_Column=c(8:ncol(ESCC_Jupiter06_Mutation)),
    Set_Comparison_A=c("Mutant"),
    Set_Comparison_B=c("Wild-type")
  )

ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub<-
  subset(ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub,
         Pc>=0.05&Items%in%DriverGene_Sub$Hugo_Symbol)

ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$Status[
  ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$PFS.Inter.Pval<=0.1&
    ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$PFS.Inter.k>0
]<-"Sig-Risk"
ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$Status[
  ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$PFS.Inter.Pval<=0.1&
    ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$PFS.Inter.k<0
]<-"Sig-Protective"
ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$Status[
  ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$PFS.Inter.Pval>0.1&
    ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$PFS.Inter.k>0
]<-"NonSig-Risk"
ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$Status[
  ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$PFS.Inter.Pval>0.1&
    ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$PFS.Inter.k<0
]<-"NonSig-Protective"

Plot<-
  ggplot(ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub,
         aes(y=(-log10(PFS.Inter.Pval)),x=PFS.Inter.k,
             color=Status))+
  geom_point(aes(size=Pc))+
  scale_color_manual(values=c("#FF8484","#9191FF",
                              "#203864","#840000"))+
  theme_classic()+
  theme(axis.text = element_text(size=15,color="black"),
        axis.title = element_text(size=18),
        legend.position = "top",
        legend.title = element_blank())+
  geom_text_repel(data=ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub[ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$PFS.Inter.Pval<0.1,],
                  label=ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub[ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$PFS.Inter.Pval<0.1,
                                                                          1],
                  size=5)+
  geom_hline(yintercept = 1,
             linetype="dashed")+
  ylab("-Log10(PFS-Inter-Pvalue)")+
  xlab("PFS-InteractionK");Plot
pdf("./Figure/Figure4/Mutation_Event_Screening_PFS(Figure4a).pdf",
    width = 6,height = 6)
print(Plot)
dev.off()

ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$Status[
  round(ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$OS.Inter.Pval,digits = 3)<=0.1&
    ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$OS.Inter.k>0
]<-"Sig-Risk"
ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$Status[
  round(ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$OS.Inter.Pval,digits = 3)<=0.1&
    ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$OS.Inter.k<0
]<-"Sig-Protective"
ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$Status[
  round(ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$OS.Inter.Pval,digits = 3)>0.1&
    ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$OS.Inter.k>0
]<-"NonSig-Risk"
ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$Status[
  round(ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$OS.Inter.Pval,digits = 3)>0.1&
    ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$OS.Inter.k<0
]<-"NonSig-Protective"

Plot<-
  ggplot(ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub,
         aes(y=(-log10(OS.Inter.Pval)),x=OS.Inter.k,
             color=Status))+
  geom_point(aes(size=Pc))+
  scale_color_manual(values=c("#FF8484","#9191FF",
                              "#840000","#203864"))+
  theme_classic()+
  theme(axis.text = element_text(size=15,color="black"),
        axis.title = element_text(size=18),
        legend.position = "top",
        legend.title = element_blank())+
  geom_text_repel(data=ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub[ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$Status%in%c("Sig-Protective",
                                                                                                                                       "Sig-Risk"),],
                  label=ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub[ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub$Status%in%c("Sig-Protective",
                                                                                                                                        "Sig-Risk"),
                                                                          1],
                  size=5)+
  geom_hline(yintercept = 1,
             linetype="dashed")+
  ylab("-Log10(OS-Inter-Pvalue)")+
  xlab("OS-InteractionK");Plot
pdf("./Figure/Figure4/Mutation_Event_Screening_OS(Figure4a).pdf",
    width = 6,height = 6)
print(Plot)
dev.off()

write.csv(ESCC_Jupiter06_Mutation_ScreeningResult_DriverSub[c(1:14)],
          "./Figure/Figure4/ESCC_Jupiter06_DriverGeneMut_ScreeningResult.csv",
          row.names = FALSE)

##Figure4b and Supplementary----

ESCC_JS001_Mutation_PIK3CA_TET2<-
  readRDS("./SourceData/Alteration/ESCC_JS001_Mutation_PIK3CA_TET2.RDS")%>%
  merge(ESCC_Jupiter06_Clinical,
        .,by="SUBJID")

InputGene<-"TET2"
{
  InputData_Select<-
    ESCC_JS001_Mutation_PIK3CA_TET2[c("SUBJID","PFS","PFS.status.","OS","OS.status.",
                                      "Grouping",InputGene)]%>%unique()
  InputData_Select$Grouping<-
    factor(InputData_Select$Grouping,
           levels = c("Chemotherapy+Placebo","Chemotherapy+Toripalimab"))
  names(InputData_Select)[7]<-"Item"
  InputData_Select$Item[
    InputData_Select$Item%in%c("Wildtype")
  ]<-"Control"
  InputData_Select$Item[
    InputData_Select$Item%in%c("Mut_Onco")
  ]<-"Pos"
  InputData_Select<-
    subset(InputData_Select,Item%in%c("Control","Pos"))
  InputData_Select$Item<-
    factor(InputData_Select$Item,
           levels = c("Control","Pos"))
  
  SurvFit<-
    survfit(Surv(`PFS`,`PFS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Pos"))
  PFS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  PFS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`PFS`,`PFS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Pos"))%>%
    summary()
  Plot1<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Pos"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median PFS:\n",
                             PFS_JS001_Median," Months vs ",PFS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  SurvFit<-
    survfit(Surv(`PFS`,`PFS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Control"))
  PFS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  PFS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`PFS`,`PFS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Control"))%>%
    summary()
  Plot2<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Control"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median PFS:\n",
                             PFS_JS001_Median," Months vs ",PFS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  SurvFit<-
    survfit(Surv(`OS`,`OS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Pos"))
  OS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  OS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`OS`,`OS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Pos"))%>%
    summary()
  Plot3<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Pos"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median OS:\n",
                             OS_JS001_Median," Months vs ",OS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  SurvFit<-
    survfit(Surv(`OS`,`OS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Control"))
  OS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  OS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`OS`,`OS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Control"))%>%
    summary()
  Plot4<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Control"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median OS:\n",
                             OS_JS001_Median," Months vs ",OS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  TotalPlot<-
    Plot1$plot+ylab("Progression-free Survival")+xlab("Time (Months)")+
    Plot2$plot+ylab("Progression-free Survival")+xlab("Time (Months)")+
    Plot3$plot+ylab("Overall Survival")+xlab("Time (Months)")+
    Plot4$plot+ylab("Overall Survival")+xlab("Time (Months)")
  
  print(TotalPlot)
}
pdf("./Figure/Figure4/TET2_Curve.pdf",
    width = 10,height = 8)
print(TotalPlot)
dev.off()

InputGene<-"PIK3CA"
{
  InputData_Select<-
    ESCC_JS001_Mutation_PIK3CA_TET2[c("SUBJID","PFS","PFS.status.","OS","OS.status.",
                                      "Grouping",InputGene)]%>%unique()
  InputData_Select$Grouping<-
    factor(InputData_Select$Grouping,
           levels = c("Chemotherapy+Placebo","Chemotherapy+Toripalimab"))
  names(InputData_Select)[7]<-"Item"
  InputData_Select$Item[
    InputData_Select$Item%in%c("Wildtype")
  ]<-"Control"
  InputData_Select$Item[
    InputData_Select$Item%in%c("Mut_Onco","Mut-VUS")
  ]<-"Pos"
  InputData_Select<-
    subset(InputData_Select,Item%in%c("Control","Pos"))
  InputData_Select$Item<-
    factor(InputData_Select$Item,
           levels = c("Control","Pos"))
  
  SurvFit<-
    survfit(Surv(`PFS`,`PFS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Pos"))
  PFS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  PFS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`PFS`,`PFS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Pos"))%>%
    summary()
  Plot1<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Pos"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median PFS:\n",
                             PFS_JS001_Median," Months vs ",PFS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  
  SurvFit<-
    survfit(Surv(`PFS`,`PFS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Control"))
  PFS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  PFS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`PFS`,`PFS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Control"))%>%
    summary()
  Plot2<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Control"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median PFS:\n",
                             PFS_JS001_Median," Months vs ",PFS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  SurvFit<-
    survfit(Surv(`OS`,`OS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Pos"))
  OS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  OS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`OS`,`OS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Pos"))%>%
    summary()
  Plot3<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Pos"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median OS:\n",
                             OS_JS001_Median," Months vs ",OS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  SurvFit<-
    survfit(Surv(`OS`,`OS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Control"))
  OS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  OS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`OS`,`OS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Control"))%>%
    summary()
  Plot4<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Control"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median OS:\n",
                             OS_JS001_Median," Months vs ",OS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  TotalPlot<-
    Plot1$plot+ylab("Progression-free Survival")+xlab("Time (Months)")+
    Plot2$plot+ylab("Progression-free Survival")+xlab("Time (Months)")+
    Plot3$plot+ylab("Overall Survival")+xlab("Time (Months)")+
    Plot4$plot+ylab("Overall Survival")+xlab("Time (Months)")
  
  print(TotalPlot)
}
pdf("./Figure/Figure4/PIK3CA_Curve.pdf",
    width = 10,height = 8)
print(TotalPlot)
dev.off()

##Figure4c----


###Lesion Level----

# ESCC_Jupiter06_CNV_Amp<-
#   ESCC_Jupiter06_Clinical_SupplementaryTable2[
#     c(1:7,grep("Amp",names(ESCC_Jupiter06_Clinical_SupplementaryTable2)))
#   ]
# 
# ESCC_JS001_LesionLevelCNV_SequenzaAmp_ScreeningResult<-
#   Interation_Screening_SV_BV(
#     ESCC_Jupiter06_CNV_Amp,
#     Set_Column=c(8:ncol(ESCC_Jupiter06_CNV_Amp)),
#     Set_Comparison_A=c("High-level amplification"),
#     Set_Comparison_B=c("No amplification")
#   )
# 
# ESCC_Jupiter06_CNV_Del<-
#   ESCC_Jupiter06_Clinical_SupplementaryTable2[
#     c(1:7,grep("Del",names(ESCC_Jupiter06_Clinical_SupplementaryTable2)))
#   ]
# ESCC_JS001_LesionLevelCNV_SequenzaDel_ScreeningResult<-
#   Interation_Screening_SV_BV(
#     ESCC_Jupiter06_CNV_Del,
#     Set_Column=c(8:ncol(ESCC_Jupiter06_CNV_Del)),
#     Set_Comparison_A=c("Deep deletion"),
#     Set_Comparison_B=c("No deletion")
#   )
# ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult<-
#   rbind(ESCC_JS001_LesionLevelCNV_SequenzaAmp_ScreeningResult,
#         ESCC_JS001_LesionLevelCNV_SequenzaDel_ScreeningResult)
# saveRDS(ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult,
#         "./SourceData/Alteration/ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult.RDS")

ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult<-
  readRDS("./SourceData/Alteration/ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult.RDS")
write.csv(ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult[c(1:14)],
          "./Figure/Figure4/ESCC_Jupiter06_LesionLevelCNA_ScreeningResult.csv",
          row.names = FALSE)

ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$Status[
  ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$PFS.Inter.Pval<=0.1&
    ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$PFS.Inter.k>0
]<-"Sig-Risk"
ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$Status[
  ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$PFS.Inter.Pval<=0.1&
    ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$PFS.Inter.k<0
]<-"Sig-Protective"
ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$Status[
  ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$PFS.Inter.Pval>0.1&
    ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$PFS.Inter.k>0
]<-"NonSig-Risk"
ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$Status[
  ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$PFS.Inter.Pval>0.1&
    ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$PFS.Inter.k<0
]<-"NonSig-Protective"

Plot<-
  ggplot(ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult,
         aes(y=(-log10(PFS.Inter.Pval)),x=PFS.Inter.k,
             color=Status))+
  geom_point(aes(size=Pc))+
  scale_color_manual(values=c("#FF8484","#9191FF",
                              "#840000","#203864"))+
  theme_classic()+
  theme(axis.text = element_text(size=15,color="black"),
        axis.title = element_text(size=18),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=12))+
  geom_hline(yintercept = 1,
             linetype="dashed")+
  ylab("-Log10(PFS-Inter-Pvalue)")+
  xlab("PFS-InteractionK");Plot

pdf("./Figure/Figure4/LesionCNV_Screening_PFS(Figure4c).pdf",
    width = 6,height = 6)
print(Plot)
dev.off()

ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$Status[
  ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$OS.Inter.Pval<0.1&
    ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$OS.Inter.k>0
]<-"Sig-Risk"
ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$Status[
  ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$OS.Inter.Pval<0.1&
    ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$OS.Inter.k<0
]<-"Sig-Protective"
ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$Status[
  ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$OS.Inter.Pval>=0.1&
    ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$OS.Inter.k>0
]<-"NonSig-Risk"
ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$Status[
  ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$OS.Inter.Pval>=0.1&
    ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$OS.Inter.k<0
]<-"NonSig-Protective"

Plot<-
  ggplot(ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult,
         aes(y=(-log10(OS.Inter.Pval)),x=OS.Inter.k,
             color=Status))+
  geom_point(aes(size=Pc))+
  scale_color_manual(values=c("#FF8484","#9191FF",
                              "#840000","#203864"))+
  theme_classic()+
  theme(axis.text = element_text(size=15,color="black"),
        axis.title = element_text(size=18),
        legend.position = "top",
        legend.title = element_blank())+
  geom_hline(yintercept = 1,
             linetype="dashed")+
  geom_text_repel(data=ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult[ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$Status%in%c("Sig-Protective",
                                                                                                                                       "Sig-Risk"),],
                  label=ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult[ESCC_JS001_LesionLevelCNV_Sequenza_ScreeningResult$Status%in%c("Sig-Protective",
                                                                                                                                        "Sig-Risk"),
                                                                          1],
                  size=5)+
  ylab("-Log10(OS-Inter-Pvalue)")+
  xlab("OS-InteractionK");Plot

pdf("./Figure/Figure4/LesionCNV_Screening_OS(Figure4c).pdf",
    width = 6,height = 6)
print(Plot)
dev.off()

##Figure4d and Supplementary----

InputGene<-"22q11.21-Amp"
{
  InputData_Select<-
    ESCC_Jupiter06_CNV_Amp[c("SUBJID","PFS","PFS.status.","OS","OS.status.",
                             "Grouping",InputGene)]%>%unique()
  InputData_Select$Grouping<-
    factor(InputData_Select$Grouping,
           levels = c("Chemotherapy+Placebo","Chemotherapy+Toripalimab"))
  names(InputData_Select)[7]<-"Item"
  InputData_Select$Item[
    InputData_Select$Item%in%c("No amplification")
  ]<-"Control"
  InputData_Select$Item[
    InputData_Select$Item%in%c("High-level amplification")
  ]<-"Pos"
  InputData_Select<-
    subset(InputData_Select,Item%in%c("Control","Pos"))
  InputData_Select$Item<-
    factor(InputData_Select$Item,
           levels = c("Control","Pos"))
  
  SurvFit<-
    survfit(Surv(`PFS`,`PFS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Pos"))
  PFS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  PFS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`PFS`,`PFS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Pos"))%>%
    summary()
  Plot1<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Pos"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median PFS:\n",
                             PFS_JS001_Median," Months vs ",PFS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  SurvFit<-
    survfit(Surv(`PFS`,`PFS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Control"))
  PFS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  PFS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`PFS`,`PFS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Control"))%>%
    summary()
  Plot2<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Control"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median PFS:\n",
                             PFS_JS001_Median," Months vs ",PFS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  SurvFit<-
    survfit(Surv(`OS`,`OS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Pos"))
  OS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  OS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`OS`,`OS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Pos"))%>%
    summary()
  Plot3<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Pos"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median OS:\n",
                             OS_JS001_Median," Months vs ",OS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  SurvFit<-
    survfit(Surv(`OS`,`OS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Control"))
  OS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  OS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`OS`,`OS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Control"))%>%
    summary()
  Plot4<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Control"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median OS:\n",
                             OS_JS001_Median," Months vs ",OS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  TotalPlot<-
    Plot1$plot+ylab("Progression-free Survival")+xlab("Time (Months)")+
    Plot2$plot+ylab("Progression-free Survival")+xlab("Time (Months)")+
    Plot3$plot+ylab("Overall Survival")+xlab("Time (Months)")+
    Plot4$plot+ylab("Overall Survival")+xlab("Time (Months)")
  
  print(TotalPlot)
}
pdf("./Figure/Figure4/22q11.21Amp_Curve.pdf",
    width = 10,height = 8)
print(TotalPlot)
dev.off()

InputGene<-"1q21.3-Amp"
{
  InputData_Select<-
    ESCC_Jupiter06_CNV_Amp[c("SUBJID","PFS","PFS.status.","OS","OS.status.",
                             "Grouping",InputGene)]%>%unique()
  InputData_Select$Grouping<-
    factor(InputData_Select$Grouping,
           levels = c("Chemotherapy+Placebo","Chemotherapy+Toripalimab"))
  names(InputData_Select)[7]<-"Item"
  InputData_Select$Item[
    InputData_Select$Item%in%c("No amplification")
  ]<-"Control"
  InputData_Select$Item[
    InputData_Select$Item%in%c("High-level amplification")
  ]<-"Pos"
  InputData_Select<-
    subset(InputData_Select,Item%in%c("Control","Pos"))
  InputData_Select$Item<-
    factor(InputData_Select$Item,
           levels = c("Control","Pos"))
  
  SurvFit<-
    survfit(Surv(`PFS`,`PFS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Pos"))
  PFS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  PFS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`PFS`,`PFS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Pos"))%>%
    summary()
  Plot1<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Pos"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median PFS:\n",
                             PFS_JS001_Median," Months vs ",PFS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  SurvFit<-
    survfit(Surv(`PFS`,`PFS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Control"))
  PFS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  PFS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`PFS`,`PFS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Control"))%>%
    summary()
  Plot2<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Control"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median PFS:\n",
                             PFS_JS001_Median," Months vs ",PFS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  SurvFit<-
    survfit(Surv(`OS`,`OS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Pos"))
  OS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  OS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`OS`,`OS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Pos"))%>%
    summary()
  Plot3<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Pos"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median OS:\n",
                             OS_JS001_Median," Months vs ",OS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  SurvFit<-
    survfit(Surv(`OS`,`OS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Control"))
  OS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  OS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`OS`,`OS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Control"))%>%
    summary()
  Plot4<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Control"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median OS:\n",
                             OS_JS001_Median," Months vs ",OS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  TotalPlot<-
    Plot1$plot+ylab("Progression-free Survival")+xlab("Time (Months)")+
    Plot2$plot+ylab("Progression-free Survival")+xlab("Time (Months)")+
    Plot3$plot+ylab("Overall Survival")+xlab("Time (Months)")+
    Plot4$plot+ylab("Overall Survival")+xlab("Time (Months)")
  
  print(TotalPlot)
}
pdf("./Figure/Figure4/1q21.3Amp_Curve.pdf",
    width = 10,height = 8)
print(TotalPlot)
dev.off()

##Figure4e----
###Lesion and Immune in TCGA-ESCC

TCGA_Immune_CibersortAbs_CNV<-
  readRDS("./SourceData/Amp&Immune/TCGA_Immune_CibersortAbs_22q11.21Amp.RDS")
CibersortAbs_CNV_Cor<-
  corr.test(TCGA_Immune_CibersortAbs_CNV[c(4:25)],
            TCGA_Immune_CibersortAbs_CNV[c(26)],
            adjust ="none")
CibersortAbs_CNV_Cor$r<-
  CibersortAbs_CNV_Cor$r[
    order(CibersortAbs_CNV_Cor$r[,1]),
  ]
Plot<-
  ggcorrplot(CibersortAbs_CNV_Cor$r%>%as.data.frame(),
             method = "circle",
             colors = c("darkblue","white","darkred"));Plot
pdf("./Figure/Figure4/22q11.21_Amp_Immune(Figure4e).pdf",
    width = 12,height = 4)
print(Plot)
dev.off()

###Lesion and Expression in TCGA-ESCC

Result_Frame_Cor<-
  readRDS("./SourceData/Amp&Immune/TCGA_ESCC_22q11.21genes_ExprCorCNV.RDS")
Result_Frame_Cor$Sig<-
  ifelse(Result_Frame_Cor$Pvalue<0.05,
         "Sig","Non-Sig")
Result_Frame_Cor$CorStatus<-
  ifelse(Result_Frame_Cor$Cor>=0.4,
         "Cor>0.4","Cor<0.4")

Plot<-
  ggplot(subset(Result_Frame_Cor,!is.na(Result_Frame_Cor$Pvalue)),
         aes(x=Cor,y=(-log10(Pvalue)),
             color=Sig))+
  geom_point()+
  geom_hline(yintercept = (-log10(0.05)),
             linetype="dashed")+
  ylab("-log10(P-value)")+
  xlab("Correlation")+
  theme_classic()+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=15,
                                 colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size=15))+
  scale_color_manual(values=c("#0072B5","#BC3C29"));Plot
ggsave(plot = Plot,
       filename = "./Figure/Figure4/TCGA_ESCC_22q11.21genes_ExprCorCNV.pdf",
       width = 6,height = 4.5)

TCGA_ESCC_Expression_FPKM_Select_t<-
  readRDS("./SourceData/Amp&Immune/TCGA_ESCC_22q11.21genesExpression_FPKM.RDS")
TCGA_ESCC_Expression22q11.21_FPKMcor<-
  cor(TCGA_ESCC_Expression_FPKM_Select_t)
Plot<-
  ggcorrplot(TCGA_ESCC_Expression22q11.21_FPKMcor[
    !is.na(TCGA_ESCC_Expression22q11.21_FPKMcor[,1]),
    !is.na(TCGA_ESCC_Expression22q11.21_FPKMcor[1,])
  ],
  ggtheme = theme(axis.text = element_text(size=3)))+
  theme(axis.text = element_text(size=3))
ggsave(filename = "./Figure/Figure4/TCGA_ESCC_Expression22q11.21_FPKMcor.pdf",
       plot = Plot,
       width = 16,height = 16)

###Expression and Immune for CLTCL1

####CLTCL1 in TCGA cohort

TCGA_Immune_CibersortAbs_FPKM<-
  readRDS("./SourceData/Amp&Immune/TCGA_Immune_CibersortAbs_CLTCL1fpkm.RDS")

CibersortAbs_FPKM_Cor<-
  corr.test(TCGA_Immune_CibersortAbs_FPKM[c(4:25)],
            TCGA_Immune_CibersortAbs_FPKM[c(26)],
            adjust ="none")
CibersortAbs_FPKM_Cor$r<-
  CibersortAbs_FPKM_Cor$r[
    order(CibersortAbs_FPKM_Cor$r[,1]),
  ]

cor.test(TCGA_Immune_CibersortAbs_FPKM$CLTCL1,
         TCGA_Immune_CibersortAbs_FPKM$T.cells.CD8,
         method = "spearman")
Plot<-
  ggplot(TCGA_Immune_CibersortAbs_FPKM,
         aes(x=CLTCL1,y=T.cells.CD8))+
  geom_point()+
  geom_smooth(method=lm)+
  xlab(paste("CLTCL-Log2(FPKM+1)"))+
  ylab(paste("T.cells.CD8 Infiltration"))+
  theme_bw()+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=15,color = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(size = 20))+
  annotate(geom = "text",x=1.8,y=0.075,
           label="Spearman Rho = -0.431\nP-value<0.001");Plot

ggsave(plot = Plot,
       filename = "./Figure/Figure4/CLTCL1_ExprCorCD8.pdf",
       width = 6.5,height = 6)

####CLTCL1 in NC2017 cohort(PMID: 28548104)

ESCC_NC2017_CibersortAbs_CLTCL1<-
  readRDS("./SourceData/Amp&Immune/ESCC_NC2017_CibersortAbs_CLTCL1.RDS")
cor.test(ESCC_NC2017_CibersortAbs_CLTCL1[,"CLTCL1"],
         ESCC_NC2017_CibersortAbs_CLTCL1$T.cells.CD8)
Plot<-
  ggplot(ESCC_NC2017_CibersortAbs_CLTCL1,
         aes(x=CLTCL1,y=T.cells.CD8))+
  geom_point()+
  geom_smooth(method=lm)+
  xlab(paste("CLTCL-Log2(TPM+1)"))+
  ylab(paste("T.cells.CD8 Infiltration"))+
  theme_bw()+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=15,color = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(size = 20))+
  annotate(geom = "text",x=6.5,y=0.4,
           label="Spearman Rho = -0.241\nP-value=0.019");Plot
ggsave(plot = Plot,
       filename = "./Figure/Figure4/CLTCL1_ExprCorCD8_NC2017.pdf",
       width = 6.5,height = 6)

##Figure4-Oncogenic-Event-Combination----

ESCC_JS001_RiskMolEventCombine<-
  merge(ESCC_Jupiter06_Mutation[c("SUBJID","PIK3CA","TET2")],
        ESCC_Jupiter06_CNV_Amp[c("SUBJID","22q11.21-Amp","1q21.3-Amp")],
        by="SUBJID")

ESCC_JS001_RiskMolEventCombine$PIK3CA<-
  ifelse(ESCC_JS001_RiskMolEventCombine$PIK3CA!="Wild-type",
         1,0)
ESCC_JS001_RiskMolEventCombine$TET2<-
  ifelse(ESCC_JS001_RiskMolEventCombine$TET2!="Wild-type",
         1,0)
ESCC_JS001_RiskMolEventCombine$`1q21.3-Amp`<-
  ifelse(ESCC_JS001_RiskMolEventCombine$`1q21.3-Amp`=="High-level amplification",
         1,0)
ESCC_JS001_RiskMolEventCombine$`22q11.21-Amp`<-
  ifelse(ESCC_JS001_RiskMolEventCombine$`22q11.21-Amp`=="High-level amplification",
         1,0)
ESCC_JS001_RiskMolEventCombine$EventCount<-
  apply(ESCC_JS001_RiskMolEventCombine[c(2:5)],1,function(x){
    sum(x)
  })

EventPropTable<-
  table(ESCC_JS001_RiskMolEventCombine$EventCount)%>%
  prop.table()%>%
  as.data.frame()
EventPropTable$Var1<-as.numeric(EventPropTable$Var1%>%as.character())
EventPropTable$Status<-
  ifelse(EventPropTable[,1]>0,"Positive","Negative")
EventPropTable$Var1<-
  as.factor(EventPropTable$Var1)

###Frequency

Plot<-
  ggplot(EventPropTable,
         aes(x=Status,fill=Var1,y=Freq))+
  geom_bar(stat = "identity",position = "stack")+
  theme_classic()+
  theme(axis.text = element_text(size=15,color = "black"),
        axis.title = element_text(size=18),
        legend.position = "top",
        legend.text = element_text(size=15),
        legend.title = element_text(size=18))+
  xlab("Event Status")+
  ylab("Frequency")+
  scale_fill_nejm();Plot
pdf("./Figure/Figure4/EventStatus_Frequency.pdf",
    width = 4,height = 8)
print(Plot)
dev.off()

###Survival Differentiation

ESCC_JS001_RiskMolEventCombine$EventStatus<-
  ifelse(
    ESCC_JS001_RiskMolEventCombine$EventCount>0,
    "Positive","Negative"
  )

ESCC_JS001_RiskMolEventCombine<-
  merge(ESCC_JS001_RiskMolEventCombine,
        ESCC_Jupiter06_Clinical,
        by="SUBJID")
RiskMolEvent_HR<-
  SubgroupHR_Plot(ESCC_JS001_RiskMolEventCombine,
                  Set_Column=which(names(ESCC_JS001_RiskMolEventCombine)=="EventStatus"),
                  Set_Comparison_Group=c("Positive","Negative"))

pdf("./Figure/Figure4/OncogenicEvent_forrest_PFS.pdf",
    width = 10,height = 8)
print(RiskMolEvent_HR$PFS)
dev.off()
pdf("./Figure/Figure4/OncogenicEvent_forrest_OS.pdf",
    width = 10,height = 8)
print(RiskMolEvent_HR$OS)
dev.off()

InputGene<-"EventStatus"
{
  InputData_Select<-
    ESCC_JS001_RiskMolEventCombine[c("SUBJID","PFS","PFS.status.","OS","OS.status.",
                                     "Grouping",InputGene)]%>%unique()
  InputData_Select$Grouping<-
    factor(InputData_Select$Grouping,
           levels = c("Chemotherapy+Placebo","Chemotherapy+Toripalimab"))
  names(InputData_Select)[7]<-"Item"
  InputData_Select$Item[
    InputData_Select$Item%in%c("Negative")
  ]<-"Control"
  InputData_Select$Item[
    InputData_Select$Item%in%c("Positive")
  ]<-"Pos"
  InputData_Select<-
    subset(InputData_Select,Item%in%c("Control","Pos"))
  InputData_Select$Item<-
    factor(InputData_Select$Item,
           levels = c("Control","Pos"))
  
  SurvFit<-
    survfit(Surv(`PFS`,`PFS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Pos"))
  PFS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  PFS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`PFS`,`PFS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Pos"))%>%
    summary()
  Plot1<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Pos"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median PFS:\n",
                             PFS_JS001_Median," Months vs ",PFS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  SurvFit<-
    survfit(Surv(`PFS`,`PFS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Control"))
  PFS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  PFS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`PFS`,`PFS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Control"))%>%
    summary()
  Plot2<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Control"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median PFS:\n",
                             PFS_JS001_Median," Months vs ",PFS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  SurvFit<-
    survfit(Surv(`OS`,`OS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Pos"))
  OS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  OS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`OS`,`OS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Pos"))%>%
    summary()
  Plot3<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Pos"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median OS:\n",
                             OS_JS001_Median," Months vs ",OS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  SurvFit<-
    survfit(Surv(`OS`,`OS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="Control"))
  OS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  OS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`OS`,`OS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="Control"))%>%
    summary()
  Plot4<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="Control"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median OS:\n",
                             OS_JS001_Median," Months vs ",OS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  TotalPlot<-
    Plot1$plot+ylab("Progression-free Survival")+xlab("Time (Months)")+
    Plot2$plot+ylab("Progression-free Survival")+xlab("Time (Months)")+
    Plot3$plot+ylab("Overall Survival")+xlab("Time (Months)")+
    Plot4$plot+ylab("Overall Survival")+xlab("Time (Months)")
  
  print(TotalPlot)
}
pdf("./Figure/Figure4/OncogenicEvent_Curve.pdf",
    width = 10,height = 8)
print(TotalPlot)
dev.off()

##Figure4-Oncogenic-Event-Combination-Validation----

Miao_Melanoma_Cohort_30150660_Events<-
  readRDS("./SourceData/Mono_IO_Cohort/Miao30150660_Melanoma_EventsCombine_Clinical.RDS")
Braun_RCC_Cohort_32472114_Clinical_Event<-
  readRDS("./SourceData/Mono_IO_Cohort/Braun32472114_RCC_EventsCombine_Clinical.RDS")
GC_integrated_Cohort_Clinical_Event<-
  readRDS("./SourceData/Mono_IO_Cohort/GC_IO_CombinatorialBiomarker.RDS")%>%
  .[c("WES_ID","EBV_Status","MSI_Status",
      "ORR","CNV","Mutation","Events")]

ResponseRate_Genomic_Combine_Freq<-
  rbind(subset(GC_integrated_Cohort_Clinical_Event,EBV_Status!="Positive"&MSI_Status!="MSI"&ORR!="UK")[c("ORR","Events")]%>%
          mutate(CancerType="GC-MSS&EBVneg"),
        subset(Miao_Melanoma_Cohort_30150660_Events,CANCER_TYPE=="Melanoma"&RECIST!="X")[c("ORR","Events")]%>%
          mutate(ORR=ifelse(.$ORR=="CR/PR","Response","Non-Response"))%>%
          mutate(CancerType="Melanoma"),
        subset(Braun_RCC_Cohort_32472114_Clinical_Event,Arm=="NIVOLUMAB"&Cohort=="CM-025"&
                 ORR!="NE")[c("BOR","Events")]%>%mutate(ORR=ifelse(.$BOR=="R","Response","Non-Response"))%>%
          .[c("ORR","Events")]%>%
          mutate(CancerType="RCC"))

table(ResponseRate_Genomic_Combine_Freq[c("ORR","Events","CancerType")])%>%
  mantelhaen.test()

RCC_Genomic_Combine<-
  table(subset(Braun_RCC_Cohort_32472114_Clinical_Event,Arm=="NIVOLUMAB"&Cohort=="CM-025"&
                 ORR!="NE")[c("BOR","Events")])%>%
  prop.table(margin = 2)%>%
  as.data.frame()%>%
  mutate(ORR=ifelse(.$BOR=="R","Response","Non-Response"))%>%
  mutate(CancerType="RCC")
Mel_Genomic_Combine<-
  table(subset(Miao_Melanoma_Cohort_30150660_Events,CANCER_TYPE=="Melanoma"&RECIST!="X")[c("ORR","Events")])%>%
  prop.table(margin = 2)%>%
  as.data.frame()%>%
  mutate(ORR=ifelse(.$ORR=="CR/PR","Response","Non-Response"))%>%
  mutate(CancerType="Melanoma")
GC_MSS_Genomic_Combine<-
  table(subset(GC_integrated_Cohort_Clinical_Event,EBV_Status!="Positive"&MSI_Status!="MSI"&ORR!="UK")[c("ORR","Events")])%>%
  prop.table(margin = 2)%>%
  as.data.frame()%>%
  mutate(CancerType="GC-MSS&EBVneg")

ResponseRate_Genomic_Combine<-
  rbind(RCC_Genomic_Combine[c("CancerType","Events","ORR","Freq")],
        GC_MSS_Genomic_Combine[c("CancerType","Events","ORR","Freq")],
        Mel_Genomic_Combine[c("CancerType","Events","ORR","Freq")])%>%
  subset(.,ORR=="Response")
ResponseRate_Genomic_Combine$Proportion<-
  round(ResponseRate_Genomic_Combine$Freq,
        digits = 2)
Plot<-
  ggplot(ResponseRate_Genomic_Combine,aes(x=CancerType,y=Freq,fill=Events))+
  geom_bar(stat = "identity",position = "dodge")+
  # geom_text(aes(label=Proportion, y=Proportion+0.01), position=position_dodge(0.9), vjust=0,size=4.5)+
  # scale_y_continuous(labels=scale::percent)+
  ylab("ORR")+
  theme_classic()+
  theme(axis.text = element_text(size = 12,color="black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        legend.position = "top",
        legend.title = element_text(size=15),
        legend.text = element_text(size=12))+
  scale_fill_jco();Plot
ggsave(plot=Plot,
       filename = "./Figure/Figure4/3CancerType_OncogenicEvent_Validation.pdf",
       width = 6,height = 5)

# ---------------------------------------------------------------
## Scripts to Reproduce the EGIC----
# ---------------------------------------------------------------

dir.create("./Figure/Figure5/")

ESCC_Jupiter06_EGIC<-
  merge(ESCC_Jupiter06_Immunogenicity[c("SUBJID","RiskCombine_Status")],
        ESCC_JS001_RiskMolEventCombine[c("SUBJID","EventStatus")],
        by="SUBJID")%>%
  merge(.,
        ESCC_Jupiter06_Clinical,
        by="SUBJID")
ESCC_Jupiter06_EGIC$FinalStatus<-
  case_when(ESCC_Jupiter06_EGIC$EventStatus=="Positive"&
              ESCC_Jupiter06_EGIC$RiskCombine_Status=="immune feature-favorable"~"EGIC2",
            ESCC_Jupiter06_EGIC$EventStatus=="Negative"&
              ESCC_Jupiter06_EGIC$RiskCombine_Status=="immune feature-favorable"~"EGIC1",
            ESCC_Jupiter06_EGIC$EventStatus=="Positive"&
              ESCC_Jupiter06_EGIC$RiskCombine_Status=="immune feature-unfavorable"~"EGIC3",
            ESCC_Jupiter06_EGIC$EventStatus=="Negative"&
              ESCC_Jupiter06_EGIC$RiskCombine_Status=="immune feature-unfavorable"~"EGIC2")

##Figure5a----

ESCC_Jupiter06_EGIC<-
  ESCC_Jupiter06_EGIC[
    order(ESCC_Jupiter06_EGIC$EventStatus,
          ESCC_Jupiter06_EGIC$RiskCombine_Status),
  ]
# 
# ESCC_Jupiter06_EGIC$Treatment_Strategy<-
#   ifelse(ESCC_Jupiter06_EGIC$EventStatus=="Positive",
#          "Chemo+Anti-PD-1+Targeted Therapy",
#          "Chemo+Anti-PD-1")

circlize_plot<-function(){
  circos.par(gap.degree=90)
  # {
  #   col_RiskEvent = structure(c("#77A6DF","#112D4E"),
  #                             names = unique(ESCC_Jupiter06_EGIC$Treatment_Strategy%>%
  #                                              as.character()))
  #   circos.heatmap(ESCC_Jupiter06_EGIC$Treatment_Strategy%>%as.character(),
  #                  col = col_RiskEvent, track.height = 0.1)
  # }
  
  ##EGIC
  {
    col_EGIC = structure(c("#3A9679","#FF9A00","#3EC1D3"),
                         names = unique(ESCC_Jupiter06_EGIC$FinalStatus%>%
                                          as.character()))
    circos.heatmap(ESCC_Jupiter06_EGIC$FinalStatus%>%as.character(),
                   col = col_EGIC, track.height = 0.2) 
  }
  
  {
    col_ImmuneFeatures = structure(c("#B20600","#FE938F"),
                                   names = unique(ESCC_Jupiter06_EGIC$RiskCombine_Status%>%
                                                    as.character()))
    circos.heatmap(ESCC_Jupiter06_EGIC$RiskCombine_Status%>%as.character(),
                   col = col_ImmuneFeatures, track.height = 0.1)
  }
  
  {
    col_RiskEvent = structure(c("#77A6DF","#112D4E"),
                              names = unique(ESCC_Jupiter06_EGIC$EventStatus%>%
                                               as.character()))
    circos.heatmap(ESCC_Jupiter06_EGIC$EventStatus%>%as.character(),
                   col = col_RiskEvent, track.height = 0.1)
  }
  
  circos.clear()
}
col_EGIC = structure(c("#3A9679","#FF9A00","#3EC1D3"),
                     names = unique(ESCC_Jupiter06_EGIC$FinalStatus%>%
                                      as.character()))
col_ImmuneFeatures = structure(c("#B20600","#FE938F"),
                               names = unique(ESCC_Jupiter06_EGIC$RiskCombine_Status%>%
                                                as.character()))
col_RiskEvent = structure(c("#77A6DF","#112D4E"),
                          names = unique(ESCC_Jupiter06_EGIC$EventStatus%>%
                                           as.character()))

lgd_EGIC = Legend(title = "EGIC", at = names(col_EGIC), 
                  legend_gp = gpar(fill = col_EGIC))
lgd_ImmuneFeatures = Legend(title = "Immunogenicity Features", at = names(col_ImmuneFeatures), 
                            legend_gp = gpar(fill = col_ImmuneFeatures))
lgd_RiskEvent = Legend(title = "Risk Events", at = names(col_RiskEvent), 
                       legend_gp = gpar(fill = col_RiskEvent))
circlize_plot()


h = dev.size()[2]
lgd_list = packLegend(lgd_EGIC,lgd_ImmuneFeatures,lgd_RiskEvent,
                      max_height = unit(0.9*h, "inch"))
draw(lgd_list, just = "left")

dev.off()

##Figure5b-c and Supplementary----

###Curve

InputGene="FinalStatus"
{
  InputData_Select<-
    ESCC_Jupiter06_EGIC[c("SUBJID","PFS","PFS.status.","OS","OS.status.",
                          "Grouping",InputGene)]%>%unique()
  InputData_Select$Grouping<-
    factor(InputData_Select$Grouping,
           levels = c("Chemotherapy+Placebo","Chemotherapy+Toripalimab"))
  names(InputData_Select)[7]<-"Item"
  
  SurvFit<-
    survfit(Surv(`PFS`,`PFS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="EGIC1"))
  PFS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  PFS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`PFS`,`PFS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="EGIC1"))%>%
    summary()
  Plot1<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="EGIC1"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median PFS:\n",
                             PFS_JS001_Median," Months vs ",PFS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  
  SurvFit<-
    survfit(Surv(`PFS`,`PFS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="EGIC2"))
  PFS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  PFS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`PFS`,`PFS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="EGIC2"))%>%
    summary()
  Plot2<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="EGIC2"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median PFS:\n",
                             PFS_JS001_Median," Months vs ",PFS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  SurvFit<-
    survfit(Surv(`PFS`,`PFS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="EGIC3"))
  PFS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  PFS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`PFS`,`PFS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="EGIC3"))%>%
    summary()
  Plot3<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="EGIC3"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median PFS:\n",
                             PFS_JS001_Median," Months vs ",PFS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  TotalPlot<-
    Plot1$plot+ylab("Progression-free Survival")+xlab("Time (Months)")+
    Plot2$plot+ylab("Progression-free Survival")+xlab("Time (Months)")+
    Plot3$plot+ylab("Progression-free Survival")+xlab("Time (Months)")
  
  print(TotalPlot)
}
pdf("./Figure/Figure5/EGIC_Curve_PFS.pdf",
    width = 20,height = 6)
print(TotalPlot)
dev.off()

InputGene="FinalStatus"
{
  InputData_Select<-
    ESCC_Jupiter06_EGIC[c("SUBJID","PFS","PFS.status.","OS","OS.status.",
                          "Grouping",InputGene)]%>%unique()
  InputData_Select$Grouping<-
    factor(InputData_Select$Grouping,
           levels = c("Chemotherapy+Placebo","Chemotherapy+Toripalimab"))
  names(InputData_Select)[7]<-"Item"
  
  SurvFit<-
    survfit(Surv(`OS`,`OS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="EGIC1"))
  OS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  OS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`OS`,`OS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="EGIC1"))%>%
    summary()
  Plot1<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="EGIC1"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median OS:\n",
                             OS_JS001_Median," Months vs ",OS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  
  SurvFit<-
    survfit(Surv(`OS`,`OS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="EGIC2"))
  OS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  OS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`OS`,`OS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="EGIC2"))%>%
    summary()
  Plot2<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="EGIC2"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median OS:\n",
                             OS_JS001_Median," Months vs ",OS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  SurvFit<-
    survfit(Surv(`OS`,`OS.status.`)~Grouping,
            data = subset(InputData_Select,Item=="EGIC3"))
  OS_Placebo_Median<-round(surv_median(SurvFit)[1,2],digits = 1)
  OS_JS001_Median<-round(surv_median(SurvFit)[2,2],digits = 1)
  CoxResult<-
    coxph(Surv(`OS`,`OS.status.`)~Grouping,
          data = subset(InputData_Select,Item=="EGIC3"))%>%
    summary()
  Plot3<-
    ggsurvplot(SurvFit,
               subset(InputData_Select,Item=="EGIC3"),
               legend=c(0.8,0.9),
               legend.title="",
               legend.labs=c("Chemo+Placebo","Chemo+JS001"),
               pval = paste0("Median OS:\n",
                             OS_JS001_Median," Months vs ",OS_Placebo_Median," Months",
                             "\nPvalue=",CoxResult$coefficients[1,5]%>%round(digits = 2),
                             "\nHR=",CoxResult$coefficients[1,2]%>%round(digits = 2)),
               palette = c("#2D5662","#F79019"),
               break.x.by=3)
  
  TotalPlot<-
    Plot1$plot+ylab("Overall Survival")+xlab("Time (Months)")+
    Plot2$plot+ylab("Overall Survival")+xlab("Time (Months)")+
    Plot3$plot+ylab("Overall Survival")+xlab("Time (Months)")
  
  print(TotalPlot)
}
pdf("./Figure/Figure5/EGIC_Curve_OS.pdf",
    width = 20,height = 6)
print(TotalPlot)
dev.off()

###Forrest plot

CombineFinalStatus_HR<-
  SubgroupHR_Plot(ESCC_Jupiter06_EGIC,
                  Set_Column=which(names(ESCC_Jupiter06_EGIC)=="FinalStatus"),
                  Set_Comparison_Group=c("EGIC1",
                                         "EGIC2",
                                         "EGIC3"))

pdf("./Figure/Figure5/MolecularEventPlusCI_forrest_PFS_EGIC.pdf",
    width = 10,height = 8)
CombineFinalStatus_HR$PFS
dev.off()
pdf("./Figure/Figure5/MolecularEventPlusCI_forrest_OS_EGIC.pdf",
    width = 10,height = 8)
CombineFinalStatus_HR$OS
dev.off()

##Random Selection(Supplementary)----

# Ransom_Selection_Result70<-
#   data.frame(
#     Seed=1:1000,
#     EGIC1_PFS_HR=NA,
#     EGIC2_PFS_HR=NA,
#     EGIC3_PFS_HR=NA,
#     EGIC1_OS_HR=NA,
#     EGIC2_OS_HR=NA,
#     EGIC3_OS_HR=NA
#   )
# Random_Selection_Result_PFS<-
#   matrix(NA,nrow = 1000,ncol = 9)%>%
#   as.data.frame()
# names(Random_Selection_Result_PFS)<-
#   c("Seed",seq(0.2,0.9,by=0.1))
# Random_Selection_Result_OS<-Random_Selection_Result_PFS
# 
# Random_Selection_Result_Total<-
#   list()
# Random_Selection_Result_Total<-
#   list()
# 
# for(j in 2:9){
#   for(i in 1:1000){
#     set.seed(i)
#     ESCC_JS001_RiskMolEventCombine_CIscore_Select70<-
#       ESCC_JS001_RiskMolEventCombine_CIscore%>%
#       group_by(FinalStatus)%>%
#       sample_frac(size = j/10)%>%
#       as.data.frame()
#     
#     ESCC_JS001_RiskMolEventCombine_CIscore_Select70$Grouping<-
#       factor(ESCC_JS001_RiskMolEventCombine_CIscore_Select70$Grouping,
#              levels = c("Placebo","JS001"))
#     
#     PFS_Summary_IO<-coxph(Surv(`PFS(day)`,`PFS.status.`)~Grouping,
#                           data = subset(ESCC_JS001_RiskMolEventCombine_CIscore_Select70,
#                                         FinalStatus=="EGIC1"))%>%summary() 
#     OS_Summary_IO<-coxph(Surv(`OS(day)`,`OS.status.`)~Grouping,
#                          data = subset(ESCC_JS001_RiskMolEventCombine_CIscore_Select70,
#                                        FinalStatus=="EGIC1"))%>%summary() 
#     Ransom_Selection_Result70$EGIC1_PFS_HR[i]<-
#       PFS_Summary_IO$coefficients[1,2]
#     Ransom_Selection_Result70$EGIC1_OS_HR[i]<-
#       OS_Summary_IO$coefficients[1,2]
#     
#     PFS_Summary_IO<-coxph(Surv(`PFS(day)`,`PFS.status.`)~Grouping,
#                           data = subset(ESCC_JS001_RiskMolEventCombine_CIscore_Select70,
#                                         FinalStatus=="EGIC2"))%>%summary() 
#     OS_Summary_IO<-coxph(Surv(`OS(day)`,`OS.status.`)~Grouping,
#                          data = subset(ESCC_JS001_RiskMolEventCombine_CIscore_Select70,
#                                        FinalStatus=="EGIC2"))%>%summary() 
#     Ransom_Selection_Result70$EGIC2_PFS_HR[i]<-
#       PFS_Summary_IO$coefficients[1,2]
#     Ransom_Selection_Result70$EGIC2_OS_HR[i]<-
#       OS_Summary_IO$coefficients[1,2]
#     
#     PFS_Summary_IO<-coxph(Surv(`PFS(day)`,`PFS.status.`)~Grouping,
#                           data = subset(ESCC_JS001_RiskMolEventCombine_CIscore_Select70,
#                                         FinalStatus=="EGIC3"))%>%summary() 
#     OS_Summary_IO<-coxph(Surv(`OS(day)`,`OS.status.`)~Grouping,
#                          data = subset(ESCC_JS001_RiskMolEventCombine_CIscore_Select70,
#                                        FinalStatus=="EGIC3"))%>%summary()
#     Ransom_Selection_Result70$EGIC3_PFS_HR[i]<-
#       PFS_Summary_IO$coefficients[1,2]
#     Ransom_Selection_Result70$EGIC3_OS_HR[i]<-
#       OS_Summary_IO$coefficients[1,2]
#   }
#   Random_Selection_Result_Total[[as.character(j)]]<-Ransom_Selection_Result70
# }

###PFS

# for(i in names(Random_Selection_Result_Total)[-1]){
#   if(i==3){
#     Random_Selection_Result_Total_PFS<-
#       Random_Selection_Result_Total[[i]][,c(1:4)]
#   }else{
#     Random_Selection_Result_Total_PFS<-
#       cbind(Random_Selection_Result_Total_PFS,
#             Random_Selection_Result_Total[[i]][,c(2:4)])
#   }
# }
# names(Random_Selection_Result_Total_PFS)[2:ncol(Random_Selection_Result_Total_PFS)]<-
#   paste(names(Random_Selection_Result_Total_PFS)[2:ncol(Random_Selection_Result_Total_PFS)],
#         rep(names(Random_Selection_Result_Total)[-1],each=3),sep = "_")
# 
# Random_Selection_Result_Total_PFS<-
#   reshape2::melt(Random_Selection_Result_Total_PFS,
#                  id="Seed")
# 
# Random_Selection_Result_Total_PFS$EGIC<-
#   gsub("_.*","",
#        Random_Selection_Result_Total_PFS$variable)
# Random_Selection_Result_Total_PFS$Fraction<-
#   gsub(".*_","",
#        Random_Selection_Result_Total_PFS$variable)
# Random_Selection_Result_Total_PFS$Fraction<-
#   paste("Fraction =",
#         as.numeric(Random_Selection_Result_Total_PFS$Fraction)/10)
# 
# saveRDS(Random_Selection_Result_Total_PFS,
#         "./SourceData/Random_Selection_EGIC_PFS_HR.RDS")

Random_Selection_Result_Total_PFS<-
  readRDS("./SourceData/InternalValidation/Random_Selection_EGIC_PFS_HR.RDS")
Plot<-
  ggplot(Random_Selection_Result_Total_PFS,
         aes(x=EGIC,fill=EGIC,y=value))+geom_boxplot()+
  facet_grid(~Fraction)+
  scale_fill_manual(values=c("#399579",
                             "#FF9A00",
                             "#3DC1D2"))+
  ylab("PFS - Hazard Ratio")+
  geom_hline(yintercept = 1,linetype="dashed")+
  theme(axis.text = element_text(color="black",
                                 size=15),
        axis.text.x = element_text(angle = 30,hjust = 1),
        axis.title = element_text(size=18),
        axis.title.x = element_blank(),
        strip.text = element_text(size=15),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=15));Plot
ggsave(plot =  Plot,
       filename = "./Figure/Figure5/RandomSelection_EGIC_PFSHR.pdf",
       height = 6,width = 12)

###OS

# for(i in names(Random_Selection_Result_Total)[-1]){
#   if(i==3){
#     Random_Selection_Result_Total_OS<-
#       Random_Selection_Result_Total[[i]][,c(1,5:7)]
#   }else{
#     Random_Selection_Result_Total_OS<-
#       cbind(Random_Selection_Result_Total_OS,
#             Random_Selection_Result_Total[[i]][,c(5:7)])
#   }
# }
# names(Random_Selection_Result_Total_OS)[2:ncol(Random_Selection_Result_Total_OS)]<-
#   paste(names(Random_Selection_Result_Total_OS)[2:ncol(Random_Selection_Result_Total_OS)],
#         rep(names(Random_Selection_Result_Total)[-1],each=3),sep = "_")
# 
# Random_Selection_Result_Total_OS<-
#   reshape2::melt(Random_Selection_Result_Total_OS,
#                  id="Seed")
# 
# Random_Selection_Result_Total_OS$EGIC<-
#   gsub("_.*","",
#        Random_Selection_Result_Total_OS$variable)
# Random_Selection_Result_Total_OS$Fraction<-
#   gsub(".*_","",
#        Random_Selection_Result_Total_OS$variable)
# 
# Random_Selection_Result_Total_OS$Fraction<-
#   paste("Fraction =",
#         as.numeric(Random_Selection_Result_Total_OS$Fraction)/10)
# saveRDS(Random_Selection_Result_Total_OS,
#         "./SourceData/Random_Selection_EGIC_OS_HR.RDS")

Random_Selection_Result_Total_OS<-
  readRDS("./SourceData/InternalValidation/Random_Selection_EGIC_OS_HR.RDS")
Plot<-
  ggplot(Random_Selection_Result_Total_OS,
         aes(x=EGIC,fill=EGIC,y=value))+geom_boxplot()+
  facet_grid(~Fraction)+
  scale_fill_manual(values=c("#399579",
                             "#FF9A00",
                             "#3DC1D2"))+
  ylab("OS - Hazard Ratio")+
  geom_hline(yintercept = 1,linetype="dashed")+
  theme(axis.text = element_text(color="black",
                                 size=15),
        axis.text.x = element_text(angle = 30,hjust = 1),
        axis.title = element_text(size=18),
        axis.title.x = element_blank(),
        strip.text = element_text(size=15),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=15));Plot
ggsave(plot =  Plot,
       filename = "./Figure/Figure5/RandomSelection_EGIC_OSHR.pdf",
       height = 6,width = 12)



##GC mono-cohort for validation----

GC_IO_Clinical_EGICstrategy<-
  readRDS("./SourceData/Mono_IO_Cohort/GC_IO_CombinatorialBiomarker.RDS")%>%
  .[c("WES_ID","EBV_Status","MSI_Status",
      "ORR","EGIC","EGIC_binary")]%>%
  subset(.,!is.na(EGIC))

GC_IO_Clinical_EGIC_ORR_EBVneg_Table<-
  table(subset(GC_IO_Clinical_EGICstrategy,ORR!="UK"&EBV_Status!="Positive")[c("ORR","EGIC")])
GC_IO_Clinical_EGIC_ORR_EBVneg<-
  table(subset(GC_IO_Clinical_EGICstrategy,ORR!="UK"&EBV_Status!="Positive")[c("ORR","EGIC")])%>%
  prop.table(margin = 2)%>%as.data.frame()%>%subset(ORR=="Response")
GC_IO_Clinical_EGIC_ORR_EBVneg_TestResult<-
  table(subset(GC_IO_Clinical_EGICstrategy,ORR!="UK"&EBV_Status!="Positive")[c("ORR","EGIC_binary")])%>%
  fisher.test()

###EBVneg

Plot2<-
  ggplot(GC_IO_Clinical_EGIC_ORR_EBVneg,
         aes(x=EGIC,y=Freq,fill=EGIC))+
  geom_bar(stat = "identity",position = "dodge",width = 0.7)+
  ylab("ORR")+
  scale_fill_manual(values=c("#46927D","#F69230","#39BAC6"))+
  theme_classic()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title = element_text(size=15),
        axis.text = element_text(size=12,color = "black"),
        axis.title.x = element_blank())+
  annotate(geom = "text",x="EGIC1",y=0.25,
           label=paste0("Fisher-Exact Test\nEGIC1/2 vs EGIC3\nP-value=",
                        round(GC_IO_Clinical_EGIC_ORR_EBVneg_TestResult$p.value,digits = 3)))+
  annotate(geom = "text",x="EGIC1",y=GC_IO_Clinical_EGIC_ORR_EBVneg$Freq[1],
           label=paste0(round(GC_IO_Clinical_EGIC_ORR_EBVneg$Freq[1],digits = 3)%>%scales::percent(),
                        "\n(",
                        GC_IO_Clinical_EGIC_ORR_EBVneg_Table[2,1],"/",
                        GC_IO_Clinical_EGIC_ORR_EBVneg_Table[2,1]+GC_IO_Clinical_EGIC_ORR_EBVneg_Table[1,1],
                        ")"))+
  annotate(geom = "text",x="EGIC2",y=GC_IO_Clinical_EGIC_ORR_EBVneg$Freq[2],
           label=paste0(round(GC_IO_Clinical_EGIC_ORR_EBVneg$Freq[2],digits = 3)%>%scales::percent(),
                        "\n(",
                        GC_IO_Clinical_EGIC_ORR_EBVneg_Table[2,2],"/",
                        GC_IO_Clinical_EGIC_ORR_EBVneg_Table[2,2]+GC_IO_Clinical_EGIC_ORR_EBVneg_Table[1,2],
                        ")"))+
  annotate(geom = "text",x="EGIC3",y=GC_IO_Clinical_EGIC_ORR_EBVneg$Freq[3],
           label=paste0(round(GC_IO_Clinical_EGIC_ORR_EBVneg$Freq[3],digits = 3)%>%scales::percent(),
                        "\n(",
                        GC_IO_Clinical_EGIC_ORR_EBVneg_Table[2,3],"/",
                        GC_IO_Clinical_EGIC_ORR_EBVneg_Table[2,3]+GC_IO_Clinical_EGIC_ORR_EBVneg_Table[1,3],
                        ")"));Plot2
ggsave(plot =  Plot2,
       filename = "./Figure/Figure5/GCmonoIO_EBVneg_EGIC_ORR.pdf",
       height = 6,width = 12)

###EBV-neg&MSS

GC_IO_Clinical_EGIC_ORR_EBVnegMSS<-
  table(subset(GC_IO_Clinical_EGICstrategy,ORR!="UK"&EBV_Status!="Positive"&MSI_Status!="MSI")[c("ORR","EGIC")])%>%
  prop.table(margin = 2)%>%as.data.frame()%>%subset(ORR=="Response")  
GC_IO_Clinical_EGIC_ORR_EBVnegMSS_Table<-
  table(subset(GC_IO_Clinical_EGICstrategy,ORR!="UK"&EBV_Status!="Positive"&MSI_Status!="MSI")[c("ORR","EGIC")])
GC_IO_Clinical_EGIC_ORR_EBVnegMSS_TestResult<-
  table(subset(GC_IO_Clinical_EGICstrategy,ORR!="UK"&EBV_Status!="Positive"&MSI_Status!="MSI")[c("ORR","EGIC_binary")])%>%
  fisher.test()

Plot3<-
  ggplot(GC_IO_Clinical_EGIC_ORR_EBVnegMSS,
         aes(x=EGIC,y=Freq,fill=EGIC))+
  geom_bar(stat = "identity",position = "dodge",width = 0.7)+
  ylab("ORR")+
  scale_fill_manual(values=c("#46927D","#F69230","#39BAC6"))+
  theme_classic()+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title = element_text(size=15),
        axis.text = element_text(size=12,color = "black"),
        axis.title.x = element_blank())+
  annotate(geom = "text",x="EGIC1",y=0.25,
           label=paste0("Fisher-Exact Test\nEGIC1/2 vs EGIC3\nP-value=",
                        round(GC_IO_Clinical_EGIC_ORR_EBVnegMSS_TestResult$p.value,digits = 3)))+
  annotate(geom = "text",x="EGIC1",y=GC_IO_Clinical_EGIC_ORR_EBVnegMSS$Freq[1],
           label=paste0(round(GC_IO_Clinical_EGIC_ORR_EBVnegMSS$Freq[1],digits = 3)%>%scales::percent(),
                        "\n(",
                        GC_IO_Clinical_EGIC_ORR_EBVnegMSS_Table[2,1],"/",
                        GC_IO_Clinical_EGIC_ORR_EBVnegMSS_Table[2,1]+GC_IO_Clinical_EGIC_ORR_EBVnegMSS_Table[1,1],
                        ")"))+
  annotate(geom = "text",x="EGIC2",y=GC_IO_Clinical_EGIC_ORR_EBVnegMSS$Freq[2],
           label=paste0(round(GC_IO_Clinical_EGIC_ORR_EBVnegMSS$Freq[2],digits = 3)%>%scales::percent(),
                        "\n(",
                        GC_IO_Clinical_EGIC_ORR_EBVnegMSS_Table[2,2],"/",
                        GC_IO_Clinical_EGIC_ORR_EBVnegMSS_Table[2,2]+GC_IO_Clinical_EGIC_ORR_EBVnegMSS_Table[1,2],
                        ")"))+
  annotate(geom = "text",x="EGIC3",y=GC_IO_Clinical_EGIC_ORR_EBVnegMSS$Freq[3],
           label=paste0(round(GC_IO_Clinical_EGIC_ORR_EBVnegMSS$Freq[3],digits = 3)%>%scales::percent(),
                        "\n(",
                        GC_IO_Clinical_EGIC_ORR_EBVnegMSS_Table[2,3],"/",
                        GC_IO_Clinical_EGIC_ORR_EBVnegMSS_Table[2,3]+GC_IO_Clinical_EGIC_ORR_EBVnegMSS_Table[1,3],
                        ")"));Plot3
ggsave(plot =  Plot3,
       filename = "./Figure/Figure5/GCmonoIO_EBVnegMSS_EGIC_ORR.pdf",
       height = 6,width = 12)
