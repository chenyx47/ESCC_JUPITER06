#Function for Screening----

##Single Variant-Binary Value

Interation_Screening_SV_BV<-
  function(InputData,
           Set_Column,
           Set_Comparison_A,
           Set_Comparison_B){
    Result_Frame<-
      data.frame(
        Items=colnames(InputData)[Set_Column],
        Pc=NA,
        OS.HR.in.Pos=NA,
        OS.Pval.in.Pos=NA,
        OS.HR.in.Neg=NA,
        OS.Pval.in.Neg=NA,
        OS.Inter.k=NA,
        OS.Inter.Pval=NA,
        PFS.HR.in.Pos=NA,
        PFS.Pval.in.Pos=NA,
        PFS.HR.in.Neg=NA,
        PFS.Pval.in.Neg=NA,
        PFS.Inter.k=NA,
        PFS.Inter.Pval=NA,
        OS.HR.Items=NA,
        OS.Pval.Items=NA,
        PFS.HR.Items=NA,
        PFS.Pval.Items=NA,
        OS.HR.Items.JS001=NA,
        OS.Pval.Items.JS001=NA,
        PFS.HR.Items.JS001=NA,
        PFS.Pval.Items.JS001=NA
      )
    for(i in 1:length(Set_Column)){
      InputData_Select<-
        InputData[c("SUBJID","PFS","PFS.status.","OS","OS.status.",
                    "Grouping",Result_Frame$Items[i])]%>%unique()
      #Frequency
      Freq_Table<-
        table(InputData_Select[,Result_Frame$Items[i]])%>%
        prop.table()%>%as.data.frame()
      if(any(Freq_Table$Var1%in%Set_Comparison_A)){
        Result_Frame[i,2]<-
          Freq_Table[,2][Freq_Table[,1]%in%Set_Comparison_A]%>%sum()
        
        names(InputData_Select)[7]<-"Item"
        InputData_Select$Item[
          InputData_Select$Item%in%Set_Comparison_B
        ]<-"Control"
        InputData_Select$Item[
          InputData_Select$Item%in%Set_Comparison_A
        ]<-"Pos"
        InputData_Select<-
          subset(InputData_Select,Item%in%c("Control","Pos"))
        InputData_Select$Item<-
          factor(InputData_Select$Item,
                 levels = c("Control","Pos"))
        InputData_Select$Grouping<-
          factor(InputData_Select$Grouping,
                 levels = c("Chemotherapy+Placebo","Chemotherapy+Toripalimab"))
        if(Result_Frame[i,2]>0.05){
          #OS
          OS_Summary_IO<-coxph(Surv(`OS`,`OS.status.`)~Grouping,
                               data = subset(InputData_Select,Item=="Pos"))%>%summary() 
          Result_Frame[i,3]<-
            OS_Summary_IO$coefficients[1,2]
          Result_Frame[i,4]<-
            OS_Summary_IO$coefficients[1,5]
          
          OS_Summary_Placebo<-coxph(Surv(`OS`,`OS.status.`)~Grouping,
                                    data = subset(InputData_Select,Item=="Control"))%>%summary()
          Result_Frame[i,5]<-
            OS_Summary_Placebo$coefficients[1,2]
          Result_Frame[i,6]<-
            OS_Summary_Placebo$coefficients[1,5]
          
          OS_Summary<-coxph(Surv(`OS`,`OS.status.`)~Item+Grouping+Item*Grouping,
                            data = InputData_Select)%>%summary() 
          Result_Frame[i,7]<-
            OS_Summary$coefficients[3,1]
          Result_Frame[i,8]<-
            OS_Summary$coefficients[3,5]
          
          
          #PFS
          PFS_Summary_IO<-coxph(Surv(`PFS`,`PFS.status.`)~Grouping,
                                data = subset(InputData_Select,Item=="Pos"))%>%summary() 
          Result_Frame[i,9]<-
            PFS_Summary_IO$coefficients[1,2]
          Result_Frame[i,10]<-
            PFS_Summary_IO$coefficients[1,5]
          
          PFS_Summary_Placebo<-coxph(Surv(`PFS`,`PFS.status.`)~Grouping,
                                     data = subset(InputData_Select,Item=="Control"))%>%summary() 
          Result_Frame[i,11]<-
            PFS_Summary_Placebo$coefficients[1,2]
          Result_Frame[i,12]<-
            PFS_Summary_Placebo$coefficients[1,5]
          
          PFS_Summary<-coxph(Surv(`PFS`,`PFS.status.`)~Item+Grouping+Item*Grouping,
                             data = InputData_Select)%>%summary() 
          Result_Frame[i,13]<-
            PFS_Summary$coefficients[3,1]
          Result_Frame[i,14]<-
            PFS_Summary$coefficients[3,5]
          
          #All
          OS_Summary<-coxph(Surv(`OS`,`OS.status.`)~Item,
                            data = InputData_Select)%>%summary() 
          Result_Frame[i,15]<-
            OS_Summary$coefficients[1,2]
          Result_Frame[i,16]<-
            OS_Summary$coefficients[1,5]
          
          PFS_Summary<-coxph(Surv(`PFS`,`PFS.status.`)~Item,
                             data = InputData_Select)%>%summary() 
          Result_Frame[i,17]<-
            PFS_Summary$coefficients[1,2]
          Result_Frame[i,18]<-
            PFS_Summary$coefficients[1,5]
          
          #JS001
          OS_Summary<-coxph(Surv(`OS`,`OS.status.`)~Item,
                            data = subset(InputData_Select,Grouping=="Chemotherapy+Toripalimab"))%>%summary() 
          Result_Frame[i,19]<-
            OS_Summary$coefficients[1,2]
          Result_Frame[i,20]<-
            OS_Summary$coefficients[1,5]
          
          PFS_Summary<-coxph(Surv(`PFS`,`PFS.status.`)~Item,
                             data = subset(InputData_Select,Grouping=="Chemotherapy+Toripalimab"))%>%summary() 
          Result_Frame[i,21]<-
            PFS_Summary$coefficients[1,2]
          Result_Frame[i,22]<-
            PFS_Summary$coefficients[1,5]
        }
      }
    }
    return(Result_Frame)
  }

Interation_Screening_SV_BV_ForKvalue<-
  function(InputData,
           Set_Column,
           Set_Comparison_A,
           Set_Comparison_B){
    Result_Frame<-
      data.frame(
        Items=colnames(InputData)[Set_Column],
        Pc=NA,
        OS.Inter.k=NA,
        OS.Inter.Pval=NA,
        OS.Inter.k.Up=NA,
        OS.Inter.k.Down=NA,
        PFS.Inter.k=NA,
        PFS.Inter.Pval=NA,
        PFS.Inter.k.Up=NA,
        PFS.Inter.k.Down=NA
      )
    for(i in 1:length(Set_Column)){
      InputData_Select<-
        InputData[c("SUBJID","PFS","PFS.status.","OS","OS.status.",
                    "Grouping",Result_Frame$Items[i])]%>%unique()
      #Frequency
      Freq_Table<-
        table(InputData_Select[,Result_Frame$Items[i]])%>%
        prop.table()%>%as.data.frame()
      if(any(Freq_Table$Var1%in%Set_Comparison_A)){
        Result_Frame[i,2]<-
          Freq_Table[,2][Freq_Table[,1]%in%Set_Comparison_A]%>%sum()
        
        names(InputData_Select)[7]<-"Item"
        InputData_Select$Item[
          InputData_Select$Item%in%Set_Comparison_B
        ]<-"Control"
        InputData_Select$Item[
          InputData_Select$Item%in%Set_Comparison_A
        ]<-"Pos"
        InputData_Select<-
          subset(InputData_Select,Item%in%c("Control","Pos"))
        InputData_Select$Item<-
          factor(InputData_Select$Item,
                 levels = c("Control","Pos"))
        InputData_Select$Grouping<-
          factor(InputData_Select$Grouping,
                 levels = c("Chemotherapy+Placebo","Chemotherapy+Toripalimab"))
        if(Result_Frame[i,2]>0.05){
          #OS
          
          OS_Summary<-coxph(Surv(`OS`,`OS.status.`)~Item+Grouping+Item*Grouping,
                            data = InputData_Select)%>%summary() 
          Result_Frame[i,3]<-
            OS_Summary$coefficients[3,1]
          Result_Frame[i,4]<-
            OS_Summary$coefficients[3,5]
          Result_Frame[i,5]<-
            OS_Summary$conf.int[3,4]%>%log()
          Result_Frame[i,6]<-
            OS_Summary$conf.int[3,3]%>%log()
          
          #PFS
          
          PFS_Summary<-coxph(Surv(`PFS`,`PFS.status.`)~Item+Grouping+Item*Grouping,
                             data = InputData_Select)%>%summary() 
          Result_Frame[i,7]<-
            PFS_Summary$coefficients[3,1]
          Result_Frame[i,8]<-
            PFS_Summary$coefficients[3,5]
          Result_Frame[i,9]<-
            PFS_Summary$conf.int[3,4]%>%log()
          Result_Frame[i,10]<-
            PFS_Summary$conf.int[3,3]%>%log()
          
          #All
        }
      }
    }
    return(Result_Frame)
  }

##ForrestPlot----

SubgroupHR_Plot<-
  function(InputData,
           Set_Column,
           Set_Comparison_Group,
           Range=c(0.2,5)){
    ForrestPlot_Data_OS<-
      data.frame(
        Items=c(names(InputData)[Set_Column],
                Set_Comparison_Group),
        `HR (95% CI)`=NA,
        HR_exact=NA,
        Up_CI=NA,
        Low_CI=NA,
        Pvalue=NA,
        Freq=NA
      )
    ForrestPlot_Data_PFS<-
      data.frame(
        Items=c(names(InputData)[Set_Column],
                Set_Comparison_Group),
        `HR (95% CI)`=NA,
        HR_exact=NA,
        Up_CI=NA,
        Low_CI=NA,
        Pvalue=NA,
        Freq=NA
      )
    
    InputData_Select<-
      InputData[c("SUBJID","PFS","PFS.status.","OS","OS.status.",
                  "Grouping",names(InputData)[Set_Column])]%>%unique()
    names(InputData_Select)[7]<-"Item"
    InputData_Select$Grouping<-
      factor(InputData_Select$Grouping,
             levels=c("Chemotherapy+Placebo","Chemotherapy+Toripalimab"))
    ForrestPlot_Data_PFS$Freq[-1]<-
      table(InputData_Select$Item)%>%as.data.frame()%>%FeaturetoRow()%>%
      .[as.character(Set_Comparison_Group),]
    ForrestPlot_Data_OS$Freq[-1]<-
      table(InputData_Select$Item)%>%as.data.frame()%>%FeaturetoRow()%>%
      .[as.character(Set_Comparison_Group),]
    for(i in Set_Comparison_Group){
      PFS_Summary_IO<-coxph(Surv(PFS,`PFS.status.`)~Grouping,
                            data = subset(InputData_Select,Item==i))%>%summary() 
      OS_Summary_IO<-coxph(Surv(OS,`OS.status.`)~Grouping,
                           data = subset(InputData_Select,Item==i))%>%summary() 
      
      ForrestPlot_Data_PFS$HR_exact[which(ForrestPlot_Data_PFS$Items==i)]<-
        PFS_Summary_IO$coefficients[1,2]
      ForrestPlot_Data_PFS$Pvalue[which(ForrestPlot_Data_PFS$Items==i)]<-
        PFS_Summary_IO$coefficients[1,5]%>%
        round(.,digits = 3)
      ForrestPlot_Data_PFS$Up_CI[which(ForrestPlot_Data_PFS$Items==i)]<-
        PFS_Summary_IO$conf.int[1,4]
      ForrestPlot_Data_PFS$Low_CI[which(ForrestPlot_Data_PFS$Items==i)]<-
        PFS_Summary_IO$conf.int[1,3]
      
      ForrestPlot_Data_OS$HR_exact[which(ForrestPlot_Data_OS$Items==i)]<-
        OS_Summary_IO$coefficients[1,2]
      ForrestPlot_Data_OS$Pvalue[which(ForrestPlot_Data_OS$Items==i)]<-
        OS_Summary_IO$coefficients[1,5]%>%
        round(.,digits = 3)
      ForrestPlot_Data_OS$Up_CI[which(ForrestPlot_Data_OS$Items==i)]<-
        OS_Summary_IO$conf.int[1,4]
      ForrestPlot_Data_OS$Low_CI[which(ForrestPlot_Data_OS$Items==i)]<-
        OS_Summary_IO$conf.int[1,3]
    }
    ForrestPlot_Data_PFS$HR..95..CI.<-
      paste(round(ForrestPlot_Data_PFS$HR_exact,digits = 2),"(",
            round(ForrestPlot_Data_PFS$Low_CI,digits = 2)," to ",
            round(ForrestPlot_Data_PFS$Up_CI,digits = 2),")",sep = "")
    ForrestPlot_Data_OS$HR..95..CI.<-
      paste(round(ForrestPlot_Data_OS$HR_exact,digits = 2),"(",
            round(ForrestPlot_Data_OS$Low_CI,digits = 2)," to ",
            round(ForrestPlot_Data_OS$Up_CI,digits = 2),")",sep = "")
    ForrestPlot_Data_PFS[1,c(2,6,7)]<-
      c("HR(95% CI)","P-value","Freq")
    PFS_Plot<-
      forestplot(as.matrix(ForrestPlot_Data_PFS[c(1,2,7,6)]),
                 ForrestPlot_Data_PFS$HR_exact,ForrestPlot_Data_PFS$Low_CI,ForrestPlot_Data_PFS$Up_CI,
                 xlog = TRUE,
                 zero = 1, 
                 xticks = c(Range[1],1,Range[2]),
                 colgap = unit(10,"mm"),
                 graphwidth=unit(50,"mm"),
                 lineheight = unit(1,"cm"),
                 graph.pos = 4,
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
                 hrzl_lines=list("1" = gpar(lwd=1, col="black")))
    
    ForrestPlot_Data_OS[1,c(2,6,7)]<-
      c("HR(95% CI)","P-value","Freq")
    OS_Plot<-
      forestplot(as.matrix(ForrestPlot_Data_OS[c(1,2,7,6)]),
                 ForrestPlot_Data_OS$HR_exact,ForrestPlot_Data_OS$Low_CI,ForrestPlot_Data_OS$Up_CI,
                 xlog = TRUE,
                 zero = 1, 
                 xticks = c(Range[1],1,Range[2]),
                 colgap = unit(10,"mm"),
                 graphwidth=unit(50,"mm"),
                 lineheight = unit(1,"cm"),
                 graph.pos = 4,
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
                 hrzl_lines=list("1" = gpar(lwd=1, col="black")))
    return(list(Data_PFS=ForrestPlot_Data_PFS,
                Data_OS=ForrestPlot_Data_OS,
                PFS=PFS_Plot,
                OS=OS_Plot))
  }

#Plot----
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}
