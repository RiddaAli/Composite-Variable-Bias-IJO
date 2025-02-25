## SETUP - Load packages & functions -------------------------------------------
# clear deep clean memory
rm(list=ls()) 
gc(full = TRUE)  

# load required packages
list_of_packages <- c("foreign", "summarytools", "data.table","fastDummies", "ggplot2", "Rfast", "dplyr", "ggpubr")
new_packages     <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
if (packageVersion('dagitty') < "0.2.3") devtools::install_github('jtextor/dagitty/r')
for (itn in 1:length(list_of_packages)) suppressPackageStartupMessages(library(list_of_packages[itn],character.only = T))
rm(list=c("itn", "new_packages", "list_of_packages"))
Plots <- paste0(getwd(),"/Plots/")                                              # sub directory repository for plots
Data  <- paste0(getwd(),"/Data/")                                               # sub directory repository for data

# Collate all relevant model coefficients into Table
Tab_Coeffs <- function(Results, xLabel, Index) {
  Table <- NULL
  for (itn in Index) {
    df_name      <- Results[itn]                                                # get df name
    df           <- get(df_name)                                                # get df
    rownames(df) <- NULL                                                        # set row names to NULL
    colnames(df) <- CI                                                          # rename columns
    assign(df_name, df)                                                         # assign modified data frame to its name
    if (NROW(Index) == 3) Mod <- Mod[1:4]                                       # truncate focus for Table S5
    Table        <- rbind(Table, data.frame(Out = Mod, Exp = xLabel[itn], df))} # collate results into single Table 
  Table <- Table[order(match(Table$Exp, xLabel), match(Table$Out, Mod)), ]      # order Table according to Figure & DAG
  return(Table) }

# Plot function for Figure 1
PlotRTA <- function(df, Mlab, Ylab, xLab, Half = NA) {
  # rescale Sex & Ethnicity for follow-up weight
  df[(df$Exp == "Sex" | df$Exp == "Ethnicity") & df$Out == "Follow-up Weight", CI] <-
    df[(df$Exp == "Sex" | df$Exp == "Ethnicity") & df$Out == "Follow-up Weight", CI] / 10
  # text position
  TxtPos <- ifelse(exists("Pct") && Pct == 1, 2.3, 1.5)
  # full set of mods (all levels of 'Out')
  Mod <- unique(df$Out)
  # Select the desired portion of the legend based on Half argument
  if (!is.na(Half))
    if (Half == 1) { 
      Mods <- Mod[1:ceiling(length(Mod) / 2)]                                   # 1st half of the legend
    } else if (Half == 2) { 
      Mods <- Mod[(ceiling(length(Mod) / 2) + 1):length(Mod)] }                 # 2nd half of the legend
    if (is.na(Half)) Mods <- Mod                                                # all the legend
  # y-axis limits from numeric columns
  Ylim <- range(df %>% select(where(is.numeric)) %>% unlist(), na.rm = TRUE)
  if (!is.na(Half)) Ylim <- c(-1.69, 0.99)
  Title   <- ifelse(is.na(Half), 20, 18)
  Legend  <- ifelse(is.na(Half), 18, 16)
  # construct plot
  p <- ggplot(df, aes(x = factor(Exp, levels = xLab), y = Est, ymin = LCI, 
                      ymax = UCI, color = factor(Out, levels = Mod))) +
    labs(title = Mlab, x = NULL, y = Ylab, colour = NULL) +
    geom_point(position = position_dodge(width = 0.4), size = 3) +
    geom_linerange(position = position_dodge(width = 0.4), linewidth = 1) +
    geom_hline(aes(yintercept = 0), linewidth = 1) + theme_classic() +
    theme(axis.text = element_text(colour = "black", face = "bold", size = 18),
          axis.title = element_text(colour = "black", size = 20, face = "bold"),
          plot.title = element_text(size = Title, face = "bold"), 
          legend.position = "top", legend.text = element_text(size = Legend)) +
    coord_cartesian(ylim = Ylim) + scale_x_discrete(labels = xLab) +
    # show only selected legend
    scale_color_manual(values = setNames(scales::hue_pal()(length(Mod)), Mod), breaks = Mods)
  # add vertical line & annotations if xLab has 6 entries
  if (NROW(xLab) == 6) {
    p <- p + geom_vline(xintercept = 2.5, linetype = "dotted") +
      annotate(geom = "text", x = 1.2, y = TxtPos, label = "Birth variables", 
               color = "black", size = 6, hjust = 0) +
      annotate(geom = "text", x = 4.3, y = TxtPos, label = "Baseline (age 23) variables", 
               color = "black", size = 6, hjust = 0) }
  return(p) }

# write graph to disk in jpeg format
SaveJpg      <- function(plt, name) { 
  ggsave(plt, file = paste0(Plots, name,".jpg"), width = 384, height = 216, units = "mm", limitsize = FALSE) } 

## SETUP - Assign generic information & parameter values -----------------------
# Select whether relative outcomes are absolute or percentage
Pct     <- 1                                                                    # set = 1 for absolute, 100 for percentage

# Setup labels used for summary tables & graphics
Mod     <- c("Follow-up Weight", "Follow-up BMI", "Weight Change", "BMI Change", "Relative Weight Change", "Relative BMI Change")
xLabel  <- c("Ethnicity", "Sex", "Reported Height", "Baseline Weight", "Economic Status", "Malaise Score")
xLabel2 <- c("Ethnicity", "Sex", "Measured Height", "Baseline Weight", "Economic Status", "Malaise Score")
Results <- c("Eth_Results", "Sex_Results", "Ht_Results", "Wt_Results", "Econ_Results", "Mal_Results")
CI      <- c("Est", "LCI", "UCI")
Mlab    <- "Comparison of different outcome variables"
Ylab    <- "Coefficients with 95% confidence intervals"

################################################################################
## MAIN PROGRAM - Read in the data & data cleaning -----------------------------
ncds4                 <- tibble(read.csv(file = paste0(Data, "ncds4.csv"), header = TRUE, sep = ","))
names(ncds4)[1]       <- "ncdsid"
selectedCol           <- subset(ncds4, select = c("ncdsid","n622_4", "dvht23", "dvwt23", "malaise", "econstrg"))
setnames(selectedCol, old = c("ncdsid", "n622_4", "dvht23", "dvwt23", "malaise", "econstrg"), 
         new = c("ncs_id", "Sex", "BaseHt", "BaseWt", "Mal", "Econ"))           # replacing old columns names 
write.csv(selectedCol, file = paste0(Data, "selected_data_ncd_age23.csv"), row.names = FALSE)

NCD_Age23             <- tibble(read.csv(file = paste0(Data, "selected_data_ncd_age23_imputed.csv"), header = TRUE, sep = ","))
setnames(NCD_Age23, old = c("sex", "baseline_height", "baseline_weight", "Malaise_score", "Economic_status"), 
         new = c("Sex", "BaseHt", "BaseWt", "Mal", "Econ"))

ncds_age33            <- tibble(read.csv(file = paste0(Data, "ncds5cmi.csv"), header = TRUE, sep = ","))
names(ncds_age33)[1]  <- "ncdsid"
selectedCols          <- subset(ncds_age33, select = c("ncdsid","n504734", "n504731", "n504654")) 
setnames(selectedCols, old = c("ncdsid", "n504734", "n504731", "n504654"), new =  c("ncs_id", "FUpWt", "FUpHt", "Eth"))
write.csv(selectedCols, file = paste0(Data, "ncds_age33.csv"), row.names = FALSE)

merged_df             <- tibble(merge(NCD_Age23, selectedCols, by = "ncs_id"))
write.csv(merged_df, file = paste0("merged_df.csv"), row.names = FALSE)

# Data cleaning
Merged_df <- tibble(read.csv(file = paste0(Data, "merged_df.csv"), header = TRUE, sep = ","))
Merged_df <- filter(Merged_df, BaseHt != -1)                                    # remove '-1' as missing data
Merged_df <- filter(Merged_df, BaseWt != -1)                                    # ditto
Merged_df <- filter(Merged_df, Mal != -1)                                       # ditto
Merged_df <- mutate(Merged_df, Mal = ifelse(Mal == 1, 1, 0))                    # re-code as 0/1
Merged_df <- filter(Merged_df, Econ != 0)                                       # remove '0' as missing data
Merged_df <- mutate(Merged_df, Econ = ifelse(Econ >= 1 & Econ <= 2, 0, 1))      # re-code as 0/1
Merged_df <- filter(Merged_df, FUpWt >=30 & FUpWt <= 200)
Merged_df <- filter(Merged_df, FUpHt != 0 & FUpHt >= 130 & FUpHt <= 250)
Merged_df <- filter(Merged_df, Eth <= 9.0)
Merged_df <- mutate(Merged_df, Eth = ifelse(Eth == 1, 0, 1))
Merged_df <- mutate(Merged_df, Sex = ifelse(Sex == 2, 0, 1))

# why convert baseline height from meters to centimeters meters? 
Merged_df$FUpHt    <- Merged_df$FUpHt  / 100                                    # convert centimeters to meters for BMI
Merged_df$BaseBMI  <- Merged_df$BaseWt / Merged_df$BaseHt^2
Merged_df$FUpBMI   <- Merged_df$FUpWt  / Merged_df$FUpHt^2
Merged_df$ChBMI    <- Merged_df$FUpBMI - Merged_df$BaseBMI
Merged_df$ChWt     <- Merged_df$FUpWt - Merged_df$BaseWt
Merged_df$RelChBMI <- Pct*(Merged_df$FUpBMI - Merged_df$BaseBMI) / Merged_df$BaseBMI
Merged_df$RelChWt  <- Pct*(Merged_df$FUpWt  - Merged_df$BaseWt)  / Merged_df$BaseWt
Merged_df$BaseHt   <- Merged_df$BaseHt * 100                                    # convert meters to centimeters for Figure 1
Merged_df$FUpHt    <- Merged_df$FUpHt * 100                                     # convert back to meters for Figure 1
Merged_df          <- Merged_df[,-1] |> 
  dplyr::select(Eth, Sex, BaseHt:BaseWt, BaseBMI, Econ, Mal, FUpHt, FUpWt, 
                FUpBMI, ChBMI, RelChBMI, ChWt, RelChWt)                         # order variables as per the DAG
write.csv(Merged_df, paste0(Data, "merged_clean.csv"), row.names = FALSE)

# make the four binary variables factors
merged_clean           <- tibble(read.csv(file = paste0(Data, "merged_clean.csv"), header = TRUE, sep = ","))
cat_cols               <- c("Sex", "Mal", "Econ", "Eth")
merged_clean[cat_cols] <- lapply(merged_clean[cat_cols], factor)

# centre all continuous variables
numeric_cols           <- merged_clean |> dplyr::select(BaseHt:BaseBMI, FUpHt:RelChWt)
centered_numeric_cols  <- tibble(as.data.frame(scale(numeric_cols, center = TRUE, scale = FALSE)))
factor_cols            <- merged_clean |> dplyr::select(Sex, Eth, Econ, Mal)
final_df               <- tibble(cbind(factor_cols, centered_numeric_cols)) |> 
  dplyr::select(Eth, Sex, BaseHt:BaseBMI, Econ:Mal, FUpHt:RelChWt)              # order variables as per the DAG
write.csv(final_df, file = paste0(Data, "final_df.csv"), row.names = FALSE)     # save the data frame

# create data summary for Figure S2 
tmp <- final_df |> rename(Ethnicity = Eth, Baseline_Height = BaseHt, Baseline_Weight = BaseWt, 
                          Baseline_BMI = BaseBMI, Economic_Status = Econ, Malaise_Score = Mal,
                          Follow_Up_Height = FUpHt, Follow_Up_Weight = FUpWt, 
                          Follow_Up_BMI = FUpBMI, Change_Weight = ChWt, Change_BMI = ChBMI, 
                          Relative_Change_Wt = RelChWt, Relative_Change_BMI = RelChBMI)
suppressMessages(view(dfSummary(tmp), file = paste0(Plots, "Final_df.html"), footnote = NA))
suppressMessages(webshot(paste0(Plots, "Final_df.html"), file = paste0(Plots, "Final_df.jpg")))

## MAIN PROGRAM - Exploratory Data Analysis ------------------------------------
Final_df  <- tibble(read.csv(file = paste0(Data, "final_df.csv"), header = TRUE, sep = ","))

p1 <- plot_weight_Eth <- ggplot(Final_df, aes(x = factor(Eth), y = BaseWt, fill = factor(Eth))) + geom_boxplot() +
  ggtitle("Baseline weight vs ethnicity") + xlab("Ethnicity") + ylab("Baseline weight (kilograms)")+
  scale_fill_manual(values = c("deepskyblue1", "deeppink1"),labels = c("Non-white", "White"), name = "Ethnicity") +
  theme(axis.text=element_text(colour="black",size=10),axis.title=element_text(colour="black",size=12),legend.position="none") +
  scale_x_discrete(breaks = c(1,0),labels = c("Non-white", "White"), name = "Ethnicity")
SaveJpg(p1, "Weight_Ethnicity")

p2 <- plot_height_Eth <- ggplot(Final_df, aes(x = factor(Eth), y = BaseHt, fill = factor(Eth)))+ geom_boxplot() +
  ggtitle("Baseline height vs ethnicity") + xlab("Ethnicity") + ylab("Baseline height (centimetres)") +
  scale_fill_manual(values = c("deepskyblue1", "deeppink1"),labels = c("Non-white", "White"), name = "Ethnicity") +
  theme(axis.text=element_text(colour="black",size=10),axis.title=element_text(colour="black",size=12)) +
  scale_x_discrete(breaks = c(1,0), labels = c("Non-white", "White"), name = "Ethnicity")
SaveJpg(p1, "Height_Ethnicity")

plots <- ggarrange(p1, p2, ncol = 2, nrow = 1); SaveJpg(p1, "Weight_Height_Ethnicity")

################################################################################
## MAIN PROGRAM - Run Models for total average causal effects (ACE) ------------
Final_df  <- tibble(read.csv(file = paste0(Data, "final_df.csv"), header = TRUE, sep = ","))

# Ethnicity 
mod1_Eth  <- lm(FUpWt    ~ Eth + Sex, data = Final_df)  # effect of ethnicity score on follow-up weight
mod2_Eth  <- lm(FUpBMI   ~ Eth + Sex, data = Final_df)  # effect of ethnicity on BMI
mod3_Eth  <- lm(ChWt     ~ Eth + Sex, data = Final_df)  # effect of ethnicity on change in weight
mod4_Eth  <- lm(ChBMI    ~ Eth + Sex, data = Final_df)  # effect of ethnicity on change in BMI
mod5_Eth  <- lm(RelChWt  ~ Eth + Sex, data = Final_df)  # effect of ethnicity on relative change in weight
mod6_Eth  <- lm(RelChBMI ~ Eth + Sex, data = Final_df)  # effect of ethnicity on relative change in BMI

# sex 
mod1_Sex  <- lm(FUpWt    ~ Sex + Eth, data = Final_df)  # effect of sex score on follow-up weight
mod2_Sex  <- lm(FUpBMI   ~ Sex + Eth, data = Final_df)  # effect of sex on BMI
mod3_Sex  <- lm(ChWt     ~ Sex + Eth, data = Final_df)  # effect of sex on change in weight
mod4_Sex  <- lm(ChBMI    ~ Sex + Eth, data = Final_df)  # effect of sex on change in BMI
mod5_Sex  <- lm(RelChWt  ~ Sex + Eth, data = Final_df)  # effect of sex on relative change in weight
mod6_Sex  <- lm(RelChBMI ~ Sex + Eth, data = Final_df)  # effect of sex on relative change in BMI

# Self-reported height
mod1_Ht   <- lm(FUpWt    ~ BaseHt + Eth + Sex, data = Final_df)  # effect of Height score on follow-up weight
mod2_Ht   <- lm(FUpBMI   ~ BaseHt + Eth + Sex, data = Final_df)  # effect of Height on BMI
mod3_Ht   <- lm(ChWt     ~ BaseHt + Eth + Sex, data = Final_df)  # effect of Height on change in weight
mod4_Ht   <- lm(ChBMI    ~ BaseHt + Eth + Sex, data = Final_df)  # effect of Height on change in BMI
mod5_Ht   <- lm(RelChWt  ~ BaseHt + Eth + Sex, data = Final_df)  # effect of Height on relative change in weight
mod6_Ht   <- lm(RelChBMI ~ BaseHt + Eth + Sex, data = Final_df)  # effect of Height on relative change in BMI

# Weight 
mod1_Wt   <- lm(FUpWt    ~ BaseWt + BaseHt + Eth + Sex, data = Final_df)  # effect of Weight score on follow-up weight
mod2_Wt   <- lm(FUpBMI   ~ BaseWt + BaseHt + Eth + Sex, data = Final_df)  # effect of Weight on BMI
mod3_Wt   <- lm(ChWt     ~ BaseWt + BaseHt + Eth + Sex, data = Final_df)  # effect of Weight on change in weight
mod4_Wt   <- lm(ChBMI    ~ BaseWt + BaseHt + Eth + Sex, data = Final_df)  # effect of Weight on change in BMI
mod5_Wt   <- lm(RelChWt  ~ BaseWt + BaseHt + Eth + Sex, data = Final_df)  # effect of Weight on relative change in weight
mod6_Wt   <- lm(RelChBMI ~ BaseWt + BaseHt + Eth + Sex, data = Final_df)  # effect of Weight on relative change in BMI

# Economic status 
mod1_Econ <- lm(FUpWt    ~ Econ + Eth + Sex + BaseHt + BaseWt, data = Final_df)  # effect of economic status on follow-up weight
mod2_Econ <- lm(FUpBMI   ~ Econ + Eth + Sex + BaseHt + BaseWt, data = Final_df)  # effect of economic status on BMI
mod3_Econ <- lm(ChWt     ~ Econ + Eth + Sex + BaseHt + BaseWt, data = Final_df)  # effect of economic status on change in weight
mod4_Econ <- lm(ChBMI    ~ Econ + Eth + Sex + BaseHt + BaseWt, data = Final_df)  # effect of economic status on change in BMI
mod5_Econ <- lm(RelChWt  ~ Econ + Eth + Sex + BaseHt + BaseWt, data = Final_df)  # effect of economic status on relative change in weight
mod6_Econ <- lm(RelChBMI ~ Econ + Eth + Sex + BaseHt + BaseWt, data = Final_df)  # effect of economic status on relative change in BMI

# Malaise score
mod1_Mal  <- lm(FUpWt    ~ Mal + Eth + Sex + BaseHt + BaseWt + Econ, data = Final_df) # effect of malaise score on follow-up weight
mod2_Mal  <- lm(FUpBMI   ~ Mal + Eth + Sex + BaseHt + BaseWt + Econ, data = Final_df) # effect of malaise score on BMI
mod3_Mal  <- lm(ChWt     ~ Mal + Eth + Sex + BaseHt + BaseWt + Econ, data = Final_df) # effect of malaise score on change in weight
mod4_Mal  <- lm(ChBMI    ~ Mal + Eth + Sex + BaseHt + BaseWt + Econ, data = Final_df) # effect of malaise score on change in BMI
mod5_Mal  <- lm(RelChWt  ~ Mal + Eth + Sex + BaseHt + BaseWt + Econ, data = Final_df) # effect of malaise score on relative change in weight
mod6_Mal  <- lm(RelChBMI ~ Mal + Eth + Sex + BaseHt + BaseWt + Econ, data = Final_df) # effect of malaise score on relative change in BMI

#Extract coefficients with 95 confidence bounds
Eth_Results   <- rbind(data.frame(Est = coef(mod1_Eth)[2],  t(confint(mod1_Eth)[2,])),
                       data.frame(Est = coef(mod2_Eth)[2],  t(confint(mod2_Eth)[2,])),
                       data.frame(Est = coef(mod3_Eth)[2],  t(confint(mod3_Eth)[2,])),
                       data.frame(Est = coef(mod4_Eth)[2],  t(confint(mod4_Eth)[2,])),
                       data.frame(Est = coef(mod5_Eth)[2],  t(confint(mod5_Eth)[2,])),
                       data.frame(Est = coef(mod6_Eth)[2],  t(confint(mod6_Eth)[2,])))
Sex_Results   <- rbind(data.frame(Est = coef(mod1_Sex)[2],  t(confint(mod1_Sex)[2,])),
                       data.frame(Est = coef(mod2_Sex)[2],  t(confint(mod2_Sex)[2,])),
                       data.frame(Est = coef(mod3_Sex)[2],  t(confint(mod3_Sex)[2,])),
                       data.frame(Est = coef(mod4_Sex)[2],  t(confint(mod4_Sex)[2,])),
                       data.frame(Est = coef(mod5_Sex)[2],  t(confint(mod5_Sex)[2,])),
                       data.frame(Est = coef(mod6_Sex)[2],  t(confint(mod6_Sex)[2,])))
Ht_Results    <- rbind(data.frame(Est = coef(mod1_Ht)[2],   t(confint(mod1_Ht)[2,])),
                       data.frame(Est = coef(mod2_Ht)[2],   t(confint(mod2_Ht)[2,])),
                       data.frame(Est = coef(mod3_Ht)[2],   t(confint(mod3_Ht)[2,])),
                       data.frame(Est = coef(mod4_Ht)[2],   t(confint(mod4_Ht)[2,])),
                       data.frame(Est = coef(mod5_Ht)[2],   t(confint(mod5_Ht)[2,])),
                       data.frame(Est = coef(mod6_Ht)[2],   t(confint(mod6_Ht)[2,])))
Wt_Results    <- rbind(data.frame(Est = coef(mod1_Wt)[2],   t(confint(mod1_Wt)[2,])),
                       data.frame(Est = coef(mod2_Wt)[2],   t(confint(mod2_Wt)[2,])),
                       data.frame(Est = coef(mod3_Wt)[2],   t(confint(mod3_Wt)[2,])),
                       data.frame(Est = coef(mod4_Wt)[2],   t(confint(mod4_Wt)[2,])),
                       data.frame(Est = coef(mod5_Wt)[2],   t(confint(mod5_Wt)[2,])),
                       data.frame(Est = coef(mod6_Wt)[2],   t(confint(mod6_Wt)[2,])))
Econ_Results  <- rbind(data.frame(Est = coef(mod1_Econ)[2], t(confint(mod1_Econ)[2,])),
                       data.frame(Est = coef(mod2_Econ)[2], t(confint(mod2_Econ)[2,])),
                       data.frame(Est = coef(mod3_Econ)[2], t(confint(mod3_Econ)[2,])),
                       data.frame(Est = coef(mod4_Econ)[2], t(confint(mod4_Econ)[2,])),
                       data.frame(Est = coef(mod5_Econ)[2], t(confint(mod5_Econ)[2,])),
                       data.frame(Est = coef(mod6_Econ)[2], t(confint(mod6_Econ)[2,])))
Mal_Results   <- rbind(data.frame(Est = coef(mod1_Mal)[2],  t(confint(mod1_Mal)[2,])),
                       data.frame(Est = coef(mod2_Mal)[2],  t(confint(mod2_Mal)[2,])),
                       data.frame(Est = coef(mod3_Mal)[2],  t(confint(mod3_Mal)[2,])),
                       data.frame(Est = coef(mod4_Mal)[2],  t(confint(mod4_Mal)[2,])),
                       data.frame(Est = coef(mod5_Mal)[2],  t(confint(mod5_Mal)[2,])),
                       data.frame(Est = coef(mod6_Mal)[2],  t(confint(mod6_Mal)[2,])))

# Rename rows & columns of all results data frames & collate into a Table
Table <- Tab_Coeffs(Results, xLabel, 1:6)
write.table(Table, paste0(Data, "Figure 2.csv"), sep = ",", row.names = FALSE)

## MAIN PROGRAM - Generate Figure 2 --------------------------------------------
Table <- read.csv(paste0(Data, "Figure 2.csv"))
p3    <- PlotRTA(Table, Mlab, Ylab, xLabel); SaveJpg(p3, "Figure 2")

################################################################################
## MAIN PROGRAM - Sensitivity analyses for Figure S3 ---------------------------
# Analysis with follow-up measured height instead of the baseline self-reported height
Final_df  <- tibble(read.csv(file = paste0(Data, "final_df.csv"), header = TRUE, sep = ","))

# Ethnicity - same
mod1_Eth  <- lm(FUpWt    ~ Eth + Sex, data = Final_df)  # effect of ethnicity score on follow-up weight
mod2_Eth  <- lm(FUpBMI   ~ Eth + Sex, data = Final_df)  # effect of ethnicity on BMI
mod3_Eth  <- lm(ChWt     ~ Eth + Sex, data = Final_df)  # effect of ethnicity on change in weight
mod4_Eth  <- lm(ChBMI    ~ Eth + Sex, data = Final_df)  # effect of ethnicity on change in BMI
mod5_Eth  <- lm(RelChWt  ~ Eth + Sex, data = Final_df)  # effect of ethnicity on relative change in weight
mod6_Eth  <- lm(RelChBMI ~ Eth + Sex, data = Final_df)  # effect of ethnicity on relative change in BMI

# sex - same
mod1_Sex  <- lm(FUpWt    ~ Sex + Eth, data = Final_df)  # effect of sex score on follow-up weight
mod2_Sex  <- lm(FUpBMI   ~ Sex + Eth, data = Final_df)  # effect of sex on BMI
mod3_Sex  <- lm(ChWt     ~ Sex + Eth, data = Final_df)  # effect of sex on change in weight
mod4_Sex  <- lm(ChBMI    ~ Sex + Eth, data = Final_df)  # effect of sex on change in BMI
mod5_Sex  <- lm(RelChWt  ~ Sex + Eth, data = Final_df)  # effect of sex on relative change in weight
mod6_Sex  <- lm(RelChBMI ~ Sex + Eth, data = Final_df)  # effect of sex on relative change in BMI

# Height measured - different
mod1_Ht   <- lm(FUpWt    ~ FUpHt + Sex + Eth, data = Final_df)  # effect of Height score on follow-up weight
mod2_Ht   <- lm(FUpBMI   ~ FUpHt + Sex + Eth, data = Final_df)  # effect of Height on BMI
mod3_Ht   <- lm(ChWt     ~ FUpHt + Sex + Eth, data = Final_df)  # effect of Height on change in weight
mod4_Ht   <- lm(ChBMI    ~ FUpHt + Sex + Eth, data = Final_df)  # effect of Height on change in BMI
mod5_Ht   <- lm(RelChWt  ~ FUpHt + Sex + Eth, data = Final_df)  # effect of Height on relative change in weight
mod6_Ht   <- lm(RelChBMI ~ FUpHt + Sex + Eth, data = Final_df)  # effect of Height on relative change in BMI

# Weight - different
mod1_Wt   <- lm(FUpWt    ~ BaseWt + FUpHt + Sex + Eth, data = Final_df)  # effect of Weight score on follow-up weight
mod2_Wt   <- lm(FUpBMI   ~ BaseWt + FUpHt + Sex + Eth, data = Final_df)  # effect of Weight on BMI
mod3_Wt   <- lm(ChWt     ~ BaseWt + FUpHt + Sex + Eth, data = Final_df)  # effect of Weight on change in weight
mod4_Wt   <- lm(ChBMI    ~ BaseWt + FUpHt + Sex + Eth, data = Final_df)  # effect of Weight on change in BMI
mod5_Wt   <- lm(RelChWt  ~ BaseWt + FUpHt + Sex + Eth, data = Final_df)  # effect of Weight on relative change in weight
mod6_Wt   <- lm(RelChBMI ~ BaseWt + FUpHt + Sex + Eth, data = Final_df)  # effect of Weight on relative change in BMI

# Economic status - different
mod1_Econ <- lm(FUpWt    ~ Econ + Eth + Sex + FUpHt + BaseWt, data = Final_df)  # effect of economic status on follow-up weight
mod2_Econ <- lm(FUpBMI   ~ Econ + Eth + Sex + FUpHt + BaseWt, data = Final_df)  # effect of economic status on BMI
mod3_Econ <- lm(ChWt     ~ Econ + Eth + Sex + FUpHt + BaseWt, data = Final_df)  # effect of economic status on change in weight
mod4_Econ <- lm(ChBMI    ~ Econ + Eth + Sex + FUpHt + BaseWt, data = Final_df)  # effect of economic status on change in BMI
mod5_Econ <- lm(RelChWt  ~ Econ + Eth + Sex + FUpHt + BaseWt, data = Final_df)  # effect of economic status on relative change in weight
mod6_Econ <- lm(RelChBMI ~ Econ + Eth + Sex + FUpHt + BaseWt, data = Final_df)  # effect of economic status on relative change in BMI

# Malaise score - different
mod1_Mal  <- lm(FUpWt    ~ Mal + Eth + Sex + FUpHt + BaseWt + Econ, data = Final_df) # effect of malaise score on follow-up weight
mod2_Mal  <- lm(FUpBMI   ~ Mal + Eth + Sex + FUpHt + BaseWt + Econ, data = Final_df) # effect of malaise score on BMI
mod3_Mal  <- lm(ChWt     ~ Mal + Eth + Sex + FUpHt + BaseWt + Econ, data = Final_df) # effect of malaise score on change in weight
mod4_Mal  <- lm(ChBMI    ~ Mal + Eth + Sex + FUpHt + BaseWt + Econ, data = Final_df) # effect of malaise score on change in BMI
mod5_Mal  <- lm(RelChWt  ~ Mal + Eth + Sex + FUpHt + BaseWt + Econ, data = Final_df) # effect of malaise score on relative change in weight
mod6_Mal  <- lm(RelChBMI ~ Mal + Eth + Sex + FUpHt + BaseWt + Econ, data = Final_df) # effect of malaise score on relative change in BMI

# Extract coefficients with 95 confidence bounds
Eth_Results   <- rbind(data.frame(Est = coef(mod1_Eth)[2],  t(confint(mod1_Eth)[2,])),
                       data.frame(Est = coef(mod2_Eth)[2],  t(confint(mod2_Eth)[2,])),
                       data.frame(Est = coef(mod3_Eth)[2],  t(confint(mod3_Eth)[2,])),
                       data.frame(Est = coef(mod4_Eth)[2],  t(confint(mod4_Eth)[2,])),
                       data.frame(Est = coef(mod5_Eth)[2],  t(confint(mod5_Eth)[2,])),
                       data.frame(Est = coef(mod6_Eth)[2],  t(confint(mod6_Eth)[2,])))
Sex_Results   <- rbind(data.frame(Est = coef(mod1_Sex)[2],  t(confint(mod1_Sex)[2,])),
                       data.frame(Est = coef(mod2_Sex)[2],  t(confint(mod2_Sex)[2,])),
                       data.frame(Est = coef(mod3_Sex)[2],  t(confint(mod3_Sex)[2,])),
                       data.frame(Est = coef(mod4_Sex)[2],  t(confint(mod4_Sex)[2,])),
                       data.frame(Est = coef(mod5_Sex)[2],  t(confint(mod5_Sex)[2,])),
                       data.frame(Est = coef(mod6_Sex)[2],  t(confint(mod6_Sex)[2,])))
Ht_Results    <- rbind(data.frame(Est = coef(mod1_Ht)[2],   t(confint(mod1_Ht)[2,])),
                       data.frame(Est = coef(mod2_Ht)[2],   t(confint(mod2_Ht)[2,])),
                       data.frame(Est = coef(mod3_Ht)[2],   t(confint(mod3_Ht)[2,])),
                       data.frame(Est = coef(mod4_Ht)[2],   t(confint(mod4_Ht)[2,])),
                       data.frame(Est = coef(mod5_Ht)[2],   t(confint(mod5_Ht)[2,])),
                       data.frame(Est = coef(mod6_Ht)[2],   t(confint(mod6_Ht)[2,])))
Wt_Results    <- rbind(data.frame(Est = coef(mod1_Wt)[2],   t(confint(mod1_Wt)[2,])),
                       data.frame(Est = coef(mod2_Wt)[2],   t(confint(mod2_Wt)[2,])),
                       data.frame(Est = coef(mod3_Wt)[2],   t(confint(mod3_Wt)[2,])),
                       data.frame(Est = coef(mod4_Wt)[2],   t(confint(mod4_Wt)[2,])),
                       data.frame(Est = coef(mod5_Wt)[2],   t(confint(mod5_Wt)[2,])),
                       data.frame(Est = coef(mod6_Wt)[2],   t(confint(mod6_Wt)[2,])))
Econ_Results  <- rbind(data.frame(Est = coef(mod1_Econ)[2], t(confint(mod1_Econ)[2,])),
                       data.frame(Est = coef(mod2_Econ)[2], t(confint(mod2_Econ)[2,])),
                       data.frame(Est = coef(mod3_Econ)[2], t(confint(mod3_Econ)[2,])),
                       data.frame(Est = coef(mod4_Econ)[2], t(confint(mod4_Econ)[2,])),
                       data.frame(Est = coef(mod5_Econ)[2], t(confint(mod5_Econ)[2,])),
                       data.frame(Est = coef(mod6_Econ)[2], t(confint(mod6_Econ)[2,])))
Mal_Results   <- rbind(data.frame(Est = coef(mod1_Mal)[2],  t(confint(mod1_Mal)[2,])),
                       data.frame(Est = coef(mod2_Mal)[2],  t(confint(mod2_Mal)[2,])),
                       data.frame(Est = coef(mod3_Mal)[2],  t(confint(mod3_Mal)[2,])),
                       data.frame(Est = coef(mod4_Mal)[2],  t(confint(mod4_Mal)[2,])),
                       data.frame(Est = coef(mod5_Mal)[2],  t(confint(mod5_Mal)[2,])),
                       data.frame(Est = coef(mod6_Mal)[2],  t(confint(mod6_Mal)[2,])))

# Rename rows & columns of all results data frames & collate into a Table
Table  <- Tab_Coeffs(Results, xLabel, 1:6)
write.table(Table, paste0(Data, "Table S3.csv"), sep = ",", row.names = FALSE)

## MAIN PROGRAM - Generate Figure S3 -------------------------------------------
Table <- read.csv(paste0(Data, "Table S3.csv"))
p4    <- PlotRTA(Table, Mlab, Ylab, xLabel2); SaveJpg(p4, "Figure S3")

################################################################################
## MAIN PROGRAM - Sensitivity analyses for Figure S4 ---------------------------
Final_df  <- tibble(read.csv(file = paste0(Data, "final_df.csv"), header = TRUE, sep = ","))

# Economic status & malaise score ADJUSTING for baseline height & baseline weight
mod1_Econ <- lm(FUpWt    ~ Econ + Eth + Sex + BaseHt + BaseWt,        data = Final_df)  # effect of economic status on follow-up weight
mod2_Econ <- lm(FUpBMI   ~ Econ + Eth + Sex + BaseHt + BaseWt,        data = Final_df)  # effect of economic status on BMI
mod3_Econ <- lm(ChWt     ~ Econ + Eth + Sex + BaseHt + BaseWt,        data = Final_df)  # effect of economic status on change in weight
mod4_Econ <- lm(ChBMI    ~ Econ + Eth + Sex + BaseHt + BaseWt,        data = Final_df)  # effect of economic status on change in BMI
mod5_Econ <- lm(RelChWt  ~ Econ + Eth + Sex + BaseHt + BaseWt,        data = Final_df)  # effect of economic status on relative change in weight
mod6_Econ <- lm(RelChBMI ~ Econ + Eth + Sex + BaseHt + BaseWt,        data = Final_df)  # effect of economic status on relative change in BMI
mod1_Mal  <- lm(FUpWt    ~ Mal +  Eth + Sex + BaseHt + BaseWt + Econ, data = Final_df)  # effect of malaise score on follow-up weight
mod2_Mal  <- lm(FUpBMI   ~ Mal +  Eth + Sex + BaseHt + BaseWt + Econ, data = Final_df)  # effect of malaise score on BMI
mod3_Mal  <- lm(ChWt     ~ Mal +  Eth + Sex + BaseHt + BaseWt + Econ, data = Final_df)  # effect of malaise score on change in weight
mod4_Mal  <- lm(ChBMI    ~ Mal +  Eth + Sex + BaseHt + BaseWt + Econ, data = Final_df)  # effect of malaise score on change in BMI
mod5_Mal  <- lm(RelChWt  ~ Mal +  Eth + Sex + BaseHt + BaseWt + Econ, data = Final_df)  # effect of malaise score on relative change in weight
mod6_Mal  <- lm(RelChBMI ~ Mal +  Eth + Sex + BaseHt + BaseWt + Econ, data = Final_df)  # effect of malaise score on relative change in BMI

# Extract coefficients with 95 confidence bounds
Econ_Results  <- rbind(data.frame(Est = coef(mod1_Econ)[2], t(confint(mod1_Econ)[2,])),
                       data.frame(Est = coef(mod2_Econ)[2], t(confint(mod2_Econ)[2,])),
                       data.frame(Est = coef(mod3_Econ)[2], t(confint(mod3_Econ)[2,])),
                       data.frame(Est = coef(mod4_Econ)[2], t(confint(mod4_Econ)[2,])),
                       data.frame(Est = coef(mod5_Econ)[2], t(confint(mod5_Econ)[2,])),
                       data.frame(Est = coef(mod6_Econ)[2], t(confint(mod6_Econ)[2,])))
Mal_Results   <- rbind(data.frame(Est = coef(mod1_Mal)[2],  t(confint(mod1_Mal)[2,])),
                       data.frame(Est = coef(mod2_Mal)[2],  t(confint(mod2_Mal)[2,])),
                       data.frame(Est = coef(mod3_Mal)[2],  t(confint(mod3_Mal)[2,])),
                       data.frame(Est = coef(mod4_Mal)[2],  t(confint(mod4_Mal)[2,])),
                       data.frame(Est = coef(mod5_Mal)[2],  t(confint(mod5_Mal)[2,])),
                       data.frame(Est = coef(mod6_Mal)[2],  t(confint(mod6_Mal)[2,])))

Table1  <- Tab_Coeffs(Results, xLabel, 5:6)
Left    <- "Adjusted for baseline (self-reported) height and weight" 
p4a     <- PlotRTA(Table1, Left, NULL, xLabel[5:6],1)

# Economic status & malaise score on change_weight WITHOUT ADJUSTING for baseline weight or baseline BMI
mod1_Econ <- lm(FUpWt    ~ Econ + Eth + Sex,        data = Final_df)
mod2_Econ <- lm(FUpBMI   ~ Econ + Eth + Sex,        data = Final_df)
mod3_Econ <- lm(ChWt     ~ Econ + Eth + Sex,        data = Final_df)
mod4_Econ <- lm(ChBMI    ~ Econ + Eth + Sex,        data = Final_df)
mod5_Econ <- lm(RelChWt  ~ Econ + Eth + Sex,        data = Final_df)
mod6_Econ <- lm(RelChBMI ~ Econ + Eth + Sex,        data = Final_df)
mod1_Mal  <- lm(FUpWt    ~ Mal +  Eth + Sex + Econ, data = Final_df)
mod2_Mal  <- lm(FUpBMI   ~ Mal +  Eth + Sex + Econ, data = Final_df)
mod3_Mal  <- lm(ChWt     ~ Mal +  Eth + Sex + Econ, data = Final_df)
mod4_Mal  <- lm(ChBMI    ~ Mal +  Eth + Sex + Econ, data = Final_df)
mod5_Mal  <- lm(RelChWt  ~ Mal +  Eth + Sex + Econ, data = Final_df)
mod6_Mal  <- lm(RelChBMI ~ Mal +  Eth + Sex + Econ, data = Final_df)

#Extract coefficients with 95 confidence bounds
Econ_Results  <- rbind(data.frame(Est = coef(mod1_Econ)[2], t(confint(mod1_Econ)[2,])),
                       data.frame(Est = coef(mod2_Econ)[2], t(confint(mod2_Econ)[2,])),
                       data.frame(Est = coef(mod3_Econ)[2], t(confint(mod3_Econ)[2,])),
                       data.frame(Est = coef(mod4_Econ)[2], t(confint(mod4_Econ)[2,])),
                       data.frame(Est = coef(mod5_Econ)[2], t(confint(mod5_Econ)[2,])),
                       data.frame(Est = coef(mod6_Econ)[2], t(confint(mod6_Econ)[2,])))
Mal_Results   <- rbind(data.frame(Est = coef(mod1_Mal)[2],  t(confint(mod1_Mal)[2,])),
                       data.frame(Est = coef(mod2_Mal)[2],  t(confint(mod2_Mal)[2,])),
                       data.frame(Est = coef(mod3_Mal)[2],  t(confint(mod3_Mal)[2,])),
                       data.frame(Est = coef(mod4_Mal)[2],  t(confint(mod4_Mal)[2,])),
                       data.frame(Est = coef(mod5_Mal)[2],  t(confint(mod5_Mal)[2,])),
                       data.frame(Est = coef(mod6_Mal)[2],  t(confint(mod6_Mal)[2,])))

Table2  <- Tab_Coeffs(Results, xLabel, 5:6)
Right   <- "Not adjusted for baseline (self-reported) height and weight" 
p4b     <- PlotRTA(Table2, Right, NULL, xLabel[5:6], 2)
pS4     <- ggarrange(p4a, p4b, ncol = 2, nrow = 1); SaveJpg(pS4, "Figure S4")

## MAIN PROGRAM - Analyses for Table S5 ----------------------------------------
Final_df  <- tibble(read.csv(file = paste0(Data, "final_df.csv"), header = TRUE, sep = ","))

# Models adjusting for baseline height & baseline weight or adjusting for baseline BMI 
# Ethnicity 
mod1_Eth  <- lm(FUpWt    ~ Eth + Sex +           BaseHt + BaseWt, data = Final_df)
mod3_Eth  <- lm(ChWt     ~ Eth + Sex +           BaseHt + BaseWt, data = Final_df)
mod2_Eth  <- lm(FUpBMI   ~ Eth + Sex + BaseBMI,                   data = Final_df)
mod4_Eth  <- lm(ChBMI    ~ Eth + Sex + BaseBMI,                   data = Final_df)

# Economic status 
mod1_Econ <- lm(FUpWt    ~ Econ + Eth + Sex +           BaseHt + BaseWt, data = Final_df)
mod3_Econ <- lm(ChWt     ~ Econ + Eth + Sex +           BaseHt + BaseWt, data = Final_df)
mod2_Econ <- lm(FUpBMI   ~ Econ + Eth + Sex + BaseBMI,                   data = Final_df)
mod4_Econ <- lm(ChBMI    ~ Econ + Eth + Sex + BaseBMI,                   data = Final_df)

# Malaise score
mod1_Mal  <- lm(FUpWt    ~ Mal + Econ + Eth + Sex +           BaseHt + BaseWt, data = Final_df)
mod3_Mal  <- lm(ChWt     ~ Mal + Econ + Eth + Sex +           BaseHt + BaseWt, data = Final_df)
mod2_Mal  <- lm(FUpBMI   ~ Mal + Econ + Eth + Sex + BaseBMI,                   data = Final_df)
mod4_Mal  <- lm(ChBMI    ~ Mal + Econ + Eth + Sex + BaseBMI,                   data = Final_df)

# Extract coefficients with 95 confidence bounds
Eth_Results   <- rbind(data.frame(Est = coef(mod1_Eth)[2],  t(confint(mod1_Eth)[2,])),
                       data.frame(Est = coef(mod2_Eth)[2],  t(confint(mod2_Eth)[2,])),
                       data.frame(Est = coef(mod3_Eth)[2],  t(confint(mod3_Eth)[2,])),
                       data.frame(Est = coef(mod4_Eth)[2],  t(confint(mod4_Eth)[2,])))
Econ_Results  <- rbind(data.frame(Est = coef(mod1_Econ)[2], t(confint(mod1_Econ)[2,])),
                       data.frame(Est = coef(mod2_Econ)[2], t(confint(mod2_Econ)[2,])),
                       data.frame(Est = coef(mod3_Econ)[2], t(confint(mod3_Econ)[2,])),
                       data.frame(Est = coef(mod4_Econ)[2], t(confint(mod4_Econ)[2,])))
Mal_Results   <- rbind(data.frame(Est = coef(mod1_Mal)[2],  t(confint(mod1_Mal)[2,])),
                       data.frame(Est = coef(mod2_Mal)[2],  t(confint(mod2_Mal)[2,])),
                       data.frame(Est = coef(mod3_Mal)[2],  t(confint(mod3_Mal)[2,])),
                       data.frame(Est = coef(mod4_Mal)[2],  t(confint(mod4_Mal)[2,])))

# Rename rows & columns of all results data frames & collate into a Table
Table <- Tab_Coeffs(Results, xLabel, c(1,5:6))
write.table(Table, paste0(Data, "Table S5.csv"), sep = ",", row.names = FALSE)

