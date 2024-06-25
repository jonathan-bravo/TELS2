if (!require("rjson"))
  install.packages("rjson")

if (!require("lme4"))
  install.packages("lme4")

if (!require("MuMIn"))
  install.packages('MuMIn')

if (!require("dplyr"))
  install.packages('dplyr')

library(rjson)
library(lme4)
library(MuMIn)
library(dplyr)

all_sample_names = c(
  "BFV2AA","BFV2AB","BFV2AC","BFV2AMA","BFV2AMB","BFV2AMC","BFV2MA","BFV2MB","BFV2MC",
  "BFXTAA","BFXTAB","BFXTAC","BFXTAMA","BFXTAMB","BFXTAMC","BFXTMA","BFXTMB","BFXTMC",
  "HFV2AA","HFV2AB","HFV2AC","HFV2AMA","HFV2AMB","HFV2AMC","HFV2MA","HFV2MB","HFV2MC",
  "HFXTAA","HFXTAB","HFXTAC","HFXTAMA","HFXTAMB","HFXTAMC","HFXTMA","HFXTMB","HFXTMC",
  "SV2AMA","SV2AMB","SV2AMC","SV2MA","SV2MB","SV2MC","SXTAA","SXTAB","SXTAC","SXTAMB",
  "SXTAMC","SXTMA","SXTMB","SXTMC")
v2_sample_names = c(
  "BFV2AA","BFV2AB","BFV2AC","BFV2AMA","BFV2AMB","BFV2AMC","BFV2MA","BFV2MB","BFV2MC",
  "HFV2AA","HFV2AB","HFV2AC","HFV2AMA","HFV2AMB","HFV2AMC","HFV2MA","HFV2MB","HFV2MC",
  "SV2AMA","SV2AMB","SV2AMC","SV2MA","SV2MB","SV2MC")
xt_sample_names = c(
  "BFXTAA","BFXTAB","BFXTAC","BFXTAMA","BFXTAMB","BFXTAMC","BFXTMA","BFXTMB","BFXTMC",
  "HFXTAA","HFXTAB","HFXTAC","HFXTAMA","HFXTAMB","HFXTAMC","HFXTMA","HFXTMB","HFXTMC",
  "SXTAA","SXTAB","SXTAC","SXTAMB","SXTAMC","SXTMA","SXTMB","SXTMC")

# Filling data frame function
fill_data_frame <- function(df, count, sample_names, suffix){
  for (index in 1:count) {
    # Chemistry
    
    if (grepl('XT', sample_names[index]))
      df$Chemistry[index] <- 'XT'
    else
      df$Chemistry[index] <- 'V2'
    
    # Organisms
    if (substr(sample_names[index], 1, 2) == 'BF')
      df$SampleType[index] <- 'Bovine'
    else if (substr(sample_names[index], 1, 2) == 'HF')
      df$SampleType[index] <- 'Human'
    else
      df$SampleType[index] <- 'Soil'
    
    # Probe Specificity
    if (grepl('AM', sample_names[index]))
      df$Probe[index] <- 'ARG-MGE'
    else if (grepl('M', sample_names[index]))
      df$Probe[index] <- 'MGE'
    else
      df$Probe[index] <- 'ARG'
    
    # Colocalization Richness
    colocalizations_richness <- read.csv(
      paste("./TELS_output/rarefy/sequel-demultiplex.", 
            paste(
              ".ccs_", suffix,
              ".fastq_deduplicated.fastq_colocalizations_richness.csv",
              sep=''), 
            sep=df$Samples[index]), 
      header=FALSE, row.names=NULL, comment.char="#")
    df$Unique_Coloc[index] = as.integer(colocalizations_richness$V2[1])
    
    # Dedup Read Count
    stats <- read.csv(
      paste("./TELS_output/rarefy/sequel-demultiplex.", 
            paste(
              ".ccs_", suffix,
              ".fastq_deduplicated.fastq_stats.csv",
              sep=''), 
            sep=df$Samples[index]), 
      header=FALSE, row.names=NULL)
    df$Read_Count[index] <- stats$V2[2]
    
    # ARG Richness
    arg <- read.csv(
      paste("./TELS_output/rarefy/sequel-demultiplex.", 
            paste(
              ".ccs_", suffix,
              ".fastq_deduplicated.fastq_SHORT_amr_richness.csv",
              sep=''), 
            sep=df$Samples[index]), 
      row.names=NULL)
    df$ARG[index] <- as.integer(arg$Group.Richness[1])
    
    # MGE Richness
    mge <- read.csv(
      paste("./TELS_output/rarefy/sequel-demultiplex.", 
            paste(
              ".ccs_", suffix,
              ".fastq_deduplicated.fastq_SHORT_mobilome.csv",
              sep=''), 
            sep=df$Samples[index]), 
      row.names=NULL)
    df$MGE[index] <- as.integer(mge$Statistics[18])
    
    # Total Read Length
    readlength <- fromJSON(
      file = paste("./TELS_output/rarefy/sequel-demultiplex.", 
                   paste(
                     ".ccs_", suffix,
                     ".fastq_deduplicated.fastq_reads_length.json",
                     sep=''), 
                   sep=df$Samples[index]))
    df$Total_Read_Length[index] <- sum(unlist(readlength))
    
  }
  return(df)
}

# Step-wise regression
lmer.step <- function(object, steps = 1000, current_step = 0){
  if (current_step == steps)
    return(object)
  fit_change <- drop1(object)
  fit_change <- arrange(fit_change, AIC)
  smallest <- rownames(fit_change)[1]
  if (smallest == "<none>" | smallest == "1")
    return(object)
  new <- update(object, paste(".~. -" , smallest))
  return(lmer.step(new, current_step = current_step+1))
}

# Linear Mixed-Effect Regression for colocalization
lmer_coloc <- function(df, chem) {
  lmer.coloc_grand_fit <- {
    if (chem) {
      lmer(
        Unique_Coloc ~ 
          Probe + SampleType + Chemistry + Total_Read_Length
        #+ Probe:SampleType + Probe:Chemistry + Probe:Total_Read_Length
        #+ SampleType:Chemistry + SampleType:Total_Read_Length
        #+ Total_Read_Length:Chemistry
        + (1 | SampleID),
        data=df)
    }
    else {
      lmer(
        Unique_Coloc ~ 
          Probe + SampleType + Total_Read_Length
        #+ Probe:SampleType + Probe:Total_Read_Length
        #+ SampleType:Total_Read_Length
        + (1 | SampleID),
        data=df)
    }
  }
  return(lmer.step(lmer.coloc_grand_fit))
  #return(lmer.coloc_grand_fit)
}

# Linear Mixed-Effect Regression for ARG
lmer_arg <- function(df, chem) {
  lmer.arg_grand_fit <- {
    if (chem) {
      lmer(
        ARG ~ 
          Probe + SampleType + Chemistry + Total_Read_Length
        #+ Probe:SampleType + Probe:Chemistry + Probe:Total_Read_Length
        #+ SampleType:Chemistry + SampleType:Total_Read_Length
        #+ Total_Read_Length:Chemistry
        + (1 | SampleID),
        data=df)
    }
    else {
      lmer(
        ARG ~ 
          Probe + SampleType + Total_Read_Length
        #+ Probe:SampleType + Probe:Total_Read_Length
        #+ SampleType:Total_Read_Length
        + (1 | SampleID),
        data=df)
    }
  }
  return(lmer.step(lmer.arg_grand_fit))
  #return(lmer.arg_grand_fit)
}

# Linear Mixed-Effect Regression for MGE
lmer_mge <- function(df, chem) {
  lmer.mge_grand_fit <- {
    if (chem) {
      lmer(
        MGE ~ 
          Probe + SampleType + Chemistry + Total_Read_Length
        #+ Probe:SampleType + Probe:Chemistry + Probe:Total_Read_Length
        #+ SampleType:Chemistry + SampleType:Total_Read_Length
        #+ Total_Read_Length:Chemistry
        + (1 | SampleID),
        data=df)
    }
    else {
      lmer(
        MGE ~ 
          Probe + SampleType + Total_Read_Length
        #+ Probe:SampleType + Probe:Total_Read_Length
        #+ SampleType:Total_Read_Length
        + (1 | SampleID),
        data=df)
    }
  }
  return(lmer.step(lmer.mge_grand_fit))
  #return(lmer.mge_grand_fit)
}

# Linear Mixed-Effect Regression for read counts
lmer_read <- function(df, chem) {
  lmer.read_grand_fit <- {
    if (chem) {
      lmer(
        Read_Count ~ 
          Probe + SampleType + Chemistry
        #+ Probe:SampleType + Probe:Chemistry + SampleType:Chemistry
        + (1 | SampleID),
        data=df)
    }
    else {
      lmer(
        Read_Count ~ 
          Probe + SampleType
        #+ Probe:SampleType 
        + (1 | SampleID),
        data=df)
    }
  }
  return(lmer.step(lmer.read_grand_fit))
  #return(lmer.read_grand_fit)
}

plot_coloc <- function(df, lmer, chem) {
  if (chem == 'xt_') {
    plot(x = df$Unique_Coloc, 
         y = round(fitted(lmer)),
         main = "Unique Colocalization",
         xlab = "Rarefied Values",
         ylab = " ",
         xlim = c(min(round(fitted(lmer))), 4),
         ylim = c(min(round(fitted(lmer))), 4),
         las = 1)
  }
  else {
    plot(x = df$Unique_Coloc, 
         y = round(fitted(lmer)),
         main = "Unique Colocalization",
         xlab = "Rarefied Values",
         ylab = " ",
         xlim = c(min(round(fitted(lmer))), 2),
         ylim = c(min(round(fitted(lmer))), 2),
         las = 1)
  }
}

plot_arg  <- function(df, lmer, chem) {
  if (chem == 'xt_') {
    plot(x = df$ARG, 
         y = round(fitted(lmer)),
         main = "LMER Model Fitness Plot for ARG Group Richness ",
         xlab = "Rarefied Values",
         ylab = " ",
         xlim = c(min(round(fitted(lmer))), 200),
         ylim = c(min(round(fitted(lmer))), 200),
         las = 1)
  }
  else if (chem == 'v2_') {
    plot(x = df$ARG, 
         y = round(fitted(lmer)),
         main = "LMER Model Fitness Plot for ARG Group Richness ",
         xlab = "Rarefied Values",
         ylab = " ",
         xlim = c(min(round(fitted(lmer))), 60),
         ylim = c(min(round(fitted(lmer))), 60),
         las = 1)
  }
  else {
    plot(x = df$ARG, 
         y = round(fitted(lmer)),
         main = "LMER Model Fitness Plot for ARG Group Richness ",
         xlab = "ARG group richness values from rarefied samples",
         ylab = " ",
         xlim = c(min(round(fitted(lmer))), 150),
         ylim = c(min(round(fitted(lmer))), 150),
         las = 1)
  }
}

plot_mge <- function(df, lmer, chem) {
  if (chem == 'xt_') {
    plot(x = df$MGE, 
         y = round(fitted(lmer)),
         main = "LMER Model Fitness Plot for MGE Accession Richness",
         xlab = "Rarefied Values",
         ylab = " ",
         xlim = c(min(round(fitted(lmer))), 500),
         ylim = c(min(round(fitted(lmer))), 500),
         las = 1)
  }
  else if (chem == 'v2_'){
    plot(x = df$MGE, 
         y = round(fitted(lmer)),
         main = "LMER Model Fitness Plot for MGE Accession Richness",
         xlab = "Rarefied Values",
         ylab = " ",
         xlim = c(min(round(fitted(lmer))), 60),
         ylim = c(min(round(fitted(lmer))), 60),
         las = 1)
  }
  else {
    plot(x = df$MGE, 
         y = round(fitted(lmer)),
         main = "LMER Model Fitness Plot for MGE Accession Richness",
         xlab = "MGE accession richness values from rarefied samples",
         ylab = " ",
         xlim = c(min(round(fitted(lmer))), 200),
         ylim = c(min(round(fitted(lmer))), 200),
         las = 1)
  }
}

plot_read <- function(df, lmer, chem) {
  if (chem == 'xt_') {
    plot(x = df$Read_Count, 
         y = round(fitted(lmer)),
         main = "Total Reads",
         xlab = "Rarefied Values",
         ylab = " ",
         xlim = c(14900, 15328),
         ylim = c(14900, 15328),
         las = 1)
  }
  else if (chem == 'v2_') {
    plot(x = df$Read_Count, 
         y = round(fitted(lmer)),
         main = "Total Reads",
         xlab = "Rarefied Values",
         ylab = " ",
         xlim = c(5020, 5033),
         ylim = c(5020, 5033),
         las = 1)
  }
  else {
    plot(x = df$Read_Count, 
         y = round(fitted(lmer)),
         main = "Total Reads",
         xlab = "Rarefied Values",
         ylab = " ",
         xlim = c(4990, 5033),
         ylim = c(4990, 5033),
         las = 1)
  }
}

# Write LMER information
write_lmer <- function(
    lmer_function, lmer_filename, rar_filename,
    plot_function, df, chem) {
  lmer <- lmer_function(df, chem)
  lmer.summary <- summary(lmer)
  lmer.random_effects <- ranef(lmer)
  lmer.rsquared <- r.squaredGLMM(lmer)
  lmer.coefficients <- lmer.summary[["coefficients"]]
  lmer.confidence_intervals <- confint(lmer, method='Wald')
  coeff_amt = length(lmer.coefficients) / 3
  df.lmer <- data_frame(
    row.names = row.names(lmer.coefficients),
    'Estimate' = lmer.coefficients[c(1:coeff_amt)],
    'Standard Error' = lmer.coefficients[c((coeff_amt+1):(2*coeff_amt))],
    't-Value' = lmer.coefficients[c((2*coeff_amt+1):(coeff_amt*3))],
    '2.5%' = lmer.confidence_intervals[c(3:(2+coeff_amt))],
    '97.5%' = lmer.confidence_intervals[c((5+coeff_amt):(2*(coeff_amt+2)))])
  write.csv(df.lmer, file=paste(
    output_folder, rar_filename, lmer_filename, 'coefficients.csv'))
  write.csv(lmer.random_effects, file=paste(
    output_folder, rar_filename, lmer_filename, 'rand_effects.csv'))
  write.csv(lmer.rsquared, file=paste(
    output_folder, rar_filename, lmer_filename, 'rsquared.csv'))
  
  png(file=paste(output_folder, rar_filename, lmer_filename, 'model_fit.png'),
      width=600, height=600)
  op <- par(mar = c(5,6,4,2) + 0.1)
  plot_function(df, lmer, rar_filename)
  title(ylab = "Fitted values from LMER model", cex.lab = 1, line = 4.5)
  par(op)
  dev.off()
  data = data.frame(row.names = paste(rar_filename, lmer_filename),
                    'AIC' = extractAIC(lmer)[2],
                    'edf' = extractAIC(lmer)[1])
  write.table(data, file=paste(output_folder, 'AIC.csv'), append=TRUE)
  return(lmer.summary[["AICtab"]][["REML"]])
}

v2_count = 24
v2_df=data.frame(
  Samples=c(
    "BFV2AA","BFV2AB","BFV2AC","BFV2AMA","BFV2AMB","BFV2AMC","BFV2MA","BFV2MB","BFV2MC",
    "HFV2AA","HFV2AB","HFV2AC","HFV2AMA","HFV2AMB","HFV2AMC","HFV2MA","HFV2MB","HFV2MC",
    "SV2AMA","SV2AMB","SV2AMC","SV2MA","SV2MB","SV2MC"),
  SampleID=c(
    "BFV2A","BFV2A","BFV2A","BFV2AM","BFV2AM","BFV2AM","BFV2M","BFV2M","BFV2M",
    "HFV2A","HFV2A","HFV2A","HFV2AM","HFV2AM","HFV2AM","HFV2M","HFV2M","HFV2M",
    "SV2AM","SV2AM","SV2AM","SV2M","SV2M","SV2M"),
  Chemistry = integer(v2_count),
  SampleType = integer(v2_count),
  Probe = integer(v2_count),
  Read_Count = integer(v2_count),
  Total_Read_Length = integer(v2_count),
  Unique_Coloc = integer(v2_count),
  ARG = integer(v2_count),
  MGE = integer(v2_count)
)
v2_df <- fill_data_frame(v2_df, v2_count, v2_sample_names, 'chem')

xt_count = 26
xt_df=data.frame(
  Samples=c(
    "BFXTAA","BFXTAB","BFXTAC","BFXTAMA","BFXTAMB","BFXTAMC","BFXTMA","BFXTMB","BFXTMC",
    "HFXTAA","HFXTAB","HFXTAC","HFXTAMA","HFXTAMB","HFXTAMC","HFXTMA","HFXTMB","HFXTMC",
    "SXTAA","SXTAB","SXTAC","SXTAMB","SXTAMC","SXTMA","SXTMB","SXTMC"),
  SampleID=c(
    "BFXTA","BFXTA","BFXTA","BFXTAM","BFXTAM","BFXTAM","BFXTM","BFXTM","BFXTM",
    "HFXTA","HFXTA","HFXTA","HFXTAM","HFXTAM","HFXTAM","HFXTM","HFXTM","HFXTM",
    "SXTA","SXTA","SXTA","SXTAM","SXTAM","SXTM","SXTM","SXTM"),
  Chemistry = integer(xt_count),
  SampleType = integer(xt_count),
  Probe = integer(xt_count),
  Read_Count = integer(xt_count),
  Total_Read_Length = integer(xt_count),
  Unique_Coloc = integer(xt_count),
  ARG = integer(xt_count),
  MGE = integer(xt_count)
)
xt_df <- fill_data_frame(xt_df, xt_count, xt_sample_names, 'chem')

all_count = 50
all_df=data.frame(
  Samples=c(
    "BFV2AA","BFV2AB","BFV2AC","BFV2AMA","BFV2AMB","BFV2AMC","BFV2MA","BFV2MB","BFV2MC",
    "BFXTAA","BFXTAB","BFXTAC","BFXTAMA","BFXTAMB","BFXTAMC","BFXTMA","BFXTMB","BFXTMC",
    "HFV2AA","HFV2AB","HFV2AC","HFV2AMA","HFV2AMB","HFV2AMC","HFV2MA","HFV2MB","HFV2MC",
    "HFXTAA","HFXTAB","HFXTAC","HFXTAMA","HFXTAMB","HFXTAMC","HFXTMA","HFXTMB","HFXTMC",
    "SV2AMA","SV2AMB","SV2AMC","SV2MA","SV2MB","SV2MC","SXTAA","SXTAB","SXTAC","SXTAMB",
    "SXTAMC","SXTMA","SXTMB","SXTMC"),
  SampleID=c(
    "BFV2A","BFV2A","BFV2A","BFV2AM","BFV2AM","BFV2AM","BFV2M","BFV2M","BFV2M",
    "BFXTA","BFXTA","BFXTA","BFXTAM","BFXTAM","BFXTAM","BFXTM","BFXTM","BFXTM",
    "HFV2A","HFV2A","HFV2A","HFV2AM","HFV2AM","HFV2AM","HFV2M","HFV2M","HFV2M",
    "HFXTA","HFXTA","HFXTA","HFXTAM","HFXTAM","HFXTAM","HFXTM","HFXTM","HFXTM",
    "SV2AM","SV2AM","SV2AM","SV2M","SV2M","SV2M","SXTA","SXTA","SXTA","SXTAM",
    "SXTAM","SXTM","SXTM","SXTM"),
  Chemistry = integer(all_count),
  SampleType = integer(all_count),
  Probe = integer(all_count),
  Read_Count = integer(all_count),
  Total_Read_Length = integer(all_count),
  Unique_Coloc = integer(all_count),
  ARG = integer(all_count),
  MGE = integer(all_count)
)
all_df <- fill_data_frame(all_df, all_count, all_sample_names, 'all')

output_folder <- './output/regression_results/actual/step_rarefy/'

if (file.exists(output_folder) == FALSE)
  dir.create(output_folder)

lmer_function_list <- c(lmer_coloc, lmer_arg, lmer_mge, lmer_read)
plot_function_list <- c(plot_coloc, plot_arg, plot_mge, plot_read)
lmer_filename_list <- c("coloc_", "arg_", "mge_", "read_counts_")
rar_filename_list <- c("all_", "v2_", "xt_")
df_list <- list(all_df, v2_df, xt_df)

summary <- list()
name <- list()

for (i in 1:4) {
  for (j in 1:3) {
    if (i == 1 && j == 2) {
      next
    }
    summary <- append(summary, write_lmer(
      lmer_function_list[[i]], lmer_filename_list[i], rar_filename_list[j],
      plot_function_list[[i]], df_list[[j]], j < 2))
    name <- append(name, paste(rar_filename_list[j], lmer_filename_list[i]))
  }
}
df.reml <- data_frame(
  row.names = array(unlist(name)),
  'REML' = array(unlist(summary))
)

write.csv(df.reml, file=paste(output_folder, 'REML.csv'))