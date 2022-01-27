# --------------------------------------------------------------------------------------------------
# TASK  01  : SEPARATE FILE FOR EACH SCENARIO - REGRET
# TASK  01B : SEPARATE FILE FOR EACH SCENARIO - ENTROPY
# TASK  02  : DRAW GRAPHS - REGRET
# TASK  02B : DRAW GRAPHS - ENTROPY
# TASK  03  : CATEGORISE SCENARIOS
# TASK  04  : REGRET AND ENTROPY PER SCENARIO CATEGORY
# TASK  05  : FIX CARDINALITY, WHICH METHOD IS BEST, WHAT IS THE REGRET AND ENTROPY
# TASK  06  : FOR REGRET < 10, WHAT IS THE SMALLEST CARDINALITY FOR EACH METHOD PER SCENARIO
# TASK  07  : AVERAGE REGRET OVER ALL SCENARIOS
# TASK  08  : QUADRANTS
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
# TASK  01  : SEPARATE FILE FOR EACH SCENARIO
# --------------------------------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(readr)

# Read MIP file
filename <- "Data_Input/MIP.xls"
sheets_mip <- readxl::excel_sheets(filename)
x_mip <- lapply(sheets_mip, function(X) readxl::read_excel(filename, sheet = X))
names(x_mip) <- sheets_mip
x_mip_a <- average_over_epsilon(x_mip)


# Read TOPK file
filename <- "Data_Input/TopK.xls"
sheets_topk <- readxl::excel_sheets(filename)
x_topk <- lapply(sheets_topk, function(X) readxl::read_excel(filename, sheet = X))
names(x_topk) <- sheets_topk
x_topk_a <- average_over_epsilon(x_topk)
identical(names(x_mip), names(x_topk))
xmip_names <- names(x_mip)


# Read BP file
filename <- "Data_Input/BP.xls"
sheets_bp <- readxl::excel_sheets(filename)
x_bp <- lapply(sheets_bp, function(X) readxl::read_excel(filename, sheet = X))
names(x_bp) <- sheets_bp
x_bp_a <- average_over_epsilon(x_bp) # _a for averaged over epsilon
identical(names(x_mip), names(x_bp))
x_bp_ar <- x_bp_a[xmip_names]


# Read RND file
filename <- "Data_Input/RND.xls"
sheets_rnd <- readxl::excel_sheets(filename)
x_rnd <- lapply(sheets_rnd, function(X) readxl::read_excel(filename, sheet = X))
names(x_rnd) <- sheets_rnd
x_rnd_a <- average_over_epsilon(x_rnd) # _a for averaged over epsilon
x_rnd_ar <- x_rnd_a[xmip_names] # _r for reordered


# Read SB file
filename <- "Data_Input/SB.xls"
sheets_sb <- readxl::excel_sheets(filename)
x_sb <- lapply(sheets_sb, function(X) readxl::read_excel(filename, sheet = X))
names(x_sb) <- sheets_sb
x_sb_a <- average_over_epsilon(x_sb) # _a for averaged over epsilon
x_sb_ar <- x_sb_a[xmip_names] # _r for reordered


# Read SF file
filename <- "Data_Input/SF.xls"
sheets_sf <- readxl::excel_sheets(filename)
x_sf <- lapply(sheets_sf, function(X) readxl::read_excel(filename, sheet = X))
names(x_sf) <- sheets_sf
x_sf_a <- average_over_epsilon(x_sf) # _a for averaged over epsilon
x_sf_ar <- x_sf_a[xmip_names] # _r for reordered


# Read ICARUS file
filename <- "Data_Input/ICARUS.xls"
sheets_icarus <- readxl::excel_sheets(filename)
x_icarus <- lapply(sheets_icarus, function(X) readxl::read_excel(filename, sheet = X))
names(x_icarus) <- sheets_icarus
x_icarus_a <- average_over_epsilon(x_icarus) # _a for averaged over epsilon
x_icarus_ar <- x_icarus_a[xmip_names] # _r for reordered


for(jj in 1:length(xmip_names)){
  if(!is.null(x_mip_a[[jj]])){
    mip_sheet <- x_mip_a[[jj]] %>% select(Cardinality, CV_mean_regret, CV_std_regret) %>% rename("MIP_mean_regret" = "CV_mean_regret", "MIP_std_regret" = "CV_std_regret")
  }else{
    mip_sheet <- tibble(Cardinality = numeric(), MIP_mean_regret = numeric(), MIP_std_regret = numeric())
  }

  if(!is.null(x_topk_a[[jj]])){
    topk_sheet <- x_topk_a[[jj]] %>% select(Cardinality, CV_mean_regret, CV_std_regret) %>% rename("TOPK_mean_regret" = "CV_mean_regret", "TOPK_std_regret" = "CV_std_regret")
  }else{
    topk_sheet <- tibble(Cardinality = numeric(), TOPK_mean_regret = numeric(), TOPK_std_regret = numeric())
  }

  if(!is.null(x_bp_ar[[jj]])){
    bp_sheet <- x_bp_ar[[jj]] %>% select(Cardinality, CV_mean_regret, CV_std_regret) %>% rename("BP_mean_regret" = "CV_mean_regret", "BP_std_regret" = "CV_std_regret")
  }else{
    bp_sheet <- tibble(Cardinality = numeric(), BP_mean_regret = numeric(), BP_std_regret = numeric())
  }

  if(!is.null(x_rnd_ar[[jj]])){
    rnd_sheet <- x_rnd_ar[[jj]] %>% select(Cardinality, CV_mean_regret, CV_std_regret) %>% rename("RND_mean_regret" = "CV_mean_regret", "RND_std_regret" = "CV_std_regret")
  }else{
    rnd_sheet <- tibble(Cardinality = numeric(), RND_mean_regret = numeric(), RND_std_regret = numeric())
  }

  if(!is.null(x_sb_ar[[jj]])){
    sb_sheet <- x_sb_ar[[jj]] %>% select(Cardinality, CV_mean_regret, CV_std_regret) %>% rename("SB_mean_regret" = "CV_mean_regret", "SB_std_regret" = "CV_std_regret")
  }else{
    sb_sheet <- tibble(Cardinality = numeric(), SB_mean_regret = numeric(), SB_std_regret = numeric())
  }

  if(!is.null(x_sf_ar[[jj]])){
    sf_sheet <- x_sf_ar[[jj]] %>% select(Cardinality, CV_mean_regret, CV_std_regret) %>% rename("SF_mean_regret" = "CV_mean_regret", "SF_std_regret" = "CV_std_regret")
  }else{
    sf_sheet <- tibble(Cardinality = numeric(), SF_mean_regret = numeric(), SF_std_regret = numeric())
  }

  if(!is.null(x_icarus_ar[[jj]])){
    icarus_sheet <- x_icarus_ar[[jj]] %>% select(Cardinality, CV_mean_regret, CV_std_regret) %>% rename("ICARUS_mean_regret" = "CV_mean_regret", "ICARUS_std_regret" = "CV_std_regret")
  }else{
    icarus_sheet <- tibble(Cardinality = numeric(), ICARUS_mean_regret = numeric(), ICARUS_std_regret = numeric())
  }

  full_sheet <- full_join(mip_sheet, topk_sheet) %>% full_join(bp_sheet) %>% full_join(rnd_sheet) %>% full_join(sb_sheet) %>% full_join(sf_sheet) %>% full_join(icarus_sheet)
  fname <- paste(xmip_names[[jj]], ".csv", sep="")
  write_csv(full_sheet, paste("Data_Output/Regret_by_Scenario/", fname, sep=""))
}

# --------------------------------------------------------------------------------------------------
# TASK  01B : SEPARATE FILE FOR EACH SCENARIO - ENTROPY
# --------------------------------------------------------------------------------------------------
for(jj in 1:length(xmip_names)){
  if(!is.null(x_mip_a[[jj]])){
    mip_sheet <- x_mip_a[[jj]] %>% select(Cardinality, CV_mean_entropy, CV_std_entropy) %>% rename("MIP_mean_entropy" = "CV_mean_entropy", "MIP_std_entropy" = "CV_std_entropy")
  }else{
    mip_sheet <- tibble(Cardinality = numeric(), MIP_mean_entropy = numeric(), MIP_std_entropy = numeric())
  }

  if(!is.null(x_topk_a[[jj]])){
    topk_sheet <- x_topk_a[[jj]] %>% select(Cardinality, CV_mean_entropy, CV_std_entropy) %>% rename("TOPK_mean_entropy" = "CV_mean_entropy", "TOPK_std_entropy" = "CV_std_entropy")
  }else{
    topk_sheet <- tibble(Cardinality = numeric(), TOPK_mean_entropy = numeric(), TOPK_std_entropy = numeric())
  }

  if(!is.null(x_bp_ar[[jj]])){
    bp_sheet <- x_bp_ar[[jj]] %>% select(Cardinality, CV_mean_entropy, CV_std_entropy) %>% rename("BP_mean_entropy" = "CV_mean_entropy", "BP_std_entropy" = "CV_std_entropy")
  }else{
    bp_sheet <- tibble(Cardinality = numeric(), BP_mean_entropy = numeric(), BP_std_entropy = numeric())
  }

  if(!is.null(x_rnd_ar[[jj]])){
    rnd_sheet <- x_rnd_ar[[jj]] %>%select(Cardinality, CV_mean_entropy, CV_std_entropy) %>% rename("RND_mean_entropy" = "CV_mean_entropy", "RND_std_entropy" = "CV_std_entropy")
  }else{
    rnd_sheet <- tibble(Cardinality = numeric(), RND_mean_entropy = numeric(), RND_std_entropy = numeric())
  }

  if(!is.null(x_sb_ar[[jj]])){
    sb_sheet <- x_sb_ar[[jj]] %>% select(Cardinality, CV_mean_entropy, CV_std_entropy) %>% rename("SB_mean_entropy" = "CV_mean_entropy", "SB_std_entropy" = "CV_std_entropy")
  }else{
    sb_sheet <- tibble(Cardinality = numeric(), SB_mean_entropy = numeric(), SB_std_entropy = numeric())
  }

  if(!is.null(x_sf_ar[[jj]])){
    sf_sheet <- x_sf_ar[[jj]] %>% select(Cardinality, CV_mean_entropy, CV_std_entropy) %>% rename("SF_mean_entropy" = "CV_mean_entropy", "SF_std_entropy" = "CV_std_entropy")
  }else{
    sf_sheet <- tibble(Cardinality = numeric(), SF_mean_entropy = numeric(), SF_std_entropy = numeric())
  }

  if(!is.null(x_icarus_ar[[jj]])){
    icarus_sheet <- x_icarus_ar[[jj]] %>% select(Cardinality, CV_mean_entropy, CV_std_entropy) %>% rename("ICARUS_mean_entropy" = "CV_mean_entropy", "ICARUS_std_entropy" = "CV_std_entropy")
  }else{
    icarus_sheet <- tibble(Cardinality = numeric(), ICARUS_mean_entropy = numeric(), ICARUS_std_entropy = numeric())
  }

  full_sheet <- full_join(mip_sheet, topk_sheet) %>% full_join(bp_sheet) %>% full_join(rnd_sheet) %>% full_join(sb_sheet) %>% full_join(sf_sheet) %>% full_join(icarus_sheet)
  fname <- paste(xmip_names[[jj]], ".csv", sep="")
  write_csv(full_sheet, paste("Data_Output/Entropy_by_Scenario/", fname, sep=""))
}

# --------------------------------------------------------------------------------------------------
# TASK  02  : DRAW GRAPHS - REGRET
# --------------------------------------------------------------------------------------------------
folder <- "Data_Output/Regret_by_Scenario/"
files_list <- list.files(folder)
regret_auc <-  data.frame(MIP = numeric(), TOPK = numeric(), BP = numeric(), RND = numeric(), SB = numeric(), SF = numeric(), ICARUS = numeric())
for(kk in 1:length(files_list)){
  fname <- files_list[kk]
  datori <- read_csv(paste(folder, fname, sep=""))
  dat <- normalize_for_AUC(datori, opt = 2)
  mip_auc <- pracma::trapz(pull(dat, Cardinality), pull(dat, MIP_mean_regret ))
  topk_auc <- pracma::trapz(pull(dat, Cardinality), pull(dat, TOPK_mean_regret ))
  bp_auc <- pracma::trapz(pull(dat, Cardinality), pull(dat, BP_mean_regret ))
  rnd_auc <- pracma::trapz(pull(dat, Cardinality), pull(dat, RND_mean_regret  ))
  sb_auc <- pracma::trapz(pull(dat, Cardinality), pull(dat, SB_mean_regret  ))
  sf_auc <- pracma::trapz(pull(dat, Cardinality), pull(dat, SF_mean_regret  ))
  icarus_auc <- pracma::trapz(pull(dat, Cardinality), pull(dat, ICARUS_mean_regret  ))
  regret_auc[kk, ] <- c(mip_auc, topk_auc, bp_auc, rnd_auc, sb_auc, sf_auc, icarus_auc)
}
sc_names <- substring(files_list, 1, (nchar(files_list) - 4))
regret_auc2 <- regret_auc %>% bind_cols(scenario = sc_names) %>% relocate(scenario)
write_csv(regret_auc2, "Data_Output/Summaries/Regret_AUC.csv")

inds <- grep("ALGO", regret_auc2$scenario)
regret_auc2 <- regret_auc2[-inds, ]

# Graphs
regret_auc3 <- regret_auc2 %>% bind_cols(sc_num = 1:dim(regret_auc2)[1]) %>% relocate(sc_num)
regret_auc3 <- regret_auc3[ ,-c(5, 9)]

library(tidyr)
library(ggplot2)

regret_lng <- pivot_longer(regret_auc3, 3:dim(regret_auc3)[2])
colnames(regret_lng)[3:4] <- c("Meta", "AUC")

one <- paste(1:43)
tt <- unique(regret_lng$scenario)

xlabels <- unique(regret_lng$scenario)
write.table(matrix(as.character(xlabels),nrow=1), sep=",",
            row.names=FALSE, col.names=FALSE)
ones <- 1:33
write.table(matrix(as.character(ones),nrow=1), sep=",",
            row.names=FALSE, col.names=FALSE)
xlabels <- c("ASP-POTASSCO","BNSL-2016","CPMP-2015","CSP-2010","CSP-Minizinc-Obj-2016","CSP-Minizinc-Time-2016","CSP-MZN-2013","GLUHACK-18","GRAPHS-2015","MAXSAT-PMS-2016","MAXSAT-WPMS-2016","MAXSAT12-PMS","MAXSAT15-PMS-INDU","MAXSAT19-UCMS","MIP-2016","OPENML-WEKA-2017","PROTEUS-2014","QBF-2011","QBF-2014","QBF-2016","SAT03-16_INDU","SAT11-HAND","SAT11-INDU","SAT11-RAND","SAT12-ALL","SAT12-HAND","SAT12-INDU","SAT12-RAND","SAT15-INDU","SAT18-EXP","SAT20-MAIN","TSP-LION2015","TTP-2016")

thelabs <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33")


ggplot(regret_lng, aes(x = factor(sc_num), y = AUC, color = Meta, shape = Meta)) + geom_point()  +
 theme_bw()  + scale_x_discrete(labels= xlabels) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Scenarios") + ylab("Regret AUC")

ggplot(regret_lng, aes(x = sc_num, y = AUC, color = Meta, shape = Meta)) + geom_point() + facet_wrap(~Meta) + theme_bw() + ylab("Regret AUC")


# --------------------------------------------------------------------------------------------------
# TASK  02B : DRAW GRAPHS - ENTROPY
# --------------------------------------------------------------------------------------------------
folder <- "Data_Output/Entropy_by_Scenario/"
files_list <- list.files(folder)
entropy_auc <-  data.frame(MIP = numeric(), TOPK = numeric(), BP = numeric(), RND = numeric(), SB = numeric(), SF = numeric(), ICARUS = numeric())
for(kk in 1:length(files_list)){
  fname <- files_list[kk]
  datori <- read_csv(paste(folder, fname, sep=""))
  dat <- normalize_for_AUC(datori, opt = 1)
  dat <- dat[-1, ]
  mip_auc <- pracma::trapz(pull(dat, Cardinality), pull(dat, MIP_mean_entropy ))
  topk_auc <- pracma::trapz(pull(dat, Cardinality), pull(dat, TOPK_mean_entropy ))
  bp_auc <- pracma::trapz(pull(dat, Cardinality), pull(dat, BP_mean_entropy ))
  rnd_auc <- pracma::trapz(pull(dat, Cardinality), pull(dat, RND_mean_entropy  ))
  sb_auc <- pracma::trapz(pull(dat, Cardinality), pull(dat, SB_mean_entropy  ))
  sf_auc <- pracma::trapz(pull(dat, Cardinality), pull(dat, SF_mean_entropy  ))
  icarus_auc <- pracma::trapz(pull(dat, Cardinality), pull(dat, ICARUS_mean_entropy  ))
  entropy_auc[kk, ] <- c(mip_auc, topk_auc, bp_auc, rnd_auc, sb_auc, sf_auc, icarus_auc)
}
sc_names <- substring(files_list, 1, (nchar(files_list) - 4))
entropy_auc2 <- entropy_auc %>% bind_cols(scenario = sc_names) %>% relocate(scenario)
write_csv(entropy_auc2, "Data_Output/Summaries/Entropy_AUC.csv")


inds <- grep("ALGO", entropy_auc2$scenario)
entropy_auc2 <- entropy_auc2[-inds, ]

# Graphs
entropy_auc3 <- entropy_auc2 %>% bind_cols(sc_num = 1:dim(entropy_auc2)[1]) %>% relocate(sc_num)
entropy_auc3 <- entropy_auc3[ ,-c(5, 9)]

library(tidyr)
library(ggplot2)

entropy_lng <- pivot_longer(entropy_auc3, 3:dim(entropy_auc3)[2])
colnames(entropy_lng)[3:4] <- c("Meta", "AUC")

one <- paste(1:43)
tt <- unique(entropy_lng$scenario)

xlabels <- unique(entropy_lng$scenario)
write.table(matrix(as.character(xlabels),nrow=1), sep=",",
            row.names=FALSE, col.names=FALSE)
# ones <- 1:33
# write.table(matrix(as.character(ones),nrow=1), sep=",",
#             row.names=FALSE, col.names=FALSE)
xlabels <- c("ASP-POTASSCO","BNSL-2016","CPMP-2015","CSP-2010","CSP-Minizinc-Obj-2016","CSP-Minizinc-Time-2016","CSP-MZN-2013","GLUHACK-18","GRAPHS-2015","MAXSAT-PMS-2016","MAXSAT-WPMS-2016","MAXSAT12-PMS","MAXSAT15-PMS-INDU","MAXSAT19-UCMS","MIP-2016","OPENML-WEKA-2017","PROTEUS-2014","QBF-2011","QBF-2014","QBF-2016","SAT03-16_INDU","SAT11-HAND","SAT11-INDU","SAT11-RAND","SAT12-ALL","SAT12-HAND","SAT12-INDU","SAT12-RAND","SAT15-INDU","SAT18-EXP","SAT20-MAIN","TSP-LION2015","TTP-2016")


ggplot(entropy_lng, aes(x = factor(sc_num), y = AUC, color = Meta, shape = Meta)) + geom_point()  +
  theme_bw()  + scale_x_discrete(labels= xlabels) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Scenarios") + ylab("Entropy AUC")

ggplot(entropy_lng, aes(x = sc_num, y = AUC, color = Meta, shape = Meta)) + geom_point() + facet_wrap(~Meta) + theme_bw() + ylab("Entropy AUC")


# --------------------------------------------------------------------------------------------------
# TASK  03  : CATEGORISE SCENARIOS
# --------------------------------------------------------------------------------------------------
library(stringr)
scenarios <- unique(entropy_lng$scenario)
stubs <- rep(" ", length(scenarios))
for(i in 1:length(scenarios)){
  pos <- str_locate_all(scenarios[i], "-")[[1]][1]
  stubs[i] <- str_sub(scenarios[i], 1, (pos-1))
}

df <- tibble(scenario = scenarios, short = stubs)
# write_csv(df, "Data_Output/Scenario_Stubs.csv")


# --------------------------------------------------------------------------------------------------
# TASK  04  : REGRET AND ENTROPY PER SCENARIO CATEGORY
# --------------------------------------------------------------------------------------------------
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)

df <- read_csv("Data_Output/Scenario_Stubs.csv")
df2 <- df[ ,c(1,4)]

regret_auc <- read_csv("Data_Output/Summaries/Regret_AUC.csv")
regret_auc <- regret_auc[ ,-c(4,8)]
regret_auc2 <- full_join(regret_auc, df2)
colnames(regret_auc2)[dim(regret_auc2)[2]] <- "category"
inds <- grep("ALGO", regret_auc2$scenario)
regret_auc2 <- regret_auc2[-inds, ]
regret_lng <- pivot_longer(regret_auc2, cols = 2:6)
colnames(regret_lng)[3:4] <- c("Meta", "Regret")

ggplot(regret_lng, aes(x = category , y = Regret, color = category)) + geom_point()  +   theme_bw()  +  ylab("Regret AUC") + facet_wrap(~Meta)  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + xlab("Scenario Category")


# entropy
entropy_auc <- read_csv("Data_Output/Summaries/Entropy_AUC.csv")
entropy_auc <- entropy_auc[ ,-c(4,8)]
entropy_auc2 <- full_join(entropy_auc, df2)
colnames(entropy_auc2)[dim(entropy_auc2)[2]] <- "category"
inds <- grep("ALGO", entropy_auc2$scenario)
entropy_auc2 <- entropy_auc2[-inds, ]
entropy_lng <- pivot_longer(entropy_auc2, cols = 2:6)
colnames(entropy_lng)[3:4] <- c("Meta", "Entropy")


ggplot(entropy_lng, aes(x = category , y = Entropy, color = category)) + geom_point()  +   theme_bw()  +  ylab("Entropy AUC") + facet_wrap(~Meta)  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + xlab("Scenario Category")


entropy_regret <- full_join(entropy_lng, regret_lng)
entropy_regret <- entropy_regret %>% mutate(invRegret = 1/Regret)
ggplot(entropy_regret, aes(x = Entropy, y =Regret, color = category)) + geom_point() + facet_wrap(~Meta) + geom_vline(xintercept = max(entropy_regret$Entropy, na.rm = TRUE)/2, linetype = "dotted") + geom_hline(yintercept = max(entropy_regret$Regret, na.rm = TRUE)/2, linetype = "dotted") + xlab("Entropy AUC") + ylab("Regret AUC")

# --------------------------------------------------------------------------------------------------
# TASK  05  : FIX CARDINALITY, WHICH METHOD IS BEST, WHAT IS THE REGRET AND ENTROPY
# --------------------------------------------------------------------------------------------------

folder1 <- "Data_Output/Regret_by_Scenario/"
folder2 <- "Data_Output/Entropy_by_Scenario/"
files_list <- list.files(folder1)

## -------------------------------------
# CHANGE OPTIONS FROM 1 TO 4
## -------------------------------------
option <- 4
# opt 1 for 25% cardinality
# opt 2 for 50% cardinality
# opt 3 for 75% cardinality
# opt 4 for 90% cardinality

regret_fix_cardinality <-  data.frame(Scenario = character(), Cardinality = numeric(), Best_regret =  numeric(), Best_method_regret = character(), entropy = numeric())


for(kk in 1:length(files_list)){
  fname1 <- files_list[kk]
  datori <- read_csv(paste(folder1, fname1, sep=""))
  # dat <- normalize_for_AUC(datori, opt = 2)


  colinds <- grep("mean", colnames(datori))
  dat <- datori[ ,c(1, colinds)]

  if(option == 1){
    ind1 <- ceiling(dim(dat)[1]/4)
  }else if(option == 2){
    ind1 <- floor(dim(dat)[1]/2)
  }else if(option == 3){
    ind1 <- floor(dim(dat)[1]/4*3)
  }else{
    ind1 <- floor(dim(dat)[1]/10*9)
  }


  row1 <- dat[ind1, -1]
  best_regret1 <- min(row1, na.rm = TRUE)
  method1 <- names(which.min(row1))
  pos <- regexpr("_", method1)
  method_short <- substring(method1, 1, (pos-1))

  datentropy <- read_csv(paste(folder2, fname1, sep=""))
  entropy_method <- paste(substring(method1,1, (nchar(method1)-6)), "entropy", sep="")
  entropy <- datentropy %>% filter(Cardinality == ind1) %>% select(entropy_method)

  pos2 <- regexpr(".csv", fname1)
  scenario <- substring(fname1, 1, (pos2 -1))

  regret_fix_cardinality[kk, 1] <- scenario
  regret_fix_cardinality[kk, 2] <- ind1
  regret_fix_cardinality[kk, 3] <- best_regret1
  regret_fix_cardinality[kk, 4] <- method_short
  regret_fix_cardinality[kk, 5] <- entropy

}

inds <- grep("ALGO", regret_fix_cardinality$Scenario)
regret_fix_cardinality <- regret_fix_cardinality[-inds, ]


if(option == 1){
  write_csv(regret_fix_cardinality, "Data_Output/Summaries/Best_Regret_At_25_Percent_Cardinality.csv")
}else if(option == 2){
  write_csv(regret_fix_cardinality, "Data_Output/Summaries/Best_Regret_At_50_Percent_Cardinality.csv")
}else if(option == 3){
  write_csv(regret_fix_cardinality, "Data_Output/Summaries/Best_Regret_At_75_Percent_Cardinality.csv")
}else{
  write_csv(regret_fix_cardinality, "Data_Output/Summaries/Best_Regret_At_90_Percent_Cardinality.csv")
}


# EXAMINE TABLE OF BEST METHODS
option <- 1

if(option == 1){
  regret_fix_cardinality <- read_csv("Data_Output/Summaries/Best_Regret_At_25_Percent_Cardinality.csv")
}else if(option == 2){
  regret_fix_cardinality <- read_csv("Data_Output/Summaries/Best_Regret_At_50_Percent_Cardinality.csv")
}else if(option == 3){
  regret_fix_cardinality <- read_csv("Data_Output/Summaries/Best_Regret_At_75_Percent_Cardinality.csv")
}else{
  regret_fix_cardinality <- read_csv("Data_Output/Summaries/Best_Regret_At_90_Percent_Cardinality.csv")
}

table(regret_fix_cardinality$Best_method_regret)


# --------------------------------------------------------------------------------------------------
# TASK  06  : FOR REGRET < 10, WHAT IS THE SMALLEST CARDINALITY FOR EACH METHOD PER SCENARIO
# --------------------------------------------------------------------------------------------------
folder1 <- "Data_Output/Regret_by_Scenario/"
folder2 <- "Data_Output/Entropy_by_Scenario/"
files_list <- list.files(folder1)

# regret_less10_cardinality <- data.frame(scenario = character(), MIP = numeric(), TOPK = numeric(), BP = numeric(), RND = numeric(), SB = numeric() )
for(kk in 1:length(files_list)){  #
  fname1 <- files_list[kk]
  datori <- read_csv(paste(folder1, fname1, sep=""))

  colinds <- grep("mean", colnames(datori))
  dat <- datori[ ,c(1, colinds)]

  out <- apply(dat[ ,-1], 2, function(x) min(which(x < 10)))

  pos2 <- regexpr(".csv", fname1)
  scenario <- substring(fname1, 1, (pos2 -1))

  out <- as.data.frame(out)
  out2 <- t(out)

  # ---------------------------
  # Get entropy of these instances
  datori <- read_csv(paste(folder2, fname1, sep=""))
  entropy_colnames <- paste(substring(colnames(out2), 1, nchar(colnames(out2))-6), "entropy", sep="")
  entropy <- rep(0, 7)
  for(ll in 1:7){
    entropy[ll] <- unlist(datori[out2[ 1, ll] ,entropy_colnames[ll]])
  }
  entropy <- data.frame(entropy)
  entropy <- t(entropy)
  colnames(entropy) <- entropy_colnames
  temp <- cbind.data.frame(out2, entropy)

  if(kk ==1){
    regret_less10_cardinality <- as.data.frame(temp) %>% mutate(scenario = scenario, full_dim = dim(dat)[1] )
  }else{
    temp <- as.data.frame(temp) %>% mutate(scenario = scenario, full_dim = dim(dat)[1])
    regret_less10_cardinality <- bind_rows(regret_less10_cardinality, temp)
  }
}

inds <- grep("ALGO", regret_less10_cardinality$scenario)
regret_less10_cardinality <- regret_less10_cardinality[-inds, ]

write_csv(regret_less10_cardinality, "Data_Output/Cardinality_for_Regret_less_than_10.csv")


# --------------------------------------------------------------------------------------------------
# TASK  07  : AVERAGE REGRET OVER ALL SCENARIOS
# --------------------------------------------------------------------------------------------------
folder1 <- "Data_Output/Regret_by_Scenario/"
files_list <- list.files(folder1)

for(kk in 1:length(files_list)){  #
  fname1 <- files_list[kk]
  datori <- read_csv(paste(folder1, fname1, sep=""))
  colinds <- grep("mean", colnames(datori))
  dat <- datori[ ,c(1, colinds)]

  dat <- normalize_for_AUC(dat, opt = 3)
  cols_names <- colnames(dat)[-1]
  regobj <- regexpr("_", cols_names)
  colnames(dat)[-1]  <- substring(cols_names, 1, (regobj-1))
  dfl <- pivot_longer(dat, cols = 2:dim(dat)[2])
  colnames(dfl)[2:3] <- c("Meta", "Regretlog10")

  pos2 <- regexpr(".csv", fname1)
  scenario <- substring(fname1, 1, (pos2 -1))
  temp <- dfl %>% mutate(scenario = scenario)

  if(kk == 1){
    df_all <- temp
  }else{
    df_all <- bind_rows(df_all, temp)
  }
}


ggplot(df_all, aes(x = Cardinality, y= Regretlog10)) + geom_point(aes(color = Meta)) + facet_wrap(~Meta) + geom_smooth()

write_csv(df_all, "Data_Output/Summaries/Log10Regret_over_Percentage_Cardinality_for_Scenarios.csv")

# --------------------------------------------------------------------------------------------------
# TASK  08  : QUADRANTS
# --------------------------------------------------------------------------------------------------

df <- data.frame(x = c(0.15, 0.6, 0.15, 0.6), y = c(1, 1, 3.5, 3.5), descr = c("Similar benchmarks, low regret", "Best quadrant", "Worst quadrant", "Diverse benchmarks, high regret"))
ggplot(df, aes(x = x, y = y, label = descr )) + geom_point(size = 0.00001) + xlim(0, 0.8) + ylim(0, 4.5) + geom_hline(yintercept = 2.2, linetype = "dotted") + geom_vline(xintercept = 0.4, linetype = "dotted") + geom_text() + xlab("Entropy") + ylab("Regret")
