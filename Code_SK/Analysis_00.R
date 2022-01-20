# --------------------------------------------------------------------------------------------------
# TASK  01  : SEPARATE FILE FOR EACH SCENARIO - REGRET
# TASK  01B : SEPARATE FILE FOR EACH SCENARIO - ENTROPY
# TASK  02  : DRAW GRAPHS - REGRET
# TASK  02B : DRAW GRAPHS - ENTROPY
# TASK  03  : CATEGORISE SCENARIOS
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

