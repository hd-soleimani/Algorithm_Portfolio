# --------------------------------------------------------------------------------------------------
# TASK  01  : SEPARATE FILE FOR EACH SCENARIO
# TASK  02  : DRAW GRAPHS
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
# TASK  01  : SEPARATE FILE FOR EACH SCENARIO
# --------------------------------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(readr)

# Read MIP file
filename <- "Data_Input/Results_OPT_211121_Final.xls"
sheets_mip <- readxl::excel_sheets(filename)
x_mip <- lapply(sheets_mip, function(X) readxl::read_excel(filename, sheet = X))
names(x_mip) <- sheets_mip
for (df_name in names(x_mip)) {
  x_mip[[df_name]] <- rename(x_mip[[df_name]],
                             Cardinality  = K, CV_mean_regret = CV_prediction_mean,  CV_mean_entropy = Shanon_eps_0,  CV_std_regret = CV_prediction_STD)
}


# Read TOPK file
filename <- "Data_Input/Results_TopK_211119_Final_eps0.xls"
sheets_topk <- readxl::excel_sheets(filename)
x_topk <- lapply(sheets_topk, function(X) readxl::read_excel(filename, sheet = X))
names(x_topk) <- sheets_topk
for (df_name in names(x_topk)) {
  x_topk[[df_name]] <- rename(x_topk[[df_name]],
                              Cardinality  = K, CV_mean_regret = CV_score_mean,  CV_mean_entropy = Shanon_eps_0,  CV_std_regret = CV_score_std)
}
identical(names(x_mip), names(x_topk))


# Get the sheet to tally
regobj <- regexpr(".csv", names(x_mip))
xmip_names <- substring(names(x_mip), 1, (regobj[1:43]-1))


# Read BP file
filename <- "Data_Input/BP.xls"
sheets_bp <- readxl::excel_sheets(filename)
x_bp <- lapply(sheets_bp, function(X) readxl::read_excel(filename, sheet = X))
names(x_bp) <- sheets_bp
x_bp_a <- average_over_epsilon(x_bp) # _a for averaged over epsilon
x_bp_ar <- x_bp_a[xmip_names] # _r for reordered


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
  if(!is.null(x_mip[[jj]])){
    mip_sheet <- x_mip[[jj]] %>% select(Cardinality, CV_mean_regret, CV_std_regret) %>% rename("MIP_mean_regret" = "CV_mean_regret", "MIP_std_regret" = "CV_std_regret")
  }else{
    mip_sheet <- tibble(Cardinality = numeric(), MIP_mean_regret = numeric(), MIP_std_regret = numeric())
  }

  if(!is.null(x_topk[[jj]])){
    topk_sheet <- x_topk[[jj]] %>% select(Cardinality, CV_mean_regret, CV_std_regret) %>% rename("TOPK_mean_regret" = "CV_mean_regret", "TOPK_std_regret" = "CV_std_regret")
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
  write_csv(full_sheet, paste("Data_Output/Regret/", fname, sep=""))
}

# tt2 <- x_sf[xmip_names]
# setdiff(names(x_mip), names(x_icarus))
# xx <- x[[38]]
# pracma::trapz(pull(xx, K), pull(xx, Regret_train))

# --------------------------------------------------------------------------------------------------
# TASK  02  : DRAW GRAPHS
# --------------------------------------------------------------------------------------------------
folder <- "Data_Output/Regret/"
files_list <- list.files(folder)
regret_auc <-  data.frame(MIP = numeric(), TOPK = numeric(), BP = numeric(), RND = numeric(), SB = numeric(), SF = numeric(), ICARUS = numeric())
for(kk in 1:length(files_list)){
  fname <- files_list[kk]
  dat <- read_csv(paste(folder, fname, sep=""))
  mip_auc <- pracma::trapz(pull(dat, Cardinality), pull(dat, MIP_mean_regret ))
  topk_auc <- pracma::trapz(pull(dat, Cardinality), pull(dat, TOPK_mean_regret ))
  bp_auc <- pracma::trapz(pull(dat, Cardinality), pull(dat, BP_mean_regret ))
  rnd_auc <- pracma::trapz(pull(dat, Cardinality), pull(dat, RND_mean_regret  ))
  sb_auc <- pracma::trapz(pull(dat, Cardinality), pull(dat, SB_mean_regret  ))
  sf_auc <- pracma::trapz(pull(dat, Cardinality), pull(dat, SF_mean_regret  ))
  icarus_auc <- pracma::trapz(pull(dat, Cardinality), pull(dat, ICARUS_mean_regret  ))
  regret_auc[kk, ] <- c(mip_auc, topk_auc, bp_auc, rnd_auc, sb_auc, sf_auc, icarus_auc)
}
write_csv(regret_auc, "Data_Output/Summaries/Regreat_AUC.csv")
