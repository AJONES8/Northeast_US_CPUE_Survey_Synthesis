# =========================================================================
# CPUE ANALYSIS PIPELINE: TRACKING, ASYMMETRY, & HYPERSTABILITY
# =========================================================================
library(tidyverse)
library(patchwork)
library(ggthemes) 
library(ggrepel)
# Load required package for Deming regression
library(deming)
# Load required package for AR(1) generalized least squares
library(nlme)


# =========================================================================
# 0. Functions outside of packages
# =========================================================================

#Function to 
get_ar1_pval <- function(x, y) {
  if (length(x) < 5) return(NA_real_)
  tryCatch({
    mod <- gls(y ~ x, correlation = corAR1(form = ~ 1))
    summary(mod)$tTable[2, "p-value"]
  }, error = function(e) {
    NA_real_ # Returns NA if the AR(1) model fails to converge
  })
}

# Safe helper function to extract AR-adjusted p-value using the Chelton method
# v1 with gls
get_ar1_pval <- function(x, y) {
  if (length(x) < 5) return(NA_real_)
  tryCatch({
    mod <- gls(y ~ x, correlation = corAR1(form = ~ 1))
    summary(mod)$tTable[2, "p-value"]
  }, error = function(e) {
    NA_real_ # Returns NA if the AR(1) model fails to converge
  })
}

#v2 Modified per Pyper and Peterman (1998)
get_chelton_pval <- function(x, y, method = "bartlett") {
  # Retain only complete pairs
  valid <- complete.cases(x, y)
  x <- x[valid]
  y <- y[valid]
  N <- length(x)
  
  # Require minimum length for meaningful autocorrelation analysis
  if (N < 5) return(NA_real_)
  
  # Set maximum lag (k) to N/5 as recommended for fisheries time series
  k <- floor(N / 5)
  if (k < 1) k <- 1
  
  # Extract autocorrelations (index 2 corresponds to lag-1)
  rx <- acf(x, lag.max = k, plot = FALSE, na.action = na.pass)$acf[2:(k+1)]
  ry <- acf(y, lag.max = k, plot = FALSE, na.action = na.pass)$acf[2:(k+1)]
  
  # Apply the selected weighting function
  lags <- 1:k
  if (method == "chelton") {
    W <- rep(1, k) # Uniform weighting (Chelton 1984)
  } else {
    W <- (N - lags) / N # Linearly decreasing weighting (Bartlett 1946)
  }
  
  # Calculate effective sample size (N_eff)
  sum_lags <- sum(W * rx * ry)
  N_eff <- N / (1 + 2 * sum_lags)
  
  # Constrain N_eff bounds (Pyper & Peterman 1998 modification)
  N_eff <- min(N, max(3, N_eff)) 
  
  # Calculate standard Pearson correlation
  r_val <- cor(x, y, method = "pearson")
  
  # Calculate adjusted significance
  df_adj <- N_eff - 2
  t_adj <- r_val * sqrt(df_adj / (1 - r_val^2))
  p_adj <- 2 * (1 - pt(abs(t_adj), df = df_adj))
  
  return(p_adj)
}

# =========================================================================
# 1. METADATA & COLOR CONFIGURATION (Z-SCORED DATA)
# =========================================================================
# Used for Tracking (Tau) and Asymmetry analysis
parent_meta <- list(
  list(file = "atlantic_cod_all_indices_combined.csv", name = "Cod (EGOM)", cpue = "LPUE_EGOM_zscored", survey = c("NMFS_Fall_EGOM_zscored", "NMFS_Spring_EGOM_zscored")),
  list(file = "atlantic_cod_all_indices_combined.csv", name = "Cod (WGOM)", cpue = "LPUE_WGOM_Combined_zscored", survey = c("NMFS_Fall_WGOM_zscored", "NMFS_Spring_WGOM_zscored")),
  list(file = "atlantic_cod_all_indices_combined.csv", name = "Cod (GB)", cpue = "LPUE_GB_zscored", survey = c("NMFS_Fall_GB_zscored", "NMFS_Spring_GB_zscored")),
  list(file = "atlantic_cod_all_indices_combined.csv", name = "Cod (SNE)", cpue = "LPUE_SNE_Commercial_zscored", survey = c("NMFS_Fall_SNE_zscored", "NMFS_Spring_SNE_zscored")),
  list(file = "monkfish_final.csv", name = "Monkfish (N)", cpue = c('Memo_N_Gillnet_Centered', 'Memo_N_Trawl_Centered', 'Memo_N_Dredge_Centered', 'Rich_N_Trawl_Logbook_Centered', 'Rich_N_Trawl_Observer_Centered'), survey = c('NMFS_Fall_North_Survey', 'NMFS_Spring_North_Survey')),
  list(file = "monkfish_final.csv", name = "Monkfish (S)", cpue = c('Memo_S_Gillnet_Centered', 'Memo_S_Trawl_Centered', 'Memo_S_Dredge_Centered', 'Rich_S_Gillnet_Logbook_Centered', 'Rich_S_Gillnet_Observer_Centered'), survey = c('NMFS_Fall_South_Survey', 'NMFS_Spring_South_Survey')),
  list(file = "yellowtail_flounder.csv", name = "Yellowtail (CC/GOM)", cpue = "CCGOM_LPUE_ZScaled", survey = "CCGOM_SpringSurvey_ZScaled"),
  list(file = "yellowtail_flounder.csv", name = "Yellowtail (GB)", cpue = "GB_LPUE_ZScaled", survey = "GB_SpringSurvey_ZScaled"),
  list(file = "yellowtail_flounder.csv", name = "Yellowtail (SNE/MA)", cpue = "SNEMA_LPUE_ZScaled", survey = "SNEMA_SpringSurvey_ZScaled"),
  list(file = "black_sea_bass_fresh_zscored_dataset.csv", name = "Black Sea Bass (N)", cpue = c('North_GAM_CPUE_zscore', 'North_Nominal_CPUE_zscore', 'North_Rec_CPUE_zscore'), survey = 'North_Survey_CPUE_zscore'),
  list(file = "black_sea_bass_fresh_zscored_dataset.csv", name = "Black Sea Bass (S)", cpue = c('South_GAM_CPUE_zscore', 'South_Nominal_CPUE_zscore', 'South_Rec_CPUE_zscore'), survey = 'South_Survey_CPUE_zscore'),
  list(file = "american_plaice_fully_updated.csv", name = "American Plaice", cpue = c('WP9_Terceiro_Dealer_LPUE', 'WP7_Hankowsky_Logbook_LPUE', 'WP8_Jones_SpatioTemporal_CPUE', 'WP8_Jones_SpatialOnly_CPUE'), survey = c('NMFS_Spring_Survey_Est', 'NMFS_Fall_Survey_Est')),
  list(file = "bluefish_indices_trim.csv", name = "Bluefish", cpue = c('Z_MRIP_Guild', 'Z_MRIP_Directed'), survey = c('Z_Albatross', 'Z_Bigelow', 'Z_NEAMAP')),
  list(file = "gom_haddock_all_centered_indices_trim.csv", name = "Haddock (GOM)", cpue = c('Z_Logbook_LPUE', 'Z_StudyFleet_CPUE'), survey = c('Z_NMFS_Fall_Bio', 'Z_NMFS_Spring_Bio')),
  list(file = "scup_centered_indices_trim.csv", name = "Scup", cpue = 'Z_MRIP_CPUE_Points', survey = c('Z_Fall_ALB', 'Z_Spring_ALB', 'Z_Winter', 'Z_Spring_BIG', 'Z_Fall_BIG')),
  list(file = "illex_squid_lpue_hand.csv", name = "Shortfin Squid", cpue = c('Dealer_WetBoats', 'Observer_WetBoats', 'StudyFleet_WetBoats'), survey = 'NEFSC_Fall'),
  list(file = "Atlantic_halibut_values.csv", name = "Halibut", cpue = 'Standardized_CPUE_Est', survey = c('NEFSC_Fall_Survey_Est', 'NEFSC_Spring_Survey_Est')),
  list(file = "golden_tilefish_final_indices.csv", name = "Golden Tilefish", cpue = "GLM_LPUE", survey = c("Survey_Abundance_zscore", "Survey_Biomass_zscore"))
)

# Extract and Alphabetize Stock Names
all_stocks <- sort(unique(unlist(lapply(parent_meta, function(x) x$name))))

# Lock Tableau 20 Colors Alphabetically
tableau_colors <- ggthemes::tableau_color_pal("Tableau 20")(length(all_stocks))
fixed_stock_colors <- setNames(tableau_colors, all_stocks)

# =========================================================================
# 2. DATA PROCESSING: DIRECTIONAL TRACKING (TAU, PEARSON & VOLATILITY)
# =========================================================================
print("Processing tracking data...")
pairs_list <- list(); stock_list <- list(); long_list <- list()

for (m in parent_meta) {
  if (!file.exists(m$file)) next
  dat <- read_csv(m$file, show_col_types = FALSE)
  if ("YEAR" %in% names(dat)) dat <- dat %>% rename(Year = YEAR)
  
  c_cols <- intersect(m$cpue, names(dat)); s_cols <- intersect(m$survey, names(dat))
  
  if(length(c_cols) > 0 & length(s_cols) > 0) {
    dat_clean <- dat %>%
      mutate(Mean_C = rowMeans(dplyr::select(., dplyr::all_of(c_cols)), na.rm = TRUE),
             Mean_S = rowMeans(dplyr::select(., dplyr::all_of(s_cols)), na.rm = TRUE),
             Stock = m$name)
    
    stock_df <- dat_clean %>% dplyr::select(Year, Mean_C, Mean_S, Stock) %>% drop_na()
    if(nrow(stock_df) > 4) {
      stock_list[[m$name]] <- stock_df %>% 
        dplyr::group_by(Stock) %>%
        dplyr::summarize(
          Tau = cor(Mean_C, Mean_S, method = "kendall"),
          Pearson_r = cor(Mean_C, Mean_S, method = "pearson"),
          Pearson_p_raw = cor.test(Mean_C, Mean_S, method = "pearson")$p.value,
          # Using the v2 of this for now
          Pearson_p_adj = get_chelton_pval(Mean_C, Mean_S),
          N_CPUE = length(c_cols), 
          .groups = "drop"
        )
    }
    
    long_list[[m$name]] <- dat_clean %>% 
      dplyr::select(Year, Stock, Mean_CPUE = Mean_C, Mean_Survey = Mean_S, dplyr::all_of(c_cols), dplyr::all_of(s_cols)) %>%
      pivot_longer(-c(Year, Stock, Mean_CPUE, Mean_Survey), names_to = "Index", values_to = "Value") %>%
      mutate(Category = if_else(Index %in% c_cols, "CPUE/LPUE", "Survey"), Line_Type = "Raw") %>%
      drop_na(Value)
  }
  
  for (cc in c_cols) {
    for (ss in s_cols) {
      pair <- dat %>% dplyr::select(dplyr::all_of(cc), dplyr::all_of(ss)) %>% drop_na()
      if (nrow(pair) > 4) {
        
        # Calculate base correlations safely
        r_val <- cor(pair[[1]], pair[[2]], method = "pearson")
        p_raw <- cor.test(pair[[1]], pair[[2]], method = "pearson")$p.value
        
        pairs_list[[paste0(m$name, cc, ss)]] <- tibble(
          Stock = m$name, 
          CPUE_Index = cc, 
          Survey_Index = ss,
          Kendall_tau = cor(pair[[1]], pair[[2]], method = "kendall"),
          Pearson_r = r_val,
          Pearson_p_raw = p_raw,
          Pearson_p_adj = get_ar1_pval(pair[[1]], pair[[2]]),
          Survey_RMSSD = sqrt(mean(diff(pair[[2]])^2)), 
          Index_Type = if_else(str_detect(cc, "LPUE|Logbook"), "LPUE", "CPUE"),
          Comparison_Label = paste0(m$name, ": ", cc, " vs\n", ss)
        )
      }
    }
  }
}
df_pairs <- bind_rows(pairs_list) %>% mutate(Stock = factor(Stock, levels = all_stocks))
df_stock <- bind_rows(stock_list) %>% mutate(Stock = factor(Stock, levels = all_stocks))
df_long <- bind_rows(long_list) %>% mutate(Stock = factor(Stock, levels = all_stocks))

# =========================================================================
# 3. DATA PROCESSING: ASYMMETRY & BOOTSTRAPPING
# =========================================================================
print("Running Asymmetry Bootstrap simulations...")
asym_list <- list(); boot_ci_list <- list()
set.seed(42) # For reproducibility
n_iterations <- 1000

for (m in parent_meta) {
  if (!file.exists(m$file)) next
  dat <- read_csv(m$file, show_col_types = FALSE) 
  if ("YEAR" %in% names(dat)) dat <- dat %>% rename(Year = YEAR)
  dat <- dat %>% arrange(Year)
  
  c_cols <- intersect(m$cpue, names(dat)); s_cols <- intersect(m$survey, names(dat))
  
  if(length(c_cols) > 0 & length(s_cols) > 0) {
    dat_clean <- dat %>%
      mutate(Mean_C = rowMeans(dplyr::select(., dplyr::all_of(c_cols)), na.rm = TRUE),
             Mean_S = rowMeans(dplyr::select(., dplyr::all_of(s_cols)), na.rm = TRUE)) %>%
      dplyr::select(Year, Mean_C, Mean_S) %>% drop_na() %>%
      mutate(
        Delta_C = Mean_C - lag(Mean_C),
        Delta_S = Mean_S - lag(Mean_S),
        Survey_Direction = if_else(Delta_S > 0, "Survey Increased", "Survey Decreased"),
        Matched_Raw = sign(Delta_C) == sign(Delta_S)
      ) %>% drop_na(Matched_Raw)
    
    if(nrow(dat_clean) < 3) next
    
    # Base summary
    asym_list[[m$name]] <- dat_clean %>%
      dplyr::group_by(Survey_Direction) %>%
      dplyr::summarize(Concordance_Rate = mean(Matched_Raw, na.rm = TRUE), N_Years = n(), .groups = 'drop') %>%
      mutate(Stock = m$name)
    
    # Bootstrap
    boot_scores <- numeric(n_iterations)
    for (i in 1:n_iterations) {
      boot_samp <- dat_clean %>% slice_sample(n = nrow(dat_clean), replace = TRUE)
      boot_summary <- boot_samp %>%
        dplyr::group_by(Survey_Direction) %>%
        dplyr::summarize(Rate = mean(Matched_Raw, na.rm = TRUE), .groups = 'drop')
      
      rate_inc <- if("Survey Increased" %in% boot_summary$Survey_Direction) boot_summary$Rate[boot_summary$Survey_Direction == "Survey Increased"] else NA
      rate_dec <- if("Survey Decreased" %in% boot_summary$Survey_Direction) boot_summary$Rate[boot_summary$Survey_Direction == "Survey Decreased"] else NA
      boot_scores[i] <- rate_inc - rate_dec
    }
    
    boot_ci_list[[m$name]] <- tibble(
      Stock = m$name,
      CI_Lower = quantile(boot_scores, 0.025, na.rm = TRUE),
      CI_Upper = quantile(boot_scores, 0.975, na.rm = TRUE)
    )
  }
}

df_asym <- bind_rows(asym_list) %>% drop_na()
asym_wide <- df_asym %>%
  pivot_wider(id_cols = Stock, names_from = Survey_Direction, values_from = Concordance_Rate) %>%
  mutate(Asymmetry_Score = `Survey Increased` - `Survey Decreased`)

asym_wide_ci <- asym_wide %>% left_join(bind_rows(boot_ci_list), by = "Stock")

# =========================================================================
# 4. DATA PROCESSING: HARLEY ET AL. (2001) DEMING REGRESSION (ALL PAIRS)
# =========================================================================
print("Calculating Deming Regressions for ALL Pairs (Geometric Mean Scaling)...")

# --- Inject raw metadata specifically for Deming regressions ---
parent_meta_raw <- list(
  list(file = "atlantic_cod_all_indices_combined.csv", name = "Cod (EGOM)", cpue = "LPUE_EGOM_raw", survey = c("Survey_Fall_EGOM_raw", "Survey_Spring_EGOM_raw")),
  list(file = "atlantic_cod_all_indices_combined.csv", name = "Cod (WGOM)", cpue = "LPUE_WGOM_raw", survey = c("Survey_Fall_WGOM_raw", "Survey_Spring_WGOM_raw")),
  list(file = "atlantic_cod_all_indices_combined.csv", name = "Cod (GB)", cpue = "LPUE_GB_raw", survey = c("Survey_Fall_GB_raw", "Survey_Spring_GB_raw")),
  list(file = "atlantic_cod_all_indices_combined.csv", name = "Cod (SNE)", cpue = "LPUE_SNE_raw", survey = c("Survey_Fall_SNE_raw", "Survey_Spring_SNE_raw")),
  list(file = "monkfish_final.csv", name = "Monkfish (N)", cpue = c('Jones_North_Gillnet', 'Jones_North_Trawl', 'Jones_North_Dredge', 'Richardson_N_Trawl_Logbook_Raw', 'Richardson_N_Trawl_Observer_Raw'), survey = c('Survey_North_Fall_Raw', 'Survey_North_Spring_Raw')),
  list(file = "monkfish_final.csv", name = "Monkfish (S)", cpue = c('Jones_South_Gillnet', 'Jones_South_Trawl', 'Jones_South_Dredge', 'Richardson_S_Gillnet_Logbook_Raw', 'Richardson_S_Gillnet_Observer_Raw'), survey = c('Survey_South_Fall_Raw', 'Survey_South_Spring_Raw')),
  list(file = "yellowtail_flounder.csv", name = "Yellowtail (CC/GOM)", cpue = "CCGOM_LPUE_Raw", survey = "CCGOM_Biomass_Spring_Index"),
  list(file = "yellowtail_flounder.csv", name = "Yellowtail (GB)", cpue = "GB_LPUE_Raw", survey = "GB_Biomass_Spring_Index"),
  list(file = "yellowtail_flounder.csv", name = "Yellowtail (SNE/MA)", cpue = "SNEMA_LPUE_Raw", survey = "SNEMA_Biomass_Spring_Index"),
  list(file = "black_sea_bass_fresh_zscored_dataset.csv", name = "Black Sea Bass (N)", cpue = c('North_GAM_CPUE', 'North_Nominal_CPUE', 'North_Rec_CPUE'), survey = 'North_Survey_CPUE'),
  list(file = "black_sea_bass_fresh_zscored_dataset.csv", name = "Black Sea Bass (S)", cpue = c('South_GAM_CPUE', 'South_Nominal_CPUE', 'South_Rec_CPUE'), survey = 'South_Survey_CPUE'),
  list(file = "american_plaice_fully_updated.csv", name = "American Plaice", cpue = c('Terceiro_Raw_Index', 'Hankowsky_Raw_Index', 'Jones_SpatioTemporal_Raw', 'Jones_SpatialOnly_Raw'), survey = c('Spring_Wt_Tow', 'Fall_Wt_Tow')),
  list(file = "bluefish_indices_trim.csv", name = "Bluefish", cpue = c('MRIP_Guild_CPUE', 'MRIP_Directed_CPUE'), survey = c('NEFSC_Albatross_Fall', 'NEFSC_Bigelow_Fall', 'NEAMAP_Fall')),
  list(file = "gom_haddock_all_centered_indices_trim.csv", name = "Haddock (GOM)", cpue = c('GOM_Logbook_LPUE', 'GOM_StudyFleet_CPUE'), survey = c('NMFS_Fall_Bio', 'NMFS_Spring_Bio')),
  list(file = "scup_centered_indices_trim.csv", name = "Scup", cpue = 'MRIP_CPUE_Points', survey = c('Fall_ALB', 'Spring_ALB', 'Winter', 'Spring_BIG', 'Fall_BIG')),
  list(file = "illex_squid_lpue_hand.csv", name = "Shortfin Squid", cpue = c('Dealer_WetBoats', 'Observer_WetBoats', 'StudyFleet_WetBoats'), survey = 'NEFSC_Fall'),
  list(file = "Atlantic_halibut_values.csv", name = "Halibut", cpue = 'Standardized_CPUE_Est', survey = c('NEFSC_Fall_Survey_Est', 'NEFSC_Spring_Survey_Est')),
  list(file = "golden_tilefish_final_indices.csv", name = "Golden Tilefish", cpue = "GLM_LPUE", survey = c("Survey_Abundance_index", "Survey_Biomass_index"))
)

harley_list <- list(); beta_list <- list()

for (m in parent_meta_raw) { 
  if (!file.exists(m$file)) next
  dat <- read_csv(m$file, show_col_types = FALSE) 
  if ("YEAR" %in% names(dat)) dat <- dat %>% rename(Year = YEAR)
  
  c_cols <- intersect(m$cpue, names(dat)); s_cols <- intersect(m$survey, names(dat))
  
  # Loop through every possible pair combination
  for (cc in c_cols) {
    for (ss in s_cols) {
      dat_pair <- dat %>% dplyr::select(Year, C_Val = dplyr::all_of(cc), S_Val = dplyr::all_of(ss)) %>% drop_na()
      if(nrow(dat_pair) < 3) next
      
      dat_clean <- dat_pair %>%
        mutate(
          # 1. Scale to Geometric Mean
          C_Scaled = C_Val / exp(mean(log(C_Val))),
          S_Scaled = S_Val / exp(mean(log(S_Val))),
          
          # 2. Add safety constant
          Pos_C = C_Scaled + 0.001,
          Pos_S = S_Scaled + 0.001,
          
          # 3. Log transform
          Log_C = log(Pos_C),
          Log_S = log(Pos_S),
          
          Stock = m$name,
          CPUE_Index = cc,
          Survey_Index = ss,
          Pair_Name = paste0(cc, "_vs_", ss)
        )
      
      harley_list[[paste0(m$name, cc, ss)]] <- dat_clean
      
      fit_deming <- tryCatch({ deming(Log_C ~ Log_S + 0, data = dat_clean) }, error = function(e) NULL)
      
      if(!is.null(fit_deming)) {
        beta_deming <- coef(fit_deming)[2]
        q_deming <- exp(coef(fit_deming)[1])
        
        # Override negative slopes
        if (beta_deming <= 0) {
          beta_deming <- 0.0001 
          q_deming <- exp(mean(dat_clean$Log_C, na.rm = TRUE)) 
        }
      } else {
        beta_deming <- NA; q_deming <- NA
      }
      
      beta_list[[paste0(m$name, cc, ss)]] <- tibble(
        Stock = m$name, CPUE_Index = cc, Survey_Index = ss, Pair_Name = paste0(cc, "_vs_", ss),
        Beta_Deming = beta_deming, Q_Deming = q_deming
      )
    }
  }
}
df_harley <- bind_rows(harley_list) %>% mutate(Stock = factor(Stock, levels = all_stocks))
df_beta <- bind_rows(beta_list) %>% mutate(Stock = factor(Stock, levels = all_stocks)) %>% drop_na()

# GENERATE DEMING CURVES FOR EVERY SINGLE PAIR (FITTED TO LOCAL EXTENT)
pred_list_full <- list()
for(i in 1:nrow(df_beta)) {
  st <- as.character(df_beta$Stock[i])
  pair_nm <- df_beta$Pair_Name[i]
  b_dem <- df_beta$Beta_Deming[i]
  q_dem <- df_beta$Q_Deming[i]
  
  # Calculate max x-value for this specific stock so facets remain bounded
  local_max_s <- max(df_harley$Pos_S[df_harley$Stock == st], na.rm = TRUE) * 1.05
  x_seq_local <- seq(0, local_max_s, length.out = 100) 
  
  pred_list_full[[i]] <- tibble(Stock = st, Pair_Name = pair_nm, S_Ratio = x_seq_local, Full_Deming = q_dem * (x_seq_local^b_dem))
}
df_curve_full <- bind_rows(pred_list_full) %>% mutate(Stock = factor(Stock, levels = all_stocks))

# Create median label per stock for Plot 9 (using explicit dplyr to prevent package conflicts)
df_beta_labels <- df_beta %>% 
  dplyr::group_by(Stock) %>% 
  dplyr::summarize(Median_Beta = median(Beta_Deming, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(Label = paste0("Median \u03B2: ", sprintf("%.2f", Median_Beta)))

# =========================================================================
# 5. GENERATE ALL PLOTS
# =========================================================================
print("Generating Plots...")

# --- PLOT 1: Histogram ---
p1 <- ggplot(df_pairs, aes(x = Kendall_tau)) +
  geom_histogram(bins = 20, fill = "orange4", color = "white", alpha = 0.6) +
  geom_vline(xintercept = 0, color = "#E15759", linetype = "dashed", linewidth = 1.2) +
  geom_vline(xintercept = median(df_pairs$Kendall_tau, na.rm = TRUE), color = "#4E79A7", linetype = "dotdash", linewidth = 1.2) +
  labs(title = expression("Distribution of Kendall's " * tau * " (Concordance)"), x = expression("Kendall's " * tau), y = "Frequency") +
  theme_bw(base_size = 14) + theme(plot.title = element_text(face = "bold"))

# --- PLOT 2: Top 20 & Bottom 20 ---
bar_df <- df_pairs %>% arrange(desc(Kendall_tau))
p_top <- ggplot(head(bar_df, 20), aes(x = Kendall_tau, y = reorder(Comparison_Label, Kendall_tau), fill = Index_Type)) +
  geom_col() + scale_fill_manual(values = c("CPUE" = "orange4", "LPUE" = "orange")) +
  coord_cartesian(xlim = c(0, 0.8)) + labs(title = "20 Indices with highest correspondence", x = expression("Kendall's " * tau), y = NULL) +
  theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold"), legend.position = "none")
p_bot <- ggplot(tail(bar_df, 20), aes(x = Kendall_tau, y = reorder(Comparison_Label, Kendall_tau), fill = Index_Type)) +
  geom_col() + scale_fill_manual(values = c("CPUE" = "orange4", "LPUE" = "orange")) +
  coord_cartesian(xlim = c(-0.8, 0.3)) + labs(title = "20 Indices with lowest correspondence", x = expression("Kendall's " * tau), y = NULL) +
  theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(hjust = 0))
p2 <- p_bot + p_top + plot_layout(guides = "collect")

# --- PLOT 3: Volatility vs Concordance ---
p3 <- ggplot(df_pairs, aes(x = Survey_RMSSD, y = Kendall_tau, color = Stock)) +
  geom_hline(yintercept = 0, color = "gray50", linetype = "solid") +
  geom_point(size = 4, alpha = 0.7) +
  geom_smooth(aes(group = 1), method = "lm", color = "black", linetype = "dashed", se = FALSE, linewidth = 1) +
  scale_color_manual(values = fixed_stock_colors) +
  labs(title = "Does High Survey Volatility Impact Concordance?", subtitle = "Analysis of all individual CPUE-to-Survey pairs", x = "Survey Volatility (RMSSD)", y = expression("Kendall's " * tau)) +
  theme_bw(base_size = 14) + theme(plot.title = element_text(face = "bold"))

# --- PLOT 4: Volatility vs Concordance (Aggregated by Stock) ---
agg_df <- df_pairs %>% 
  dplyr::group_by(Stock) %>%
  dplyr::summarize(Mean_Tau = mean(Kendall_tau), SD_Tau = sd(Kendall_tau), Mean_RMSSD = mean(Survey_RMSSD), SD_RMSSD = sd(Survey_RMSSD), .groups = "drop") %>%
  dplyr::mutate(SD_Tau = replace_na(SD_Tau, 0), SD_RMSSD = replace_na(SD_RMSSD, 0))
p4 <- ggplot(agg_df, aes(x = Mean_RMSSD, y = Mean_Tau, color = Stock)) +
  geom_hline(yintercept = 0, color = "gray50") +
  geom_errorbar(aes(ymin = Mean_Tau - SD_Tau, ymax = Mean_Tau + SD_Tau), width = 0.02, linewidth = 0.8, alpha = 0.6) +
  geom_errorbarh(aes(xmin = Mean_RMSSD - SD_RMSSD, xmax = Mean_RMSSD + SD_RMSSD), height = 0.02, linewidth = 0.8, alpha = 0.6) +
  geom_point(size = 5, alpha = 0.9) +
  geom_smooth(aes(group = 1), method = "lm", color = "black", linetype = "dashed", se = FALSE, linewidth = 1) +
  scale_color_manual(values = fixed_stock_colors) +
  labs(x = "Mean Survey Volatility (RMSSD)", y = expression("Mean Kendall's " * tau)) +
  theme_bw(base_size = 14) + theme(plot.title = element_text(face = "bold"))

# --- PLOT 5: Portfolio Effect ---
p5 <- ggplot(df_stock, aes(x = factor(N_CPUE), y = Tau)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", alpha = 0.6) +
  geom_boxplot(aes(group = N_CPUE), width = 0.4, outlier.shape = NA, alpha = 0.1) +
  geom_jitter(aes(color = Stock), size = 4, width = 0.15, alpha = 0.8) +
  scale_color_manual(values = fixed_stock_colors) +
  labs(x = "Number of CPUE/LPUE Indices Averaged Together", y = expression("Concordance of the Means (Kendall's " * tau * ")")) +
  theme_bw(base_size = 14) + theme(plot.title = element_text(face = "bold"))

# --- PLOT 6: Time-Series Trends Grid ---
p6 <- ggplot(df_long) +
  geom_hline(yintercept = 0, color = "gray50", linetype = "dashed", alpha = 0.6) +
  geom_line(aes(x = Year, y = Value, group = Index, color = Category), alpha = 0.35, linewidth = 0.6) +
  geom_line(aes(x = Year, y = Mean_CPUE, color = "CPUE/LPUE"), linewidth = 1.2) +
  geom_line(aes(x = Year, y = Mean_Survey, color = "Survey"), linewidth = 1.2) +
  scale_color_manual(values = c("CPUE/LPUE" = "darkorange", "Survey" = "#008080")) +
  facet_wrap(~Stock, ncol = 5, scales = "free") +
  labs(#title = "CPUE vs. Fishery-Independent Survey Trends",
       x = "Year", y = "Standardized Index Value", color = "Index Category") +
  theme_bw(base_size = 14) + 
  theme(
    plot.title = element_text(face = "bold"), 
    strip.text = element_text(face = "bold"), 
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"  # <--- This moves the legend
  )
# --- PLOT 7: Asymmetry Error Bars ---
p7 <- ggplot(asym_wide_ci, aes(x = Asymmetry_Score, y = reorder(Stock, Asymmetry_Score))) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1.2) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.3, color = "black", linewidth = 1.2) +
  geom_point(aes(fill = Asymmetry_Score > 0), color = "black", shape = 21, size = 6, stroke = 1.5) +
  scale_fill_manual(values = c("TRUE" = "#E15759", "FALSE" = "#4E79A7")) +
  scale_x_continuous(labels = scales::percent_format()) +
  labs(x = "Difference in Tracking Concordance (% Match Ups vs % Match Downs)", y = NULL) +
  theme_bw(base_size = 14) + theme(legend.position = "none", plot.title = element_text(face = "bold"))

# --- PLOT 8: Conceptual Quadrant (WGOM Cod) ---
ex_dat <- read_csv("atlantic_cod_all_indices_combined.csv", show_col_types = FALSE)
ex_clean <- ex_dat %>%
  mutate(Mean_C = LPUE_WGOM_Combined_zscored, Mean_S = rowMeans(dplyr::select(., NMFS_Fall_WGOM_zscored, NMFS_Spring_WGOM_zscored), na.rm = TRUE)) %>%
  dplyr::select(Year, Mean_C, Mean_S) %>% drop_na() %>%
  mutate(Delta_C = Mean_C - lag(Mean_C), Delta_S = Mean_S - lag(Mean_S),
         Survey_Direction = if_else(Delta_S > 0, "Survey Increased", "Survey Decreased"),
         Matched = sign(Delta_C) == sign(Delta_S)) %>% drop_na(Survey_Direction)

p8 <- ggplot(ex_clean, aes(x = Delta_S, y = Delta_C)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "#E15759", alpha = 0.15) + 
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = "#4E79A7", alpha = 0.15) + 
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill = "gray90", alpha = 0.5) + 
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = "gray90", alpha = 0.5) + 
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1) + geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
  geom_point(aes(fill = Survey_Direction, shape = Matched), size = 6, color = "black", stroke = 1.2) +
  scale_shape_manual(values = c("TRUE" = 21, "FALSE" = 24)) + 
  scale_fill_manual(values = c("Survey Increased" = "#E15759", "Survey Decreased" = "#4E79A7")) +
  annotate("text", x = 1.5, y = 1.5, label = "CPUE correctly\ntracked UP", color = "#E15759", fontface = "bold", size = 5) +
  annotate("text", x = -1.5, y = -1.5, label = "CPUE correctly\ntracked DOWN", color = "#4E79A7", fontface = "bold", size = 5) +
  annotate("text", x = 1.5, y = -1.5, label = "MISMATCH\n(Survey Up, CPUE Down)", color = "gray40", fontface = "italic") +
  annotate("text", x = -1.5, y = 1.5, label = "MISMATCH\n(Survey Down, CPUE Up)", color = "gray40", fontface = "italic") +
  labs(title = "WGOM Cod: Year-over-Year Delta Tracking", x = "Change in Survey Index", y = "Change in Commercial LPUE") +
  coord_fixed(xlim = c(-max(abs(ex_clean$Delta_S),na.rm=T)-0.5, max(abs(ex_clean$Delta_S),na.rm=T)+0.5), 
              ylim = c(-max(abs(ex_clean$Delta_C),na.rm=T)-0.5, max(abs(ex_clean$Delta_C),na.rm=T)+0.5)) +
  theme_bw(base_size = 14) + theme(plot.title = element_text(face = "bold"), guides = guide_legend(override.aes = list(fill = "gray50"))) 

# --- PLOT 9: Deming Curves Faceted (All Pairs) ---
p9 <- ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  geom_point(data = df_harley, aes(x = Pos_S, y = Pos_C, color = Stock), alpha = 0.3, size = 1.5) +
  geom_line(data = df_curve_full, aes(x = S_Ratio, y = Full_Deming, group = Pair_Name, color = Stock), linewidth = 0.8, alpha = 0.6) +
  geom_text(data = df_beta_labels, aes(x = -Inf, y = Inf, label = Label), hjust = -0.1, vjust = 1.2, size = 4, fontface = "bold", color = "gray20") +
  scale_color_manual(values = fixed_stock_colors) +
  facet_wrap(~Stock, scales = "free", ncol = 5) +
  labs(x = "Survey Abundance (Geometric Mean Scaled)", y = "Commercial CPUE (Geometric Mean Scaled)") +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  theme_bw(base_size = 12) + 
  theme(plot.title = element_text(face = "bold"), legend.position = "none", 
        strip.background = element_rect(fill = "gray90"), strip.text = element_text(face = "bold"))

# --- PLOT 10: Appendix - Kendall vs Pearson ---
# Prepare the dataset by categorizing significance based on the AR(1) adjusted p-value
df_pairs_plot10 <- df_pairs %>%
  filter(!is.na(Pearson_r) & !is.na(Kendall_tau)) %>%
  mutate(
    Significant = factor(
      if_else(Pearson_p_adj < 0.05, "Significant (p < 0.05)", "Not Significant"),
      levels = c("Significant (p < 0.05)", "Not Significant")
    )
  )

p10 <- ggplot(df_pairs_plot10, aes(x = Kendall_tau, y = Pearson_r)) +
  # Reference lines
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  
  # Plot points with conditional fill for open/closed circles
  geom_point(aes(fill = Significant), shape = 21, size = 3.5, color = "black", stroke = 0.8, alpha = 0.8) +
  
  # Styling
  scale_fill_manual(values = c("Significant (p < 0.05)" = "#4E79A7", "Not Significant" = "white")) +
  labs(
    title = "Comparison of Pearson (r) vs. Kendall's Rank (\u03C4) Correlations",
    subtitle = "Solid points indicate significance using AR(1) adjusted degrees of freedom",
    x = expression("Kendall's " * tau),
    y = "Pearson Correlation (r)",
    fill = "AR(1) Adjusted Significance"
  ) +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

# --- PLOT 11: Beta Variation (Boxplot + Jitter) ---
p11 <- ggplot(df_beta, aes(x = reorder(Stock, Beta_Deming, FUN = median), y = Beta_Deming)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1.2, alpha = 0.8) +
  geom_boxplot(outlier.shape = NA, fill = "gray90", alpha = 0.5, width = 0.5) +
  geom_jitter(aes(color = Stock), width = 0.15, size = 3, alpha = 0.8) +
  scale_color_manual(values = fixed_stock_colors) +
  coord_flip() + 
  labs(title = "Variation in Proportionality (\u03B2) Across Fleet-Survey Pairs",
       subtitle = "Values near 1.0 (Red Line) indicate 1:1 proportional tracking",
       x = NULL, y = expression("Estimated Deming Slope (" * beta * ")")) +
  theme_bw(base_size = 14) + theme(plot.title = element_text(face = "bold"), legend.position = "none")


# =========================================================================
# 6. SAVE ALL PLOT OUTPUTS 
# =========================================================================
print("Saving outputs...")

ggsave("01_kendall_tau_histogram.png", p1, width = 8, height = 6, dpi = 300)
ggsave("02_top_bottom_20_tracking.png", p2, width = 16, height = 8, dpi = 300)
#ggsave("03_tau_vs_rmssd_all_points.png", p3, width = 10, height = 6, dpi = 300)
ggsave("04_tau_vs_rmssd_aggregated.png", p4, width = 10, height = 6, dpi = 300)
ggsave("05_n_cpue_vs_tracking.png", p5, width = 10, height = 7, dpi = 300)
ggsave("06_stock_specific_trends.png", p6, width = 22, height = 12, dpi = 300)
ggsave("07_Asymmetry_Styled_ErrorBars.png", p7, width = 12, height = 8, dpi = 300)
#ggsave("08_WGOM_Cod_Quadrant.png", p8, width = 10, height = 10, dpi = 300)
ggsave("09_Harley_Faceted_Pairs.png", p9, width = 14, height = 12, dpi = 300)
ggsave("10_Appendix_Kendall_vs_Pearson.png", p10, width = 8, height = 8, dpi = 300)
ggsave("11_Harley_Beta_Boxplots.png", p11, width = 12, height = 8, dpi = 300)


# =========================================================================
# 7. PRINT CONSOLE STATISTICS & GENERATE APPENDIX CSVs
# =========================================================================

cat("\n=========================================\n")
cat("      BETA STATISTICS (\u03B2) REPORT\n")
cat("=========================================\n")

# Mean/Median excluding artificial 0.0001 placeholders
mean_beta <- mean(df_beta$Beta_Deming[df_beta$Beta_Deming > 0.001], na.rm = TRUE)
median_beta <- median(df_beta$Beta_Deming[df_beta$Beta_Deming > 0.001], na.rm = TRUE)

cat("Overall Mean Beta:  ", round(mean_beta, 2), "\n")
cat("Overall Median Beta:", round(median_beta, 2), "\n\n")

# Dunn et al. (2000) Plausible Window Filtering
beta_filtered <- df_beta$Beta_Deming[df_beta$Beta_Deming >= 0.5 & df_beta$Beta_Deming <= 2.0]
mean_beta_filtered <- mean(beta_filtered, na.rm = TRUE)
median_beta_filtered <- median(beta_filtered, na.rm = TRUE)

cat("Filtered Mean Beta (0.5 - 2.0):  ", round(mean_beta_filtered, 2), "\n")
cat("Filtered Median Beta (0.5 - 2.0):", round(median_beta_filtered, 2), "\n\n")

cat("Total fleet-survey pairs:    ", nrow(df_beta), "\n")
cat("Pairs within 0.5-2.0 window: ", length(beta_filtered), "\n")
cat("Percent within window:       ", round((length(beta_filtered) / nrow(df_beta)) * 100, 1), "%\n")


cat("\n=========================================\n")
cat("      KENDALL'S TAU TABLES GENERATION\n")
cat("=========================================\n")

# 1. ENSEMBLE Tau (Stock-level averages)
ensemble_tau_report <- df_stock %>%
  dplyr::select(Stock, N_CPUE, Tau) %>%
  dplyr::mutate(Tau = round(Tau, 3)) %>%
  dplyr::rename(
    `Number of CPUE Indices` = N_CPUE,
    `Ensemble Kendall's Tau` = Tau
  ) %>%
  dplyr::arrange(desc(`Ensemble Kendall's Tau`))

write_csv(ensemble_tau_report, "Report_Ensemble_Kendalls_Tau.csv")
cat("-> Saved 'Report_Ensemble_Kendalls_Tau.csv'\n")

# 2. PAIRWISE Tau (Individual Fleet vs Survey)
pairwise_tau_report <- df_pairs %>%
  dplyr::select(Stock, CPUE_Index, Survey_Index, Kendall_tau) %>%
  dplyr::mutate(Kendall_tau = round(Kendall_tau, 3)) %>%
  dplyr::rename(
    `Commercial CPUE Series` = CPUE_Index,
    `Scientific Survey Series` = Survey_Index,
    `Pairwise Kendall's Tau` = Kendall_tau
  ) %>%
  dplyr::arrange(Stock, desc(`Pairwise Kendall's Tau`))

write_csv(pairwise_tau_report, "Report_Pairwise_Kendalls_Tau.csv")
cat("-> Saved 'Report_Pairwise_Kendalls_Tau.csv'\n")

print("Analysis pipeline successfully finished!")