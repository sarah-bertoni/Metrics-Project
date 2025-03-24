
#' Internal function of did_multiplegt_stat to handle unbalancedness
#' @param df df
#' @returns A matrix
#' @importFrom stats reshape
#' @noRd
balanced_map <- function(df) {
  df <- df[c("ID_XX", "T_XX", "tsfilled_XX")]
  df$H_t_XX <- as.numeric(1 - df$tsfilled_XX)
  df$tsfilled_XX <- NULL
  df <- df[order(df$ID_XX, df$T_XX), ]
  df <- df %>% group_by(.data$ID_XX) %>%
    mutate(H_t_m_1_XX = lag(.data$H_t_XX)) %>% ungroup()
  df <- subset(df, !is.na(df$H_t_m_1_XX))
  df$H_t <- as.numeric(df$H_t_XX == 1 & df$H_t_m_1_XX == 1)
  df$H_t_XX <- df$H_t_m_1_XX <- NULL
  times <- levels(factor(df$T_XX))
  df <- stats::reshape(data = as.data.frame(df), idvar = "ID_XX", timevar = "T_XX", 
                       direction = "wide", v.names = "H_t") 
  colnames(df) <- c("ID_XX", sapply(times, function(x) return(paste0("H_",x))))
  rownames(df) <- 1:nrow(df)
  return(df)
}
#' Internal function for by_fd option graph
#' @param obj obj
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @returns The did_multiplegt_stat object + a ggplot object
#' @noRd
by_fd_graph <- function(obj) {
  suppressWarnings({
    pe_set <- as.data.frame(matrix(NA, ncol = 6, nrow = 0))
    names(pe_set) <- c("model", "pe", "lb", "ub", "nswitchers")
    models <- c("aoss", "waoss", "ivwaoss")
    for (i in 1:3) {
      if (models[i] %in% obj$args$estimator) {
        pe_set_temp <- as.data.frame(matrix(NA, ncol = 6, nrow = length(obj$by_levels)))
        for (j in 1:length(obj$by_levels)) {
          subobj <- obj[[paste0("results_by_",j)]]
          pe_set_temp[j,2] <- subobj$table[subobj$pairs * (i-1) + 1, 1]
          pe_set_temp[j,3] <- subobj$table[subobj$pairs * (i-1) + 1, 3]
          pe_set_temp[j,4] <- subobj$table[subobj$pairs * (i-1) + 1, 4]
          pe_set_temp[j,5] <- obj$switchers_df$N_partition_XX[j]
          pe_set_temp[j,6] <- obj$switchers_df$Med_delta_pre_XX[j]
        }
        names(pe_set_temp) <- c("model", "pe", "lb", "ub", "nswitchers", "median")
        pe_set_temp$model <- models[i]
        pe_set_temp$lbin <- obj$quantiles[2, 1:ncol(obj$quantiles) - 1]
        pe_set_temp$ubin <- obj$quantiles[2, 2:ncol(obj$quantiles)]
        pe_set_temp$lbin_cdf <- obj$quantiles[1, 1:ncol(obj$quantiles) - 1] * 100
        pe_set_temp$ubin_cdf <- obj$quantiles[1, 2:ncol(obj$quantiles)] * 100
        pe_set_temp$id <- 1:length(obj$by_levels)
        pe_set_temp$include <- ifelse(pe_set_temp$id == 1, "[", "(")
        
        pe_set <- rbind(pe_set, pe_set_temp)
      }
    }
    for (c in 1:ncol(pe_set)) {
      pe_set[,c] <- ifelse(is.nan(pe_set[,c]), NA, pe_set[,c])
    }
    var_gr <- ifelse("ivwaoss" %in% pe_set$model, "Z", "D")
    pe_set$colname <- sprintf("%.2f\n(%.0f%%-%.0f%%)\n%s%.2f,%.2f]\nN=%.0f\n", pe_set$median, pe_set$lbin_cdf, pe_set$ubin_cdf,pe_set_temp$include, pe_set_temp$lbin, pe_set_temp$ubin, pe_set$nswitchers)
    
    ticks <- c()
    labels <- c()
    ticks_size <- c()
    for (j in 1:(nrow(pe_set)/length(levels(as.factor(pe_set$model))))) {
      ticks <- c(ticks, pe_set$median[j])
      labels <- c(labels, pe_set$colname[j])
    }
    ticks <- c(ticks, pe_set$ubin[nrow(pe_set)])
    labels <- c(labels, "")
    ticks_size <- c(ticks_size, 1)
    
    by_graph_tot <- NULL
    by_graph_1_XX <- NULL
    by_graph_2_XX <- NULL
    font <- 18/length(obj$args$estimator)
    tot_lim <- c(0,0)
    for (j in 1:length(obj$args$estimator))  {
      pe_set_temp <- subset(pe_set, pe_set$model == obj$args$estimator[j])
      
      by_graph <- ggplot(data = pe_set_temp, aes(x = .data$median, y = .data$pe)) + 
        geom_point(size = 4) + geom_line(size = 0.2, linetype = "dashed") +
        geom_errorbar(aes(ymin = .data$lb, ymax = .data$ub, fill = .data$model), width = 0) +
        xlab(sprintf("|\U0394%s| - %s", var_gr, toupper(obj$args$estimator[j]))) + ylab("") +
        scale_x_continuous(breaks= ticks, labels = labels) +
        theme(plot.title = element_text(hjust = 0.5, size = 2*font), axis.ticks.x = element_line(), axis.line.x = element_line(color="black", size = 0.5), axis.ticks.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_line(color = "black", size = 0.2), axis.text = element_text(size = font), axis.title = element_text(size = font)) + geom_hline(yintercept = 0, color = "black", size = 0.5)
      
      tot_lim[1] <- ifelse(layer_scales(by_graph)$y$range$range[1] > tot_lim[1], tot_lim[1], layer_scales(by_graph)$y$range$range[1])
      tot_lim[2] <- ifelse(layer_scales(by_graph)$y$range$range[2] > tot_lim[2], layer_scales(by_graph)$y$range$range[2], tot_lim[2])
      assign(paste0("by_graph_",j,"_XX"), by_graph)
      pe_set_temp <- by_graph <- NULL
    }
    for (j in 1:length(obj$args$estimator))  {
      assign(paste0("by_graph_",j,"_XX"), get(paste0("by_graph_",j,"_XX")) + ylim(tot_lim))
    }
    if (length(obj$args$estimator) == 1) {
      by_graph_tot <- by_graph_1_XX
    } else {
      by_graph_tot <- plot_grid(by_graph_1_XX, by_graph_2_XX, nrow = 1)
    }
  })
  return(by_graph_tot)
}

#' Internal function for by option graph
#' @param obj obj
#' @import ggplot2
#' @returns The did_multiplegt_stat object + a ggplot object
#' @noRd
by_graph <- function(
    obj
){
  
  by_str <- obj$args$by[1]
  if (length(obj$args$by) > 1) {
    for (j in 2:length(obj$args$by)) {
      by_str <- paste0(by_str,",", obj$args$by[j])
    }      
  }
  
  pe_set <- as.data.frame(matrix(NA, ncol = 5, nrow = 0))
  names(pe_set) <- c("model", "pe", "lb", "ub")
  models <- c("aoss", "waoss", "ivwaoss")
  x_count <- 1
  n_ci <- length(obj$by_levels)
  v_off <- 0.6
  offset <- -(v_off/2) + (v_off/(n_ci -1)) * (1:n_ci -1)
  
  gr_col <- c()
  for (i in 1:3) {
    if (models[i] %in% obj$args$estimator) {
      pe_set_temp <- as.data.frame(matrix(NA, ncol = 5, nrow = length(obj$by_levels)))
      for (j in 1:length(obj$by_levels)) {
        subobj <- obj[[paste0("results_by_",j)]]
        pe_set_temp[j,2] <- subobj$table[subobj$pairs * (i-1) + 1, 1]
        pe_set_temp[j,3] <- subobj$table[subobj$pairs * (i-1) + 1, 3]
        pe_set_temp[j,4] <- subobj$table[subobj$pairs * (i-1) + 1, 4]
        pe_set_temp[j,5] <- sprintf("(%s)", obj$by_levels[j])
      }
      names(pe_set_temp) <- c("model", "pe", "lb", "ub", "id")
      pe_set_temp$model <- models[i]
      pe_set_temp$x_pos <- x_count - offset
      x_count <- x_count + 1            
      pe_set <- rbind(pe_set, pe_set_temp)
    }
  }
  pe_set <- subset(pe_set, !is.na(pe_set$pe) & !is.nan(pe_set$pe))
  by_graph <- ggplot(data = pe_set, aes(x = .data$x_pos, y = .data$pe, group = 1)) + 
    geom_point(size = 2, aes(color = factor(.data$model))) + 
    geom_errorbar(aes(ymin = .data$lb, ymax = .data$ub, color = factor(.data$model)), width = 0.1) +
    scale_x_continuous(breaks= pe_set$x_pos, labels = pe_set$id) + ggtitle(sprintf("did_multiplegt_stat results by %s", by_str)) + ylab(" ") + xlab(" ") +
    theme(plot.title = element_text(hjust = 0.5), axis.ticks.y = element_blank(), legend.position = "right",
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(), panel.grid.major.y = element_line(color = "black", size = 0.2)) + geom_hline(yintercept = 0, color = "black", size = 0.2) + labs(color = "Estimators")
  return(by_graph)
}


#' Main interface for did_multiplegt_stat
#' @importFrom haven read_dta
#' @md
#' @description Estimation of Heterogeneity-robust Difference-in-Differences Estimators, with a Binary, Discrete, or Continuous Treatment or Instrument, in Designs with Stayers.
#' @param df (data.frame) A dataframe object.
#' @param Y (char) Outcome variable.
#' @param ID (char) Identifier of the unit of analysis.
#' @param Time (char) Time variable.
#' @param D (char) Treatment variable.
#' @param Z (char) Instrumental variable. This option is only required when the IV-related estimator (the so-called ivwaoss) is requested.
#' @param estimator (char vector) Estimator(s) to be computed. The allowed arguments are: (1) "aoss", i.e the Average Of Switchers’ Slopes which is the average, across switchers, of the effect on their period-(t) outcome of moving their treatment from its period-(t-1) to its period-(t) value, scaled by the difference between these two values. (2) "waoss" which corresponds to a weighted version of "aoss" where slopes receive a weight proportional to switchers’ absolute treatment change from period-(t-1) to period-(t). (3) "ivwaoss" which generalizes "waoss" to the instrumental-variable case, and is equal to the reduced-form "waoss" effect of the instrument on the outcome, divided by the first-stage "waoss" effect of the instrument on the treatment. If this option is not specified: by default, the command estimates both "aoss" and "waoss" if the instrumental-variable Z is not specified, or only "ivwaoss" otherwise. 
#' @param estimation_method (char) This option allows to specify which estimation method to use when estimating the waoss or the ivwaoss. The allowed arguments are "ra" (regression adjustment-based approach), "ps" (propensity-based approach), "dr" (double robust-based approach). By default, a doubly-robust estimator is used.
#' @param order (int) when the exact_match option is not specified, this option specifies the polynomial order to be used in the OLS regressions of \eqn{Y_t-Y_{t-1}} on a polynomial in \eqn{D_{t-1}} and/or in the logistic regressions of an indicator for \eqn{(t-1)}-to-\eqn{t} switchers on a polynomial in \eqn{D_{t-1}}. By default, a polynomial of order 1 is used.
#' @param switchers (char)  if the argument \code{up} is inputted, the command estimates the AOSS, WAOSS, or IV-WAOSS for switchers-up only, i.e for units whose treatment (or instrument) increases from period \eqn{(t-1)} to \eqn{t}. If the argument \code{down} is inputted, the command estimates the AOSS, WAOSS, or IV-WAOSS for switchers-down only, i.e. for units whose treament (or instrument) decreases from period \eqn{(t-1)} to \eqn{t}. By default, the command estimates those parameters for all switchers.
#' @param placebo (logical) This option allows to estimate the placebos versions of the estimators requested in the estimator option.
#' @param aoss_vs_waoss (logical) When this option is specified, the command performs and displays the test of the equality between the aoss and  the waoss. Note that the use of this option requires specifying in the estimator option both aoss and waoss.
#' @param other_treatments (character, len \eqn{\geq 1}) This option allows controlling for other treatments that may also change over the panel.
#' @param exact_match (logical) With this option, the DID estimators computed by the command compare the outcome evolution of switchers and stayers with the same period-\eqn{(t-1)} treatment (or instrument) value. This option can only be used when the treatment (or instrument) is binary or discrete: with a continuously distributed treatment (or instrument), one cannot find switchers and stayers with the exact same period-\eqn{(t-1)} treatment (or instrument). With a discrete treatment taking a large number of values, specifying this option may be undesirable: then, there may only be few switchers that can be matched to a stayer with the exact same period-\eqn{(t-1)} treatment, thus restricting the estimation sample.
#' @param noextrapolation (logical) when this option is specified, the command only keeps switchers whose period-\eqn{(t-1)} treatment (or instrument) is between the minimum and the maximum values of the period-\eqn{(t-1)} treatment (or instrument) of stayers.
#' @param by (character) runs the program by each level of varname specified. Only time-invariant variables are allowed.
#' @param by_fd (numeric integer) This option can be used if one wants to assess the heterogeneity of the effect according to the absolute value of the changes in the treatment. For example, if \code{by_fd = 5} is specified, the command will split the switchers into 5 groups delimited by the 4 quantiles of the distribution of \eqn{|\Delta D_t|} (or \eqn{|\Delta Z_t|}), and computes the models of each sample.
#' @param disaggregate (logical)  when this option is specified, the command shows the estimated AOSS, WAOSS, or IV-WAOSS effects for each pair of consecutive time periods, on top of the effects aggregated across all time periods. By default, the command only shows effects aggregated across all time periods.
#' @param cluster cluster
#' @details
#' # Overview
#' ## Data and design
#' The command uses panel data at the \eqn{(\text{ID},T)} level to estimate heterogeneity-robust DID estimators, 
#' with a binary, discrete, or continuous treatment (or instrument). The command can be used in designs where there is at least one 
#' pair of consecutive time periods between which the treatment of some units, the switchers, changes, while the treatment of some
#' other units, the stayers, does not change. 
#' 
#' ## Target parameters
#' The command can estimate the Average Of Switchers' Slopes (AOSS) and the Weighted Average Of Switchers' Slopes (WAOSS) parameters
#' introduced in de Chaisemartin et al ([2022](https://ssrn.com/abstract=4011782)).
#' The AOSS is the average, across switchers, of the slope: 
#' \deqn{\dfrac{Y_t(D_t)-Y_t(D_{t-1})}{D_t-D_{t-1}}} 
#' that is, the effect on their period-t outcome of moving their period-t
#' treatment from its period-(t-1) to its period-t value, scaled by the difference between these two values. The WAOSS is a weighted 
#' average of switchers' slopes, where slopes receive a weight proportional to \eqn{|D_t-D_{t-1}|},
#' switchers’ absolute treatment change from period-\eqn{(t-1)} to period-\eqn{t}. The variance of the WAOSS estimator is often smaller 
#' than that of the AOSS estimator, especially when there are switchers that experience a small treatment change. 
#' The WAOSS estimator is also amenable to doubly-robust estimation, unlike the AOSS estimator.
#' 
#' ## Assumptions
#' When the data has more than two time periods, the command assumes a static model: units' outcome at period \eqn{t} only depends on their period-\eqn{t}
#' treatment, not on their lagged treatments. See the [did_multiplegt_dyn](https://cran.r-project.org/web/packages/DIDmultiplegtDYN/index.html) 
#' command for estimators allowing for dynamic effects. 
#' The command also makes a parallel trends assumption: the counterfactual outcome evolution
#' switchers would have experienced if their treatment had not changed is assumed to be equal to
#' the outcome evolution of stayers with the same baseline treatment. To test that assumption, the command can compute placebo estimators comparing
#' the outcome evolution of switchers and stayers with the same baseline treatment before switchers' treatment changes. 

#' ## Estimators, when the exact_match option is specified
#' With a binary or discrete treatment, if the \code{exact_match} option is specified, the estimators computed by the command compare the 
#' outcome evolution of switchers and stayers 
#' with the same period-\eqn{(t-1)} treatment. Then, the WAOSS estimator computed by \code{did_multiplegt_stat()} 
#' is numerically equivalent to the DID_M estimator proposed by de Chaisemartin and D'Haultfoeuille 
#' ([2020a](https://aeaweb.org/articles?id=10.1257/aer.20181169)) and computed by the 
#' [did_multiplegt](https://cran.r-project.org/web/packages/DIDmultiplegt/index.html) command. \code{did_multiplegt_stat()} uses an analytic formula
#' to compute the estimator's variance, while \code{did_multiplegt()} uses the bootstrap. 
#' Thus, the run time of \code{did_multiplegt_stat()} is typically much lower. 
#' The \code{exact_match} option can only be specified when the treatment is binary or discrete:
#' with a continuously distributed treatment, one cannot find switchers and stayers 
#' with the exact same period-\eqn{(t-1)} treatment. 
#' With a discrete treatment taking a large number of values, specifying this option may be undesirable: 
#' then, there may only be few switchers that can be
#' matched to a stayer with the exact same period-\code{(t-1)} treatment, thus restricting the estimation sample.

#' ## Estimators, when the exact_match option is not specified
#' When the \code{exact_match} option is not specified, the command can use a regression adjustment to recover switchers' 
#' counterfactual outcome evolution: for all \eqn{t}, it runs an OLS regression of \eqn{Y_t-Y_{t-1}} 
#' on a polynomial in \eqn{D_{t-1}} in the sample of \eqn{(t-1)}-to-\eqn{t} stayers, and uses that regression
#' to predict switchers' counterfactual outcome evolution. Alternatively, when it estimates the WAOSS, the command can also use propensity-score 
#' reweighting to recover switchers' counterfactual outcome evolution. First, for all \eqn{t} it estimates a logistic regression of an indicator 
#' for \eqn{(t-1)}-to-\eqn{t} switchers on a polynomial in \eqn{D_{t-1}}, to predict units' probability of being a switcher. 
#' Then, it computes a weighted average of stayers' outcome evolution, upweighting stayers with a large probability of being switchers, 
#' and downweighting stayers with a low probability of being switchers. Finally, when it estimates the WAOSS, the command can also combine 
#' regression-adjustment and propensity-score reweighting, thus yielding a doubly-robust estimator.

#' ## Instrumental-variable case
#' There may be instances where the parallel-trends assumption fails, but one
#' has at hand an instrument satisfying a similar parallel-trends assumption. For instance, one may
#' be interested in estimating the price-elasticity of a good's consumption, but prices respond to
#' supply and demand shocks, and the counterfactual consumption evolution of units experiencing and not experiencing a price 
#' change may therefore not be the same. On the other hand, taxes
#' may not respond to supply and demand shocks and may satisfy a parallel-trends assumption. In such cases, the command
#' can compute the IV-WAOSS estimator introduced in de Chaisemartin et al ([2022](https://ssrn.com/abstract=4011782))
#' The IV-WAOSS estimator is equal to the WAOSS estimator of the instrument's reduced-form effect on the outcome, divided by the
#' WAOSS estimator of the instrument's first-stage effect on the treatment.

#' @section FAQ:
#' TBD
#' @section Comparison with Stata command:
#' Stata "logit" and R "glm" functions handle binary prediction with slightly different conventions. These discrepancies are usually negligible, but they may add up to detectable (yet small) differences in the final estimates. The next code blocks showcase an instance where the logit and glm predictions differ. We estimate a logistic regression of a binary variable D on an order 2 polynomial of a continuous variable X. The binary variable D takes value 1 only for D = 38.4. Due to these sample features, the logit command fails to converge. Both commands yield non missing predictions from their respective regression outputs. However, the predicted values at rows 8 and 9 are strictly different, with Stata reporting way larger predictions than R.
#' ## Stata
#' 
#' global rep "https://raw.githubusercontent.com/chaisemartinPackages"
#' 
#' use "$rep/ApplicationData/main/Tests/logit_tests.dta", clear
#' 
#' cap logit D X X_sq, asis
#' 
#' predict D_hat, pr asif
#' 
#' browse
#' 
#' ## R
#' 
#' library(haven)
#' 
#' library(stats)
#' 
#' rep <- "https://raw.githubusercontent.com/chaisemartinPackages"
#' 
#' data <- haven::read_dta(paste0(rep,"/ApplicationData/main/Tests/logit_tests.dta"))
#' 
#' model <- glm(D ~ X + X_sq, data = data, family = binomial(link = "logit"))
#' 
#' data$D_hat <- predict(model, newdata = data, type="response")
#' 
#' View(data)
#' 
#' @section References:
#' de Chaisemartin, C, D'Haultfoeuille, X, Pasquier, F, Vazquez‐Bare, G (2022). [Difference-in-Differences for Continuous Treatments and Instruments with Stayers](https://ssrn.com/abstract=4011782)
#' 
#' de Chaisemartin, C, D'Haultfoeuille, X (2020a) [Two-Way Fixed Effects Estimators with Heterogeneous Treatment Effects](https://cran.r-project.org/web/packages/DIDmultiplegt/index.html)
#' 
#' de Chaisemartin, C, D'Haultfoeuille, X (2020b) [Two-way fixed effects regressions with several treatments.](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3751060)
#' 
#' Li, S, Linn, J, Muehlegger, E (2014) [Gasoline Taxes and Consumer Behavior](https://www.aeaweb.org/articles?id=10.1257/pol.6.4.302)
#' 
#' @examples
#' # In the following example, we use data from Li et al. (2014). 
#' # The dataset can be downloaded from GitHub:
#' repo <- "chaisemartinPackages/ApplicationData/main" 
#' file <- "data_gazoline.dta"
#' url <- paste("https://raw.githubusercontent.com", repo, file, sep = "/")
#' gazoline <-  haven::read_dta(url)

#' # Estimating the effect of gasoline taxes on gasoline consumption and prices
#' summary(did_multiplegt_stat(
#'     df = gazoline,
#'     Y = "lngca",
#'     ID = "id",
#'     Time = "year",
#'     D = "tau",
#'     order = 1,
#'     estimator = "waoss",
#'     estimation_method = "dr",
#'     noextrapolation = TRUE
#' ))
#' @returns A list with two sublists. The first sublist includes the arguments used in the command. The second sublist includes the results from the estimation. Regardless of the options, the output in object$results$table will be a (3 x object$results$pairs) x 6 matrix, where only the requested output (that is, the rows corresponding to the estimators requested) will be non-missing. The list is given a did_multiplegt_stat class to trigger custom method dispatches by the print and summary functions. 
#' @export
did_multiplegt_stat <- function(
    df,
    Y,
    ID,
    Time,
    D,
    Z = NULL,
    estimator = NULL,
    estimation_method = NULL,
    order = 1,
    noextrapolation = FALSE,
    placebo = FALSE,
    switchers = NULL,
    disaggregate = FALSE,
    aoss_vs_waoss = FALSE,
    exact_match = FALSE,
    by = NULL,
    by_fd = NULL,
    other_treatments = NULL,
    cluster = NULL,
    twfe_suboptions= full_sample
) {
  
  params <- as.list(match.call())[-1]
  ## For now, the weight, cluster and by_fd options will be shut down until further theoretical results about the appropriate way to perform weighting and clustering while aggregating the IFs
  weight <- NULL
  
  args <- list()  
  for (v in names(formals(did_multiplegt_stat))) {
    if (!is.null(get(v))) {
      if (v == "df" & !inherits(get(v), "data.frame")) {
        stop(sprintf("Syntax error in %s option. Dataframe required.",v))
      } else if (v %in% c("Y", "ID", "Time", "D", "Z", "estimation_method", "switchers", "cluster")) {
        if (!(inherits(get(v), "character") & length(get(v)) == 1)) {
          stop(sprintf("Syntax error in %s option. The option requires a single string.",v))
        }
      } else if (v == "estimator") {
        if (!inherits(get(v), "character")) {
          stop(sprintf("Syntax error in %s option. Character vector required.",v))
        }
      } else if (v %in% c("noextrapolation", "placebo", "disaggregate", "aoss_vs_waoss", "exact_match")) {
        if (!inherits(get(v), "logical")) {
          stop(sprintf("Syntax error in %s option. Logical required.",v))
        }
      } else if (v %in% c("order", "by_fd")) {
        if (!(inherits(get(v), "numeric") & get(v) %% 1 == 0)) {
          stop(sprintf("Syntax error in %s option. Integer required.",v))
        }
      } else if (v == "by") {
        if (!(inherits(get(v), "character"))) {
          stop(sprintf("Syntax error in %s option. Character array required.",v))
        }
      }
    }
    if (v != "df") {
      args[[v]] <- get(v)
    } else {
      args$df <- params$df
    }
  }
  params <- NULL
  
  if (is.null(estimator) & is.null(Z)) {
    estimator <-  c("aoss", "waoss")
  } else if (is.null(estimator) & !is.null(Z) ) {
    estimator <- "ivwaoss"
  }
  
  if (!is.null(estimator) & length(intersect(estimator, c("aoss","waoss","ivwaoss"))) != length(estimator)) {
    stop("Syntax error in estimator option: only aoss, waoss and ivwaoss allowed.")
  }
  
  if (!is.null(switchers)) {
    if (!(switchers %in% c("up", "down"))) {
      stop("Switchers could be either NULL, up or down")          
    }
  }
  
  if (is.null(estimation_method)) {
    if (isFALSE(exact_match)) {
      estimation_method <- "dr"
    } else {
      estimation_method <- "ra"
    }
  }
  
  if (length(estimator) == 1) {
    if (estimator == "aoss") {
      estimation_method <- "ra"
    }
  }
  
  if (isTRUE(exact_match)) {
    if (estimation_method != "ra") {
      stop("The exact_match option is only compatible with regression adjustment estimation method")
    }
    if (isTRUE(noextrapolation)) {
      message("When exact_match and noextrapolation are both specified, the command will only consider the option exact_match.")
      noextrapolation <- FALSE
    }
    if (order != 1) {
      stop("The order option is not compatible with exact_match.")
    } else {
      order <- 1
    }
  }
  
  if (!(estimation_method %in% c("ra", "dr", "ps"))) {
    stop("Syntax error in estimation_method option.")
  }
  if (length(estimator) == 1) {
    if (estimation_method %in% c("dr","ps")  & estimator == "aoss") {
      stop("The propensity score-based approach is only available for the waoss and the ivwaoss.")
    }
  }
  
  if ("ivwaoss" %in% estimator & sum(c("aoss", "waoss") %in% estimator)) {
    stop("The estimation of AOSS or WAOSS cannot be combined with the estimation of IV-WAOSS (see helpfile).")
  }
  
  if (isTRUE(aoss_vs_waoss) & sum(c("aoss","waoss") %in% estimator) != 2) {
    stop("To test the equility between AOSS and WAOSS you must specify aoss and waoss in the estimator option.")
  }
  
  if ("ivwaoss" %in% estimator & is.null(Z)) {
    stop("To compute the ivwaoss you must specify the IV variable.")
  }
  
  did_multiplegt_stat <- list(args)
  names(did_multiplegt_stat) <- c("args")
  
  if (!is.null(by) | !is.null(by_fd)) {
    if (!is.null(by)) {
      df$by_total <- ""
      iter <- 1
      for (v in by) {
        if (!by_check(df, ID, v)) {
          stop("The ID variable should be nested within the by variable.")
        } 
        if (iter == 1) {
          df$by_total <- sprintf("%s", df[[v]])
        } else {
          df$by_total <- sprintf("%s,%s", df$by_total, df[[v]])
        }
        iter <- iter + 1
      }
      by_levels <- levels(factor(df$by_total))
      by_str <- by[1]
      if (length(by) > 1) {
        for (j in 2:length(by)) {
          by_str <- paste0(by_str,",", by[j])
        }      
      }
      
      did_multiplegt_stat <- append(did_multiplegt_stat, list(by_levels))
      
    } else if (!is.null(by_fd)) {
      if (100 %% by_fd != 0) {
        stop("Syntax error in by option. When the by option is specified with an integer, the argument must be divisible by 100 to allow for an integer subsetting of the quantiles.")
      }
      q_levels <- c(0) 
      for (k in 1:by_fd) {
        q_levels <- c(q_levels, q_levels[length(q_levels)] + 1/by_fd)
      }
      by_set <- did_multiplegt_stat_quantiles(df = df, ID = ID, Time = Time, D = D, Z = Z,by_opt = by_fd, quantiles = q_levels)
      df <- by_set$df
      val_quantiles <- by_set$val_quantiles
      quantiles <- by_set$quantiles
      switch_df <- by_set$switch_df
      quantiles_plot <- by_set$quantiles_plot
      by_set <- NULL
      by_levels <- levels(factor(subset(df, df$partition_XX != 0)$partition_XX))
      quantiles_mat <- as.matrix(rbind(quantiles, val_quantiles))
      did_multiplegt_stat <- append(did_multiplegt_stat, list(quantiles_mat))
      names(did_multiplegt_stat)[length(did_multiplegt_stat)] <- "quantiles"      
      did_multiplegt_stat <- append(did_multiplegt_stat, list(switch_df))
      names(did_multiplegt_stat)[length(did_multiplegt_stat)] <- "switchers_df"             
      did_multiplegt_stat <- append(did_multiplegt_stat, list(quantiles_plot))
      names(did_multiplegt_stat)[length(did_multiplegt_stat)] <- "cdf_plot"             
      quantiles_mat <- switch_df <- NULL
      
      if (length(by_levels) != by_fd) {
        message(sprintf("Point mass > %.0f%% detected. %.0f bin(s) collapsed.", 100/by_fd, by_fd - length(by_levels)))
      }
      
      did_multiplegt_stat <- append(did_multiplegt_stat, list(levels(factor(subset(df, df$partition_XX != 0)$partition_XX))))
    }
    names(did_multiplegt_stat)[length(did_multiplegt_stat)] <- "by_levels"
  } else {
    by_levels <- "_no_by"
  } 
  did_multiplegt_stat$args$estimator <- estimator
  
  df_main <- df
  obj_name <- "results"
  by_fd_opt <- NULL
  
  for (by_lev in 1:length(by_levels)) {
    if (by_levels[by_lev] != "_no_by" & !is.null(by)) {
      df_main <- subset(df, df$by_total == by_levels[by_lev])
      obj_name <- paste0("results_by_", by_lev)
      message(sprintf("Running did_multiplegt_stat with %s = %s", by_str, by_levels[by_lev]))
    } else if (by_levels[by_lev] != "_no_by" & !is.null(by_fd)) {
      obj_name <- paste0("results_by_", by_lev)
      diff_var <- ifelse("ivwaoss" %in% estimator, "Z", "D")
      sep <- ifelse(by_lev == 1, "[", "(")
      message(sprintf("Running did_multiplegt_stat with switchers s.t. \U0394%s \U2208 %s%.3f,%.3f] <%.0f%%-%.0f%% quantiles>.", diff_var, sep, val_quantiles[by_lev], val_quantiles[by_lev + 1], quantiles[by_lev] * 100, quantiles[by_lev+1]*100))
      if (val_quantiles[by_lev] == val_quantiles[by_lev + 1])  {
        warning(sprintf("(%.0f%%, %0.f%%) quantile bin dropped: upper and lower bounds are equal.", quantiles[by_lev] * 100, quantiles[by_lev+1]*100))
      }
      by_fd_opt <- by_levels[by_lev]
    }
    results <- did_multiplegt_stat_main(df = df_main, Y = Y, ID = ID, Time = Time, D = D, Z = Z, estimator = estimator, estimation_method = estimation_method, order = order, noextrapolation = noextrapolation, placebo = placebo, weight = weight, switchers = switchers, disaggregate = disaggregate, aoss_vs_waoss = aoss_vs_waoss, exact_match = exact_match, cluster = cluster, by_fd_opt = by_fd_opt, other_treatments = other_treatments)
    did_multiplegt_stat <- append(did_multiplegt_stat, list(results))
    names(did_multiplegt_stat)[length(did_multiplegt_stat)] <- obj_name
  }
  
  if (!is.null(did_multiplegt_stat$args[["by"]])) {
    did_multiplegt_stat <- append(did_multiplegt_stat, list(by_graph(obj = did_multiplegt_stat)))
    names(did_multiplegt_stat)[length(did_multiplegt_stat)] <- "by_graph"
  }
  if (!is.null(did_multiplegt_stat$quantiles)) {
    did_multiplegt_stat <- append(did_multiplegt_stat, list(by_fd_graph(obj = did_multiplegt_stat)))
    names(did_multiplegt_stat)[length(did_multiplegt_stat)] <- "by_fd_graph"
  }
  
  class(did_multiplegt_stat) <- "did_multiplegt_stat"
  return(did_multiplegt_stat)
}
#' Internal function of did_multiplegt_stat for the computation of displayed results.
#' @param df df
#' @param Y Y
#' @param ID ID
#' @param Time Time
#' @param D D
#' @param Z Z
#' @param estimator estimator
#' @param estimation_method estimation_method
#' @param order order
#' @param noextrapolation noextrapolation
#' @param placebo placebo
#' @param switchers switchers
#' @param disaggregate disaggregate
#' @param aoss_vs_waoss aoss_vs_waoss
#' @param exact_match exact_match
#' @param weight weight
#' @param cluster cluster
#' @param by_fd_opt by_fd_opt
#' @param other_treatments other_treatments
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang := 
#' @importFrom rlang .data
#' @importFrom plm pdata.frame make.pbalanced
#' @importFrom stats sd pnorm
#' @returns A list with the raw output to be displayed.
#' @noRd
did_multiplegt_stat_main <- function(
    df,
    Y,
    ID,
    Time,
    D,
    Z,
    estimator,
    estimation_method,
    order,
    noextrapolation,
    placebo,
    switchers,
    disaggregate,
    aoss_vs_waoss,
    exact_match,
    weight,
    cluster,
    by_fd_opt,
    other_treatments
) {
  suppressWarnings({
    # Preallocation of scalars
    aoss_XX <- NULL
    waoss_XX <- NULL
    ivwaoss_XX <- NULL
    
    # Patching the estimator option
    for (v in c("aoss", "waoss", "ivwaoss")) {
      assign(paste0(v,"_XX"), as.numeric(v %in% estimator))
    }
    
    # Layer 1: keep only variables of interest, as to speed up what follows
    varlist <- c()
    for (v in c(Y, ID, Time, D, Z, weight, cluster, other_treatments)) {
      if (!is.null(v)) {
        if (!(v %in% varlist)) {
          varlist <- c(varlist, v)
        }
      }
    }
    if (!is.null(df$partition_XX)) {
      varlist <- c(varlist, "partition_XX")
    }
    df <- df[varlist]
    df_base <- list(Y = Y, ID = ID, T = Time, D = D, Z = Z, weight = weight, cluster = cluster)
    for (i in 1:length(df_base)) {
      if (!is.null(df_base[[i]])) {
        col <- as.character(df_base[[i]])
        df[[paste0(names(df_base)[i],"_XX")]] <- df[[col]]
      }
    }
    df_base <- NULL
    
    if (!is.null(cluster)) {
      if (cluster == ID) {
        cluster <- NULL
        message("The cluster option should be different from (and coarser than) the ID variable. The command will ignore the cluster option.")
      } else {
        df <- df %>% group_by(.data$ID_XX) %>% 
          mutate(cluster_sd_XX = sd(.data$cluster_XX, na.rm = TRUE))
        if (max(df$cluster_sd_XX, na.rm = TRUE) > 0) {
          stop("The ID variable should be nested within the clustering variable.")
        } else {
          n_clus_XX <- length(unique(df$cluster_XX))
        }
      }
    }
    
    df$to_drop_XX <- (is.na(df$T_XX) | is.na(df$D_XX) | is.na(df$ID_XX))
    IV_req_XX <- 0
    if (ivwaoss_XX == 1) {
      df$to_drop_XX <- (is.na(df$Z_XX) | df$to_drop_XX)
      IV_req_XX <- 1
    }
    df <- subset(df, df$to_drop_XX == 0)
    
    # Layer 2: balancing the panel and then keeping again only variables of interest
    df$tsfilled_XX <- 0
    df <- pdata.frame(df, index = c("ID_XX", "T_XX")) 
    df <- make.pbalanced(df, balance.type = "fill")
    df$tsfilled_XX <- is.na(df$tsfilled_XX)
    df$T_temp_XX <- NULL; df$T_temp_XX <- df$T_XX; df$T_XX <- NULL; 
    df <- df %>% group_by(.data$T_temp_XX) %>% mutate(T_XX = cur_group_id())
    df$T_temp_XX <- NULL   
    
    if (is.null(weight)) {
      df$weight_XX <- 1
      df$weight_c_XX <- 1
    } else {
      df$weight_XX <- ifelse(is.na(df$weight_XX), 0, df$weight_XX)
    }
    if (!is.null(cluster)) {
      df <- df %>% group_by(.data$cluster_XX, .data$T_XX) %>%
        mutate(weight_c_XX = sum(.data$weight_XX, na.rm = TRUE))
    }
    
    # Further useful steps prior to the estimation
    IDs_XX <- as.data.frame(unique(factor(df$ID_XX)))
    names(IDs_XX) <- "ID_XX"
    if (!is.null(cluster)) {
      cluster_df <- df %>% group_by(.data$ID_XX) %>% 
        summarise(cluster_XX = mean(.data$cluster_XX)) %>% ungroup()
      IDs_XX <- merge(IDs_XX, cluster_df, by = "ID_XX")
      cluster_df <- NULL
    }
    
    max_T_XX <- max(df$T_XX, na.rm = TRUE)
    scalars <- list(
      PS_sum_XX = 0,
      delta_1_1_XX = 0,
      E_abs_delta_D_sum_XX = 0,
      delta_2_1_XX = 0,
      denom_delta_IV_sum_XX = 0,
      delta_3_1_XX = 0,
      N_Switchers_1_1_XX = 0,
      N_Stayers_1_1_XX = 0,
      N_Switchers_2_1_XX = 0,
      N_Stayers_2_1_XX = 0,
      
      N_Switchers_3_1_XX = 0,
      N_Stayers_3_1_XX = 0,
      denom_delta_IV_sum_XX = 0,
      N_drop_total_XX = 0,
      N_drop_total_C_XX = 0,
      IV_req_XX = IV_req_XX
    )
    if (isTRUE(placebo)) {
      scalars <- c(scalars,
                   PS_sum_pl_XX = 0,
                   delta_1_1_pl_XX = 0,
                   E_abs_delta_D_sum_pl_XX = 0,
                   delta_2_1_pl_XX = 0,
                   denom_delta_IV_sum_pl_XX = 0,
                   delta_3_1_pl_XX = 0,
                   N_Switchers_1_1_pl_XX = 0,
                   N_Stayers_1_1_pl_XX = 0,
                   N_Switchers_2_1_pl_XX = 0,
                   N_Stayers_2_1_pl_XX = 0,
                   N_Switchers_3_1_pl_XX = 0,
                   N_Stayers_3_1_pl_XX = 0,
                   denom_delta_IV_sum_pl_XX = 0)
    }
    
    ## Computing E(H_t)
    balanced_df <- balanced_map(df)
    
    for (p in 2:max_T_XX) {
      
      est_out <- did_multiplegt_stat_pairwise(df = df, Y = "Y_ID", ID = "ID_XX", Time = "T_XX", D = "D_XX", Z = "Z_XX", estimator = estimator, order = order, noextrapolation = noextrapolation, weight = "weight_XX", switchers = switchers, pairwise = p, aoss = aoss_XX, waoss = waoss_XX, ivwaoss = ivwaoss_XX, estimation_method = estimation_method, scalars = scalars, placebo = FALSE, exact_match = exact_match, cluster = cluster, by_fd_opt = by_fd_opt, other_treatments = other_treatments)
      
      IDs_XX <- merge(IDs_XX, est_out$to_add, by = "ID_XX", all = TRUE) 
      IDs_XX <- IDs_XX[order(IDs_XX$ID_XX), ]
      scalars <- est_out$scalars;
      est_out <- NULL;
      
      if (aoss_XX == 1) {
        scalars$delta_1_1_XX <- scalars$delta_1_1_XX + 
          scalars[[paste0("P_",p,"_XX")]] * scalars[[paste0("delta_1_",p,"_XX")]]
        
        if (scalars[[paste0("N_Stayers_1_",p,"_XX")]] > 1 & !is.na(scalars[[paste0("N_Stayers_1_",p,"_XX")]]))  {
          scalars$N_Switchers_1_1_XX <- scalars$N_Switchers_1_1_XX + scalars[[paste0("N_Switchers_1_",p,"_XX")]]
        }
        if (scalars[[paste0("N_Switchers_1_",p,"_XX")]] > 0 & !is.na(scalars[[paste0("N_Switchers_1_",p,"_XX")]]))  {
          scalars$N_Stayers_1_1_XX <- scalars$N_Stayers_1_1_XX + scalars[[paste0("N_Stayers_1_",p,"_XX")]]
        }
      }
      
      if (waoss_XX == 1) {
        scalars$delta_2_1_XX <- scalars$delta_2_1_XX + 
          scalars[[paste0("E_abs_delta_D_",p,"_XX")]] * scalars[[paste0("delta_2_",p,"_XX")]]
        
        if (scalars[[paste0("N_Stayers_2_",p,"_XX")]] > 1 & !is.na(scalars[[paste0("N_Stayers_2_",p,"_XX")]]))  {
          scalars$N_Switchers_2_1_XX <- scalars$N_Switchers_2_1_XX + scalars[[paste0("N_Switchers_2_",p,"_XX")]]
        }
        if (scalars[[paste0("N_Switchers_2_",p,"_XX")]] > 0 & !is.na(scalars[[paste0("N_Switchers_2_",p,"_XX")]]))  {
          scalars$N_Stayers_2_1_XX <- scalars$N_Stayers_2_1_XX + scalars[[paste0("N_Stayers_2_",p,"_XX")]]
        }
      }
      
      if (ivwaoss_XX == 1) {
        scalars$delta_3_1_XX <- scalars$delta_3_1_XX + 
          scalars[[paste0("denom_delta_IV_",p,"_XX")]] * scalars[[paste0("delta_3_",p,"_XX")]]
        
        if (scalars[[paste0("N_Stayers_3_",p,"_XX")]] > 1 & !is.na(scalars[[paste0("N_Stayers_3_",p,"_XX")]]))  {
          scalars$N_Switchers_3_1_XX <- scalars$N_Switchers_3_1_XX + scalars[[paste0("N_Switchers_3_",p,"_XX")]]
        }
        if (scalars[[paste0("N_Switchers_3_",p,"_XX")]] > 0 & !is.na(scalars[[paste0("N_Switchers_3_",p,"_XX")]]))  {
          scalars$N_Stayers_3_1_XX <- scalars$N_Stayers_3_1_XX + scalars[[paste0("N_Stayers_3_",p,"_XX")]]
        }
      }
    }
    
    if (isTRUE(placebo)) {
      for (p in 3:max_T_XX) {
        
        est_out <- did_multiplegt_stat_pairwise(df = df, Y = "Y_ID", ID = "ID_XX", Time = "T_XX", D = "D_XX", Z = "Z_XX", estimator = estimator, order = order, noextrapolation = noextrapolation, weight = "weight_XX", switchers = switchers, pairwise = p, aoss = aoss_XX, waoss = waoss_XX, ivwaoss = ivwaoss_XX, estimation_method = estimation_method, scalars = scalars, placebo = TRUE, exact_match = exact_match, cluster = cluster, by_fd_opt = by_fd_opt, other_treatments = other_treatments)
        
        if (!is.null(est_out$to_add)) {
          IDs_XX <- merge(IDs_XX, est_out$to_add, by = "ID_XX", all = TRUE) 
        }
        IDs_XX <- IDs_XX[order(IDs_XX$ID_XX), ]
        scalars <- est_out$scalars;
        est_out <- NULL;
        
        if (aoss_XX == 1) {
          scalars$delta_1_1_pl_XX <- scalars$delta_1_1_pl_XX + 
            scalars[[paste0("P_",p,"_pl_XX")]] * scalars[[paste0("delta_1_",p,"_pl_XX")]]
          
          
          if (scalars[[paste0("N_Stayers_1_",p,"_pl_XX")]] > 1 & !is.na(paste0("N_Stayers_1_",p,"_pl_XX")))  {
            scalars$N_Switchers_1_1_pl_XX <- scalars$N_Switchers_1_1_pl_XX + scalars[[paste0("N_Switchers_1_",p,"_pl_XX")]]
          }
          if (scalars[[paste0("N_Switchers_1_",p,"_pl_XX")]] > 0 & !is.na(paste0("N_Switchres_1_",p,"_pl_XX")))  {
            scalars$N_Stayers_1_1_pl_XX <- scalars$N_Stayers_1_1_pl_XX + scalars[[paste0("N_Stayers_1_",p,"_pl_XX")]]
          }
        }
        
        if (waoss_XX == 1) {
          scalars$delta_2_1_pl_XX <- scalars$delta_2_1_pl_XX + 
            scalars[[paste0("E_abs_delta_D_",p,"_pl_XX")]] * scalars[[paste0("delta_2_",p,"_pl_XX")]]
          
          if (scalars[[paste0("N_Stayers_2_",p,"_pl_XX")]] > 1 & !is.na(paste0("N_Stayers_2_",p,"_pl_XX")))  {
            scalars$N_Switchers_2_1_pl_XX <- scalars$N_Switchers_2_1_pl_XX + scalars[[paste0("N_Switchers_2_",p,"_pl_XX")]]
          }
          if (scalars[[paste0("N_Switchers_2_",p,"_pl_XX")]] > 0 & !is.na(paste0("N_Switchers_2_",p,"_pl_XX")))  {
            scalars$N_Stayers_2_1_pl_XX <- scalars$N_Stayers_2_1_pl_XX + scalars[[paste0("N_Stayers_2_",p,"_pl_XX")]]
          }
        }
        
        if (ivwaoss_XX == 1) {
          scalars$delta_3_1_pl_XX <- scalars$delta_3_1_pl_XX + 
            scalars[[paste0("denom_delta_IV_",p,"_pl_XX")]] * scalars[[paste0("delta_3_",p,"_pl_XX")]]
          
          if (scalars[[paste0("N_Stayers_3_",p,"_pl_XX")]] > 1 & !is.na(paste0("N_Stayers_3_",p,"_pl_XX")))  {
            scalars$N_Switchers_3_1_pl_XX <- scalars$N_Switchers_3_1_pl_XX + scalars[[paste0("N_Switchers_3_",p,"_pl_XX")]]
          }
          if (scalars[[paste0("N_Switchers_3_",p,"_pl_XX")]] > 0 & !is.na(paste0("N_Switchers_3_",p,"_pl_XX")))  {
            scalars$N_Stayers_3_1_pl_XX <- scalars$N_Stayers_3_1_pl_XX + scalars[[paste0("N_Stayers_3_",p,"_pl_XX")]]
          }
        }
      }
    }
    
    # Compute the aggregated estimators
    if (aoss_XX == 1) {
      scalars$delta_1_1_XX <- scalars$delta_1_1_XX / scalars$PS_sum_XX
      if (isTRUE(placebo)) {
        scalars$delta_1_1_pl_XX <- scalars$delta_1_1_pl_XX / scalars$PS_sum_pl_XX
      }
    }
    
    if (waoss_XX == 1) {
      scalars$delta_2_1_XX <- scalars$delta_2_1_XX / scalars$E_abs_delta_D_sum_XX
      if (isTRUE(placebo)) {
        scalars$delta_2_1_pl_XX <- scalars$delta_2_1_pl_XX / scalars$E_abs_delta_D_sum_pl_XX
      }
    }
    
    if (ivwaoss_XX == 1) {
      scalars$delta_3_1_XX <- scalars$delta_3_1_XX / scalars$denom_delta_IV_sum_XX
      if (isTRUE(placebo)) {
        scalars$delta_3_1_pl_XX <- scalars$delta_3_1_pl_XX / scalars$denom_delta_IV_sum_pl_XX
      }
    }
    
    # Compute the influence functions
    for (i in 1:3) {
      for (pl in c("","_pl")) {
        IDs_XX[[paste0("Phi_",i,pl,"_XX")]] <- 0
      }
    }
    counter_XX <- 0
    
    IDs_XX <- IDs_XX[order(IDs_XX$ID_XX), ]
    for (p in 2:max_T_XX) {
      if (aoss_XX == 1 & scalars[[paste0("non_missing_",p,"_XX")]] == 1) {
        IDs_XX[[paste0("Phi_1_",p,"_XX")]] <- (scalars[[paste0("P_",p,"_XX")]]*IDs_XX[[paste0("Phi_1_",p,"_XX")]] + (scalars[[paste0("delta_1_",p,"_XX")]] - scalars$delta_1_1_XX) * (IDs_XX[[paste0("S_",p,"_XX")]] - scalars[[paste0("P_",p,"_XX")]])) / scalars$PS_sum_XX
        
        IDs_XX$Phi_1_XX <- ifelse(is.na(IDs_XX[[paste0("Phi_1_",p,"_XX")]]), IDs_XX$Phi_1_XX, IDs_XX$Phi_1_XX + IDs_XX[[paste0("Phi_1_",p,"_XX")]])
        
        if (isTRUE(placebo) & p > 2) {
          if (scalars[[paste0("non_missing_",p,"_pl_XX")]] == 1) {
            IDs_XX[[paste0("Phi_1_",p,"_pl_XX")]] <- (scalars[[paste0("P_",p,"_pl_XX")]]*IDs_XX[[paste0("Phi_1_",p,"_pl_XX")]] + (scalars[[paste0("delta_1_",p,"_pl_XX")]] - scalars$delta_1_1_pl_XX) * (IDs_XX[[paste0("S_",p,"_pl_XX")]] - scalars[[paste0("P_",p,"_pl_XX")]])) / scalars$PS_sum_pl_XX
            
            IDs_XX$Phi_1_pl_XX <- ifelse(is.na(IDs_XX[[paste0("Phi_1_",p,"_pl_XX")]]), IDs_XX$Phi_1_pl_XX, IDs_XX$Phi_1_pl_XX + IDs_XX[[paste0("Phi_1_",p,"_pl_XX")]])
          }
        }
      }
      if (waoss_XX == 1 & scalars[[paste0("non_missing_",p,"_XX")]] == 1) {
        IDs_XX[[paste0("Phi_2_",p,"_XX")]] <- (scalars[[paste0("E_abs_delta_D_",p,"_XX")]]*IDs_XX[[paste0("Phi_2_",p,"_XX")]] + (scalars[[paste0("delta_2_",p,"_XX")]] - scalars$delta_2_1_XX) * (IDs_XX[[paste0("abs_delta_D_",p,"_XX")]] - scalars[[paste0("E_abs_delta_D_",p,"_XX")]])) / scalars$E_abs_delta_D_sum_XX
        
        IDs_XX$Phi_2_XX <- ifelse(is.na(IDs_XX[[paste0("Phi_2_",p,"_XX")]]), IDs_XX$Phi_2_XX, IDs_XX$Phi_2_XX + IDs_XX[[paste0("Phi_2_",p,"_XX")]])
        
        if (isTRUE(placebo) & p > 2) {
          if (scalars[[paste0("non_missing_",p,"_pl_XX")]] == 1) {
            IDs_XX[[paste0("Phi_2_",p,"_pl_XX")]] <- (scalars[[paste0("E_abs_delta_D_",p,"_pl_XX")]]*IDs_XX[[paste0("Phi_2_",p,"_pl_XX")]] + (scalars[[paste0("delta_2_",p,"_pl_XX")]] - scalars$delta_2_1_pl_XX) * (IDs_XX[[paste0("abs_delta_D_",p,"_pl_XX")]] - scalars[[paste0("E_abs_delta_D_",p,"_pl_XX")]])) / scalars$E_abs_delta_D_sum_pl_XX
            
            
            IDs_XX$Phi_2_pl_XX <- ifelse(is.na(IDs_XX[[paste0("Phi_2_",p,"_pl_XX")]]), IDs_XX$Phi_2_pl_XX, IDs_XX$Phi_2_pl_XX + IDs_XX[[paste0("Phi_2_",p,"_pl_XX")]])
          }
        }
      }
      
      if (ivwaoss_XX == 1 & scalars[[paste0("non_missing_",p,"_XX")]] == 1) {
        
        IDs_XX[[paste0("Phi_3_",p,"_XX")]] <- ((scalars[[paste0("denom_delta_IV_",p,"_XX")]] * IDs_XX[[paste0("Phi_3_",p,"_XX")]]) + (scalars[[paste0("delta_3_",p,"_XX")]] - scalars$delta_3_1_XX) * (IDs_XX[[paste0("inner_sum_IV_denom_",p,"_XX")]] - scalars[[paste0("denom_delta_IV_",p,"_XX")]])) / scalars$denom_delta_IV_sum_XX
        
        IDs_XX$Phi_3_XX <- ifelse(is.na(IDs_XX[[paste0("Phi_3_",p,"_XX")]]), IDs_XX$Phi_3_XX, IDs_XX$Phi_3_XX + IDs_XX[[paste0("Phi_3_",p,"_XX")]])
        
        if (isTRUE(placebo) & p > 2) {
          if (scalars[[paste0("non_missing_",p,"_pl_XX")]] == 1) {
            IDs_XX[[paste0("Phi_3_",p,"_pl_XX")]] <- (scalars[[paste0("denom_delta_IV_",p,"_pl_XX")]] * IDs_XX[[paste0("Phi_3_",p,"_pl_XX")]] + (scalars[[paste0("delta_3_",p,"_pl_XX")]] - scalars$delta_3_1_pl_XX) * (IDs_XX[[paste0("inner_sum_IV_denom_",p,"_pl_XX")]] - scalars[[paste0("denom_delta_IV_",p,"_pl_XX")]])) / scalars$denom_delta_IV_sum_pl_XX
            
            IDs_XX$Phi_3_pl_XX <- ifelse(is.na(IDs_XX[[paste0("Phi_3_",p,"_pl_XX")]]), IDs_XX$Phi_3_pl_XX, IDs_XX$Phi_3_pl_XX + IDs_XX[[paste0("Phi_3_",p,"_pl_XX")]])
          }
        }
      }
      
      counter_XX <- counter_XX + 1
    }
    
    if (!is.null(cluster)) {
      df <- df %>% group_by(.data$ID_XX) %>%
        mutate(gr_id = row_number()) %>% ungroup()
      df$id_temp <- as.numeric(df$gr_id == 1)
      df <- df %>% group_by(.data$cluster_XX) %>%
        mutate(N_c_XX = sum(.data$id_temp, na.rm = TRUE)) %>% ungroup()
      N_bar_c_XX <-  mean(df$N_c_XX, na.rm = TRUE)
      df$id_temp <- df$N_c_XX <- df$gr_id <- NULL
    }
    
    if (aoss_XX == 1) {
      scalars$mean_IF1 <- ifelse(counter_XX == 0, NA, mean(IDs_XX$Phi_1_XX, na.rm = TRUE))
      if (!is.null(cluster)) {
        IDs_XX <- IDs_XX %>% group_by(.data$cluster_XX) %>% 
          mutate(Phi_1_c_XX = sum(.data$Phi_1_XX, na.rm = TRUE)) %>%
          mutate(first_obs_by_clus = row_number() == 1) %>% ungroup()
        IDs_XX$Phi_1_c_XX <- ifelse(IDs_XX$first_obs_by_clus == 1, IDs_XX$Phi_1_c_XX, NA) / N_bar_c_XX
        n_obs <- nrow(subset(IDs_XX, !is.na(IDs_XX$Phi_1_c_XX)))
        scalars$sd_delta_1_1_XX <- ifelse(counter_XX == 0, NA, sd(IDs_XX$Phi_1_c_XX, na.rm = TRUE)/ sqrt(n_obs))
        df$first_obs_by_clus <- NULL
      } else {
        n_obs <- nrow(subset(IDs_XX, !is.na(IDs_XX$Phi_1_XX)))
        scalars$sd_delta_1_1_XX <- ifelse(counter_XX == 0, NA, sd(IDs_XX$Phi_1_XX, na.rm = TRUE)/ sqrt(n_obs))
      }
      scalars$LB_1_1_XX <-  scalars$delta_1_1_XX - 1.96 * scalars$sd_delta_1_1_XX
      scalars$UB_1_1_XX <-  scalars$delta_1_1_XX + 1.96 * scalars$sd_delta_1_1_XX            
      
      if (isTRUE(placebo)) {
        scalars$mean_IF1_pl <- ifelse(counter_XX == 0, NA, mean(IDs_XX$Phi_1_pl_XX, na.rm = TRUE))
        if (!is.null(cluster)) {
          IDs_XX <- IDs_XX %>% group_by(.data$cluster_XX) %>% 
            mutate(Phi_1_pl_c_XX = sum(.data$Phi_1_pl_XX, na.rm = TRUE)) %>%
            mutate(first_obs_by_clus = row_number() == 1) %>% ungroup()
          IDs_XX$Phi_1_pl_c_XX <- ifelse(IDs_XX$first_obs_by_clus == 1, IDs_XX$Phi_1_pl_c_XX, NA) / N_bar_c_XX
          n_obs_pl <- nrow(subset(IDs_XX, !is.na(IDs_XX$Phi_1_pl_c_XX)))
          scalars$sd_delta_1_1_pl_XX <- ifelse(counter_XX == 0, NA, sd(IDs_XX$Phi_1_pl_c_XX, na.rm = TRUE)/ sqrt(n_obs_pl))
          df$first_obs_by_clus <- NULL
        } else {
          n_obs_pl <- nrow(subset(IDs_XX, !is.na(IDs_XX$Phi_1_pl_XX)))
          scalars$sd_delta_1_1_pl_XX <- ifelse(counter_XX == 0, NA, sd(IDs_XX$Phi_1_pl_XX, na.rm = TRUE)/ sqrt(n_obs_pl))
        }
        scalars$LB_1_1_pl_XX <-  scalars$delta_1_1_pl_XX - 1.96 * scalars$sd_delta_1_1_pl_XX
        scalars$UB_1_1_pl_XX <-  scalars$delta_1_1_pl_XX + 1.96 * scalars$sd_delta_1_1_pl_XX            
      }
    }
    
    if (waoss_XX == 1) {
      scalars$mean_IF1 <- ifelse(counter_XX == 0, NA, mean(IDs_XX$Phi_2_XX, na.rm = TRUE))
      if (!is.null(cluster)) {
        IDs_XX <- IDs_XX %>% group_by(.data$cluster_XX) %>% 
          mutate(Phi_2_c_XX = sum(.data$Phi_2_XX, na.rm = TRUE)) %>%
          mutate(first_obs_by_clus = row_number() == 1) %>% ungroup()
        IDs_XX$Phi_2_c_XX <- ifelse(IDs_XX$first_obs_by_clus == 1, IDs_XX$Phi_2_c_XX, NA) / N_bar_c_XX
        n_obs <- nrow(subset(IDs_XX, !is.na(IDs_XX$Phi_2_c_XX)))
        scalars$sd_delta_2_1_XX <- ifelse(counter_XX == 0, NA, sd(IDs_XX$Phi_2_c_XX, na.rm = TRUE)/ sqrt(n_obs))
        df$first_obs_by_clus <- NULL
      } else {
        n_obs <- nrow(subset(IDs_XX, !is.na(IDs_XX$Phi_2_XX)))
        scalars$sd_delta_2_1_XX <- ifelse(counter_XX == 0, NA, sd(IDs_XX$Phi_2_XX, na.rm = TRUE)/ sqrt(n_obs))
      }
      scalars$LB_2_1_XX <-  scalars$delta_2_1_XX - 1.96 * scalars$sd_delta_2_1_XX
      scalars$UB_2_1_XX <-  scalars$delta_2_1_XX + 1.96 * scalars$sd_delta_2_1_XX            
      
      if (isTRUE(placebo)) {
        n_obs_pl <- nrow(subset(IDs_XX, !is.na(IDs_XX$Phi_2_pl_XX)))
        scalars$mean_IF1_pl <- ifelse(counter_XX == 0, NA, mean(IDs_XX$Phi_2_pl_XX, na.rm = TRUE))
        if (!is.null(cluster)) {
          IDs_XX <- IDs_XX %>% group_by(.data$cluster_XX) %>% 
            mutate(Phi_2_pl_c_XX = sum(.data$Phi_2_pl_XX, na.rm = TRUE)) %>%
            mutate(first_obs_by_clus = row_number() == 1) %>% ungroup()
          IDs_XX$Phi_2_pl_c_XX <- ifelse(IDs_XX$first_obs_by_clus == 1, IDs_XX$Phi_2_pl_c_XX, NA) / N_bar_c_XX
          n_obs_pl <- nrow(subset(IDs_XX, !is.na(IDs_XX$Phi_2_pl_c_XX)))
          scalars$sd_delta_2_1_pl_XX <- ifelse(counter_XX == 0, NA, sd(IDs_XX$Phi_2_pl_c_XX, na.rm = TRUE)/ sqrt(n_obs_pl))
          df$first_obs_by_clus <- NULL
        } else {
          n_obs_pl <- nrow(subset(IDs_XX, !is.na(IDs_XX$Phi_2_pl_XX)))
          scalars$sd_delta_2_1_pl_XX <- ifelse(counter_XX == 0, NA, sd(IDs_XX$Phi_2_pl_XX, na.rm = TRUE)/ sqrt(n_obs_pl))
        }
        scalars$LB_2_1_pl_XX <-  scalars$delta_2_1_pl_XX - 1.96 * scalars$sd_delta_2_1_pl_XX
        scalars$UB_2_1_pl_XX <-  scalars$delta_2_1_pl_XX + 1.96 * scalars$sd_delta_2_1_pl_XX            
      }
    }
    
    if (ivwaoss_XX == 1) {
      scalars$mean_IF3 <- ifelse(counter_XX == 0, NA, mean(IDs_XX$Phi_3_XX, na.rm = TRUE))
      if (!is.null(cluster)) {
        IDs_XX <- IDs_XX %>% group_by(.data$cluster_XX) %>% 
          mutate(Phi_3_c_XX = sum(.data$Phi_3_XX, na.rm = TRUE)) %>%
          mutate(first_obs_by_clus = row_number() == 1) %>% ungroup()
        IDs_XX$Phi_3_c_XX <- ifelse(IDs_XX$first_obs_by_clus == 1, IDs_XX$Phi_3_c_XX, NA) / N_bar_c_XX
        n_obs <- nrow(subset(IDs_XX, !is.na(IDs_XX$Phi_3_c_XX)))
        scalars$sd_delta_3_1_XX <- ifelse(counter_XX == 0, NA, sd(IDs_XX$Phi_3_c_XX, na.rm = TRUE)/ sqrt(n_obs))
        df$first_obs_by_clus <- NULL
      } else {
        n_obs <- nrow(subset(IDs_XX, !is.na(IDs_XX$Phi_3_XX)))
        scalars$sd_delta_3_1_XX <- ifelse(counter_XX == 0, NA, sd(IDs_XX$Phi_3_XX, na.rm = TRUE)/ sqrt(n_obs))
      }
      scalars$LB_3_1_XX <-  scalars$delta_3_1_XX - 1.96 * scalars$sd_delta_3_1_XX
      scalars$UB_3_1_XX <-  scalars$delta_3_1_XX + 1.96 * scalars$sd_delta_3_1_XX            
      
      if (isTRUE(placebo)) {
        n_obs_pl <- nrow(subset(IDs_XX, !is.na(IDs_XX$Phi_3_pl_XX)))
        scalars$mean_IF3_pl <- ifelse(counter_XX == 0, NA, mean(IDs_XX$Phi_3_pl_XX, na.rm = TRUE))
        if (!is.null(cluster)) {
          IDs_XX <- IDs_XX %>% group_by(.data$cluster_XX) %>% 
            mutate(Phi_3_pl_c_XX = sum(.data$Phi_3_pl_XX, na.rm = TRUE)) %>%
            mutate(first_obs_by_clus = row_number() == 1) %>% ungroup()
          IDs_XX$Phi_3_pl_c_XX <- ifelse(IDs_XX$first_obs_by_clus == 1, IDs_XX$Phi_3_pl_c_XX, NA) / N_bar_c_XX
          n_obs_pl <- nrow(subset(IDs_XX, !is.na(IDs_XX$Phi_3_pl_c_XX)))
          scalars$sd_delta_3_1_pl_XX <- ifelse(counter_XX == 0, NA, sd(IDs_XX$Phi_3_pl_c_XX, na.rm = TRUE)/ sqrt(n_obs_pl))
          df$first_obs_by_clus <- NULL
        } else {
          n_obs_pl <- nrow(subset(IDs_XX, !is.na(IDs_XX$Phi_3_pl_XX)))
          scalars$sd_delta_3_1_pl_XX <- ifelse(counter_XX == 0, NA, sd(IDs_XX$Phi_3_pl_XX, na.rm = TRUE)/ sqrt(n_obs_pl))
        }
        scalars$LB_3_1_pl_XX <-  scalars$delta_3_1_pl_XX - 1.96 * scalars$sd_delta_3_1_pl_XX
        scalars$UB_3_1_pl_XX <-  scalars$delta_3_1_pl_XX + 1.96 * scalars$sd_delta_3_1_pl_XX            
      }
    }
    
    
    # AOSS vs WAOSS
    if (isTRUE(aoss_vs_waoss)) {
      diff_delta_1_2_XX <- scalars$delta_1_1_XX - scalars$delta_2_1_XX
      IDs_XX$diff_Phi_1_2_XX <- IDs_XX$Phi_1_XX - IDs_XX$Phi_2_XX
      if (!is.null(cluster)) {
        IDs_XX <- IDs_XX %>% group_by(.data$cluster_XX) %>% 
          mutate(diff_Phi_1_2_c_XX = sum(.data$diff_Phi_1_2_XX, na.rm = TRUE)) %>%
          mutate(first_obs_by_clus = row_number() == 1) %>% ungroup()
        IDs_XX$diff_Phi_1_2_c_XX <- ifelse(IDs_XX$first_obs_by_clus == 1, IDs_XX$diff_Phi_1_2_c_XX, NA)
        sd_diff_delta_1_2_XX <- sd(IDs_XX$diff_Phi_1_2_c_XX, na.rm = TRUE)
        df$first_obs_by_clus <- NULL
      } else {
        sd_diff_delta_1_2_XX <- sd(IDs_XX$diff_Phi_1_2_XX, na.rm = TRUE)
      }
      
      t_diff_delta_1_2_XX <- diff_delta_1_2_XX * sqrt(nrow(IDs_XX)) / sd_diff_delta_1_2_XX
      p_diff_delta_1_2_XX <- 2 * (1 - pnorm(abs(t_diff_delta_1_2_XX)))
      LB_diff_delta_1_2_XX <- diff_delta_1_2_XX - 1.96 * sd_diff_delta_1_2_XX / sqrt(nrow(IDs_XX))
      UB_diff_delta_1_2_XX <- diff_delta_1_2_XX + 1.96 * sd_diff_delta_1_2_XX / sqrt(nrow(IDs_XX))
      
      t_mat <- matrix(0, nrow = 1, ncol = 6)
      t_i <- 1
      
      for (v in c("", "sd_", "LB_", "UB_", "t_", "p_")) {
        t_mat[1,t_i] <- get(paste0(v,"diff_delta_1_2_XX"))
        t_i <- t_i + 1
      }
      rownames(t_mat) <- c("Diff.")
      colnames(t_mat) <- c("Estimate", "SE", "LB CI", "UB CI", "t stat.", "pval.")
    }
    
    
    # Returning the results #
    
    ## Message for quasi stayers ##
    if (aoss_XX == 1 & waoss_XX == 1) {
      if (!is.na(scalars$delta_1_1_XX) & !is.na(scalars$delta_2_1_XX)) {
        if (scalars$delta_1_1_XX / scalars$delta_2_1_XX > 10) {
          message("You might have quasi-stayers in your data. The aoss estimand is likely to be biased.")
        }
      }
    }
    
    if (isTRUE(noextrapolation) | isTRUE(exact_match)) {
      if (scalars$N_drop_total_XX > 0) {
        message(sprintf("%.0f switchers are dropped out of the estimation because their baseline treatments do not belong to the support of stayers' baseline treatments.", scalars$N_drop_total_XX))
      }
      if (isTRUE(exact_match) & scalars$N_drop_total_C_XX > 0) {
        message(sprintf("%.0f stayers are dropped out of the estimation because their baseline treatments do not belong to the support of switchers' baseline treatments.", scalars$N_drop_total_C_XX))
      }
    }
    
    IDs_XX <- NULL
    estims <- c("aoss", "waoss", "ivwaoss")
    
    ret_mat_XX <- matrix(NA, nrow = 3*max_T_XX, ncol = 6)
    if (isTRUE(placebo)) {
      ret_mat_pl_XX <- matrix(NA, nrow = 3*max_T_XX, ncol = 6)
    }
    rown <- c()
    for (j in 1:length(estims)) {
      for (p in 1:max_T_XX) {
        if (get(paste0(estims[j],"_XX")) == 1) {
          if (((is.na(scalars[[paste0("N_Stayers_",j,"_",p,"_XX")]]) & is.na(scalars[[paste0("N_Switchers_",j,"_",p,"_XX")]])) | scalars[[paste0("N_Stayers_",j,"_",p,"_XX")]] < 2 | scalars[[paste0("N_Switchers_",j,"_",p,"_XX")]] == 0) & p != 1) {
            scalars[[paste0("delta_",j,"_",p,"_XX")]] <- NA
          }
          ret_mat_XX[(j-1)*max_T_XX + p, 1] <- scalars[[paste0("delta_",j,"_",p,"_XX")]]
          ret_mat_XX[(j-1)*max_T_XX + p, 2] <- scalars[[paste0("sd_delta_",j,"_",p,"_XX")]]
          ret_mat_XX[(j-1)*max_T_XX + p, 3] <- scalars[[paste0("LB_",j,"_",p,"_XX")]]
          ret_mat_XX[(j-1)*max_T_XX + p, 4] <- scalars[[paste0("UB_",j,"_",p,"_XX")]]
          ret_mat_XX[(j-1)*max_T_XX + p, 5] <- scalars[[paste0("N_Switchers_",j,"_",p,"_XX")]]
          ret_mat_XX[(j-1)*max_T_XX + p, 6] <- scalars[[paste0("N_Stayers_",j,"_",p,"_XX")]]   
          
          if (isTRUE(placebo) & p != 2) {
            if (((is.na(scalars[[paste0("N_Stayers_",j,"_",p,"_pl_XX")]]) & is.na(scalars[[paste0("N_Switchers_",j,"_",p,"_pl_XX")]])) | scalars[[paste0("N_Stayers_",j,"_",p,"_pl_XX")]] < 2 | scalars[[paste0("N_Switchers_",j,"_",p,"_pl_XX")]] == 0) & p != 1) {
              scalars[[paste0("delta_",j,"_",p,"_pl_XX")]] <- NA
            }
            ret_mat_pl_XX[(j-1)*max_T_XX + p, 1] <- scalars[[paste0("delta_",j,"_",p,"_pl_XX")]]
            ret_mat_pl_XX[(j-1)*max_T_XX + p, 2] <- scalars[[paste0("sd_delta_",j,"_",p,"_pl_XX")]]
            ret_mat_pl_XX[(j-1)*max_T_XX + p, 3] <- scalars[[paste0("LB_",j,"_",p,"_pl_XX")]]
            ret_mat_pl_XX[(j-1)*max_T_XX + p, 4] <- scalars[[paste0("UB_",j,"_",p,"_pl_XX")]]
            ret_mat_pl_XX[(j-1)*max_T_XX + p, 5] <- scalars[[paste0("N_Switchers_",j,"_",p,"_pl_XX")]]
            ret_mat_pl_XX[(j-1)*max_T_XX + p, 6] <- scalars[[paste0("N_Stayers_",j,"_",p,"_pl_XX")]]   
          }
        }
        
        if (p == 1) {
          rown <- c(rown, toupper(estims[j]))
        } else {
          rown <- c(rown, paste0(estims[j],"_",p))
        }
      }
    }
    rownames(ret_mat_XX) <- rown
    colnames(ret_mat_XX) <- c("Estimate", "SE", "LB CI", "UB CI", "Switchers", "Stayers")
    
    out <- list(table = ret_mat_XX, pairs = max_T_XX)
    if (isTRUE(placebo)) {
      rownames(ret_mat_pl_XX) <- rown
      colnames(ret_mat_pl_XX) <- c("Estimate", "SE", "LB CI", "UB CI", "Switchers", "Stayers")
      out <- c(out, list(ret_mat_pl_XX))
      names(out)[length(out)] <- "table_placebo"
    }
    if (isTRUE(aoss_vs_waoss)) {
      out <- c(out, list(t_mat))
      names(out)[length(out)] <- "aoss_vs_waoss"
    }
    if (!is.null(cluster)) {
      out <- c(out, list(n_clus_XX))
      names(out)[length(out)] <- "n_clusters"
    }
  })
  return(out)
}

#' Internal function to generate quantile subsets based on delta D or delta Z
#' @param df df
#' @param ID ID
#' @param Time Time
#' @param D D
#' @param Z Z
#' @param by_opt by_opt
#' @param quantiles quantiles
#' @importFrom plm pdata.frame make.pbalanced
#' @importFrom stats quantile median
#' @import ggplot2
#' @returns A list with the df object and the relevant quantiles
#' @noRd
did_multiplegt_stat_quantiles <- function(
    df,
    ID,
    Time,
    D,  
    Z,
    by_opt,
    quantiles
) {
  
  df <- subset(df, !(is.na(df[[ID]]) | is.na(df[[Time]]) | is.na(df[[D]]) ))
  if (!is.null(Z)) {
    df <- subset(df, !is.na(df[[Z]]))
  }
  
  ## To be added: case with multiple observations per cell
  
  ## Balance the panel
  df <- pdata.frame(df, index = c(ID, Time)) 
  df <- make.pbalanced(df, balance.type = "fill")
  
  if (is.null(Z)) {
    df$delta_pre_XX <- abs(diff(df[[D]]))
    gr_title <- paste0("\U0394", "D")
  } else {
    df$delta_pre_XX <- abs(diff(df[[Z]]))
    gr_title <- paste0("\U0394", "Z")
  }
  
  df$switchers_dummy_XX <- df$delta_pre_XX != 0
  df <- df %>% group_by(.data[[Time]]) %>%
    mutate(switchers_N_XX = sum(.data$switchers_dummy_XX, na.rm = TRUE)) %>% 
    mutate(stayers_N_XX = sum(1-.data$switchers_dummy_XX, na.rm = TRUE)) %>% 
    ungroup()
  df$in_aggregation_XX <- df$switchers_N_XX > 0 & df$stayers_N_XX > 1
  df$switchers_dummy_XX <- df$switchers_N_XX <- df$stayers_N_XX <- NULL
  
  df_switch <- subset(df, !is.na(df$delta_pre_XX) & df$delta_pre_XX != 0 & df$in_aggregation_XX == 1)
  N_switchers_plot <- nrow(df_switch)
  df_switch$unit_XX <- 1
  df_switch <- df_switch %>% group_by(.data$delta_pre_XX) %>% 
    mutate(tot_delta_XX = sum(.data$unit_XX, na.rm = TRUE))
  df_switch$unit_XX <- NULL; 
  df_switch <- df_switch[c("delta_pre_XX", "tot_delta_XX")]
  df_switch <- df_switch %>% group_by(.data$delta_pre_XX) %>% summarise(tot_delta_XX = mean(.data$tot_delta_XX, na.rm = TRUE))
  df_switch$tot_delta_XX <- df_switch$tot_delta_XX/sum(df_switch$tot_delta_XX, na.rm = TRUE)
  df_switch <- df_switch[order(df_switch$delta_pre_XX), ]
  df_switch$cdf <- cumsum(df_switch$tot_delta_XX)
  df_switch$partition_XX <- by_opt
  cut_off <- c()
  quantiles_temp <- c(0)
  for (j in 2:length(quantiles)) {
    df_switch$partition_XX <- ifelse(df_switch$cdf >= quantiles[j-1] & df_switch$cdf < quantiles[j], j - 1, df_switch$partition_XX)
    if (nrow(subset(df_switch, df_switch$partition_XX == j-1)) > 0) {
      cut_off <- c(cut_off, min(subset(df_switch, df_switch$partition_XX == j-1)$delta_pre_XX,na.rm = TRUE))
      quantiles_temp <- c(quantiles_temp, max(subset(df_switch, df_switch$partition_XX == j-1)$cdf,na.rm = TRUE))
    }
  }
  cut_off <- c(cut_off, max(df_switch$delta_pre_XX, na.rm = TRUE))
  
  quantiles_plot <- ggplot(data = df_switch, aes(x = .data$delta_pre_XX, y = .data$cdf)) + geom_line(size = 0.5) + scale_x_continuous(breaks= cut_off, labels = sprintf("%.2f", cut_off)) + theme(plot.title = element_text(hjust = 0.5), panel.grid.minor = element_blank()) + ylab("CDF") + xlab(gr_title) + ggtitle(sprintf("Empirical distribution of %s", gr_title)) + labs(caption = sprintf("N = %.0f. Quantiles bins cutoffs reported as x axis ticks.", N_switchers_plot))
  quantiles <- quantiles_temp
  df_switch <- quantiles_temp <- NULL
  
  
  df$switchers_XX <- df$delta_pre_XX != 0 & !is.na(df$delta_pre_XX) & df$in_aggregation_XX == 1
  df$partition_XX <- ifelse(df$switchers_XX, 1, 0)
  for (p in 2:length(cut_off)) {
    df$partition_XX <- ifelse(df$switchers_XX & (df$delta_pre_XX > cut_off[p-1] & df$delta_pre_XX <= cut_off[p]), p-1, df$partition_XX)
  }
  df$partition_XX <- as.numeric(df$partition_XX)
  df$partition_XX <- ifelse(df$in_aggregation_XX, df$partition_XX, NA)
  names(cut_off) <- c()
  class(df) <- "data.frame"
  df$it_XX <- 1
  switch_df <- df %>% filter(.data$partition_XX != 0) %>% 
    group_by(.data$partition_XX) %>%
    summarise(N_partition_XX = sum(.data$it_XX, na.rm = TRUE), Med_delta_pre_XX = median(.data$delta_pre_XX, na.rm = TRUE)) %>% ungroup()
  df$it_XX <- switch_df$it_XX <- NULL
  ret <- list(df = df, val_quantiles = cut_off, quantiles = quantiles, switch_df = switch_df, quantiles_plot = quantiles_plot)
  df$delta_pre_XX <- df$switchers_XX <-  NULL     
  return(ret)
}

#' Internal function for estimation of pairwise DiD between consecutive time periods.
#' @param df df
#' @param Y Y
#' @param ID ID
#' @param Time Time
#' @param D D
#' @param Z Z
#' @param estimator estimator
#' @param order order
#' @param noextrapolation noextrapolation
#' @param weight weight
#' @param switchers switchers
#' @param pairwise pairwise
#' @param aoss aoss
#' @param waoss waoss
#' @param ivwaoss ivwaoss
#' @param estimation_method estimation_method
#' @param scalars scalars
#' @param placebo placebo
#' @param exact_match exact_match
#' @param cluster cluster
#' @param by_fd_opt by_fd_opt
#' @param other_treatments other_treatments
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang := 
#' @importFrom rlang .data
#' @importFrom plm pdata.frame make.pbalanced
#' @importFrom stats as.formula lm sd
#' @returns A list with two elements. The first element is a sublist of scalars that are updated throughout the run of the command and across each pair of consecutive periods. The second element is a dataframe with the group-specific variables (among which, the influence function) that will be aggregated in the main program.
#' @noRd
did_multiplegt_stat_pairwise <- function(
    df,
    Y,
    ID,
    Time,
    D,
    Z,
    estimator,
    order,
    noextrapolation,
    weight,
    switchers,
    pairwise,
    IDs,
    aoss,
    waoss,
    ivwaoss,
    estimation_method,
    scalars,
    placebo,
    exact_match,
    cluster,
    by_fd_opt,
    other_treatments
) {
  # Preallocation of scalars
  IV_req_XX <- NULL
  PS_0_XX <- NULL
  
  for (v in names(scalars)) {
    assign(v, scalars[[v]])
  }
  
  ## Placebo ##
  if (isTRUE(placebo)) {
    df <- subset(df, df$T_XX %in% c(pairwise-2,pairwise-1,pairwise))
    pl <- "_pl"
    
  } else {
    df <- subset(df, df$T_XX %in% c(pairwise-1,pairwise))
    pl <- ""
  }
  
  ## Start of the program
  df <- df %>% group_by(.data$T_XX) %>% 
    mutate(tsfilled_min_XX = min(.data$tsfilled_XX, na.rm = TRUE)) 
  gap_XX <- max(df$tsfilled_min_XX, na.rm = TRUE)
  
  df <- df %>% group_by(.data$T_XX) %>% mutate(Tbis_XX = cur_group_id()) 
  df$T_XX <- df$Tbis_XX
  df$Tbis_XX <- NULL
  
  df <- pdata.frame(df, index = c("ID_XX", "T_XX")) 
  df$ID_XX <- as.numeric(as.character(df$ID_XX))
  df$T_XX <- as.numeric(as.character(df$T_XX))
  df <- df[order(df$ID_XX, df$T_XX), ]
  
  df$delta_Y_XX <- diff(df$Y_XX)
  df$delta_D_XX <- diff(df$D_XX)
  if (ivwaoss == 1) {
    df$delta_Z_XX <- diff(df$Z_XX)        
  }
  if (!is.null(other_treatments)) {
    for (v in other_treatments) {
      df[[paste0("fd_",v,"_temp_XX")]] <- diff(df[[v]])
    }
  }
  
  if (!is.null(other_treatments)) {
    for (v in other_treatments) {
      df <- df %>% group_by(.data$ID_XX) %>% 
        mutate(!!paste0("fd_",v,"_XX") := sum(.data[[paste0("fd_",v,"_temp_XX")]], na.rm = TRUE)) %>%ungroup()
      df[[paste0("fd_",v,"_temp_XX")]] <- NULL
    }
  }
  
  if (!is.null(df$partition_XX)) {
    df$partition_lead_XX <- lead(df$partition_XX)
    df$partition_XX <- NULL
  }
  
  if (isTRUE(placebo))  {
    df$delta_temp <- ifelse(df$T_XX == 2, df$delta_Y_XX, NA)
    df <- df %>% group_by(.data$ID_XX) %>% 
      mutate(delta_temp2 = mean(.data$delta_temp, na.rm = TRUE))
    df$delta_Y_XX <- NULL; df$delta_Y_XX <- df$delta_temp2; 
    df$delta_temp <- NULL; df$delta_temp2 <- NULL;
  } else {
    # Generate deltaY_t = Y_t - Y_(t-1) and put it at the same level
    df <- df %>% group_by(.data$ID_XX) %>% 
      mutate(delta_Y_temp_XX = mean(.data$delta_Y_XX, na.rm = TRUE)) %>% ungroup()
    df$delta_Y_XX <- df$delta_Y_temp_XX
    df$delta_Y_temp_XX <- NULL
  }
  
  if (isTRUE(placebo) & (aoss == 1 | waoss == 1)) {
    # Units s.t. D_{t-2} = D_{t-1}
    df$inSamplePlacebo_temp_XX <- (df$delta_D_XX == 0 & df$T_XX == 2)  
    df$inSamplePlacebo_temp_XX <- ifelse(is.na(df$delta_D_XX), NA,df$inSamplePlacebo_temp_XX)
    df <- df %>% group_by(.data$ID_XX) %>%
      mutate(inSamplePlacebo_XX = max(.data$inSamplePlacebo_temp_XX, na.rm = TRUE))
    
    # Only keep Units such that D_{t-2} = D_{t-1}
    #df <- subset(df, df$inSamplePlacebo_XX == 1)
    # we do not need that line since we've already computed y_{t-2} - y_{t-1}, and selected units such that d_{t-2} = d_{t-1} // and evrything that follows is the same as the computation of the effects:)
    df <- subset(df, df$T_XX != 1)
    # We need the DeltaD_t only and we will take the mean after to keep the same value for all the dates
    df$delta_D_XX <- ifelse(df$T_XX != 3, NA,  df$delta_D_XX)
  } 
  
  if (isTRUE(placebo) & ivwaoss == 1) {
    df$inSamplePlacebo_IV_temp_XX <- (df$delta_Z_XX == 0 & df$T_XX == 2)  
    df <- df %>% group_by(.data$ID_XX) %>%
      mutate(inSamplePlacebo_XX = max(.data$inSamplePlacebo_IV_temp_XX, na.rm = TRUE))
    
    # Only keep Units such that D_{t-2} = D_{t-1}
    #df <- subset(df, df$inSamplePlacebo_IV_XX == 1)
    # we do not need that line since we've already computed y_{t-2} - y_{t-1}, and selected units such that z_{t-2} = z_{t-1} // and evrything that follows is the same as the computation of the effects:)
    df <- subset(df, df$T_XX != 1)
    # We need the Delta_Z_t only and we will take the mean after to keep the same value for all the dates
    df$delta_Z_XX <- ifelse(df$T_XX != 3, NA,  df$delta_Z_XX)
  }
  if (nrow(df) == 0) {
    # Since there are no obs, we exit the program
    for (v in names(scalars)) {
      scalars[[v]] <- get(v)
    }
    
    if (aoss == 1) {
      scalars[[paste0("P_",pairwise,pl,"_XX")]] <- 0
    }
    if (waoss == 1) {
      scalars[[paste0("E_abs_delta_D_",pairwise,pl,"_XX")]] <- 0
    }
    
    if (ivwaoss == 1) {
      scalars[[paste0("denom_delta_IV_",pairwise,pl,"_XX")]] <- 0
    }
    
    scalars[[paste0("non_missing_",pairwise,pl,"_XX")]] <- 0
    for (v in c("Switchers", "Stayers")) {
      for (n in 1:3) {
        scalars[[paste0("N_",v,"_",n,"_",pairwise,pl,"_XX")]] <- 0
      }
    }
    
    estims <- c("aoss", "waoss", "ivwaoss")
    indices <- c() 
    for (j in 1:length(estims)) {
      if (get(estims[j]) == 1) {
        indices <- c(indices, j)
      }
    }
    
    for (i in indices) {
      scalars[[paste0("delta_",i,"_",pairwise,pl,"_XX")]] <- 0
      scalars[[paste0("sd_delta_",i,"_",pairwise,pl,"_XX")]] <- NA
      scalars[[paste0("LB_",i,"_",pairwise,pl,"_XX")]] <- NA
      scalars[[paste0("UB_",i,"_",pairwise,pl,"_XX")]] <- NA
    }
    
    out_res <- list(scalars = scalars, to_add = NULL)
    names(out_res) <- c("scalars", "to_add")
    return(out_res)
  }
  
  # Generate deltaD_t = D_t - D_(t-1) and put it at the same level
  df <- df %>% group_by(.data$ID_XX) %>% 
    mutate(delta_D_temp_XX = mean(.data$delta_D_XX, na.rm = TRUE)) %>% ungroup()
  df$delta_D_XX <- df$delta_D_temp_XX
  df$delta_D_temp_XX <- NULL
  
  if (ivwaoss == 1) {
    df <- df %>% group_by(.data$ID_XX) %>% 
      mutate(delta_Z_temp_XX = mean(.data$delta_Z_XX, na.rm = TRUE)) %>% ungroup()
    df$delta_Z_XX <- df$delta_Z_temp_XX
    df$delta_Z_temp_XX <- NULL
    
    df$SI_XX <- (df$delta_Z_XX > 0) - (df$delta_Z_XX < 0) 
    # This is equivalent to sgn(delta_X)
    df$Z1_XX <- NULL; df$Z1_XX <- df$Z_XX; df$Z_XX <- NULL;
  }
  
  df[[paste0("used_in_",pairwise,"_XX")]] <- as.numeric(!is.na(df$delta_Y_XX) & !is.na(df$delta_D_XX))
  if (ivwaoss == 1) {
    df[[paste0("used_in_IV_",pairwise,"_XX")]] <- as.numeric(df[[paste0("used_in_",pairwise,"_XX")]] == 1 & !is.na(df$delta_Z_XX))
    df <- subset(df, df[[paste0("used_in_IV_", pairwise, "_XX")]] == 1)
    
  } else {
    #df <- subset(df, df[[paste0("used_in_", pairwise, "_XX")]] == 1)
  }
  # Generate Switcher : S = 1 if switcher-up, -1 if switcher-down, 0 if stayer
  df$S_XX <- (df$delta_D_XX > 0) - (df$delta_D_XX < 0)
  
  if (waoss == 1 | aoss == 1) {
    df$abs_delta_D_XX <- df$S_XX * df$delta_D_XX
    if (!is.null(switchers)) {
      df <- subset(df, !(df$S_XX == (switchers == "down") - (switchers == "up")))
    }
  }
  if (ivwaoss == 1) {
    if (!is.null(switchers)) {
      df <- subset(df, !(df$SI_XX == (switchers == "down") - (switchers == "up")))
    }
    df$abs_delta_Z_XX <- df$SI_XX * df$delta_Z_XX
  }
  
  # We have all the variable we need at the first year so we can drop the 'second' year line
  df <- subset(df, df$T_XX != max(df$T_XX, na.rm = TRUE))
  df$D1_XX <- df$D_XX; df$D_XX <- NULL;
  
  df$Ht_XX <- as.numeric(!is.na(df$delta_D_XX) & !is.na(df$delta_Y_XX))
  df$S_XX <- ifelse(df$Ht_XX == 0, NA, df$S_XX)
  if (ivwaoss == 1) {
    df$Ht_XX <- as.numeric(df$Ht_XX == 1 & !is.na(df$delta_Z_XX))
    df$SI_XX <- ifelse(df$Ht_XX == 0, NA, df$SI_XX)
  }
  
  if (!is.null(by_fd_opt)) {
    df <- subset(df, df$partition_lead_XX == 0 | df$partition_lead_XX == by_fd_opt)
  }
  
  ## Imbalanced panel adjustments ##
  vars_to_set_missing <- c("S_XX", "delta_D_XX", "delta_Y_XX", "D1_XX")
  if (aoss == 1 | waoss == 1) {
    vars_to_set_missing <- c(vars_to_set_missing, "abs_delta_D_XX")
  } else {
    vars_to_set_missing <- c(vars_to_set_missing, "Z1_XX", "SI_XX")
  }
  
  if (isTRUE(placebo)) {
    for (v in vars_to_set_missing) {
      df[[v]] <- ifelse(df$inSamplePlacebo_XX == 0, NA, df[[v]])
    }
    df$Ht_XX <- ifelse(df$inSamplePlacebo_XX == 0, 0, df$Ht_XX)
  }
  if (!is.null(other_treatments)) {
    for (ot in other_treatments) {
      for (v in vars_to_set_missing) {
        df[[v]] <- ifelse(df[[paste0("fd_",ot,"_XX")]] != 0, NA, df[[v]])
      }
      df$Ht_XX <- ifelse(df[[paste0("fd_",ot,"_XX")]] != 0, 0, df$Ht_XX)
    }
  }
  
  if (isTRUE(noextrapolation)) {
    if (aoss == 1 | waoss == 1) {
      assign(paste0("max_D1",pl,"_XX"), max(df$D1_XX[df$S_XX == 0], na.rm = TRUE))
      assign(paste0("min_D1",pl,"_XX"), min(df$D1_XX[df$S_XX == 0], na.rm = TRUE))
      df$outofBounds_XX <- (df$D1_XX < get(paste0("min_D1",pl,"_XX")) |
                              df$D1_XX > get(paste0("max_D1",pl,"_XX"))) 
      assign(paste0("N_drop_",pairwise,pl,"_XX"), sum(df$outofBounds_XX, na.rm = TRUE))
      df <- subset(df, df$outofBounds_XX != 1)
      
      if (get(paste0("N_drop_",pairwise,pl,"_XX")) > 0 & isFALSE(placebo) & gap_XX == 0 & get(paste0("N_drop_", pairwise,pl,"_XX")) < nrow(df)-1) {
        N_drop_total_XX <- N_drop_total_XX + get(paste0("N_drop_",pairwise,pl,"_XX"))
        
      }           
    }
    if (ivwaoss == 1) {
      assign(paste0("max_Z1",pl,"_XX"), max(df$Z1_XX[df$SI_XX == 0], na.rm = TRUE))
      assign(paste0("min_Z1",pl,"_XX"), min(df$Z1_XX[df$SI_XX == 0], na.rm = TRUE))
      df$outofBoundsIV_XX <- (df$Z1_XX < get(paste0("min_Z1",pl,"_XX")) |
                                df$Z1_XX > get(paste0("max_Z1",pl,"_XX"))) 
      assign(paste0("N_IVdrop_",pairwise,pl,"_XX"), sum(df$outofBoundsIV_XX, na.rm = TRUE))
      df <- subset(df, df$outofBoundsIV_XX != 1)
      
      if (get(paste0("N_IVdrop_",pairwise,pl,"_XX")) > 0 & isFALSE(placebo) & gap_XX == 0 & get(paste0("N_IVdrop_",pairwise,pl,"_XX")) < nrow(df)-1) {
        N_drop_total_XX <- N_drop_total_XX + get(paste0("N_IVdrop_",pairwise,pl,"_XX"))
      }           
    }
  }
  if (isTRUE(exact_match)) {
    if (aoss == 1 | waoss == 1) {
      df$all_treat_XX <- df[c("D1_XX", other_treatments)]
      df <- df %>% group_by(.data$all_treat_XX) %>% 
        mutate(has_match_min_XX = min(.data$abs_delta_D_XX[!is.na(.data$S_XX)], na.rm = TRUE)) %>%
        mutate(has_match_max_XX = max(.data$abs_delta_D_XX[!is.na(.data$S_XX)], na.rm = TRUE)) %>% ungroup()
      df$s_has_match_XX <- ifelse(!is.na(df$S_XX), as.numeric(df$has_match_min_XX == 0), -1)
      df$s_has_match_XX <- ifelse(df$S_XX == 0, -1, df$s_has_match_XX)
      df$c_has_match_XX <- ifelse(!is.na(df$S_XX), as.numeric(df$has_match_max_XX > 0), -1)
      df$c_has_match_XX <- ifelse(df$S_XX != 0 & !is.na(df$S_XX), -1, df$c_has_match_XX)
    } else if (ivwaoss == 1) {
      df$all_treat_XX <- df[c("Z1_XX", other_treatments)]
      df <- df %>% group_by(.data$all_treat_XX) %>% 
        mutate(has_match_min_XX = min(.data$abs_delta_Z_XX[!is.na(.data$SI_XX)], na.rm = TRUE)) %>%
        mutate(has_match_max_XX = max(.data$abs_delta_Z_XX[!is.na(.data$SI_XX)], na.rm = TRUE)) %>% ungroup()
      df$s_has_match_XX <- ifelse(!is.na(df$SI_XX), as.numeric(df$has_match_min_XX == 0, na.rm = TRUE), -1)
      df$s_has_match_XX <- ifelse(df$SI_XX == 0, -1, df$s_has_match_XX)
      df$c_has_match_XX <- ifelse(!is.na(df$SI_XX), as.numeric(df$has_match_max_XX > 0, na.rm = TRUE), -1)
      df$c_has_match_XX <- ifelse(df$SI_XX != 0 & !is.na(df$S_XX), -1, df$c_has_match_XX)
    }
    df$all_treat_XX <- NULL
    df$has_match_min_XX <- df$has_match_max_XX <- NULL
    
    assign(paste0("N_drop_",pairwise,pl,"_XX"), nrow(subset(df, df$s_has_match_XX == 0)))
    assign(paste0("N_drop_",pairwise,pl,"_C_XX"), nrow(subset(df, df$c_has_match_XX == 0)))
    if (get(paste0("N_drop_",pairwise,pl,"_XX")) > 0 & get(paste0("N_drop_",pairwise,pl,"_XX")) != nrow(df) & gap_XX == 0) {
      N_drop_total_XX <- N_drop_total_XX + get(paste0("N_drop_",pairwise,pl,"_XX")) 
    }
    if (get(paste0("N_drop_",pairwise,pl,"_C_XX")) > 0 & get(paste0("N_drop_",pairwise,pl,"_C_XX")) != nrow(df) & gap_XX == 0) {
      N_drop_total_C_XX <- N_drop_total_C_XX + get(paste0("N_drop_",pairwise,pl,"_C_XX")) 
    }
    for (v in vars_to_set_missing) {
      df[[v]] <- ifelse((df$s_has_match_XX == 0 | df$c_has_match_XX == 0), NA, df[[v]])
    }
    df$Ht_XX <- ifelse((df$s_has_match_XX == 0 | df$c_has_match_XX == 0), 0, df$Ht_XX)
    df$c_has_match_XX <- df$s_has_match_XX <- NULL
    
    order <- length(unique(df$D1_XX))
  }
  assign(paste0("W",pl,"_XX"), sum(df$weight_XX, na.rm = TRUE))
  assign(paste0("N",pl,"_XX"), nrow(df))
  
  # Panel with gaps (using tsfilled_XX) and cases where we have only switchers or only stayers (using count)
  if (waoss == 1 | aoss == 1) {
    assign(paste0("N_Switchers",pl,"_XX"), nrow(subset(df, df$S_XX != 0 & !is.na(df$S_XX))))
    assign(paste0("N_Stayers",pl,"_XX"),nrow(subset(df, df$S_XX == 0)))
  }
  if (ivwaoss == 1) {
    assign(paste0("N_Switchers_IV",pl,"_XX"), nrow(subset(df, df$SI_XX != 0 & !is.na(df$SI_XX))))
    assign(paste0("N_Stayers_IV",pl,"_XX"),nrow(subset(df, df$SI_XX == 0)))
  }
  
  factor_temp <- FALSE
  vars_pol_XX <- c()
  for (pol_level in 1:order) {
    df[[paste0("D1_",pol_level,"_XX")]] <- df$D1_XX^pol_level
    vars_pol_XX <- c(vars_pol_XX, paste0("D1_",pol_level,"_XX"))
  }
  reg_pol_XX <- ""
  for (level in 1:length(vars_pol_XX)) {
    if (level > 1) {
      reg_pol_XX <- paste0(reg_pol_XX," + ")
    }
    reg_pol_XX <- paste0(reg_pol_XX,paste0("D1_",level,"_XX"))
  } 
  if (!is.null(other_treatments)) {
    if (isTRUE(factor_temp)) {
      df$D1_1_XXFACT <- as.factor(df$D1_1_XX) 
      for (v in other_treatments) {
        df[[paste0(v,"FACT")]] <- as.factor(df[[v]])
        other_treatments[other_treatments == v] <- paste0(v,"FACT")
      }
      interact <- "D1_1_XXFACT"
      for (v in other_treatments) {
        interact <- paste0(interact," * ",v)
      }
      reg_pol_XX <- paste0(reg_pol_XX," + ", interact)
      saturated_reg <- power_set(c("D1_1_XXFACT", other_treatments))
      vars_pol_XX <- union(vars_pol_XX, saturated_reg)
    } else {
      vars_pol_XX <- power_set(c("D1_1_XX", other_treatments))
      interact <- "D1_1_XX"
      for (v in other_treatments) {
        interact <- paste0(interact," * ",v)                
      }
      reg_pol_XX <- paste0(reg_pol_XX," + ", interact)
    }
  }
  
  if (ivwaoss == 1) {
    IV_vars_pol_XX <- c()
    for (pol_level in 1:order) {
      df[[paste0("Z1_",pol_level,"_XX")]] <- df$Z1_XX^pol_level
      IV_vars_pol_XX <- c(IV_vars_pol_XX, paste0("Z1_",pol_level,"_XX"))
    }
    IV_reg_pol_XX <- ""
    for (level in 1:length(IV_vars_pol_XX)) {
      if (level > 1) {
        IV_reg_pol_XX <- paste0(IV_reg_pol_XX," + ")
      }
      IV_reg_pol_XX <- paste0(IV_reg_pol_XX,paste0("Z1_",level,"_XX"))
    } 
    if (!is.null(other_treatments)) {
      interact <- "Z1_1_XX"
      for (v in other_treatments) {
        interact <- paste0(interact," * ",v)
      }
      IV_reg_pol_XX <- gsub("Z1_1_XX", interact, IV_reg_pol_XX)
    }
  }
  
  df$S_bis_XX <- ifelse(!is.na(df$S_XX), as.numeric(df$S_XX!= 0), NA)
  
  # Feasibility conditions:
  feasible_est <- FALSE
  if (aoss == 1 | waoss == 1) {
    feasible_est <- (gap_XX == 0 & get(paste0("N_Switchers",pl,"_XX")) > 0 & get(paste0("N_Stayers",pl,"_XX")) > 1)
  } else if (ivwaoss == 1) {
    feasible_est <- (gap_XX == 0 & get(paste0("N_Switchers_IV",pl,"_XX")) > 0 & get(paste0("N_Stayers_IV",pl,"_XX")) > 1)
  }
  
  #fact_reg <- !is.null(other_treatments)
  fact_reg <- FALSE
  assign(paste0("P_Ht_",pairwise,pl,"_XX"), Mean("Ht_XX", df))
  
  if (!is.null(cluster)) {
    df <- df %>% group_by(.data$ID_XX) %>%
      mutate(gr_id = row_number()) %>% ungroup()
    df$id_temp <- as.numeric(df$gr_id == 1)
    df <- df %>% group_by(.data$cluster_XX) %>%
      mutate(N_c_XX = sum(.data$id_temp, na.rm = TRUE)) %>% ungroup()
    assign(paste0("N_bar_c_",pairwise,pl,"_XX"), mean(df$N_c_XX, na.rm = TRUE))
    df$id_temp <- df$N_c_XX <- df$gr_id <- NULL
  }
  # Start of feasible estimation
  if (feasible_est) {       
    
    if (waoss == 1 | aoss == 1) {
      
      df0 <- subset(df, df$S_XX == 0)
      # \hat{E}(deltaY|D1, S=0)
      model <- lm(as.formula(paste("delta_Y_XX", adj_fact(reg_pol_XX, df0, fact_reg), sep = "~")), data = df0, weights = df0$weight_XX, na.rm = TRUE)
      df <- lpredict(df, "mean_pred_XX", model, vars_pol_XX, factor = fact_reg)
      
      df$inner_sum_delta_1_2_XX <- df$delta_Y_XX -  df$mean_pred_XX
      
      df$S0_XX <- 1 - df$S_bis_XX
      
      if (isFALSE(exact_match)) {
        ## 1. Estimation of P(S = 0|D_1) 
        model <- stata_logit(as.formula(paste("S0_XX",adj_fact(reg_pol_XX,df,fact_reg),sep="~")), df)
        df <- lpredict(df,"PS_0_D_1_XX", model, vars_pol_XX, prob = TRUE, factor = fact_reg)
      } else {
        ## Estimation of 1-E(S|D1) = P(S = 0|D_1) in exact_match case
        model <- lm(as.formula(paste("S_bis_XX",adj_fact(reg_pol_XX,df,fact_reg),sep="~")),data=df,weights= df$weight_XX)
        df <- lpredict(df, "ES_bis_XX_D_1", model, vars_pol_XX, factor = fact_reg)
        
        ## Estimation of \hat{E}(S+-S-|D1) for both \Phi_2:
        model <- lm(as.formula(paste("S_XX",adj_fact(reg_pol_XX,df,fact_reg),sep="~")),data=df,weights= df$weight_XX)
        df <- lpredict(df, "ES_XX_D_1", model, vars_pol_XX, factor = fact_reg)
      }
      assign(paste0("PS_0",pl,"_XX"), Mean("S0_XX", df))
    }
    
    ####################################################### AOSS ##########
    
    if (aoss == 1) {
      # 0) Compute P_t = P(S_t = 1) = E(S_t) for the aggregation afterward
      assign(paste0("P_",pairwise,pl,"_XX"), Mean("S_bis_XX", df) * get(paste0("P_Ht_",pairwise,pl,"_XX")))
      assign(paste0("PS_sum",pl,"_XX"), get(paste0("PS_sum",pl,"_XX")) + get(paste0("P_",pairwise,pl,"_XX")))
      assign(paste0("ES",pl,"_XX"), Mean("S_bis_XX", df))
      
      # 1) Compute \hat{delta}_1
      df$inner_sum_delta_1_XX <- df$inner_sum_delta_1_2_XX / df$delta_D_XX
      df$inner_sum_delta_1_XX <- ifelse(df$delta_D_XX == 0, NA, df$inner_sum_delta_1_XX)
      assign(paste0("delta_1_",pairwise,pl,"_XX"), Mean("inner_sum_delta_1_XX",df))
      
      # 2) Compute the variance of \hat{delta}_1
      df$S_over_delta_D_XX <- df$S_bis_XX / df$delta_D_XX
      df$S_over_delta_D_XX <- ifelse(df$S_bis_XX == 0, 0, df$S_over_delta_D_XX)
      
      modelvar <- lm(as.formula(paste("S_over_delta_D_XX",adj_fact(reg_pol_XX,df,fact_reg), sep = "~")), data = df, weights = df$weight_XX)
      df <- lpredict(df, "mean_S_over_delta_D_XX", modelvar, vars_pol_XX, factor = fact_reg)
      
      # i. estimation of \hat{E}(S/deltaD|D1)
      if (isFALSE(exact_match)) {
        df[[paste0("Phi_1_",pairwise,pl,"_XX")]] <- (df$S_over_delta_D_XX - df$mean_S_over_delta_D_XX * ((1 - df$S_bis_XX)/(df$PS_0_D_1_XX))) * df$inner_sum_delta_1_2_XX
      } else {
        df[[paste0("Phi_1_",pairwise,pl,"_XX")]] <- (df$S_over_delta_D_XX - df$mean_S_over_delta_D_XX * ((1 - df$S_bis_XX)/(1-df$ES_bis_XX_D_1))) * df$inner_sum_delta_1_2_XX
      }
      
      df[[paste0("Phi_1_",pairwise,pl,"_XX")]] <- (df[[paste0("Phi_1_",pairwise,pl,"_XX")]] - (get(paste0("delta_1_", pairwise,pl,"_XX")) * df$S_bis_XX)) / (get(paste0("ES",pl,"_XX")) * get(paste0("P_Ht_",pairwise,pl,"_XX")))
      df[[paste0("Phi_1_",pairwise,pl,"_XX")]] <- ifelse(df$Ht_XX == 0, 0, df[[paste0("Phi_1_",pairwise,pl,"_XX")]])
      
      assign(paste0("mean_IF_1_",pairwise,pl),
             Mean(paste0("Phi_1_",pairwise,pl,"_XX"), df))
      
      if (!is.null(cluster)) {
        df <- df %>% group_by(.data$cluster_XX) %>% 
          mutate(!!paste0("Phi_1_",pairwise,pl,"_c_XX") := sum(.data[[paste0("Phi_1_", pairwise,pl,"_XX")]], na.rm = TRUE)) %>%
          mutate(first_obs_by_clus = row_number() == 1) %>% ungroup()
        df[[paste0("Phi_1_",pairwise,pl,"_c_XX")]] <- ifelse(df$first_obs_by_clus == 1, df[[paste0("Phi_1_",pairwise,pl,"_c_XX")]], NA) / get(paste0("N_bar_c_",pairwise,pl,"_XX"))
        nobs_c_XX <- wSum(subset(df, !is.na(df[[paste0("Phi_1_",pairwise,pl,"_c_XX")]])), w = "weight_c_XX")
        assign(paste0("sd_delta_1_",pairwise,pl,"_XX"),
               Sd(paste0("Phi_1_",pairwise,pl,"_c_XX"), df, w = "weight_c_XX")/sqrt(nobs_c_XX))
        nobs_c_XX <- df$first_obs_by_clus <- NULL
      } else {
        assign(paste0("sd_delta_1_",pairwise,pl,"_XX"),
               sd(df[[paste0("Phi_1_",pairwise,pl,"_XX")]], na.rm = TRUE)/sqrt(wSum(df)))
      }
      assign(paste0("LB_1_",pairwise,pl,"_XX"),
             get(paste0("delta_1_",pairwise,pl,"_XX")) - 1.96 * get(paste0("sd_delta_1_", pairwise,pl,"_XX")))
      assign(paste0("UB_1_",pairwise,pl,"_XX"),
             get(paste0("delta_1_",pairwise,pl,"_XX")) + 1.96 * get(paste0("sd_delta_1_", pairwise,pl,"_XX")))
      
      df[[paste0("S_",pairwise,pl,"_XX")]] <- df$S_bis_XX
      df[[paste0("S_",pairwise,pl,"_XX")]] <- ifelse(df$Ht_XX == 0, df$Ht_XX, df[[paste0("S_",pairwise,pl,"_XX")]] )
    }
    
    ###################################################### WAOSS ##########
    
    if (waoss == 1) {
      assign(paste0("E_abs_delta_D",pl,"_XX"), Mean("abs_delta_D_XX", df))
      assign(paste0("E_abs_delta_D_",pairwise,pl,"_XX"), get(paste0("E_abs_delta_D",pl,"_XX")) * get(paste0("P_Ht_",pairwise,pl,"_XX")))
      assign(paste0("E_abs_delta_D_sum",pl,"_XX"), get(paste0("E_abs_delta_D_sum",pl,"_XX")) + get(paste0("E_abs_delta_D_",pairwise,pl,"_XX")))
      
      for (suffix in c("Minus", "Plus")) {
        df$Ster_XX <- NULL
        df$Ster_XX <- df$S_XX == as.numeric((suffix == "Plus") - (suffix == "Minus"))
        
        ## Computing the contribution weights ##
        df$prod_sgn_delta_D_delta_D_XX <- df$S_XX * df$delta_D_XX
        sum_prod_sgn_delta_D_delta_D_XX <- Sum("prod_sgn_delta_D_delta_D_XX", subset(df,df$Ster_XX == 1))
        assign(paste0("w_",suffix,"_",pairwise,pl,"_XX"), sum_prod_sgn_delta_D_delta_D_XX/get(paste0("N",pl,"_XX")))
        assign(paste0("denom_delta_2_",suffix,"_",pairwise,pl,"_XX"), Sum("delta_D_XX", subset(df, df$Ster_XX == 1)))
        
        
        if (estimation_method == "ra") {
          if (get(paste0("denom_delta_2_",suffix,"_",pairwise,pl,"_XX")) == 0) {
            assign(paste0("denom_delta_2_",suffix,"_",pairwise,pl,"_XX"), 1)
            # in case it is zero set it to 1 to avoid dividing by 0, in that case the numerator is also equal to 0	
          }
          assign(paste0("num_delta_2_",suffix,"_",pairwise,pl,"_XX"),
                 Sum("inner_sum_delta_1_2_XX", subset(df, df$Ster_XX == 1)))
          
          assign(paste0("delta_2_",suffix,"_",pairwise,pl,"_XX"), 
                 get(paste0("num_delta_2_",suffix,"_",pairwise,pl,"_XX")) /
                   get(paste0("denom_delta_2_",suffix,"_",pairwise,pl,"_XX")))
        } 
        
        assign(paste0("nb_Switchers_",suffix,pl,"_XX"), nrow(subset(df, df$Ster_XX ==1)))
        
        assign(paste0("PS_",suffix,"1",pl,"_XX"), get(paste0("nb_Switchers_",suffix,pl,"_XX"))/get(paste0("N",pl,"_XX")))
        
        if (isFALSE(exact_match)) {
          if (get(paste0("PS_",suffix,"1",pl,"_XX")) == 0) {
            # The regression is performed iff there is at least one switcher up/down.
            assign(paste0("delta_2_",suffix,"_",pairwise,pl,"_XX"), 0)
            df[[paste0("PS_1_",suffix,"_D_1_XX")]] <- 0
          } else {
            model <- stata_logit(as.formula(paste("Ster_XX",reg_pol_XX,sep="~")), df)
            df <- lpredict(df,paste0("PS_1_",suffix,"_D_1_XX"),model, vars_pol_XX, prob = TRUE, factor = fact_reg)
            
            
            if (estimation_method == "ps") {
              df[[paste0("delta_Y_P_",suffix,"_XX")]] <- df$delta_Y_XX * (df[[paste0("PS_1_",suffix,"_D_1_XX")]]/df$PS_0_D_1_XX) * (get(paste0("PS_0",pl,"_XX")) / get(paste0("PS_",suffix,"1",pl,"_XX")))
              
              assign(paste0("mean_delta_Y_P_",suffix,pl,"_XX"), 
                     Mean(paste0("delta_Y_P_",suffix,"_XX"),subset(df,df$S_XX == 0)))
              assign(paste0("mean_delta_Y",pl,"_XX"), Mean("delta_Y_XX", subset(df, df$Ster_XX == 1)))
              assign(paste0("mean_delta_D",pl,"_XX"), Mean("delta_D_XX", subset(df, df$Ster_XX == 1)))
              assign(paste0("delta_2_",suffix,"_",pairwise,pl,"_XX"), 
                     (get(paste0("mean_delta_Y",pl,"_XX")) - get(paste0("mean_delta_Y_P_",suffix,pl,"_XX"))) / get(paste0("mean_delta_D",pl,"_XX")))
            }
          }                
        }
      }
      
      if (estimation_method == "ra" | estimation_method == "ps") {
        
        ## Computing the final weights ##
        assign(paste0("W_Plus_",pairwise,pl,"_XX"), get(paste0("w_Plus_",pairwise,pl,"_XX")) /
                 (get(paste0("w_Plus_",pairwise,pl,"_XX")) + get(paste0("w_Minus_",pairwise,pl,"_XX"))))
      }
      
      if (isFALSE(exact_match)) {
        df$dr_delta_Y_XX <- (df$S_XX - ((df$PS_1_Plus_D_1_XX - df$PS_1_Minus_D_1_XX)/df$PS_0_D_1_XX) * (1 - df$S_bis_XX)) * df$inner_sum_delta_1_2_XX
        assign(paste0("denom_dr_delta_2",pl,"_XX"), Sum("dr_delta_Y_XX", df))
      }
      
      if (estimation_method == "ra" | estimation_method == "ps") {
        assign(paste0("delta_2_",pairwise,pl,"_XX"), 
               get(paste0("W_Plus_",pairwise,pl,"_XX"))*get(paste0("delta_2_Plus_",pairwise,pl,"_XX")) + (1-get(paste0("W_Plus_",pairwise,pl,"_XX")))*get(paste0("delta_2_Minus_",pairwise,pl,"_XX")))
      } else if (estimation_method == "dr") {
        sum_abs_delta_D_XX <- Sum("abs_delta_D_XX",df)
        assign(paste0("delta_2_",pairwise,pl,"_XX"), get(paste0("denom_dr_delta_2",pl,"_XX")) / sum_abs_delta_D_XX)
      }
      
      # Computing the variance
      if (isFALSE(exact_match)) {
        df[[paste0("Phi_2_",pairwise,pl,"_XX")]] <- (df$dr_delta_Y_XX - get(paste0("delta_2_", pairwise,pl,"_XX")) * df$abs_delta_D_XX)
      } else {
        df[[paste0("Phi_2_",pairwise,pl,"_XX")]] <- ((df$S_XX - df$ES_XX_D_1 * ((1 - df$S_bis_XX)/(1 - df$ES_bis_XX_D_1))) * df$inner_sum_delta_1_2_XX - get(paste0("delta_2_", pairwise,pl,"_XX")) * df$abs_delta_D_XX)
      }
      df[[paste0("Phi_2_",pairwise,pl,"_XX")]] <- df[[paste0("Phi_2_",pairwise,pl,"_XX")]]/(get(paste0("P_Ht_",pairwise,pl,"_XX"))*get(paste0("E_abs_delta_D",pl,"_XX")))
      df[[paste0("Phi_2_",pairwise,pl,"_XX")]]  <- ifelse(df$Ht_XX == 0, df$Ht_XX, df[[paste0("Phi_2_",pairwise,pl,"_XX")]])
      
      assign(paste0("mean_IF_2_",pairwise,pl),
             Mean(paste0("Phi_2_",pairwise,pl,"_XX"), df))
      
      if (!is.null(cluster)) {
        df <- df %>% group_by(.data$cluster_XX) %>% 
          mutate(!!paste0("Phi_2_",pairwise,pl,"_c_XX") := sum(.data[[paste0("Phi_2_", pairwise,pl,"_XX")]], na.rm = TRUE)) %>%
          mutate(first_obs_by_clus = row_number() == 1) %>% ungroup()
        df[[paste0("Phi_2_",pairwise,pl,"_c_XX")]] <- ifelse(df$first_obs_by_clus == 1, df[[paste0("Phi_2_",pairwise,pl,"_c_XX")]], NA) / get(paste0("N_bar_c_",pairwise, pl,"_XX"))
        nobs_c_XX <- wSum(subset(df, !is.na(df[[paste0("Phi_2_",pairwise,pl,"_c_XX")]])), w = "weight_c_XX")
        assign(paste0("sd_delta_2_",pairwise,pl,"_XX"),
               Sd(paste0("Phi_2_",pairwise,pl,"_c_XX"), df, w = "weight_c_XX")/sqrt(nobs_c_XX))
        nobs_c_XX <- df$first_obs_by_clus <- NULL
      } else {
        assign(paste0("sd_delta_2_",pairwise,pl,"_XX"),
               Sd(paste0("Phi_2_",pairwise,pl,"_XX"), df)/sqrt(wSum(df)))
      }
      assign(paste0("LB_2_",pairwise,pl,"_XX"), get(paste0("delta_2_",pairwise,pl,"_XX")) - 1.96 * get(paste0("sd_delta_2_", pairwise,pl,"_XX")))
      assign(paste0("UB_2_",pairwise,pl,"_XX"), get(paste0("delta_2_",pairwise,pl,"_XX")) + 1.96 * get(paste0("sd_delta_2_", pairwise,pl,"_XX")))
      
      df[[paste0("abs_delta_D_",pairwise,pl,"_XX")]] <- ifelse(df$Ht_XX == 0, 0, df$abs_delta_D_XX)
    }
    
    ##################################################### ivwaoss ##########
    if (ivwaoss == 1) {
      assign(paste0("E_abs_delta_Z",pl,"_XX"), Mean("abs_delta_Z_XX", df))
      
      df$SI_bis_XX <- (df$SI_XX != 0 & !is.na(df$SI_XX))
      df$SI_Plus_XX <- (df$SI_XX == 1)
      df$SI_Minus_XX <- (df$SI_XX == -1)
      
      # Preliminaries: logit regression
      df$S_IV_0_XX <- 1 - df$SI_bis_XX
      if (isFALSE(exact_match)) {
        model <- stata_logit(as.formula(paste("S_IV_0_XX",IV_reg_pol_XX,sep="~")), df)
        df <- lpredict(df,"PS_IV_0_Z_1_XX",model, IV_vars_pol_XX, prob = TRUE, factor = fact_reg)
      } else {
        model <- lm(as.formula(paste("SI_bis_XX",IV_reg_pol_XX,sep="~")), data = df, weights = df$weight_XX)
        df <- lpredict(df, "ES_I_bis_XX_Z_1", model, IV_vars_pol_XX, factor = fact_reg)
        model <- lm(as.formula(paste("SI_XX",IV_reg_pol_XX,sep="~")), data = df, weights = df$weight_XX)
        df <- lpredict(df, "ES_I_XX_Z_1", model, IV_vars_pol_XX, factor = fact_reg)
      }
      assign(paste0("PS_IV_0",pl,"_XX"), Mean("S_IV_0_XX", df))
      
      for (suffix in c("Minus", "Plus")) {
        assign(paste0("nb_Switchers_I_",suffix,pl,"_XX"), 
               nrow(subset(df, df[[paste0("SI_",suffix,"_XX")]] == 1)))
        assign(paste0("PS_I_",suffix,"_1",pl,"_XX"), get(paste0("nb_Switchers_I_",suffix,pl,"_XX")) / get(paste0("N",pl,"_XX")))
        
        if (get(paste0("PS_I_",suffix,"_1",pl,"_XX")) == 0) {
          df[[paste0("PS_I_",suffix,"_1_Z_1_XX")]] <- 0
        } else {
          if (isFALSE(exact_match)) {
            model <- stata_logit(as.formula(paste(paste0("SI_",suffix,"_XX"),IV_reg_pol_XX,sep="~")), df)
            df <- lpredict(df,paste0("PS_I_",suffix,"_1_Z_1_XX"),model, IV_vars_pol_XX, prob = TRUE, factor = fact_reg)
          }
        }
      }
      
      # i. Estimation of  \hat{E}(deltaY|Z1, SI=0)
      df_temp <- subset(df, df$SI_XX == 0)
      model <- lm(as.formula(paste("delta_Y_XX",IV_reg_pol_XX, sep = "~")), 
                  data = df_temp, weights = df_temp$weight_XX)
      df <- lpredict(df,"mean_delta_Y_pred_IV_XX", model, IV_vars_pol_XX, factor = fact_reg)   
      df$inner_sum_IV_num_XX <- df$delta_Y_XX - df$mean_delta_Y_pred_IV_XX
      
      # i. Estimation of  \hat{E}(deltaD|Z1, SI=0)
      model <- lm(as.formula(paste("delta_D_XX",IV_reg_pol_XX, sep = "~")), 
                  data = df_temp, weights = df_temp$weight_XX)
      df <- lpredict(df,"mean_delta_D_pred_IV_XX", model, IV_vars_pol_XX, factor = fact_reg)   
      df$inner_sum_IV_denom_XX <- df$delta_D_XX - df$mean_delta_D_pred_IV_XX
      df_temp <- NULL
      
      if (estimation_method == "ra") {
        for (v in c("num","denom")) {
          df[[paste0("inner_sum_IV_",v,"_XX")]] <- df[[paste0("inner_sum_IV_",v,"_XX")]] * df$SI_XX
          assign(paste0(v,"_delta_IV_",pairwise,pl,"_XX"), 
                 Mean(paste0("inner_sum_IV_",v,"_XX"), df))
        }
      }
      
      if (estimation_method == "ps") {
        # Numerator
        df$delta_Y_P_IV_XX <- df$delta_Y_XX * ((df$PS_I_Plus_1_Z_1_XX - df$PS_I_Minus_1_Z_1_XX)/df$PS_IV_0_Z_1_XX) * get(paste0("PS_IV_0",pl,"_XX"))
        assign(paste0("mean_delta_Y_P_IV",pl,"_XX"), Mean("delta_Y_P_IV_XX", subset(df, df$SI_bis_XX == 0)))
        df$prod_sgn_delta_Z_delta_Y_XX <- df$SI_XX * df$delta_Y_XX
        assign(paste0("mean_sgn_delta_Z_delta_Y",pl,"_XX"), Mean("prod_sgn_delta_Z_delta_Y_XX", df))
        assign(paste0("num_delta_IV_",pairwise,pl,"_XX"), 
               get(paste0("mean_sgn_delta_Z_delta_Y",pl,"_XX")) - get(paste0("mean_delta_Y_P_IV",pl,"_XX")))
        
        # Denominator
        df$delta_D_P_IV_XX <- df$delta_D_XX * ((df$PS_I_Plus_1_Z_1_XX - df$PS_I_Minus_1_Z_1_XX)/df$PS_IV_0_Z_1_XX) * get(paste0("PS_IV_0",pl,"_XX"))
        assign(paste0("mean_delta_D_P_IV",pl,"_XX"), Mean("delta_D_P_IV_XX", subset(df, df$SI_bis_XX == 0)))
        df$prod_sgn_delta_Z_delta_D_XX <- df$SI_XX * df$delta_D_XX
        assign(paste0("mean_sgn_delta_Z_delta_D",pl,"_XX"), Mean("prod_sgn_delta_Z_delta_D_XX", df))
        assign(paste0("denom_delta_IV_",pairwise,pl,"_XX"), 
               get(paste0("mean_sgn_delta_Z_delta_D",pl,"_XX")) - get(paste0("mean_delta_D_P_IV",pl,"_XX")))
      }
      
      if (estimation_method == "dr") {
        df$dr_IV_delta_Y_XX <- (df$SI_XX - ((df$PS_I_Plus_1_Z_1_XX - df$PS_I_Minus_1_Z_1_XX) / df$PS_IV_0_Z_1_XX) * (1 - df$SI_bis_XX)) * df$inner_sum_IV_num_XX
        assign(paste0("num_delta_IV_",pairwise,pl,"_XX"), Mean("dr_IV_delta_Y_XX", df))
        
        df$dr_IV_delta_D_XX <- (df$SI_XX - ((df$PS_I_Plus_1_Z_1_XX - df$PS_I_Minus_1_Z_1_XX) / df$PS_IV_0_Z_1_XX) * (1 - df$SI_bis_XX)) * df$inner_sum_IV_denom_XX
        assign(paste0("denom_delta_IV_",pairwise,pl,"_XX"), Mean("dr_IV_delta_D_XX", df))
      }
      
      assign(paste0("delta_3_",pairwise,pl,"_XX"), 
             get(paste0("num_delta_IV_",pairwise,pl,"_XX"))/get(paste0("denom_delta_IV_",pairwise,pl,"_XX")))
      
      assign(paste0("denom_delta_IV_sum",pl,"_XX"), get(paste0("denom_delta_IV_sum",pl,"_XX")) + get(paste0("denom_delta_IV_",pairwise,pl,"_XX")))
      
      assign(paste0("delta_Y",pl,"_XX"), Mean("inner_sum_IV_num_XX", df)) 
      df_temp <- subset(df, df$SI_XX == 0)
      model <- lm(as.formula(paste("delta_Y_XX", reg_pol_XX, sep = "~")), 
                  data = df_temp, weights = df_temp$weight_XX)
      df <- lpredict(df,"mean_pred_Y_IV_XX", model, vars_pol_XX, factor = fact_reg)   
      if (isFALSE(exact_match)) {
        df$Phi_Y_XX <- ((df$SI_XX - (df$PS_I_Plus_1_Z_1_XX - df$PS_I_Minus_1_Z_1_XX) * (1 - df$SI_bis_XX) / df$PS_IV_0_Z_1_XX) * (df$delta_Y_XX - df$mean_pred_Y_IV_XX) - get(paste0("delta_Y",pl,"_XX")) * df$abs_delta_Z_XX) / get(paste0("E_abs_delta_Z",pl,"_XX"))
      } else {
        df$Phi_Y_XX <- ((df$SI_XX - (df$ES_I_XX_Z_1) * ((1 - df$SI_bis_XX) / (1- df$ES_I_bis_XX_Z_1))) * (df$delta_Y_XX - df$mean_pred_Y_IV_XX) - get(paste0("delta_Y",pl,"_XX")) * df$abs_delta_Z_XX) / get(paste0("E_abs_delta_Z",pl,"_XX"))
      }
      
      assign(paste0("delta_D",pl,"_XX"), Mean("inner_sum_IV_denom_XX", df)) 
      df_temp <- subset(df, df$SI_XX == 0)
      model <- lm(as.formula(paste("delta_D_XX", reg_pol_XX, sep = "~")), 
                  data = df_temp, weights = df_temp$weight_XX)
      df <- lpredict(df,"mean_pred_D_IV_XX", model, vars_pol_XX, factor = fact_reg)   
      
      if (isFALSE(exact_match)) {
        df$Phi_D_XX <- ((df$SI_XX - (df$PS_I_Plus_1_Z_1_XX - df$PS_I_Minus_1_Z_1_XX) * (1 - df$SI_bis_XX) / df$PS_IV_0_Z_1_XX) * (df$delta_D_XX - df$mean_pred_D_IV_XX) - get(paste0("delta_D",pl,"_XX")) * df$abs_delta_Z_XX) / get(paste0("E_abs_delta_Z",pl,"_XX"))
      } else {
        df$Phi_D_XX <- ((df$SI_XX - (df$ES_I_XX_Z_1) * ((1 - df$SI_bis_XX) / (1- df$ES_I_bis_XX_Z_1))) * (df$delta_D_XX - df$mean_pred_D_IV_XX) - get(paste0("delta_D",pl,"_XX")) * df$abs_delta_Z_XX) / get(paste0("E_abs_delta_Z",pl,"_XX"))
      }
      
      if (get(paste0("delta_D",pl,"_XX"))!= 0) {
        df[[paste0("Phi_3_",pairwise,pl,"_XX")]] <- (df$Phi_Y_XX - get(paste0("delta_3_",pairwise,pl,"_XX")) * df$Phi_D_XX) / get(paste0("delta_D",pl,"_XX"))
        mean_IF3 <- Mean(paste0("Phi_3_",pairwise,pl,"_XX"), df)
        
        if (!is.null(cluster)) {
          df <- df %>% group_by(.data$cluster_XX) %>% 
            mutate(!!paste0("Phi_3_",pairwise,pl,"_c_XX") := sum(.data[[paste0("Phi_3_", pairwise,pl,"_XX")]], na.rm = TRUE)) %>%
            mutate(first_obs_by_clus = row_number() == 1) %>% ungroup()
          df[[paste0("Phi_3_",pairwise,pl,"_c_XX")]] <- ifelse(df$first_obs_by_clus == 1, df[[paste0("Phi_3_",pairwise,pl,"_c_XX")]], NA) / get(paste0("N_bar_c_",pairwise,pl,"_XX"))
          nobs_c_XX <- wSum(subset(df, !is.na(df[[paste0("Phi_3_",pairwise,pl,"_c_XX")]])), w = "weight_c_XX")
          assign(paste0("sd_delta_3_",pairwise,pl,"_XX"),
                 Sd(paste0("Phi_3_",pairwise,pl,"_c_XX"), df, w = "weight_c_XX")/sqrt(nobs_c_XX))
          nobs_c_XX <- NULL
          df$first_obs_by_clus <- NULL
        } else {
          assign(paste0("sd_delta_3_",pairwise,pl,"_XX"),
                 Sd(paste0("Phi_3_",pairwise,pl,"_XX"), df)/sqrt(wSum(df)))
        }
        
        assign(paste0("LB_3_",pairwise,pl,"_XX"), get(paste0("delta_3_",pairwise,pl,"_XX")) - 1.96*get(paste0("sd_delta_3_",pairwise,pl,"_XX")))
        assign(paste0("UB_3_",pairwise,pl,"_XX"), get(paste0("delta_3_",pairwise,pl,"_XX")) + 1.96*get(paste0("sd_delta_3_",pairwise,pl,"_XX")))
      } else {
        df[[paste0("Phi_3_",pairwise,pl,"_XX")]] <- NA
        assign(paste0("sd_delta_3_",pairwise,pl,"_XX"), NA)
        assign(paste0("LB_3_",pairwise,pl,"_XX"), NA)
        assign(paste0("UB_3_",pairwise,pl,"_XX"), NA)
      }
      
      df[[paste0("inner_sum_IV_denom_",pairwise,pl,"_XX")]] <- df$inner_sum_IV_denom_XX
    }
    
    assign(paste0("non_missing_",pairwise,pl,"_XX"), 1)
  } else {
    for (i in 1:3) {
      assign(paste0("delta_",i,"_",pairwise,pl,"_XX"), 0)
      assign(paste0("sd_delta_",i,"_",pairwise,pl,"_XX"), NA)
      assign(paste0("LB_",i,"_",pairwise,pl,"_XX"), NA)
      assign(paste0("UB_",i,"_",pairwise,pl,"_XX"), NA)
      df[[paste0("Phi",i,"_",pairwise,pl,"_XX")]] <- NA
    }
    
    if (aoss == 1 | waoss == 1) {
      IVt <- ""
    } else if (ivwaoss == 1) {
      IVt <- "_IV"
    }
    
    if (gap_XX != 0) {
      assign(paste0("N_Switchers",IVt,pl,"_XX"), NA)
      assign(paste0("N_Stayers",IVt,pl,"_XX"), NA)
    }
    if (!is.na(get(paste0("N_Stayers",IVt,pl,"_XX"))) & get(paste0("N_Stayers",IVt,pl,"_XX")) < 2) {
      assign(paste0("N_Switchers",IVt,pl,"_XX"), nrow(df))
      assign(paste0("N_Stayers",IVt,pl,"_XX"), 0)
    }
    if (!is.na(get(paste0("N_Switchers",IVt,pl,"_XX"))) & get(paste0("N_Switchers",IVt,pl,"_XX")) == 0) {
      assign(paste0("N_Switchers",IVt,pl,"_XX"), 0)
      assign(paste0("N_Stayers",IVt,pl,"_XX"), nrow(df))
    }
    df[[paste0("abs_delta_D_",pairwise,pl,"_XX")]] <- NA
    df[[paste0("S_",pairwise,pl,"_XX")]] <- NA
    if (aoss == 1) {
      assign(paste0("P_",pairwise,pl,"_XX"), 0)
    }
    if (waoss == 1) {
      assign(paste0("E_abs_delta_D_",pairwise,pl,"_XX"), 0)
    }
    if (ivwaoss == 1) { 
      assign(paste0("denom_delta_IV_",pairwise,pl,"_XX"), 0)
    }
    
    assign(paste0("non_missing_",pairwise,pl,"_XX"), 0)
  }
  
  df <- df[order(df$ID_XX), ]
  to_keep <- c("ID_XX", paste0("Phi_1_",pairwise,pl,"_XX"), paste0("Phi_2_",pairwise,pl,"_XX"), paste0("Phi_3_",pairwise,pl,"_XX"), paste0("S_",pairwise,pl,"_XX"), paste0("abs_delta_D_",pairwise,pl,"_XX"), paste0("used_in_",pairwise,pl,"_XX"), paste0("inner_sum_IV_denom_",pairwise,pl,"_XX"), cluster) 
  
  df <- df %>% select(dplyr::any_of(to_keep))
  
  ## End of the program
  
  for (v in names(scalars)) {
    scalars[[v]] <- get(v)
  }
  if (aoss == 1) {
    scalars[[paste0("P_",pairwise,pl,"_XX")]] <- get(paste0("P_",pairwise,pl,"_XX"))
  }
  if (waoss == 1) {
    scalars[[paste0("E_abs_delta_D_",pairwise,pl,"_XX")]] <- get(paste0("E_abs_delta_D_",pairwise,pl,"_XX"))
  }
  
  if (ivwaoss == 1) {
    scalars[[paste0("denom_delta_IV_",pairwise,pl,"_XX")]] <- get(paste0("denom_delta_IV_",pairwise,pl,"_XX"))
  }
  
  scalars[[paste0("non_missing_",pairwise,pl,"_XX")]] <- get(paste0("non_missing_",pairwise,pl,"_XX"))
  for (v in c("Switchers", "Stayers")) {
    if (waoss == 1 | aoss == 1) {
      scalars[[paste0("N_",v,"_1_",pairwise,pl,"_XX")]] <- get(paste0("N_",v,pl,"_XX"))
      scalars[[paste0("N_",v,"_2_",pairwise,pl,"_XX")]] <- get(paste0("N_",v,pl,"_XX"))
    } else if (ivwaoss == 1) {
      scalars[[paste0("N_",v,"_3_",pairwise,pl,"_XX")]] <- get(paste0("N_",v,"_IV",pl,"_XX"))
    }
  }
  
  estims <- c("aoss", "waoss", "ivwaoss")
  indices <- c() 
  for (j in 1:length(estims)) {
    if (get(estims[j]) == 1) {
      indices <- c(indices, j)
    }
  }
  
  for (i in indices) {
    scalars[[paste0("delta_",i,"_",pairwise,pl,"_XX")]] <- get(paste0("delta_",i,"_",pairwise,pl,"_XX"))
    scalars[[paste0("sd_delta_",i,"_",pairwise,pl,"_XX")]] <- get(paste0("sd_delta_",i,"_",pairwise,pl,"_XX"))
    scalars[[paste0("LB_",i,"_",pairwise,pl,"_XX")]] <- get(paste0("LB_",i,"_",pairwise,pl,"_XX"))
    scalars[[paste0("UB_",i,"_",pairwise,pl,"_XX")]] <- get(paste0("UB_",i,"_",pairwise,pl,"_XX"))
  }
  
  out_res <- list(scalars, list(df))
  names(out_res) <- c("scalars", "to_add")
  return(out_res)
}

#' Internal function of did_multiplegt_stat that emulates Stata predict function.
#' @param df df
#' @param varname varname
#' @param model model
#' @param varlist varlist
#' @param const const
#' @param prob prob
#' @param factor factor
#' @returns The same input dataframe df with an added column of predicted values.
#' @noRd
lpredict <- function(
    df,
    varname,
    model,
    varlist,
    const = TRUE,
    prob = FALSE,
    factor = FALSE
) {
  
  sensitivity <- 10^-10
  df[[varname]] <- 0
  if (isTRUE(factor)) {
    singletons <- subset(varlist, sapply(varlist, function(x) length(grep(":", x)) == 0))
    varlist <- names(model$coefficients)[2:length(names(model$coefficients))]
    
  }
  for (v in varlist) {
    if (is.na(model$coefficients[[v]])) {
      next
    } else if (!is.na(model$coefficients[[v]]) & isTRUE(factor) & grepl("FACT",v,fixed = TRUE)) {
      var_set <- var_extract(str = v, vars = singletons)
      df[[paste0("sel",v)]] <-t(matrix(1,1,length(var_set$var)) %*% (t(as.matrix(df[var_set$var])) == var_set$val)) == length(var_set$var)
      df[[varname]] <- ifelse(df[[paste0("sel",v)]] == 1, 
                              df[[varname]] + model$coefficients[[v]], df[[varname]])
      df[[paste0("sel",v)]] <- NULL
    } else if (!is.na(model$coefficients[[v]]) & (isFALSE(factor) |(isTRUE(factor) &  !grepl("FACT",v,fixed = TRUE)))) {
      if (!grepl(":",v,fixed = TRUE)) {
        df[[varname]] <- df[[varname]] + df[[v]] * model$coefficients[[v]]
      } else  {
        df$interact_tmp <- 1
        for (var in strsplit(v,":")[[1]]) {
          df$interact_tmp <- df$interact_tmp * df[[var]]
        }
        df[[varname]] <- df[[varname]] + df$interact_tmp * model$coefficients[[v]]
        df$interact_tmp <- NULL
      }
    }
  }
  if (isTRUE(const)) {
    df[[varname]] <- df[[varname]] + model$coefficients[1]
  }
  if (isTRUE(prob)) {        
    df[[varname]] <- exp(df[[varname]]) / (1 + exp(df[[varname]]))
    df[[varname]] <- ifelse(is.nan(df[[varname]]), 1, df[[varname]])
    df[[varname]] <- ifelse(df[[varname]] < sensitivity, 0, df[[varname]])
  }
  return(df)
}

#' Internal function of did_multiplegt_stat that generates tre string handle for interactions.
#' @param df df
#' @param str str
#' @param vars vars
#' @importFrom stringr str_extract_all
#' @returns A list with the string name of the factor variable and one of the values associated to it.
#' @noRd
var_extract <- function(str, vars) {
  if (length(vars) > 25) {
    stop("Interaction limit (25) exceeded. Reduce number of other treatments.")
  }
  repl <- sapply(1:length(vars), function(x) paste0(intToUtf8(64 + x),"_XX"))
  for (i in 1:length(vars)) {
    str <- gsub(vars[i], repl[i], str)
  }
  num <- as.numeric(str_extract_all(str,"\\d+")[[1]])
  str <- paste(str_extract_all(str,"\\D+")[[1]], collapse = "")
  for (i in 1:length(vars)) {
    str <- gsub(repl[i], vars[i], str)
  }
  
  return(list(var = strsplit(str,":")[[1]], val = num))
}

power_set <- function(set, sep = ":") {
  pwdf <- data.frame(subset = "", len = 0, id = 1:(2^length(set) -1))
  for (j in 1:nrow(pwdf)) {
    temp <- set[rev(set_fill(as.bin(j), length(set))) * 1:length(set)]
    pwdf$subset[j] <- paste(temp, collapse = sep)
    pwdf$len[j] <- length(temp)
  }
  return(pwdf[order(pwdf$len, pwdf$subset), ]$subset)    
}

as.bin <- function(x) {
  bin <- rev(as.integer(intToBits(x)))
  return(bin[match(1, bin, length(bin)):length(bin)])
}

set_fill <- function(set, n) {
  if (length(set) < n) {
    set <- c(rep(0,n-length(set)), set)
  }
  return(set)
}

#' @title summary method for did_multiplegt_stat
#' @name summary.did_multiplegt_stat
#' @description A customized printed display for did_multiplegt_stat output
#' @param object A did_multiplegt_stat object
#' @param ... Undocumented
#' @returns No return, just a custom summary method for did_multiplegt_stat output.
#' @export
summary.did_multiplegt_stat <- function(object, ...) {
  estims <- list(0, 1, 2)
  names(estims) <- c("aoss", "waoss", "ivwaoss") 
  
  if (is.null(object$args$estimator) & is.null(object$args$Z)) {
    estim_list <- c("aoss","waoss")
  } else if (is.null(object$args$estimator) & !is.null(object$args$Z)) {
    estim_list <- "ivwaoss"
  } else {
    estim_list <- object$args$estimator
  }
  
  if (is.null(object$args$by)) {
    by_levs <- c("_no_by")
    by_obj <- c("results")
  } else {
    by_levs <- object$by_levels
    by_obj <- c()
    for (res in 1:length(by_levs)) {
      by_obj <- c(by_obj, paste0("results_by_", res))
    }
    
    if (!is.null(object$args$by_fd)) {
      by_name <- "quantiles"
    } else if (!is.null(object$args[["by"]])) {
      by_name <- object$args$by
    }
    by_totl <- length(by_levs)
    
    cat("\n");
    cat(noquote(strrep("#", 70)));cat("\n");
    cat(sprintf("## did_multipegt_stat by %s (%.0f levels)", by_name, by_totl));cat("\n");
    cat(noquote(strrep("#", 70)));cat("\n");
  }
  
  for (temp in 1:length(by_obj)) {
    print_obj <- object[[by_obj[temp]]]
    
    if (by_levs[temp] != "_no_by") {
      msg <- paste0(" By level: ", by_levs[temp])
      cat(noquote(strrep("#", 70 - nchar(msg))));cat(msg);cat("\n");
    }
    
    cat("\n");
    cat(noquote(strrep("-", 35)));cat("\n");
    if ("ivwaoss" %in% estim_list) {
      strdisplay("N",print_obj$table[2*print_obj$pairs+1,5] + print_obj$table[2*print_obj$pairs+1,6])
    } else {
      if ("waoss" %in% estim_list) {
        strdisplay("N",print_obj$table[print_obj$pairs+1,5] + print_obj$table[print_obj$pairs+1,6])
      } else {
        strdisplay("N",print_obj$table[1,5] + print_obj$table[1,6])
      }
    }
    methods <- list(ra = "Reg. Adjustment", dr = "Doubly Robust", ps = "Propensity Score")
    method <- ifelse(is.null(object$args$estimation_method), "dr", object$args$estimation_method)
    method <- ifelse(isTRUE(object$args$exact_match), "ra", method)
    for (m in c("waoss", "ivwaoss")) {
      if (m %in% estim_list) { 
        strdisplay(paste0(toupper(m), " Method"), methods[[method]])            
      }
    }
    if (isFALSE(object$args$exact_match)) {
      strdisplay("Polynomial Order",object$args$order)
    }
    support <- c("Exact Matching", "No Extrapolation")
    index <- 1
    for (m in c("exact_match", "noextrapolation")) {
      if (isTRUE(object$args[[m]])) {
        strdisplay("Common Support",support[index])
      }
      index <- index + 1
    }
    if (!is.null(object$args$switchers)) {
      strdisplay("Switchers", object$args$switchers)
    }
    cat(noquote(strrep("-", 35)));cat("\n");
    if (!is.null(object$args$cluster)) {
      if (object$args$cluster != object$args$ID) {
        cat(sprintf("(Std. errors adjusted for %.0f clusters in %s)\n", 
                    print_obj$n_clusters[[1]], object$args$cluster))
      }
    }
    
    for (t in names(estims)){
      if (t %in% estim_list) {
        
        cat("\n");
        cat(noquote(strrep("-", 70)));cat("\n");
        cat(strrep(" ", 20));cat(sprintf("Estimation of %s(s)", toupper(t)));cat("\n");
        cat(noquote(strrep("-", 70)));cat("\n");
        
        l_bound <- 1 + estims[[t]] * print_obj$pairs 
        u_bound <- l_bound + isTRUE(object$args$disaggregate) * (print_obj$pairs - 1)
        mat_sel <- print_obj$table[l_bound:u_bound, ]
        mat_print(mat_sel, t)
        cat("\n");
        
        if (isTRUE(object$args$placebo)) {
          
          cat("\n");
          cat(noquote(strrep("-", 70)));cat("\n");
          cat(strrep(" ", 15));cat(sprintf("Estimation of %s(s) - Placebo", toupper(t)));cat("\n");
          cat(noquote(strrep("-", 70)));cat("\n");
          
          l_bound <- 1 + estims[[t]] * print_obj$pairs 
          u_bound <- l_bound + isTRUE(object$args$disaggregate) * (print_obj$pairs - 1)
          mat_sel_placebo <- print_obj$table_placebo[l_bound:u_bound, ]
          
          mat_print(mat_sel_placebo, t)
          cat("\n");
        }
      }
    }
    
    if (isTRUE(object$args$aoss_vs_waoss)) {
      
      cat("\n");
      cat(noquote(strrep("-", 70)));cat("\n");
      cat(strrep(" ", 15));cat("Difference test: AOSS and WAOSS");cat("\n");
      cat(noquote(strrep("-", 70)));cat("\n");
      cat("H0: AOSS = WAOSS\n");
      
      tab_print(print_obj$aoss_vs_waoss)
    }
  }
}

#' @title print method for did_multiplegt_stat
#' @name print.did_multiplegt_stat
#' @description A customized printed display for did_continous output
#' @param x A did_multiplegt_stat object
#' @param ... Undocumented
#' @returns No return, just a custom summary print for did_multiplegt_stat output.
#' @export
print.did_multiplegt_stat <- function(x, ...) {
  summary(x)
}

#' Ancillary function for print/summary methods
#' @param mat mat
#' @param name name
#' @returns No return, just printing output.
#' @noRd
mat_print <- function(mat, name) {
  if (inherits(mat,"matrix")) {
    dis <- matrix(data = 0, nrow = nrow(mat), ncol = ncol(mat))
    dis[,1:4] <- sprintf("%s", format(round(mat[,1:4], 5), big.mark=",", scientific=FALSE, trim=TRUE))
    dis[,5:ncol(dis)] <- 
      sprintf("%s", format(round(mat[,5:ncol(dis)], 0), big.mark=",", scientific=FALSE, trim=TRUE))
    rownames(dis) <- rownames(mat)
    colnames(dis) <- colnames(mat)
    print(noquote(dis[, , drop = FALSE]))
  } else {
    new_mat <- t(as.matrix(mat))
    rownames(new_mat) <- toupper(name)
    mat_print(new_mat) 
  }
}

#' Ancillary function for print/summary methods
#' @param mat mat
#' @returns No return, just printing output.
#' @noRd
tab_print <- function(mat) {
  if (inherits(mat,"matrix")) {
    dis <- matrix(data = 0, nrow = nrow(mat), ncol = ncol(mat))
    dis[,1:ncol(dis)] <- sprintf("%s", format(round(mat[,1:ncol(mat)], 5), big.mark=",", scientific=FALSE, trim=TRUE))
    rownames(dis) <- rownames(mat)
    colnames(dis) <- colnames(mat)
    print(noquote(dis[, , drop = FALSE]))
  }
}

#' Ancillary function for print/summary methods
#' @param objs string object
#' @param objn numeric object
#' @returns No return, just printing output.
#' @noRd
strdisplay <- function(objs, objn) {
  ltot1 <- 16; ltot2 <- 16;
  out1 <- ifelse(nchar(objs) <= ltot1, paste0(objs,strrep(" ",ltot1 - nchar(objs))), substr(objs, 1, ltot1))
  if (inherits(objn, "character")) {
    out2 <- ifelse(nchar(objn) <= ltot2, paste0(strrep(" ",ltot2 - nchar(objn)), objn), substr(objn, 2, ltot2))
    cat(paste0(out1," = ",out2),"\n")
  } else {
    objns <- sprintf("%.0f", objn)
    strdisplay(objs, objns)
  }
}


#' @title rnames method for did_multiplegt_stat
#' @name rnames.did_multiplegt_stat
#' @description A customized rnames method for did_multiplegt_stat output
#' @param obj A did_multiplegt_stat object
#' @param ... Undocumented
#' @import rnames
#' @returns The same output as rnames.
#' @export
rnames.did_multiplegt_stat <- function(obj, ...) {
  class(obj) <- "list"
  return(rnames(obj = obj, ignore = c("by_fd_graph", "by_graph", "cdf_plot", "args")))
}

#' Internal function of did_multiplegt_stat that emulates Stata logit function.
#' @param formula formula
#' @param df df
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats binomial glm predict
#' @returns A glm object.
#' @noRd
stata_logit <- function(
    formula,
    df
) {
  suppressWarnings({
    model <- glm(formula, data = df, weights = df$weight_XX, family = binomial(link = 'logit'), maxit = 300, epsilon = 10^-8)
  })
  return(model)
}

#' Customized function for mean
#' @param object object
#' @param df df
#' @param w w
#' @importFrom stats weighted.mean
#' @returns A scalar.
#' @noRd
Mean <- function(
    object,
    df,
    w = "weight_XX"
) {
  return(weighted.mean(x = df[[object]], w = df[[w]], na.rm = TRUE))    
}

#' Customized function for sd
#' @param object object
#' @param df df
#' @param w w
#' @importFrom Hmisc wtd.var
#' @returns A scalar.
#' @noRd
Sd <- function(
    object,
    df,
    w = "weight_XX"
) {
  return(sqrt(wtd.var(x = df[[object]], weights = df[[w]], na.rm = TRUE)))
}

#' Customized function for sum_w
#' @param df df
#' @param w w
#' @returns A scalar.
#' @noRd
wSum <- function(
    df,
    w = "weight_XX"
) {
  return(sum(df[[w]], na.rm = TRUE))
}

#' Customized function for sum
#' @param object object
#' @param df df
#' @param w w
#' @returns A scalar.
#' @noRd
Sum <- function(
    object,
    df,
    w = "weight_XX"
) {
  df_nm <- subset(df, !(is.na(df[[object]] | is.na(df[[w]]))))
  return(as.numeric(t(df_nm[[object]]) %*% df_nm[[w]]))
}

#' By option consistency check 
#' @param df df
#' @param ID ID
#' @param by by 
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @returns A logical.
#' @noRd
by_check <- function(
    df,
    ID,
    by
) {
  df$temp_G <- as.numeric(factor(df[[by]]))
  df <- df %>% group_by(.data[[ID]]) %>% mutate(sd_ID = sd(.data$temp_G, na.rm = TRUE))
  return(mean(df$sd_ID) == 0)
}

adj_fact <- function(str, df, check) {
  if (isTRUE(check)) {
    base <- substr(str, 1, find_last(str,"+") + 1)
    str <- substr(str, find_last(str,"+")+2, nchar(str))
    vars <- strsplit(str, " * ", fixed = TRUE)[[1]]
    str_temp <- ""
    for (v in vars) {
      if (length(levels(factor(df[[v]]))) > 1) {
        str_temp <- paste0(str_temp,v," * ")                
      }
    }
    str <- paste0(base, substr(str_temp, 1, nchar(str_temp)-3))
  }
  return(str)
}

find_last <- function(str, char) {
  pos <- 0
  for (n in 1:nchar(str)) {
    if (substr(str,n,n) == char) {
      pos <- n
    }
  }
  return(pos)
}
