################################################################################
#' neglog
#'
#' Log transforms a vector. Handles non-positive values. Injective and strictly increasing function
#'
#' @param x vector or integer
#' @param log_base base of log function - Default is exp(1)
#'
#' @return vector of transformed values
#' @export
#'
#'
neglog = function(x, log_base = exp(1)) {
  result = sign(x)*log(abs(x)+1, base = log_base)
  return (result)
}

################################################################################

#' rem_outliers
#'
#' Shaves a certain percentage off each tail (of vector)
#'
#' @param vect vector
#' @param percentage which percentage to shave off each tail. Default = 0.01
#'
#' @return vector without outliers
#' @export
#'
rem_outliers = function(vect, percentage = 0.01) {
  quantiles = quantile(vect, probs=c(percentage, 1-percentage), na.rm = TRUE)
  return(vect[(vect <= quantiles[2]) & (vect >= quantiles[1])])
}

################################################################################
#' df_rem_outliers
#'
#' Shaves off rows where value of a given variable is an outlier
#'
#' @param df data frame to remove outliers from
#' @param col_name which column to consider outliers for (char)
#' @param percentage which percentage to remove (in each tail) - Default = 0.01
#'
#' @return data frame without outliers
#' @export
#'
df_rem_outliers = function(df, col_name, percentage = 0.01) {
  quantiles = quantile(df[,col_name], probs=c(percentage, 1-percentage), na.rm=TRUE)
  return(df[(df[,col_name] <= quantiles[2]) & (df[,col_name] >= quantiles[1]),])
}

################################################################################
#' sum_outliers
#'
#' Sums values of specified column of outliers and non-outliers of another specified column
#'
#' @param dat Data frame
#' @param outlier_col column name (char) to define which observations are considered outliers
#' @param sum_col which column (char) to sum ouliers for
#' @param percentage which percentage of each tail is considered outliers
#'
#' @return double
#' @export
#'
#'
sum_outliers = function(dat, outlier_col, sum_col, percentage = 0.1) {
  col = dat[, outlier_col]
  quantiles = quantile(col, probs=c(percentage, 1-percentage), na.rm = FALSE)
  outliers = dat[(col < quantiles[1]) | (col > quantiles[2]), ][, sum_col]
  non_outliers = dat[(col > quantiles[1]) & (col < quantiles[2]), ][, sum_col]
  return(c(outliers_sum = sum(outliers), non_outliers_sum = sum(non_outliers)))
}

################################################################################
#' normal_control
#'
#' Simple normal qqplots and histogram for a vector
#'
#' @param vect vector (numeric)
#'
#' @return None: Two plots
#' @export
#'
#'
normal_control = function(vect) {
  par(mfrow=c(1,2))
  qqnorm(vect)
  qqline(vect)
  hist(vect)
}

################################################################################
#' par_cor
#'
#' calculated partial correlation of columns of data frame with specified controls
#'
#' @param dat data frame with variables
#' @param x column name of first variable (char)
#' @param y column name of second variable (char)
#' @param controls vector of column names of control variables
#'
#' @return double
#' @export
#'
#'
par_cor = function(dat, x, y, controls=c()) {
  controls = controls[(controls != x) & (controls != y)]
  x_obs = dat[,x]
  y_obs = dat[,y]
  if (length(controls) == 0) {
    return (cor(x_obs, y_obs))
  } else if (length(controls) == 1) {
    resid_x = lm(x_obs~dat[,controls])$residuals
    resid_y = lm(y_obs~dat[,controls])$residuals
    return (cor(resid_x, resid_y))
  } else {
    controls_matrix = dat[,controls]
    resid_x = lm(x_obs~., data=controls_matrix)$residuals
    resid_y = lm(y_obs~., data=controls_matrix)$residuals
    return (cor(resid_x, resid_y))
  }
}

################################################################################

#' bs_CI
#'
#' Bootstrap is used e.g. for very skewed distributions when CLT is inaccurate. Bootstrap is not advised for large data sets, almost regardless of distribution
#' since the CIs are virtually the same (as CLT CI) in this case and bootstrap is calculation intense. Bootstrap t CI (number of resamples = 10^4)
#'
#' @param data vector
#' @param conf_level confidence level
#'
#' @return confidence interval (vector)
#' @export
#'
#'
bs_CI = function(data,conf_level=0.95){
  n = 10^4
  res = numeric(n)
  for (i in 1:n){
    x = sample(data,length(data),replace = T)
    res[i] = (mean(x)-mean(data)) / (sd(x)/sqrt(length(data)))
  }
  qs = quantile(res, c((1-conf_level)/2,1-(1-conf_level)/2))
  return (mean(data)-qs*sd(data)/sqrt(length(data)))
}


################################################################################

#'Funktionen additivitetsPlot beregner gennemsnit i hver gruppe
#' defineret ved A:B (A og B er faktorer) og afsætter disse mod A.
#' For hvert niveau af B forbindes gennemsnit.
#' Desuden angives 0.95 CI for middelværdien
#'
#' @param A faktor 1
#' @param B faktor 2
#' @param x response variable
#' @param ci type of ci (char)
#'
#' @return plots
#' @export
#'
#'
additivitetsPlot=function(A,B,x,ci=c('normal', 'bootstrap', 'none')){
  ci=match.arg(ci)
  nlevB=length(levels(B))
  NavnA <-deparse(substitute(A))
  NavnB <-deparse(substitute(B))
  interaction.plot(A,B,x,col=c(1:nlevB),fixed=TRUE,
                   xlab=NavnA,trace.label=NavnB,ylab="Average")
  me=aggregate(x,list(A:B),mean)
  sdv=aggregate(x,list(A:B),sd)
  levA=levels(A)
  nlevA=length(levA)
  n=aggregate(x,list(A:B),length)
  if (ci %in% c('normal', 'bootstrap')){
    if (ci == 'normal'){
      lower=me[,2]+sdv[,2]/sqrt(n[,2])*qnorm(0.025)
      upper=me[,2]-sdv[,2]/sqrt(n[,2])*qnorm(0.025)
    } else {
      bs_CIs=aggregate(x, list(A:B), bs_CI)
      lower=bs_CIs$x[,1]
      upper=bs_CIs$x[,2]
    }
    for (i in 1:nlevB){
      med=(c(1:nlevA)-1)*nlevB+i
      arrows(c(1:nlevA),lower[med],c(1:nlevA),upper[med],
             code=3,angle=90,length=0.05,col=i)
    }
  } else {
    for (i in 1:nlevB){
      med=(c(1:nlevA)-1)*nlevB+i
    }
  }
}


################################################################################
zip <- function(...) {
  mapply(list, ..., SIMPLIFY = FALSE)
}

enumerate <- function(...) {
  zip(ix=seq_along(..1), ...)
}


################################################################################

#' constant_friendly_t_test
#' Function returns a regular t.test except when vector is constant. When constant,
#' only "conf.int" element is returned
#'
#' @param vect vector of values to perform t-test on
#'
#' @return double
#' @export
#'
#'
constant_friendly_t_test = function(vect) {
  if (var(vect) <= 1e-15) {
    return (list(conf.int = c(mean(vect), mean(vect))))
  } else {
    return (t.test(vect))
  }
}


################################################################################


#' mark_significance
#'
#' Marks significance with an asterisk. Useful for both character and numeric
#' Always returns character
#'
#' @param x char or double
#' @param level confidence level to mark - Default = 0.05
#'
#' @return char
#' @export
#'
mark_significance = function(x, level=0.05) {
  if (as.numeric(x) <= level) {
    return (glue::glue('{x}*'))
  } else {
    return (as.character(x))
  }
}

################################################################################


#' hush
#'
#' Wrap another function call with hush() to silence output
#'
#' @param code code chunk to hush
#'
#' @return Output of original function
#' @export
#'
hush=function(code){
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}


################################################################################


#' linear_normal_model_test
#'
#' simple linear (normal univariate) model plots for model testing. Assuming dat only has columns
#' that are included in model, and every column in model is in dat
#'
#' @param model lm object
#' @param dat data frame from which model is generated
#'
#' @return plots
#' @export
#'
#'
linear_normal_model_test = function(model, dat) {
  par(mfrow=c(1, 1))
  # (assume independence of observations)

  # It is recommended to use "ggpairs(dat)" to check the following:
  # 1. Checking linear bivariate relationships between independent and dependent variables
  # 2. Checking for multicollinearity as well (-0.8 < correlations < 0.8)

  dat$resid = model$residuals
  qqnorm(dat$resid) # Normally distributed residuals
  qqline(dat$resid)
  hist(dat$resid)

  # Plotting residuals against all predictors, responses and fitted values (check for connection with variance)
  # Can also be used as a check for linear relationships, i.e. low absolute correlation with residuals means linear relationship
  # is a good approximation (not for response variable)
  sidelength = ceiling(sqrt(ncol(dat)-1))
  par(mfrow=c(sidelength, sidelength))
  for (col in colnames(dat)) {
    if (col == 'resid') {
      next
    }
    plot(
      dat[,col],
      dat$resid,
      main = col,
      xlab = col,
      ylab = 'residuals'
    )
  }
  plot(model$fitted.values, dat$resid)
}

################################################################################

#' expected_shortfall
#'
#' Calculate expected shortfall for a vector
#'
#' @param v vector (numeric)
#' @param q quantile for value-at-risk
#'
#' @return double
#' @export
#'
expected_shortfall = function(v, q=0.05) {
  q = quantile(v, 0.05)
  big_losses = v[v <= q]
  return (mean(big_losses))
}

