
### Project: Poliovirus sequencing from Matlab
### Purpose: Estimate the size of the transmission bottleneck.
### Working directory: Poliovirus_Intrahost

# ================================= Read in data, import packages ============================

library(tidyverse)
library(wesanderson)
library(bbmle)

# ======================================== Bottleneck support functions =========================================

is_wholenumber <- function(x, tol = .Machine$double.eps^0.5)
{
  abs(x - round(x)) < tol
}

dzpois <- function(x, lambda) # the probability mass function of the zero-truncated poisson distribution
{ 
  if(x <= 0 | is_wholenumber(x) == FALSE)
  {
    return(0)
  } else
  {
    (exp(-1*lambda)*lambda^x) / (factorial(x)*(1-exp(-1*lambda)))
  }
}
dzpois <- Vectorize(dzpois, vectorize.args = c("x"))

p_all <- function(p,n)
{ # probability all success - only finding the variant at frequency p in n draws
  p^n
}
p_all <- Vectorize(p_all, vectorize.args = "n")

math_fit <- function(data, Nb_max, model, threshold, acc)
{
  # this  calculation is for each position in the genome.
  stopifnot(length(unique(data$chr))==1, length(unique(data$pos))==1)
  # and each model.
  if(model == "PA")
  {
    #  In this fit we take the minority frequency to be  correct
    # and set the major frequency to 1-minority. This accounts for
    # the fact that frequencies are related by the expression :
    # minor allele + major allele + errror =1.
    #Here we make the major allele frequency = major allele + error.
    #The error is always small. if it exceeds 1% then we throw an error here.
    
    if(1 - sum(data$freq1) > 0.01)
    {
      warning("The sum of the frequencies is less than 99% make sure the assumption major freq = 1-minor freq is valid")
    }
    data$freq1[data$freq1==max(data$freq1)] <- 1 - min(data$freq1)
    
    found <- data[data$found == TRUE, ] # only alleles that were transmitted
    Nb <- 1:Nb_max 
    
    if(nrow(found) == 0 | nrow(found) > 2)
    {
      print(nrow(found))
      stop(paste0("No variant transmitted for this site or",
                  "there are more than 2 variants here"))
    } else if(nrow(found) == 1)
    { # one variant found here. All successes
      prob <- p_all(p = found$freq1, n = Nb)
      # this is a vector of probabilities for each
      # n prob[i]= the probability of only getting that
      # variant in Nb[i] (i.e. draws)
    } else if(nrow(found) == 2)
    {
      # if at least one of each allele was transmitted
      
      first_var <- p_all(p = found$freq1[1], n = Nb) # all this one
      second_var <- p_all(p = found$freq1[2], n = Nb) # all the other one
      one_each <- 1 - (first_var + second_var)
      # at least one of each -
      # This is a vector as above since R adds and subtracts the
      # elements of the vectors as expected
      prob <- one_each
    }
    
  } else if(model == "BetaBin")
    {
    if(nrow(data[data$freq1 > 0.5, ]) > 0)
    {
      data<-data[df$freq1 < 0.5, ]
      warning("The beta binomials model only uses minor alleles. Subsetting the data now.")
    }
    # 2 is the recipient
    v_r = data$freq2
    v_d = data$freq1
    gc_ul = data$gc_ul2
    threshold = 0.05
    
    Nb <- 1:Nb_max # Here are the bottlenecks
    
    prob = L.Nb.beta(v_r, v_d, Nb, gc_ul, threshold, acc)
  }
  return(tibble(Nb = Nb, prob = prob))
  
}

trans_fit <- function(data, Nb_max, model, threshold, acc, ...)
{
  
  group <- rlang::quos(...,Nb)
  original_group <- rlang::quos(...)
  probs <- data %>% dplyr::group_by(chr, pos, pair_id) %>% dplyr::do(math_fit(.data, Nb_max, model, threshold, acc))
  # For each genomic position in question
  
  #control for number of the mutations in each donor.
  counts <- data %>% dplyr::group_by(!!!original_group) %>%
    dplyr::summarise(donor_mutants = length(which(freq1 > 0 & freq1 < 0.5)))
  max_mut <- max(counts$donor_mutants)
  LL.df <- probs %>% dplyr::left_join(counts, by = "pair_id") %>%
    dplyr::group_by(!!!group) %>%
    dplyr::summarize(LL = sum(log(prob)), weighted_LL = (max_mut / unique(donor_mutants)) * LL)
  # Get the  log likelihood of the each lambda for this pair - we can sum across pairs later.
  return(LL.df)
}

dist_prob <- function(data, weight, ddist, ...)
{
  params <- rlang::enquos(...)
  ddist <- rlang::enexpr(ddist)
  
  Nb = unique(data$Nb) # probability of data
  data <- dplyr::mutate(data, prob_D = exp(LL))
  # prob of Nb
  prob_Nb <- rlang::eval_tidy(quo(Nb %>% purrr::map_dbl(.f = !!ddist,!!!params)))
  data <- data %>%
    dplyr::left_join(dplyr::tibble(Nb = Nb, prob_Nb = prob_Nb), by = "Nb") %>%
    dplyr::mutate(prob_D_and_Nb = prob_D*prob_Nb)
  
  l_by_pair <- data %>% 
    dplyr::group_by(pair_id) %>%
    dplyr::summarize(prob_D_given_dist = sum(prob_D_and_Nb), LL_D_given_dist = log(prob_D_given_dist)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(weighted_total_LL = weight$weight_factor[weight$pair_id == pair_id] * LL_D_given_dist)
  
  return(-sum(l_by_pair$weighted_total_LL))
  
}

dist_prob_wrapper<-function(ddist,params)
{
  f_string <- paste0("function(data,weight,",params,"){dist_prob(data,weight,",ddist,",", params,")}")
  eval(parse(text = f_string))
}


# ============================================ Estimate bottleneck with PA model =====================================

library(bbmle)

transmission_variants_donorPolymorphic <- read.csv("data/processed/variants_donor_polymorphic.csv", stringsAsFactors = FALSE)
transmission_variants_donorPolymorphic_run <- transmission_variants_donorPolymorphic %>% mutate(chr = "Sabin2ref")
trans_freq.comp <- transmission_variants_donorPolymorphic_run %>% mutate(pair_id = houseId, found = ifelse(freq2 == 0, FALSE, TRUE))
pa_total_fit <- trans_fit(trans_freq.comp, Nb_max = 100, model = "PA", threshold = NULL, acc = NULL, pair_id)

# Look at weighted log likelihoods
ggplot(pa_total_fit, aes(x = Nb, y = weighted_LL, color = as.factor(pair_id))) +
  geom_point()

counts <- trans_freq.comp %>% 
  group_by(pair_id) %>%
  summarize(donor_mutants = length(which(freq1 > 0 & freq1 < 0.5))) %>%
  mutate(weight_factor_kk = max(donor_mutants)/donor_mutants, weight_factor = 1)
zdpois_fit <- dist_prob_wrapper(ddist = "dzpois", params = "lambda")
dzpois_model_fit <- bbmle::mle2(minuslogl = zdpois_fit, start = list(lambda = 1), data = list(data = pa_total_fit, weight = counts))
con_int <- confint(dzpois_model_fit)
dzpois_model_fit_profile <- profile(dzpois_model_fit)
plot(dzpois_model_fit_profile)
AIC(dzpois_model_fit)
BIC(dzpois_model_fit)
mean_zpois <- function(l) l/(1 - exp(-1*l))
mean_zpois(dzpois_model_fit@coef)

# Look at fitted poisson distribution
events <- 1:15
density <- dpois(x = events, lambda = dzpois_model_fit@coef)
prob <- ppois(q = events, lambda = dzpois_model_fit@coef, lower.tail = TRUE)
df <- data.frame(events, density, prob)
ggplot(df, aes(x = factor(events), y = density)) +
  geom_point() +
  geom_line(data = df, aes(x = events, y = density)) +
  geom_line(data = df, aes(x = events, y = prob)) +
  theme_bw()

w <- 0.03
step <- 0.015
windows <- tibble(s = seq(0.02, 1 - w, by = step), end = s + w)

out <- windows %>% rowwise() %>%
  mutate(iSNV = nrow(filter(trans_freq.comp, freq1 >= s, freq1 < end)),
         transmitted = nrow(filter(trans_freq.comp, freq1 >= s, freq1 < end, found == TRUE)),
         freq = mean(trans_freq.comp$freq1[which(trans_freq.comp$freq1 >= s & trans_freq.comp$freq1 < end)]),
         prob = transmitted/iSNV,
         error_bottom = qbinom(c(0.025), iSNV, prob)/iSNV,
         error_top = qbinom(c(0.975), iSNV, prob)/iSNV,
         many = iSNV > 5)

Pt_PA <- function(x, l, max_Nb)
{
  s <- 0
  for(i in 1:max_Nb)
  {
    c <- ((1-x)^i) * (l^i)/((exp(l)-1)*factorial(i) )
    s <- s + c
  }
  return(1 - s)
}

model <- tibble(s = seq(0, 1, 0.01))
model <- mutate(model, prob = Pt_PA(s,dzpois_model_fit@coef,100),
                lower = Pt_PA(s,con_int[1],100),
                upper = Pt_PA(s,con_int[2],100))

palette <- wesanderson::wes_palette("Darjeeling1")
palette_viridis <- c("#FDE725FF", "#73D055FF", "#29AF7FFF", "#238A8DFF", "#39568CFF", "black")

# bottleneck of 10
model_10 <- tibble(s = seq(0, 1, 0.01))
model_10 <- mutate(model, prob = Pt_PA(s, 10, 100))

bottleneck.plot.cf10 <- ggplot() + 
  geom_point(data = out, aes(x = freq, y = prob, alpha = many)) +
  geom_errorbar(data = out, aes(x = freq, ymin = error_bottom, ymax = error_top, alpha = many)) +
  geom_line(data = model, aes(x = s, y = prob), color=palette[1]) +
  geom_ribbon(data = model, aes(x = s, ymin = lower, ymax = upper), alpha = 0.5, fill = palette[1]) +
  scale_alpha_manual(values = c(0,1)) +
  theme(legend.position = 'none') +
  xlab("Frequency in Donor") +
  ylab("Probability of transmission") +
  geom_point(data = trans_freq.comp, aes(x = freq1, y = as.numeric(found) + (as.numeric(found)-0.5)/10), alpha = 0.5) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  geom_line(data = model_10, aes(x = s, y = prob), color = palette[5]) +
  theme_bw() + 
  theme(legend.position = "none",  axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold")) # 8 by 5

# Save bottleneck info
PA_info <- data.frame(model = "PA", lambda = as.numeric(dzpois_model_fit@coef), lower = as.numeric(con_int[1]), upper = as.numeric(con_int[2]))
#write.csv(PA_info, "data/processed/bottleneck_estimates.csv")

# ========================================== Beta-binomial functions =================================

GetGenomeCopies <- function(gc_ul)
{
  if(gc_ul >= 4.5e4){
    x=4.5e4} else if (gc_ul >= 9e3 & gc_ul <= 4.5e4){
      x=9e3} else if (gc_ul >= 9e2 & gc_ul <= 9e3){
        x=9e2} else if (gc_ul >= 9 & gc_ul <= 9e2){
          x=9}
}
GetGenomeCopies.v <- Vectorize(GetGenomeCopies)

#' @describeIn L.Nb.beta likelihood the variant is  found in recipient

like_found.beta<-function(v_r,v_d,Nb,gc_ul,threshold,acc=accuracy_stringent)
{
  acc_gc = GetGenomeCopies.v(gc_ul)

  # If the frequency is below our threshold then the probability of being found is 0
  if(v_r<threshold){
    return(0)
  }
  # This will round the gc down to the nearest log as discussed below.
  if(v_r>=0.02 & v_r<0.05){
    sense= acc$sensitivity[which(acc$gc_ul==acc_gc & acc$freq==0.02)]
  }else if(v_r>=0.05 & v_r<0.1){
    sense= acc$sensitivity[which(acc$gc_ul==acc_gc & acc$freq==0.05)]
  }
  else {
    sense=1
  }
  prob = c()
  if(v_r<=(1-threshold)){
    for(k in 0:Nb){
      prob[k+1]=dbeta(x=v_r,shape1 = k,shape2=Nb-k)*dbinom(x=k,size=Nb,prob = v_d)*sense
    }
    prob=sum(prob)
  }else if(v_r>(1-threshold)){
    # it is fixed - whats the probability the other allele was not found
    lost_allele_freq = 1-v_d
    # The for loop over k and sum is in the like_lost.beta.uncert function
    prob = like_lost.beta.uncert(v_r=0,v_d=lost_allele_freq,Nb,gc_ul,threshold,acc=acc)
  }
  return(prob)
}

#' @describeIn L.Nb.beta likelihood the variant is not found in recipient
like_lost.beta.uncert<-function(v_r,v_d,Nb,gc_ul,threshold,acc=accuracy_stringent){ # sum over k
  stopifnot(v_r<threshold) # ensure not found in the recipient
  
  acc_gc = GetGenomeCopies.v(gc_ul)
  
  # This will round the gc down to the nearest log as discussed below.
  prob=c()
  for(k in 0:Nb){
    uncert_term=c()
    f=c(0.02,0.05,0.10)
    for(i in 1:(length(f)-1)){
      uncert=1-acc$sensitivity[which(acc$gc_ul==acc_gc & acc$freq==f[i])]
      # The prob the variant is missed because itis between f[i] and f[i+1] given
      # the sample size
      uncert_term[i]=(pbeta(q = f[i+1],shape1 = k,shape2 = Nb-k)-
                        pbeta(q = f[i],shape1 = k,shape2 = Nb-k))*
        dbinom(x=k,size=Nb,prob = v_d)*uncert
    }
    #probability the variant is below the cut off or present but missed
    prob[k+1]=pbeta(q = threshold,shape1 = k,shape2 = Nb-k)*dbinom(x=k,size=Nb,prob = v_d)+sum(uncert_term)
  }
  sum(prob)
}

#' Betabinomial Likelihood functions
#'
#' What is the probability of observing the data given a bottleneck
#' size under the beta binomial model. This is a wrapper around the two likelihood
#' functions
#' @param v_r allele frequency in the recipient
#' @param v_d allele frequency in the donor
#' @param Nb Bottleneck size (may be a vector)
#' @param gc_ul titer in recipeint sample
#' @param threshold detection level threshold
#' @param acc a data frame with accrucacy data must contain, columns
#'  freq,sensitivity, and gc_ul
#'
#' @return The likelihood of observing the data given the bottleneck size.
#' @export
L.Nb.beta<-function(v_r,v_d,Nb,gc_ul,threshold,acc=accuracy_stringent){
  if(v_r>=threshold){
    like=like_found.beta(v_r=v_r,v_d=v_d,Nb=Nb,gc_ul=gc_ul,threshold=threshold,acc=acc)
  }else if(v_r<threshold){ # Not found
    like=like_lost.beta.uncert(v_r=v_r,v_d=v_d,Nb=Nb,gc_ul=gc_ul,threshold=threshold,acc=acc)
  }
  return(like)
}
L.Nb.beta<-Vectorize(L.Nb.beta,vectorize.args = "Nb")


#' Simulate transmission Beta binomial model
#'
#' Simulate the transmission of alleles using the
#' betabinomial model. Essentially this takes
#' the frequency of the variant allele, and bottleneck
#' size and simulates the found column. It also selects a bottleneck size
#' based on lambda and a zero truncated Poisson. If simulation the whole data
#' set each pair should be run separately as to not use the same bottleneck
#'
#' @param data a data frame with with chr,pos,freq1, and  pair_id columns
#' @param lambda The lambda of the zero truncated Poisson
#' @param threshold limit of variant calling detection
#' @param acc a data frame with accuracy metrics
#'
#' @return data frame the same as data but with a simulated found column
#'
#' @examples
#' betabin_sim(small_trans,1.2,0.02,accuracy_stringent)
#' @export
betabin_sim<-function(data,lambda,threshold,acc=accuracy_stringent){
  
  pair<-unique(data$pair_id)
  if(length(pair)>1){
    warning(paste0("Running on ",length(pair)," pairs. All will have the same bottleneck."))
  }
  
  Nb<-rzpois(1,lambda)
  out<-data %>% dplyr::group_by(chr,pos,pair_id) %>%
    dplyr::do(betabin_sim_helper(.,Nb,threshold,acc))
  return(out)
  
  
  
}


#' beta bin sim helper
#'
#' Helper function to simulate beta binomial data. This takes in one loci.
#' the data frame should have 2 rows (one for each allele.)
#'  Essentially this takes
#' the frequency of the variant allele, and bottleneck
#' size and simulates the found column.
#'
#' @param data a data frame with with freq1 and pair_id columns and 2 rows
#' @param Nb the bottleneck size
#' @param threshold limit of variant calling detection
#' @param acc a data frame with accuracy metrics
#'
#' @return data frame the same as data but with a simulated found column
#'
betabin_sim_helper<-function(data,Nb,threshold,acc){
  # data refers to the 2 mutations at this position, bottlenecks-
  #  a data frame with the bottle_neck,
  #model - the name of the column in the bottlenecks that
  # contains the bottleneck size we want to use.
  if(nrow(data)!=2){
    stop(print("There are not 2 mutations at this point."))
  }
  pair<-unique(data$pair_id)
  
  if(length(pair)!=1){
    stop(print("There should only be one pair id at this point."))
  }
  
  minor<-min(data$freq1) # This is the minor allele frequency
  major<-1-min(data$freq1)
  
  # includes uncertainty -
  # with v_r = 0 or 1 these are probabilities otherwise they are densities.
  only_minor= L.Nb.beta(v_r=1,v_d=minor,Nb=Nb,gc_ul=unique(data$gc_ul2),threshold=threshold,acc)
  # includes uncertainty - The likelihood the allele was lost
  only_major = L.Nb.beta(v_r=1,v_d=major,Nb=Nb,gc_ul=unique(data$gc_ul2),threshold=threshold,acc) # includes uncertainty - The likelihood the allele was lost
  
  both = 1-(only_minor+only_major)
  
  # flip the coin
  pick = runif(1,0,1)
  data$found=F
  
  # 0  only_minor         only_major         both                        1
  # |------------------|---------------|---------------------------------|
  
  if(pick<=only_minor){
    # onlny the minor was found
    data$found[data$freq1==min(data$freq1)]=T
  }else if(pick>only_minor & pick<(only_minor+only_major)){
    # only the major was found
    data$found[data$freq1==max(data$freq1)]=T
  }else if(pick>=(only_minor+only_major)){
    # both were found
    data$found=T
  }
  return(data)
}

# ======================================== Estimate bottleneck with beta-binomial model ==========================

metadata <- read.csv("data/metadata/specimens_Sabin2Positive_CTbelow37_processedFinal_deID_runInfo_covInfo_week_final.csv", stringsAsFactors = FALSE)

# Add gc_ul
trans_freq.comp <- mutate(trans_freq.comp, gc_ul1 = metadata$S2_copiesPerUl_inTNA[match(anonSpecId1,metadata$anonSpecId)], gc_ul2 = metadata$S2_copiesPerUl_inTNA[match(anonSpecId2,metadata$anonSpecId)])

# Get accuracy metrics
accuracy_stringent <- read.csv("data/reference/accuracy_metrics_bbmodel_2.csv", stringsAsFactors = FALSE)

# Fit the model
beta_total_fit <- trans_fit(filter(trans_freq.comp, freq1 < 0.5), Nb_max = 100, model = "BetaBin", threshold = 0.05, acc = accuracy_stringent, pair_id)

# Look at weighted log likelihoods
ggplot(beta_total_fit, aes(x = Nb, y = weighted_LL, color = as.factor(pair_id))) +
  geom_point()

zdpois_fit <- dist_prob_wrapper(ddist = "dzpois", params = "lambda")

counts <- trans_freq.comp %>% 
  group_by(pair_id) %>%
  summarize(donor_mutants = length(which(freq1>0 & freq1<0.5))) %>%
  mutate(weight_factor_kk = max(donor_mutants)/donor_mutants, weight_factor = 1)

dzpois_model_fit_bb <- bbmle::mle2(minuslogl = zdpois_fit, start = list(lambda = 1), data = list(data = beta_total_fit, weight = counts))
conf_int_BB <- bbmle::confint(dzpois_model_fit_bb)
summary(dzpois_model_fit_bb)
AIC(dzpois_model_fit_bb)
BIC(dzpois_model_fit_bb)
mean_zpois(dzpois_model_fit_bb@coef)

# ============================== Get and plot individual bottlenecks =====================

straight_sum <-function(data,max_nb){
  Nb<-data$Nb[which(data$LL==max(data$LL))] # Get the  max lambda
  good_range<-subset(data,LL> (max(LL)-1.92)) # get the bottlenecks that fall in this region the 95% confidence intereval
  lower<-good_range$Nb[1]
  upper <- good_range$Nb[nrow(good_range)]
  if(length(Nb)>1){
    Nb=max(Nb)
  }
  if(Nb==max_nb){
    return(tibble(Nb=NA,lower_95=NA,
                  upper_95=NA))
  }else{
    return(tibble(Nb=Nb,lower_95=lower,
                  upper_95=upper))
  }
}

max_nb <- 100

pa_total_fit_individual <- pa_total_fit %>% group_by(pair_id) %>% do(straight_sum(., max_nb))

beta_total_fit_individual <- beta_total_fit %>% group_by(pair_id) %>% do(straight_sum(., max_nb))


PlotSingleBottleneck <- function(fit, CIs, pair)
{
  plot <- ggplot(filter(fit, pair_id == pair), aes(x = Nb, y = weighted_LL)) +
    geom_line() +
    xlab("Bottleneck size") +
    ylab("Log-Likelihood") +
    theme_bw() +
    geom_vline(xintercept = filter(CIs, pair_id == pair)$Nb, color = "red") +
    geom_vline(xintercept = filter(CIs, pair_id == pair)$lower_95, color = "red", linetype = "dotted") +
    geom_vline(xintercept = filter(CIs, pair_id == pair)$upper_95, color = "red", linetype = "dotted")
  
  return(plot)
}

pa.grid <- cowplot::plot_grid(PlotSingleBottleneck(pa_total_fit, pa_total_fit_individual, 115), 
                   PlotSingleBottleneck(pa_total_fit, pa_total_fit_individual, 171), 
                   PlotSingleBottleneck(pa_total_fit, pa_total_fit_individual, 702), 
                   PlotSingleBottleneck(pa_total_fit, pa_total_fit_individual, 927), 
                   nrow = 2, 
                   ncol = 2, 
                   labels = c("115", "171", "702", "927"), 
                   label_size = 10)

bb.grid <- cowplot::plot_grid(PlotSingleBottleneck(beta_total_fit, beta_total_fit_individual, 115), 
                              PlotSingleBottleneck(beta_total_fit, beta_total_fit_individual, 171), 
                              PlotSingleBottleneck(beta_total_fit, beta_total_fit_individual, 702), 
                              PlotSingleBottleneck(beta_total_fit, beta_total_fit_individual, 927), 
                              nrow = 2, 
                              ncol = 2, 
                              labels = c("115", "171", "702", "927"), 
                              label_size = 10)


