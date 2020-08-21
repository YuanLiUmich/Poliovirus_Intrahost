

# ====================== ROC functions ==========================

fill_in_plots <- function(x)
{
  spec <- min(x$adj.specificity[which(x$adj.specificity > 0.995)])
  sense <- max(x$adj.sensitivity[x$adj.specificity == spec]) # in case there is a step right here
  
  y <- data.frame(threshold = x$threshold[1], Id = x$Id[1], adj.specificity = 0.9950001, adj.sensitivity = sense, samp = unique(x$samp), exp.freq = x$exp.freq[1], FP = x$FP[1], TP = x$TP[1], sensitivity = x$sensitivity[1], specificity = x$specificity[1])
  
  if("gc" %in% names(x))
  {
    y$gc = x$gc[1]
  }
  return(y)
}

sum_roc <- function(x, direction)
{
  sum.df <- subset(x, category %in% c(TRUE, FALSE), select = c(category, p.val)) # get the TRUE and false variant calls and p.vals
  if(length(which(sum.df$category == TRUE)) > 0) # filter out cases where there aren't any TP found
  { 
    if(length(which(sum.df$category == FALSE)) == 0)
    {
      sum.df <- rbind(sum.df, data.frame(category = FALSE, p.val = 1))
    }  
    roc(sum.df$category ~ sum.df$p.val, plot = FALSE, CI = TRUE, direction = direction)
  }
}

roc_df <- function(roc_analysis) # get the coordinates and cut offs for all the points in the ROC object
{ 
  roc_analysis <- roc_analysis[unlist(lapply(roc_analysis, is.null)) == FALSE]
  all <- lapply(roc_analysis, coords, x = "all")
  all.long <- lapply(all, melt) 
  all.long <- lapply(all.long, function(x) { mutate(x, Id = c(0, head(as.numeric(rownames(x))%/%3, -1)))}) # Set an id column to group the threshold,specificity, and sensitivity together with the same number by integer division. I adjust with c(0, head ..., -1) since 3%/%3=1 but it should be grouped with the 0's.
  roc.ls <- lapply(all.long, function(x) dcast(x, Id ~ Var1)) # dcast into columns of threshold, specificity and sensitivity
  roc.df <- do.call(rbind, roc.ls) # combine into data frame
  
  roc.df$samp <- unlist(lapply(rownames(roc.df), function(x) strsplit(x, split = "[.]")[[1]][1])) # Add sample id column.
  return(roc.df)
}  

adjust.coords <- function(roc.df, sum.df, possible_tp, possible_vars) # adjust the sensitivity and specificity called by pROC  
{ 
  samp <- roc.df$samp[1]
  samp.df <- subset(sum.df, PercWT == samp)
  sense.factor <- length(which(samp.df$category == TRUE))/possible_tp
  TN.samp <- length(which(samp.df$category == FALSE))
  mutate(roc.df, adj.sensitivity = sensitivity*sense.factor, FP = (TN.samp - TN.samp*specificity), TP = adj.sensitivity*possible_tp, adj.specificity = (possible_vars - possible_tp - FP)/((possible_vars-possible_tp)))
}

roc_df.one <- function(roc_anal, thr) # get the coordinates and cut offs for all the points in the ROC object
{ 
  roc_analysis <- roc_anal[unlist(lapply(roc_anal, is.null)) == FALSE]
  all <- lapply(roc_analysis, coords, x = thr, input = "thr")
  #roc.ls<-lapply(all.long,function(x) dcast(x,Id~Var1)) # dcast into columns of threshold, specificity and sensitivity
  roc.df <- as.data.frame(do.call(rbind, all)) # combine into data frame
  roc.df$samp <- rownames(roc.df) # add sample id column
  
  roc.df <- mutate(roc.df, exp.freq = as.factor(as.numeric(samp)/100))
} 

miseq.roc <- function(sum.df, possible_tp, possible_vars, direction) # A function to run the roc calculations and adjustments
{ 
  roc.ls <- plyr::dlply(sum.df, ~Id, sum_roc, direction)
  
  roc.df <- roc_df(roc.ls)
  
  roc.df.adj <- ddply(roc.df, ~samp, adjust.coords, sum.df, possible_tp, possible_vars)
  
  roc.df.adj <- mutate(roc.df.adj, exp.freq = as.factor(as.numeric(samp)/100)) # new by ALV
  
  fills <- ddply(roc.df.adj, ~samp, fill_in_plots) # take out these lines to not fill in the plots
  roc.df.adj <- rbind(fills, roc.df.adj) # here too
  roc.df.adj <- roc.df.adj[order(roc.df.adj$adj.sensitivity),]
  roc.df.adj$exp.freq <- factor(roc.df.adj$exp.freq, levels = rev(levels(roc.df.adj$exp.freq)))
  return(roc.df.adj)
}

miseq.roc.table <- function(sum.df, cut.off, possible_tp, possible_vars, direction)
{
  roc.ls <- dlply(sum.df, ~Id, sum_roc, direction)
  roc.df <- roc_df.one(roc.ls, cut.off)
  roc.df.adj <- ddply(roc.df, ~samp, adjust.coords, sum.df, possible_tp, possible_vars)
  roc.df.adj <- roc.df.adj[order(roc.df.adj$exp.freq,decreasing = TRUE),]
  roc.table <- subset(roc.df.adj, select = c(exp.freq, adj.sensitivity, TP, adj.specificity, FP))
  return(roc.table)
}