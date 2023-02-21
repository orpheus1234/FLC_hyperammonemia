library(ggplot2)



mean_and_sd_for_incidence <- function(flc_df,hcc_df,dist_cur,hcc_ex,print_simulation_pdf,file_name){
  flc_df <- flc_cur
  hcc_df <- hcc_cur
  dist_cur <- 1000000
  hcc_ex <- HCCdata_ex
  zip_missing_pt_num <- nrow(hcc_ex)-(nrow(flc_df)+nrow(hcc_df))
  
  hcc_cur <- hcc_df[hcc_df$Distance<=dist_cur,]
  flc_cur <- flc_df[flc_df$Distance<=dist_cur,]
  #hcc proportion missing from this distance threshold
  hcc_to_add <- ceiling(zip_missing_pt_num*(1-(nrow(hcc_cur)+nrow(flc_cur))/(nrow(flc_df)+nrow(hcc_df))))
  
  
  #now have to actually figure out the age bins so we know how to distribute these patients
  FLCdata_ex = flc_df
  HCCdata_ex = hcc_df
  # create values for age, bins, histograms, and counts
  # age
  ageFLC_ex <- FLCdata_ex[, "age_years_dx"]
  ageHCC_ex <- HCCdata_ex[, "age_years_dx"]
  # bins
  bins_0_50 <- c(0, 10, 20, 30, 40, 50)
  # histograms
  agehist_FLC_ex <- hist(ageFLC_ex, breaks = bins_0_50)
  
  agehist_HCC_ex <- hist(ageHCC_ex, breaks = bins_0_50)
  bin_proportions <- agehist_HCC_ex$counts/sum(agehist_HCC_ex$counts)
  add_values <- hcc_to_add*bin_proportions
  # counts
  FLC_counts_ex = agehist_FLC_ex$counts
  HCC_counts_ex = agehist_HCC_ex$counts+add_values
  
  
  
  
  ####################################################################################################
  ##################Hierarchical model using mean and variance from each age range####################
  # FLC List
  # 0-1 density plot, so:
  empty_col <- rep(0, 101)
  df_post <- data.frame(zero_ten=empty_col, ten_twenty=empty_col, twenty_thirty=empty_col, thirty_forty=empty_col, forty_fifty=empty_col)
  
  p_list = seq(from = 0, to = 1, by = 0.01)
  # for iterations in all columns of df_post
  for (k in 1:ncol(df_post)) {
    # for iterations in all rows of df_post
    for (j in 1:length(p_list)) {
      # for each j,k compute binom density
      df_post[j,k] = dbinom(FLC_counts_ex[k], size = FLC_counts_ex[k]+HCC_counts_ex[k], prob = p_list[j])
    }
    # *NORMALIZATION* ?
    df_post[, k] = df_post[, k]/sum(df_post[, k])
  }
  
  df_post$x_axis = p_list
  Komodo_counts = list(c(14,54,118,259,538),
                       c(13,37,135,246,522),
                       c(16,41,124,232,528),
                       c(16,44,107,190,439))
  
  Komodo_ratio_by_year <- c(0.2680,0.2474,0.2298,0.2532)
  
  
  KC_means = rep(0,5)
  KC_std = rep(0,5)
  for (i in 1:5){
    cur_vec <- c()
    for (j in 1:3){
      cur_vec = c(cur_vec,round(Komodo_counts[[j]][i]/Komodo_ratio_by_year[j]))
    }
    KC_means[i] <- mean(cur_vec)
    KC_std[i] <- sd(cur_vec)
  }
  
  n = 50000
  zero_col = rep(0,n)
  break_list = seq(from = 0, to = 1400, by = 20)
  h_list = list()
  p_list = seq(from = 0, to = 1, by = 0.01)
  sum_list = list()
  df_rand <- data.frame(zero_ten=zero_col, ten_twenty=zero_col, twenty_thirty=zero_col, thirty_forty=zero_col, forty_fifty=zero_col)
  #cycle through age ranges first
  for (k in 1:5){
    #sample from our KOMODO patient count distribution
    pt_count_vec = round(rnorm(n,mean=KC_means[k],sd=KC_std[k]))
    
    for (i in 1:n){
      
      #sample from our binomial probability distribution
      p_cur = sample(p_list, 1, prob=df_post[,k])
      df_rand[i,k] = rbinom(1,pt_count_vec[i],p_cur)
    }
  }
  sums = rowSums(df_rand)
  if (print_simulation_pdf){
    h <- hist(sums,breaks=break_list)
    pdf(file = file_name)
    matplot(h$counts, type = c("b"),pch=1,col=4,xaxt = "n",xlab = "FLC incidence",
            ylab = "Simulation Counts",cex.axis=1.5,cex.lab=1.5,xlim = c(0, 60)) #plot
    axis(side=1, at=seq(0,60,10), labels = c('0','200','400','600','800','1000','1200'),cex.axis=1.5)
    axis(side=2, cex.axis=1.5)
    legend("topright", legend = c('Expected National Incidence'), col=4, pch=1) # optional legend
    
    dev.off()
  }
  
  output_list <- c(mean(sums),sd(sums))
}

#first create the simulation gaussian plot using all data irrespective of zip code
hcc_tot <- read.csv('HCC_Binary_file_pt_exclusions_any_time_in_pt_history_final.csv')
hcc_tot$age_years_dx <- round(hcc_tot$days_old_at_dx/365)

labels <- c('0-10','10-20','20-30','30-40','40-50','50-60','60-70','70+')
breaks <- seq(0,70,10)
breaks
hcc_tot$age_bin <- '70+'
for (i in 1:(length(breaks)-1)){
  hcc_tot$age_bin[(hcc_tot$age_years_dx>breaks[i]) & (hcc_tot$age_years_dx<=breaks[i+1])] <- labels[i]
}

exclusions_rows = c("Hepatitis", "Alc_liver_disease", "Tox_liver_disease", "Fibrosis", "sclerosis", "Other_LD", "Other_LC")

# create value for HCC data containing no positive exclusions
HCC_minus_ex <- c()
for (i in 1:nrow(hcc_tot)) {
  test <- all(hcc_tot[i, exclusions_rows] == 0)
  if (test == TRUE) { 
    HCC_minus_ex <- c(HCC_minus_ex, i)
  }
}
# dataset of HCC excluding all patients with one of the ICD10 comorbidities
HCCdata_ex <- hcc_tot[HCC_minus_ex, ]
HCCdata_ex <- HCCdata_ex[is.na(HCCdata_ex$age_years_dx)==FALSE,]
HCCdata_ex <- HCCdata_ex[HCCdata_ex$age_years_dx<=50,]






# counts
FLC_counts_all = agehist_FLC_all$counts
FLC_counts_ex = agehist_FLC_ex$counts
HCC_counts_all = agehist_HCC_all$counts
HCC_counts_ex = agehist_HCC_ex$counts

dist_cur <- 100

#splitting patients by FLC Dx
flc_cur <- HCCdata_ex[HCCdata_ex$FLC_dx==1,]
hcc_cur <- HCCdata_ex[HCCdata_ex$FLC_dx==0,]

#create the stacked bar graph by Distance

# bins
dist_thresh <- c(0,10,50,100,500,1000,10000)
distconvert_FLC <- rep(0,length(distFLC_ex))
distconvert_HCC <- rep(0,length(distHCC_ex))


dist_thresh <- c(0,10,50,100,500,1000,10000)
hcc_df <- read.csv(file='hcc_pt_with_Distance.csv')
flc_df <- read.csv(file='flc_pt_with_Distance.csv')

distFLC_ex <- flc_df[, "Distance"]
distHCC_ex <- hcc_df[, "Distance"]
# histograms
disthist_FLC_ex <- hist(distFLC_ex, breaks = dist_thresh)
disthist_HCC_ex <- hist(distHCC_ex, breaks = dist_thresh)

dist_thresh <- c(10000,2000,seq(1000,10,-50),40,30,20,10)
df_stats <- data.frame()
cur_vec <- mean_and_sd_for_incidence(flc_cur,hcc_cur,100000,HCCdata_ex,TRUE,'all_data_national_incidence_graph.pdf')
df_stats[1,'mean'] <- cur_vec[1]
df_stats[1,'sd'] <- cur_vec[2]




for (d in 2:length(dist_thresh)){
  hcc_cur <- hcc_df[hcc_df$Distance<=dist_thresh[d],]
  flc_cur <- flc_df[flc_df$Distance<=dist_thresh[d],]
  cur_vec <- mean_and_sd_for_incidence(flc_cur,hcc_cur,dist_thresh[d],HCCdata_ex,FALSE,'')
  df_stats[d,'mean'] <- cur_vec[1]
  df_stats[d,'sd'] <- cur_vec[2]
}
df_stats$loc <- log10(dist_thresh)
df_stats$lb <- df_stats$mean-df_stats$sd
df_stats$ub <- df_stats$mean+df_stats$sd
# Create the plot
pdf("incidence_plot_fixed_aspect_ratio.pdf", paper = "special", width = 12, height = 12)
ggplot(data = df_stats, aes(x = loc, y = mean)) +
  geom_line(color='blue',size=1.5) + ylim(0,850)+
  geom_linerange(aes(ymin = lb, ymax = ub),color='lightblue', size = 1)+
  scale_x_continuous(breaks=c(2,3,4),labels=c('100','1000','all'))+
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_line(color = "black"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.line.y = element_line(color = "black", size = 1))+
  xlab('Distance threshold from SF (in miles) \n for patient inclusion')+
  ylab('Mean national incidence estimate')
dev.off()