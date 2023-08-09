library(readr)
library(ggplot2)
require(dplyr)
require(tidyverse)
require(hesim)

#### START ####
## Markov model 1 ##
params_M1<-read_csv("./Markov1.csv")
trans_M1<-as.numeric(params_M1$value[params_M1$var_g1=="prob_in"])
util_vec<-as.numeric(params_M1$value[params_M1$var_g1=="util_in"])
cost_vec_1<-as.numeric(params_M1$value[params_M1$var_g1=="cost_in"])
tree_vec<-as.numeric(params_M1$value[params_M1$var_g1=="tree_in"])

params_M2<-read_csv("./Markov2.csv")
trans_M2<-as.numeric(params_M2$value[params_M2$var_g1=="prob_in"])
params_M3<-read_csv("./Markov3.csv")
trans_M3<-as.numeric(params_M3$value[params_M3$var_g1=="prob_in"])
cost_vec_2<-as.numeric(params_M3$value[params_M3$var_g1=="cost_in"])
params_M4<-read_csv("./Markov4.csv")
trans_M4<-as.numeric(params_M4$value[params_M4$var_g1=="prob_in"])

results<-matrix(data=NA, nrow=1, ncol = 2)
results<-as.data.frame(results)
colnames(results)<-c("cost","qaly")

fd<-c(rep(0.007,3),rep(0.01,5),rep(0.015,5),rep(0.023,5),rep(0.036,5),rep(0.057,5),rep(0.12,50))
init_fd<-0.007
d_cost<-c(as.numeric(params_M1$value[params_M1$var_name=="cost"]))
d_qaly<-c(as.numeric(params_M1$value[params_M1$var_name=="utility"]))

cycle_index<-c(0:70)


Markov_func<-function(trans_vec,cost_vec){
  init_ff<-trans_vec[5]
  init_fp<-trans_vec[6]
  init_pd<-trans_vec[7]
  init_pp<-trans_vec[9]
  ff<-(1-fd)*(init_ff/(init_ff+init_fp))
  fp<-(1-fd)*(init_fp/(init_ff+init_fp))
  pd<-1-((1-init_pd)/(1-init_fd))*(1-fd)
  pp<-1-pd-trans_vec[8]
  Tmat<-matrix(data=0, nrow=4, ncol=4)
  Tmat[1,]<-c(0,trans_vec[1],trans_vec[2],init_fd)
  Tmat[4,]<-c(0,0,0,1)
  
  state<-matrix(data=0,nrow = 70,ncol = 4)
  colnames(state)<-c("Firstline therapy","Progression free","Progression","Death")
  state[1,]<-c(1,0,0,0)
  
  for (i in 2:70) {
    Tmat[2,]<-c(0,ff[i-1],fp[i-1],fd[i-1])
    Tmat[3,]<-c(0,trans_vec[8],pp[i-1],pd[i-1])
    
    state[i,]<-state[(i-1),]%*%Tmat
  }
  
  ly<-as.data.frame(state)
  ly$life_years<-rowSums(ly[,1:3])
  qaly<-ly
  qaly[,1]<-ly[,1]*util_vec[1]
  qaly[,2]<-ly[,2]*util_vec[2]
  qaly[,3]<-ly[,3]*util_vec[3]
  qaly[,4]<-ly[,4]*util_vec[4]
  qaly[,5]<-rowSums(qaly[,1:4])
  colnames(qaly)[5]<-c("QALY_sum")
  
  cost<-ly
  cost[,1]<-ly[,1]*cost_vec[3]
  cost[2,2]<-ly[2,2]*cost_vec[4]
  cost[3,2]<-(ly[2,2]*ff[2]+ly[2,3]*trans_vec[8])*cost_vec[4]
  cost[2,3]<-ly[2,3]*cost_vec[6]
  cost[3,3]<-(ly[2,3]*pp[2]+ly[2,2]*fp[2])*cost_vec[6]
  for (i in 4:70) {
    f_new<-ly[i-2,3]*trans_vec[8]*ff[i-1]+ly[i-1,3]*trans_vec[8]
    p_new<-ly[i-2,2]*fp[i-2]*pp[i-1]+ly[i-1,2]*fp[i-1]
    cost[i,2]<-f_new*cost_vec[4] + (ly[i,2]-f_new)*cost_vec[5]
    cost[i,3]<-p_new*cost_vec[6] + (ly[i,3]-p_new)*cost_vec[7]
  }
  cost[,4]<-cost[,4]*0
  cost[,5]<-rowSums(cost[,1:4])
  colnames(cost)[5]<-c("cost_sum")
  
  ly_d<-ly
  ly_d[,1]<-ly[,1]*(1/(1+d_qaly)^cycle_index[1:70])
  ly_d[,2]<-ly[,2]*(1/(1+d_qaly)^cycle_index[1:70])
  ly_d[,3]<-ly[,3]*(1/(1+d_qaly)^cycle_index[1:70])
  ly_d[,5]<-rowSums(ly_d[,1:3])
  qaly_d<-qaly
  qaly_d[,1]<-qaly[,1]*(1/(1+d_qaly)^cycle_index[1:70])
  qaly_d[,2]<-qaly[,2]*(1/(1+d_qaly)^cycle_index[1:70])
  qaly_d[,3]<-qaly[,3]*(1/(1+d_qaly)^cycle_index[1:70])
  qaly_d[,5]<-rowSums(qaly_d[,1:3])
  cost_d<-cost
  cost_d[,1]<-cost[,1]*(1/(1+d_cost)^cycle_index[1:70])
  cost_d[,2]<-cost[,2]*(1/(1+d_cost)^cycle_index[1:70])
  cost_d[,3]<-cost[,3]*(1/(1+d_cost)^cycle_index[1:70])
  cost_d[,5]<-rowSums(cost_d[,1:3])
  
  ly_co<-ly_d
  ly_co[2:69,1]<-(ly_d[2:69,1]+ly_d[3:70,1])/2
  ly_co[2:69,2]<-(ly_d[2:69,2]+ly_d[3:70,2])/2
  ly_co[2:69,3]<-(ly_d[2:69,3]+ly_d[3:70,3])/2
  ly_co[,5]<-rowSums(ly_co[,1:3])
  qaly_co<-qaly_d
  qaly_co[2:69,1]<-(qaly_d[2:69,1]+qaly_d[3:70,1])/2
  qaly_co[2:69,2]<-(qaly_d[2:69,2]+qaly_d[3:70,2])/2
  qaly_co[2:69,3]<-(qaly_d[2:69,3]+qaly_d[3:70,3])/2
  qaly_co[,5]<-rowSums(qaly_co[,1:3])
  cost_co<-cost_d
  cost_co[2:69,1]<-(cost_d[2:69,1]+cost_d[3:70,1])/2
  cost_co[2:69,2]<-(cost_d[2:69,2]+cost_d[3:70,2])/2
  cost_co[2:69,3]<-(cost_d[2:69,3]+cost_d[3:70,3])/2
  cost_co[,5]<-rowSums(cost_co[,1:3])
  
  results$cost<-sum(cost_co[,5])
  results$qaly<-sum(qaly_co[,5])
  return(results)
}

Markov_func(trans_M4,cost_vec)


## ICER calculation ##
results_M1<-Markov_func(trans_M1,cost_vec_1)
results_M2<-Markov_func(trans_M2,cost_vec_1)
results_M3<-Markov_func(trans_M3,cost_vec_2)
results_M4<-Markov_func(trans_M4,cost_vec_2)

icer_func <- function(tree_vec,cost_vec,results_M1,results_M2,results_M3,results_M4){
  qaly_pet<-(1-tree_vec[2])*tree_vec[3]*results_M1$qaly + tree_vec[2]*(1-tree_vec[4])*results_M2$qaly +
    tree_vec[2]*tree_vec[4]*results_M3$qaly + (1-tree_vec[2])*(1-tree_vec[3])*results_M4$qaly
  qaly_ct<-(1-tree_vec[2])*tree_vec[5]*results_M1$qaly + tree_vec[2]*(1-tree_vec[6])*results_M2$qaly +
    tree_vec[2]*tree_vec[6]*results_M3$qaly + (1-tree_vec[2])*(1-tree_vec[5])*results_M4$qaly
  
  cost_pet<-(1-tree_vec[2])*tree_vec[3]*results_M1$cost + tree_vec[2]*(1-tree_vec[4])*results_M2$cost +
    tree_vec[2]*tree_vec[4]*results_M3$cost + (1-tree_vec[2])*(1-tree_vec[3])*results_M4$cost + cost_vec[1]*1
  cost_ct<-(1-tree_vec[2])*tree_vec[5]*results_M1$cost + tree_vec[2]*(1-tree_vec[6])*results_M2$cost +
    tree_vec[2]*tree_vec[6]*results_M3$cost + (1-tree_vec[2])*(1-tree_vec[5])*results_M4$cost + cost_vec[2]*1
  
  i_qaly<-qaly_pet-qaly_ct
  i_cost<-cost_pet-cost_ct
  icer<-i_cost/i_qaly
  
  return(c(qaly_pet,qaly_ct,cost_pet,cost_ct,i_qaly,i_cost,icer))
}

icer_base<-icer_func(tree_vec,cost_vec_1,results_M1,results_M2,results_M3,results_M4)
icer_base

## One way sensitivity analysis ##
results_dsa <-matrix(data=NA,nrow=24,ncol=3)      
results_dsa <- as.data.frame(results_dsa)      
colnames(results_dsa) <- c("param","icer_low","icer_high")      
results_dsa$param <- c('d_cost','d_effect','s_late','s_pet1','s_pet2','s_ct1','s_ct2','p_tf','p_tp','p_ff','p_fp',
                       'p_pd','p_pf','p_pp','c_pet','c_ct','c_t','c_f1','c_f2',
                       'c_p1','c_p2','u_t','u_f','u_p')
params_M1_o<-params_M1
params_M2_o<-params_M2
params_M3_o<-params_M3
params_M4_o<-params_M4

dsa_M1<-read_csv("./M1_DSA.csv")
dsa_M2<-read_csv("./M2_DSA.csv")
dsa_M3<-read_csv("./M3_DSA.csv")
dsa_M4<-read_csv("./M4_DSA.csv")



for (i in 1:24) {
  var_index <- paste("x",i,sep="")
  params_M1$value[params_M1$sens==var_index] <- dsa_M1[1,1+i]
  params_M2$value[params_M2$sens==var_index] <- dsa_M2[1,1+i]
  params_M3$value[params_M3$sens==var_index] <- dsa_M3[1,1+i]
  params_M4$value[params_M4$sens==var_index] <- dsa_M4[1,1+i]
  
  d_cost<-c(as.numeric(params_M1$value[params_M1$var_name=="cost"]))
  d_qaly<-c(as.numeric(params_M1$value[params_M1$var_name=="utility"]))
  tree_vec<-as.numeric(params_M1$value[params_M1$var_g1=="tree_in"])
  util_vec<-as.numeric(params_M1$value[params_M1$var_g1=="util_in"])
  cost_vec_1<-as.numeric(params_M1$value[params_M1$var_g1=="cost_in"])
  cost_vec_2<-as.numeric(params_M3$value[params_M3$var_g1=="cost_in"])
  trans_M1<-as.numeric(params_M1$value[params_M1$var_g1=="prob_in"])
  trans_M2<-as.numeric(params_M2$value[params_M2$var_g1=="prob_in"])
  trans_M3<-as.numeric(params_M3$value[params_M3$var_g1=="prob_in"])
  trans_M4<-as.numeric(params_M4$value[params_M4$var_g1=="prob_in"])
  
  results_dsa_M1 <- Markov_func(trans_M1,cost_vec_1)
  results_dsa_M2 <- Markov_func(trans_M2,cost_vec_1)
  results_dsa_M3 <- Markov_func(trans_M3,cost_vec_2)
  results_dsa_M4 <- Markov_func(trans_M4,cost_vec_2)
  
  results_dsa$icer_low[i]<-icer_func(tree_vec,cost_vec_1,results_dsa_M1,results_dsa_M2,results_dsa_M3,results_dsa_M4)[7]
  
  # Upper value: second row in the input_dsa 
  params_M1$value[params_M1$sens==var_index] <- dsa_M1[2,1+i]
  params_M2$value[params_M2$sens==var_index] <- dsa_M2[2,1+i]
  params_M3$value[params_M3$sens==var_index] <- dsa_M3[2,1+i]
  params_M4$value[params_M4$sens==var_index] <- dsa_M4[2,1+i]
  
  d_cost<-c(as.numeric(params_M1$value[params_M1$var_name=="cost"]))
  d_qaly<-c(as.numeric(params_M1$value[params_M1$var_name=="utility"]))
  tree_vec<-as.numeric(params_M1$value[params_M1$var_g1=="tree_in"])
  util_vec<-as.numeric(params_M1$value[params_M1$var_g1=="util_in"])
  cost_vec_1<-as.numeric(params_M1$value[params_M1$var_g1=="cost_in"])
  cost_vec_2<-as.numeric(params_M3$value[params_M3$var_g1=="cost_in"])
  trans_M1<-as.numeric(params_M1$value[params_M1$var_g1=="prob_in"])
  trans_M2<-as.numeric(params_M2$value[params_M2$var_g1=="prob_in"])
  trans_M3<-as.numeric(params_M3$value[params_M3$var_g1=="prob_in"])
  trans_M4<-as.numeric(params_M4$value[params_M4$var_g1=="prob_in"])
  
  results_dsa_M1 <- Markov_func(trans_M1,cost_vec_1)
  results_dsa_M2 <- Markov_func(trans_M2,cost_vec_1)
  results_dsa_M3 <- Markov_func(trans_M3,cost_vec_2)
  results_dsa_M4 <- Markov_func(trans_M4,cost_vec_2)
  
  results_dsa$icer_high[i]<-icer_func(tree_vec,cost_vec_1,results_dsa_M1,results_dsa_M2,results_dsa_M3,results_dsa_M4)[7]
  

  params_M1 <- params_M1_o
  params_M2 <- params_M2_o
  params_M3 <- params_M3_o
  params_M4 <- params_M4_o
  
}


base_value<-51676
results_dsa$diff <- abs(results_dsa$icer_high-results_dsa$icer_low)
df <- results_dsa
df$param<-c("Discount rate of cost","Discount rate of utility","Porportion of late-stage patients",
            "Probability early-stage patients staged correctly by PET/CT",
            "Probability late-stage patients staged correctly by PET/CT",
            "Probability early-stage patients staged correctly by CT",
            "Probability late-stage patients staged correctly by CT",
            "TP from firstline therapy to progression-free","TP from firstline therapy to progression",
            "TP from progression-free to progression-free","TP from progression-free to progression",
            "TP from progression to death","TP from progression to progression-free",
            "TP from progression to progression","Cost of PET/CT","Cost of CT","Cost of the firstline therapy",
            "Cost of surveillance for progression-free in the first 2 years",
            "Cost of surveillance for progression-free in the years after",
            "Cost of progression in the first 2 years","Cost of progression in the years after",
            "Utility values for patients receiving the firstline therapy","Utility values for progression-free patients",
            "Utility values for patients with progression")

order.df <- df %>% arrange(diff) %>%
  mutate(param=factor(x=param, levels=param)) %>%
  select(param) %>% unlist() %>% levels()

# width of columns in plot (value between 0 and 1)
width <- 0.95

# get data frame in shape for ggplot and geom_rect
df.2 <- df  %>% 
  # gather columns Lower_Bound and Upper_Bound into a single column using gather
  gather(key='type', value='output.value',icer_low:icer_high) %>%
  # just reordering columns
  select(param, type, output.value, diff) %>%
  # create the columns for geom_rect
  mutate(param=factor(param, levels=order.df),
         ymin=pmin(output.value, base_value),
         ymax=pmax(output.value, base_value),
         xmin=as.numeric(param)-width/2,
         xmax=as.numeric(param)+width/2)

# create plot
# (use scale_x_continuous to change labels in y axis to name of parameters)

ggplot() + 
  geom_rect(data = df.2, 
            aes(ymax=ymax, ymin=ymin, xmax=xmax, xmin=xmin, fill=type)) +
  theme_bw() + 
  theme(axis.title.y=element_blank(), legend.position = 'bottom',
        legend.title = element_blank()) + 
  geom_hline(yintercept = base_value) +
  geom_hline(aes(yintercept = 160000),linetype="dashed",linewidth=0.4)+
  scale_x_continuous(breaks = c(1:length(order.df)), 
                     labels = order.df) +
  scale_y_continuous(breaks = c(0,250000,500000,750000,5000000,1250000))+
  coord_flip()+
  ylab("ICER (Thai Bhat per QALY)")+
  scale_fill_discrete(breaks=c("icer_low", "icer_high"),
                      labels=c("Lower bound value", "Upper bound value"))




#### Probabilistic Sensitivity Analysis ####
psa_param_1<-read_csv("./PSA_param_1.csv")
psa_param_2<-read_csv("./PSA_param_2.csv")

## Staging Accuracy & Cost & Utility ##
input_psa_gen_1 <-as.data.frame(matrix(data=NA,nrow=5000,ncol=16))      
colnames(input_psa_gen_1) <- c("sample","s_late","s_pet1","s_pet2","s_ct1","s_ct2","c_pet","c_ct","c_t",
                             "c_f1","c_f2","c_p1","c_p2","u_t","u_f","u_p")   
input_psa_gen_1$sample<-seq(1,5000,1)
input_psa_gen_1$s_late <- rbeta(n = 5000, shape1 = psa_param_1$alpha[1], shape2 = psa_param_1$beta[1])
input_psa_gen_1$s_pet1 <- rbeta(n = 5000, shape1 = psa_param_1$alpha[2], shape2 = psa_param_1$beta[2])
input_psa_gen_1$s_pet2 <- rbeta(n = 5000, shape1 = psa_param_1$alpha[3], shape2 = psa_param_1$beta[3])
input_psa_gen_1$s_ct1 <- rbeta(n = 5000, shape1 = psa_param_1$alpha[4], shape2 = psa_param_1$beta[4])
input_psa_gen_1$s_ct2 <- rbeta(n = 5000, shape1 = psa_param_1$alpha[5], shape2 = psa_param_1$beta[5])
input_psa_gen_1$c_pet <- rgamma(n = 5000, shape = psa_param_1$alpha[6], scale = psa_param_1$beta[6])
input_psa_gen_1$c_ct <- rgamma(n = 5000, shape = psa_param_1$alpha[7], scale = psa_param_1$beta[7])
input_psa_gen_1$c_t <- rgamma(n = 5000, shape = psa_param_1$alpha[8], scale = psa_param_1$beta[8])
input_psa_gen_1$c_f1 <- rgamma(n = 5000, shape = psa_param_1$alpha[9], scale = psa_param_1$beta[9])
input_psa_gen_1$c_f2 <- rgamma(n = 5000, shape = psa_param_1$alpha[10], scale = psa_param_1$beta[10])
input_psa_gen_1$c_p1 <- rgamma(n = 5000, shape = psa_param_1$alpha[11], scale = psa_param_1$beta[11])
input_psa_gen_1$c_p2 <- rgamma(n = 5000, shape = psa_param_1$alpha[12], scale = psa_param_1$beta[12])
x_u<-runif(5000,0,1)
input_psa_gen_1$u_t <- qbeta(x_u, shape1 = psa_param_1$alpha[13], shape2 = psa_param_1$beta[13])
input_psa_gen_1$u_f <- qbeta(x_u, shape1 = psa_param_1$alpha[14], shape2 = psa_param_1$beta[14])
input_psa_gen_1$u_p <- qbeta(x_u, shape1 = psa_param_1$alpha[15], shape2 = psa_param_1$beta[15])

input_psa_gen_2<-input_psa_gen_1
input_psa_gen_2$c_t <- rgamma(n = 5000, shape = psa_param_2$alpha[8], scale = psa_param_2$beta[8])


## Transitional Probability ##
# Markov 1 #
psa_trans_M1<-read_csv("./PSA_trans_M1.csv")
input_psa_M1 <-as.data.frame(matrix(data=NA,nrow=5000,ncol=46))      
colnames(input_psa_M1) <- c("sample","p_tf","p_tp","p_td","p_ff_a","p_fp_a","p_fd_a","p_pf_a","p_pp_a","p_pd_a",
                            "p_ff_b","p_fp_b","p_fd_b","p_pf_b","p_pp_b","p_pd_b",
                            "p_ff_c","p_fp_c","p_fd_c","p_pf_c","p_pp_c","p_pd_c",
                            "p_ff_d","p_fp_d","p_fd_d","p_pf_d","p_pp_d","p_pd_d",
                            "p_ff_e","p_fp_e","p_fd_e","p_pf_e","p_pp_e","p_pd_e",
                            "p_ff_f","p_fp_f","p_fd_f","p_pf_f","p_pp_f","p_pd_f",
                            "p_ff_g","p_fp_g","p_fd_g","p_pf_g","p_pp_g","p_pd_g")   
input_psa_M1$sample<-seq(1,5000,1)
alpha_M1_a<-matrix(as.numeric(psa_trans_M1[1,2:10]),ncol = 3,byrow = TRUE)
sample_M1_a<-hesim::rdirichlet_mat(n = 5000,alpha_M1_a+1)
input_psa_M1$p_tf<-sample_M1_a[1,1,1:5000]
input_psa_M1$p_tp<-sample_M1_a[1,2,1:5000]
input_psa_M1$p_td<-sample_M1_a[1,3,1:5000]
input_psa_M1$p_ff_a<-sample_M1_a[2,1,1:5000]
input_psa_M1$p_fp_a<-sample_M1_a[2,2,1:5000]
input_psa_M1$p_fd_a<-sample_M1_a[2,3,1:5000]
input_psa_M1$p_pf_a<-sample_M1_a[3,1,1:5000]
input_psa_M1$p_pp_a<-sample_M1_a[3,2,1:5000]
input_psa_M1$p_pd_a<-sample_M1_a[3,3,1:5000]
alpha_M1_b<-matrix(as.numeric(psa_trans_M1[2,5:10]),ncol = 3,byrow = TRUE)
sample_M1_b<-hesim::rdirichlet_mat(n = 5000,alpha_M1_b+1)
input_psa_M1$p_ff_b<-sample_M1_b[1,1,1:5000]
input_psa_M1$p_fp_b<-sample_M1_b[1,2,1:5000]
input_psa_M1$p_fd_b<-sample_M1_b[1,3,1:5000]
input_psa_M1$p_pf_b<-sample_M1_b[2,1,1:5000]
input_psa_M1$p_pp_b<-sample_M1_b[2,2,1:5000]
input_psa_M1$p_pd_b<-sample_M1_b[2,3,1:5000]
alpha_M1_c<-matrix(as.numeric(psa_trans_M1[3,5:10]),ncol = 3,byrow = TRUE)
sample_M1_c<-hesim::rdirichlet_mat(n = 5000,alpha_M1_c+1)
input_psa_M1$p_ff_c<-sample_M1_c[1,1,1:5000]
input_psa_M1$p_fp_c<-sample_M1_c[1,2,1:5000]
input_psa_M1$p_fd_c<-sample_M1_c[1,3,1:5000]
input_psa_M1$p_pf_c<-sample_M1_c[2,1,1:5000]
input_psa_M1$p_pp_c<-sample_M1_c[2,2,1:5000]
input_psa_M1$p_pd_c<-sample_M1_c[2,3,1:5000]
alpha_M1_d<-matrix(as.numeric(psa_trans_M1[4,5:10]),ncol = 3,byrow = TRUE)
sample_M1_d<-hesim::rdirichlet_mat(n = 5000,alpha_M1_d+1)
input_psa_M1$p_ff_d<-sample_M1_d[1,1,1:5000]
input_psa_M1$p_fp_d<-sample_M1_d[1,2,1:5000]
input_psa_M1$p_fd_d<-sample_M1_d[1,3,1:5000]
input_psa_M1$p_pf_d<-sample_M1_d[2,1,1:5000]
input_psa_M1$p_pp_d<-sample_M1_d[2,2,1:5000]
input_psa_M1$p_pd_d<-sample_M1_d[2,3,1:5000]
alpha_M1_e<-matrix(as.numeric(psa_trans_M1[5,5:10]),ncol = 3,byrow = TRUE)
sample_M1_e<-hesim::rdirichlet_mat(n = 5000,alpha_M1_e+1)
input_psa_M1$p_ff_e<-sample_M1_e[1,1,1:5000]
input_psa_M1$p_fp_e<-sample_M1_e[1,2,1:5000]
input_psa_M1$p_fd_e<-sample_M1_e[1,3,1:5000]
input_psa_M1$p_pf_e<-sample_M1_e[2,1,1:5000]
input_psa_M1$p_pp_e<-sample_M1_e[2,2,1:5000]
input_psa_M1$p_pd_e<-sample_M1_e[2,3,1:5000]
alpha_M1_f<-matrix(as.numeric(psa_trans_M1[6,5:10]),ncol = 3,byrow = TRUE)
sample_M1_f<-hesim::rdirichlet_mat(n = 5000,alpha_M1_f+1)
input_psa_M1$p_ff_f<-sample_M1_f[1,1,1:5000]
input_psa_M1$p_fp_f<-sample_M1_f[1,2,1:5000]
input_psa_M1$p_fd_f<-sample_M1_f[1,3,1:5000]
input_psa_M1$p_pf_f<-sample_M1_f[2,1,1:5000]
input_psa_M1$p_pp_f<-sample_M1_f[2,2,1:5000]
input_psa_M1$p_pd_f<-sample_M1_f[2,3,1:5000]
alpha_M1_g<-matrix(as.numeric(psa_trans_M1[7,5:10]),ncol = 3,byrow = TRUE)
sample_M1_g<-hesim::rdirichlet_mat(n = 5000,alpha_M1_g+1)
input_psa_M1$p_ff_g<-sample_M1_g[1,1,1:5000]
input_psa_M1$p_fp_g<-sample_M1_g[1,2,1:5000]
input_psa_M1$p_fd_g<-sample_M1_g[1,3,1:5000]
input_psa_M1$p_pf_g<-sample_M1_g[2,1,1:5000]
input_psa_M1$p_pp_g<-sample_M1_g[2,2,1:5000]
input_psa_M1$p_pd_g<-sample_M1_g[2,3,1:5000]

# Markov 2 #
psa_trans_M2<-read_csv("./PSA_trans_M2.csv")
input_psa_M2 <-input_psa_M1     
alpha_M2_a<-matrix(as.numeric(psa_trans_M2[1,2:10]),ncol = 3,byrow = TRUE)
sample_M2_a<-hesim::rdirichlet_mat(n = 5000,alpha_M2_a+1)
input_psa_M2$p_tf<-sample_M2_a[1,1,1:5000]
input_psa_M2$p_tp<-sample_M2_a[1,2,1:5000]
input_psa_M2$p_td<-sample_M2_a[1,3,1:5000]
input_psa_M2$p_ff_a<-sample_M2_a[2,1,1:5000]
input_psa_M2$p_fp_a<-sample_M2_a[2,2,1:5000]
input_psa_M2$p_fd_a<-sample_M2_a[2,3,1:5000]
input_psa_M2$p_pf_a<-sample_M2_a[3,1,1:5000]
input_psa_M2$p_pp_a<-sample_M2_a[3,2,1:5000]
input_psa_M2$p_pd_a<-sample_M2_a[3,3,1:5000]
alpha_M2_b<-matrix(as.numeric(psa_trans_M2[2,5:10]),ncol = 3,byrow = TRUE)
sample_M2_b<-hesim::rdirichlet_mat(n = 5000,alpha_M2_b+1)
input_psa_M2$p_ff_b<-sample_M2_b[1,1,1:5000]
input_psa_M2$p_fp_b<-sample_M2_b[1,2,1:5000]
input_psa_M2$p_fd_b<-sample_M2_b[1,3,1:5000]
input_psa_M2$p_pf_b<-sample_M2_b[2,1,1:5000]
input_psa_M2$p_pp_b<-sample_M2_b[2,2,1:5000]
input_psa_M2$p_pd_b<-sample_M2_b[2,3,1:5000]
alpha_M2_c<-matrix(as.numeric(psa_trans_M2[3,5:10]),ncol = 3,byrow = TRUE)
sample_M2_c<-hesim::rdirichlet_mat(n = 5000,alpha_M2_c+1)
input_psa_M2$p_ff_c<-sample_M2_c[1,1,1:5000]
input_psa_M2$p_fp_c<-sample_M2_c[1,2,1:5000]
input_psa_M2$p_fd_c<-sample_M2_c[1,3,1:5000]
input_psa_M2$p_pf_c<-sample_M2_c[2,1,1:5000]
input_psa_M2$p_pp_c<-sample_M2_c[2,2,1:5000]
input_psa_M2$p_pd_c<-sample_M2_c[2,3,1:5000]
alpha_M2_d<-matrix(as.numeric(psa_trans_M2[4,5:10]),ncol = 3,byrow = TRUE)
sample_M2_d<-hesim::rdirichlet_mat(n = 5000,alpha_M2_d+1)
input_psa_M2$p_ff_d<-sample_M2_d[1,1,1:5000]
input_psa_M2$p_fp_d<-sample_M2_d[1,2,1:5000]
input_psa_M2$p_fd_d<-sample_M2_d[1,3,1:5000]
input_psa_M2$p_pf_d<-sample_M2_d[2,1,1:5000]
input_psa_M2$p_pp_d<-sample_M2_d[2,2,1:5000]
input_psa_M2$p_pd_d<-sample_M2_d[2,3,1:5000]
alpha_M2_e<-matrix(as.numeric(psa_trans_M2[5,5:10]),ncol = 3,byrow = TRUE)
sample_M2_e<-hesim::rdirichlet_mat(n = 5000,alpha_M2_e+1)
input_psa_M2$p_ff_e<-sample_M2_e[1,1,1:5000]
input_psa_M2$p_fp_e<-sample_M2_e[1,2,1:5000]
input_psa_M2$p_fd_e<-sample_M2_e[1,3,1:5000]
input_psa_M2$p_pf_e<-sample_M2_e[2,1,1:5000]
input_psa_M2$p_pp_e<-sample_M2_e[2,2,1:5000]
input_psa_M2$p_pd_e<-sample_M2_e[2,3,1:5000]
alpha_M2_f<-matrix(as.numeric(psa_trans_M2[6,5:10]),ncol = 3,byrow = TRUE)
sample_M2_f<-hesim::rdirichlet_mat(n = 5000,alpha_M2_f+1)
input_psa_M2$p_ff_f<-sample_M2_f[1,1,1:5000]
input_psa_M2$p_fp_f<-sample_M2_f[1,2,1:5000]
input_psa_M2$p_fd_f<-sample_M2_f[1,3,1:5000]
input_psa_M2$p_pf_f<-sample_M2_f[2,1,1:5000]
input_psa_M2$p_pp_f<-sample_M2_f[2,2,1:5000]
input_psa_M2$p_pd_f<-sample_M2_f[2,3,1:5000]
alpha_M2_g<-matrix(as.numeric(psa_trans_M2[7,5:10]),ncol = 3,byrow = TRUE)
sample_M2_g<-hesim::rdirichlet_mat(n = 5000,alpha_M2_g+1)
input_psa_M2$p_ff_g<-sample_M2_g[1,1,1:5000]
input_psa_M2$p_fp_g<-sample_M2_g[1,2,1:5000]
input_psa_M2$p_fd_g<-sample_M2_g[1,3,1:5000]
input_psa_M2$p_pf_g<-sample_M2_g[2,1,1:5000]
input_psa_M2$p_pp_g<-sample_M2_g[2,2,1:5000]
input_psa_M2$p_pd_g<-sample_M2_g[2,3,1:5000]

mean(input_psa_M2$p_pd_a)

# Markov 3 #
psa_trans_M3<-read_csv("./PSA_trans_M3.csv")
input_psa_M3 <-input_psa_M1     
alpha_M3_a<-matrix(as.numeric(psa_trans_M3[1,2:10]),ncol = 3,byrow = TRUE)
sample_M3_a<-hesim::rdirichlet_mat(n = 5000,alpha_M3_a+1)
input_psa_M3$p_tf<-sample_M3_a[1,1,1:5000]
input_psa_M3$p_tp<-sample_M3_a[1,2,1:5000]
input_psa_M3$p_td<-sample_M3_a[1,3,1:5000]
input_psa_M3$p_ff_a<-sample_M3_a[2,1,1:5000]
input_psa_M3$p_fp_a<-sample_M3_a[2,2,1:5000]
input_psa_M3$p_fd_a<-sample_M3_a[2,3,1:5000]
input_psa_M3$p_pf_a<-sample_M3_a[3,1,1:5000]
input_psa_M3$p_pp_a<-sample_M3_a[3,2,1:5000]
input_psa_M3$p_pd_a<-sample_M3_a[3,3,1:5000]
alpha_M3_b<-matrix(as.numeric(psa_trans_M3[2,5:10]),ncol = 3,byrow = TRUE)
sample_M3_b<-hesim::rdirichlet_mat(n = 5000,alpha_M3_b+1)
input_psa_M3$p_ff_b<-sample_M3_b[1,1,1:5000]
input_psa_M3$p_fp_b<-sample_M3_b[1,2,1:5000]
input_psa_M3$p_fd_b<-sample_M3_b[1,3,1:5000]
input_psa_M3$p_pf_b<-sample_M3_b[2,1,1:5000]
input_psa_M3$p_pp_b<-sample_M3_b[2,2,1:5000]
input_psa_M3$p_pd_b<-sample_M3_b[2,3,1:5000]
alpha_M3_c<-matrix(as.numeric(psa_trans_M3[3,5:10]),ncol = 3,byrow = TRUE)
sample_M3_c<-hesim::rdirichlet_mat(n = 5000,alpha_M3_c+1)
input_psa_M3$p_ff_c<-sample_M3_c[1,1,1:5000]
input_psa_M3$p_fp_c<-sample_M3_c[1,2,1:5000]
input_psa_M3$p_fd_c<-sample_M3_c[1,3,1:5000]
input_psa_M3$p_pf_c<-sample_M3_c[2,1,1:5000]
input_psa_M3$p_pp_c<-sample_M3_c[2,2,1:5000]
input_psa_M3$p_pd_c<-sample_M3_c[2,3,1:5000]
alpha_M3_d<-matrix(as.numeric(psa_trans_M3[4,5:10]),ncol = 3,byrow = TRUE)
sample_M3_d<-hesim::rdirichlet_mat(n = 5000,alpha_M3_d+1)
input_psa_M3$p_ff_d<-sample_M3_d[1,1,1:5000]
input_psa_M3$p_fp_d<-sample_M3_d[1,2,1:5000]
input_psa_M3$p_fd_d<-sample_M3_d[1,3,1:5000]
input_psa_M3$p_pf_d<-sample_M3_d[2,1,1:5000]
input_psa_M3$p_pp_d<-sample_M3_d[2,2,1:5000]
input_psa_M3$p_pd_d<-sample_M3_d[2,3,1:5000]
alpha_M3_e<-matrix(as.numeric(psa_trans_M3[5,5:10]),ncol = 3,byrow = TRUE)
sample_M3_e<-hesim::rdirichlet_mat(n = 5000,alpha_M3_e+1)
input_psa_M3$p_ff_e<-sample_M3_e[1,1,1:5000]
input_psa_M3$p_fp_e<-sample_M3_e[1,2,1:5000]
input_psa_M3$p_fd_e<-sample_M3_e[1,3,1:5000]
input_psa_M3$p_pf_e<-sample_M3_e[2,1,1:5000]
input_psa_M3$p_pp_e<-sample_M3_e[2,2,1:5000]
input_psa_M3$p_pd_e<-sample_M3_e[2,3,1:5000]
alpha_M3_f<-matrix(as.numeric(psa_trans_M3[6,5:10]),ncol = 3,byrow = TRUE)
sample_M3_f<-hesim::rdirichlet_mat(n = 5000,alpha_M3_f+1)
input_psa_M3$p_ff_f<-sample_M3_f[1,1,1:5000]
input_psa_M3$p_fp_f<-sample_M3_f[1,2,1:5000]
input_psa_M3$p_fd_f<-sample_M3_f[1,3,1:5000]
input_psa_M3$p_pf_f<-sample_M3_f[2,1,1:5000]
input_psa_M3$p_pp_f<-sample_M3_f[2,2,1:5000]
input_psa_M3$p_pd_f<-sample_M3_f[2,3,1:5000]
alpha_M3_g<-matrix(as.numeric(psa_trans_M3[7,5:10]),ncol = 3,byrow = TRUE)
sample_M3_g<-hesim::rdirichlet_mat(n = 5000,alpha_M3_g+1)
input_psa_M3$p_ff_g<-sample_M3_g[1,1,1:5000]
input_psa_M3$p_fp_g<-sample_M3_g[1,2,1:5000]
input_psa_M3$p_fd_g<-sample_M3_g[1,3,1:5000]
input_psa_M3$p_pf_g<-sample_M3_g[2,1,1:5000]
input_psa_M3$p_pp_g<-sample_M3_g[2,2,1:5000]
input_psa_M3$p_pd_g<-sample_M3_g[2,3,1:5000]

mean(input_psa_M3$p_ff_a)

# Markov 4 #
psa_trans_M4<-read_csv("./PSA_trans_M4.csv")
input_psa_M4 <-input_psa_M1     
alpha_M4_a<-matrix(as.numeric(psa_trans_M4[1,2:10]),ncol = 3,byrow = TRUE)
sample_M4_a<-hesim::rdirichlet_mat(n = 5000,alpha_M4_a+1)
input_psa_M4$p_tf<-sample_M4_a[1,1,1:5000]
input_psa_M4$p_tp<-sample_M4_a[1,2,1:5000]
input_psa_M4$p_td<-sample_M4_a[1,3,1:5000]
input_psa_M4$p_ff_a<-sample_M4_a[2,1,1:5000]
input_psa_M4$p_fp_a<-sample_M4_a[2,2,1:5000]
input_psa_M4$p_fd_a<-sample_M4_a[2,3,1:5000]
input_psa_M4$p_pf_a<-sample_M4_a[3,1,1:5000]
input_psa_M4$p_pp_a<-sample_M4_a[3,2,1:5000]
input_psa_M4$p_pd_a<-sample_M4_a[3,3,1:5000]
alpha_M4_b<-matrix(as.numeric(psa_trans_M4[2,5:10]),ncol = 3,byrow = TRUE)
sample_M4_b<-hesim::rdirichlet_mat(n = 5000,alpha_M4_b+1)
input_psa_M4$p_ff_b<-sample_M4_b[1,1,1:5000]
input_psa_M4$p_fp_b<-sample_M4_b[1,2,1:5000]
input_psa_M4$p_fd_b<-sample_M4_b[1,3,1:5000]
input_psa_M4$p_pf_b<-sample_M4_b[2,1,1:5000]
input_psa_M4$p_pp_b<-sample_M4_b[2,2,1:5000]
input_psa_M4$p_pd_b<-sample_M4_b[2,3,1:5000]
alpha_M4_c<-matrix(as.numeric(psa_trans_M4[3,5:10]),ncol = 3,byrow = TRUE)
sample_M4_c<-hesim::rdirichlet_mat(n = 5000,alpha_M4_c+1)
input_psa_M4$p_ff_c<-sample_M4_c[1,1,1:5000]
input_psa_M4$p_fp_c<-sample_M4_c[1,2,1:5000]
input_psa_M4$p_fd_c<-sample_M4_c[1,3,1:5000]
input_psa_M4$p_pf_c<-sample_M4_c[2,1,1:5000]
input_psa_M4$p_pp_c<-sample_M4_c[2,2,1:5000]
input_psa_M4$p_pd_c<-sample_M4_c[2,3,1:5000]
alpha_M4_d<-matrix(as.numeric(psa_trans_M4[4,5:10]),ncol = 3,byrow = TRUE)
sample_M4_d<-hesim::rdirichlet_mat(n = 5000,alpha_M4_d+1)
input_psa_M4$p_ff_d<-sample_M4_d[1,1,1:5000]
input_psa_M4$p_fp_d<-sample_M4_d[1,2,1:5000]
input_psa_M4$p_fd_d<-sample_M4_d[1,3,1:5000]
input_psa_M4$p_pf_d<-sample_M4_d[2,1,1:5000]
input_psa_M4$p_pp_d<-sample_M4_d[2,2,1:5000]
input_psa_M4$p_pd_d<-sample_M4_d[2,3,1:5000]
alpha_M4_e<-matrix(as.numeric(psa_trans_M4[5,5:10]),ncol = 3,byrow = TRUE)
sample_M4_e<-hesim::rdirichlet_mat(n = 5000,alpha_M4_e+1)
input_psa_M4$p_ff_e<-sample_M4_e[1,1,1:5000]
input_psa_M4$p_fp_e<-sample_M4_e[1,2,1:5000]
input_psa_M4$p_fd_e<-sample_M4_e[1,3,1:5000]
input_psa_M4$p_pf_e<-sample_M4_e[2,1,1:5000]
input_psa_M4$p_pp_e<-sample_M4_e[2,2,1:5000]
input_psa_M4$p_pd_e<-sample_M4_e[2,3,1:5000]
alpha_M4_f<-matrix(as.numeric(psa_trans_M4[6,5:10]),ncol = 3,byrow = TRUE)
sample_M4_f<-hesim::rdirichlet_mat(n = 5000,alpha_M4_f+1)
input_psa_M4$p_ff_f<-sample_M4_f[1,1,1:5000]
input_psa_M4$p_fp_f<-sample_M4_f[1,2,1:5000]
input_psa_M4$p_fd_f<-sample_M4_f[1,3,1:5000]
input_psa_M4$p_pf_f<-sample_M4_f[2,1,1:5000]
input_psa_M4$p_pp_f<-sample_M4_f[2,2,1:5000]
input_psa_M4$p_pd_f<-sample_M4_f[2,3,1:5000]
alpha_M4_g<-matrix(as.numeric(psa_trans_M4[7,5:10]),ncol = 3,byrow = TRUE)
sample_M4_g<-hesim::rdirichlet_mat(n = 5000,alpha_M4_g+1)
input_psa_M4$p_ff_g<-sample_M4_g[1,1,1:5000]
input_psa_M4$p_fp_g<-sample_M4_g[1,2,1:5000]
input_psa_M4$p_fd_g<-sample_M4_g[1,3,1:5000]
input_psa_M4$p_pf_g<-sample_M4_g[2,1,1:5000]
input_psa_M4$p_pp_g<-sample_M4_g[2,2,1:5000]
input_psa_M4$p_pd_g<-sample_M4_g[2,3,1:5000]

mean(input_psa_M4$p_fd_g)

## Run PSA ##
results_psa <-as.data.frame(matrix(data=NA,nrow=5000,ncol=16))      
colnames(results_psa) <- c("sample","qaly_M1","qaly_M2","qaly_M3","qaly_M4","cost_M1","cost_M2","cost_M3",
                           "cost_M4","qaly_pet","qaly_ct","cost_pet","cost_ct","i_cost","i_qaly","icer")   
results_psa$sample <- seq(1,5000,1)


psa_markov_func<-function(input_markov,input_gen){
  
  results_psa_markov<-as.data.frame(matrix(data=NA,nrow=5000,ncol=2))   
  colnames(results_psa_markov)<-c("cost","qaly")
  
  for (s in 1:5000) {
    ff<-c(rep(input_markov[s,5],3),rep(input_markov[s,11],5),rep(input_markov[s,17],5),rep(input_markov[s,23],5),
          rep(input_markov[s,29],5),rep(input_markov[s,35],5),rep(input_markov[s,41],42))
    fp<-c(rep(input_markov[s,6],3),rep(input_markov[s,12],5),rep(input_markov[s,18],5),rep(input_markov[s,24],5),
          rep(input_markov[s,30],5),rep(input_markov[s,36],5),rep(input_markov[s,42],42))
    fd<-c(rep(input_markov[s,7],3),rep(input_markov[s,13],5),rep(input_markov[s,19],5),rep(input_markov[s,25],5),
          rep(input_markov[s,31],5),rep(input_markov[s,37],5),rep(input_markov[s,43],42))
    pf<-c(rep(input_markov[s,8],3),rep(input_markov[s,14],5),rep(input_markov[s,20],5),rep(input_markov[s,26],5),
          rep(input_markov[s,32],5),rep(input_markov[s,38],5),rep(input_markov[s,44],42))
    pp<-c(rep(input_markov[s,9],3),rep(input_markov[s,15],5),rep(input_markov[s,21],5),rep(input_markov[s,27],5),
          rep(input_markov[s,33],5),rep(input_markov[s,39],5),rep(input_markov[s,45],42))
    pd<-c(rep(input_markov[s,10],3),rep(input_markov[s,16],5),rep(input_markov[s,22],5),rep(input_markov[s,28],5),
          rep(input_markov[s,34],5),rep(input_markov[s,40],5),rep(input_markov[s,46],42))
    
    state<-matrix(data=0,nrow = 70,ncol = 4)
    colnames(state)<-c("Firstline therapy","Progression free","Progression","Death")
    state[1,]<-c(1,0,0,0)
    
    Tmat<-matrix(data=0, nrow=4, ncol=4)
    Tmat[1,]<-c(0,input_markov[s,2],input_markov[s,3],input_markov[s,4])
    Tmat[4,]<-c(0,0,0,1)
    
    for (i in 2:70) {
      Tmat[2,]<-c(0,ff[i-1],fp[i-1],fd[i-1])
      Tmat[3,]<-c(0,pf[i-1],pp[i-1],pd[i-1])
      
      state[i,]<-state[(i-1),]%*%Tmat
    }
    
    ly<-as.data.frame(state)
    ly$life_years<-rowSums(ly[,1:3])
    qaly<-ly
    qaly[,1]<-ly[,1]*input_gen[s,14]
    qaly[,2]<-ly[,2]*input_gen[s,15]
    qaly[,3]<-ly[,3]*input_gen[s,16]
    qaly[,4]<-ly[,4]*0
    qaly[,5]<-rowSums(qaly[,1:4])
    colnames(qaly)[5]<-c("QALY_sum")
    
    cost<-ly
    cost[,1]<-ly[,1]*input_gen[s,9]
    cost[2,2]<-ly[2,2]*input_gen[s,10]
    cost[3,2]<-(ly[2,2]*ff[2]+ly[2,3]*pf[2])*input_gen[s,10]
    cost[2,3]<-ly[2,3]*input_gen[s,12]
    cost[3,3]<-(ly[2,3]*pp[2]+ly[2,2]*fp[2])*input_gen[s,12]
    for (i in 4:70) {
      f_new<-ly[i-2,3]*pf[i-2]*ff[i-1]+ly[i-1,3]*pf[i-1]
      p_new<-ly[i-2,2]*fp[i-2]*pp[i-1]+ly[i-1,2]*fp[i-1]
      cost[i,2]<-f_new*input_gen[s,10] + (ly[i,2]-f_new)*input_gen[s,11]
      cost[i,3]<-p_new*input_gen[s,12] + (ly[i,3]-p_new)*input_gen[s,13]
    }
    cost[,4]<-cost[,4]*0
    cost[,5]<-rowSums(cost[,1:4])
    colnames(cost)[5]<-c("cost_sum")
    
    ly_d<-ly
    ly_d[,1]<-ly[,1]*(1/(1+d_qaly)^cycle_index[1:70])
    ly_d[,2]<-ly[,2]*(1/(1+d_qaly)^cycle_index[1:70])
    ly_d[,3]<-ly[,3]*(1/(1+d_qaly)^cycle_index[1:70])
    ly_d[,5]<-rowSums(ly_d[,1:3])
    qaly_d<-qaly
    qaly_d[,1]<-qaly[,1]*(1/(1+d_qaly)^cycle_index[1:70])
    qaly_d[,2]<-qaly[,2]*(1/(1+d_qaly)^cycle_index[1:70])
    qaly_d[,3]<-qaly[,3]*(1/(1+d_qaly)^cycle_index[1:70])
    qaly_d[,5]<-rowSums(qaly_d[,1:3])
    cost_d<-cost
    cost_d[,1]<-cost[,1]*(1/(1+d_cost)^cycle_index[1:70])
    cost_d[,2]<-cost[,2]*(1/(1+d_cost)^cycle_index[1:70])
    cost_d[,3]<-cost[,3]*(1/(1+d_cost)^cycle_index[1:70])
    cost_d[,5]<-rowSums(cost_d[,1:3])
    
    ly_co<-ly_d
    ly_co[2:69,1]<-(ly_d[2:69,1]+ly_d[3:70,1])/2
    ly_co[2:69,2]<-(ly_d[2:69,2]+ly_d[3:70,2])/2
    ly_co[2:69,3]<-(ly_d[2:69,3]+ly_d[3:70,3])/2
    ly_co[,5]<-rowSums(ly_co[,1:3])
    qaly_co<-qaly_d
    qaly_co[2:69,1]<-(qaly_d[2:69,1]+qaly_d[3:70,1])/2
    qaly_co[2:69,2]<-(qaly_d[2:69,2]+qaly_d[3:70,2])/2
    qaly_co[2:69,3]<-(qaly_d[2:69,3]+qaly_d[3:70,3])/2
    qaly_co[,5]<-rowSums(qaly_co[,1:3])
    cost_co<-cost_d
    cost_co[2:69,1]<-(cost_d[2:69,1]+cost_d[3:70,1])/2
    cost_co[2:69,2]<-(cost_d[2:69,2]+cost_d[3:70,2])/2
    cost_co[2:69,3]<-(cost_d[2:69,3]+cost_d[3:70,3])/2
    cost_co[,5]<-rowSums(cost_co[,1:3])
    
    results_psa_markov$cost[s]<-sum(cost_co[,5])
    results_psa_markov$qaly[s]<-sum(qaly_co[,5])
  }
  return(results_psa_markov)
}

results_psa_M1 <- psa_markov_func(input_psa_M1,input_psa_gen_1)
results_psa_M2 <- psa_markov_func(input_psa_M2,input_psa_gen_1)
results_psa_M3 <- psa_markov_func(input_psa_M3,input_psa_gen_2)
results_psa_M4 <- psa_markov_func(input_psa_M4,input_psa_gen_2)
results_psa$qaly_M1 <- results_psa_M1$qaly
results_psa$qaly_M2 <- results_psa_M2$qaly
results_psa$qaly_M3 <- results_psa_M3$qaly
results_psa$qaly_M4 <- results_psa_M4$qaly
results_psa$cost_M1 <- results_psa_M1$cost
results_psa$cost_M2 <- results_psa_M2$cost
results_psa$cost_M3 <- results_psa_M3$cost
results_psa$cost_M4 <- results_psa_M4$cost


## QALY & Cost & ICER ##
for (s in 1:5000) {
  results_psa$cost_pet[s] <- (1-input_psa_gen_1[s,2])*input_psa_gen_1[s,3]*results_psa$cost_M1[s] + 
    input_psa_gen_1[s,2]*(1-input_psa_gen_1[s,4])*results_psa$cost_M2[s] +
    input_psa_gen_2[s,2]*input_psa_gen_2[s,4]*results_psa$cost_M3[s] + 
    (1-input_psa_gen_2[s,2])*(1-input_psa_gen_2[s,3])*results_psa$cost_M4[s] + input_psa_gen_2[s,7]*1
  
  results_psa$cost_ct[s] <- (1-input_psa_gen_1[s,2])*input_psa_gen_1[s,5]*results_psa$cost_M1[s] + 
    input_psa_gen_1[s,2]*(1-input_psa_gen_1[s,6])*results_psa$cost_M2[s] +
    input_psa_gen_2[s,2]*input_psa_gen_2[s,6]*results_psa$cost_M3[s] + 
    (1-input_psa_gen_2[s,2])*(1-input_psa_gen_2[s,5])*results_psa$cost_M4[s] + input_psa_gen_2[s,8]*1
  
  results_psa$qaly_pet[s] <- (1-input_psa_gen_1[s,2])*input_psa_gen_1[s,3]*results_psa$qaly_M1[s] + 
    input_psa_gen_1[s,2]*(1-input_psa_gen_1[s,4])*results_psa$qaly_M2[s] + 
    input_psa_gen_2[s,2]*input_psa_gen_2[s,4]*results_psa$qaly_M3[s] + 
    (1-input_psa_gen_2[s,2])*(1-input_psa_gen_2[s,3])*results_psa$qaly_M4[s]
  
  results_psa$qaly_ct[s] <- (1-input_psa_gen_1[s,2])*input_psa_gen_1[s,5]*results_psa$qaly_M1[s] + 
    input_psa_gen_1[s,2]*(1-input_psa_gen_1[s,6])*results_psa$qaly_M2[s] + 
    input_psa_gen_2[s,2]*input_psa_gen_2[s,6]*results_psa$qaly_M3[s] + 
    (1-input_psa_gen_2[s,2])*(1-input_psa_gen_2[s,5])*results_psa$qaly_M4[s]
  
  results_psa$i_cost[s] <- results_psa$cost_pet[s]-results_psa$cost_ct[s]
  results_psa$i_qaly[s] <- results_psa$qaly_pet[s]-results_psa$qaly_ct[s]
  
  results_psa$icer[s] <- results_psa$i_cost[s]/results_psa$i_qaly[s]
}


## Plot CE Plane ##
psa <- ggplot(results_psa,aes(i_qaly,i_cost)) +  geom_point(color = "grey40") +
  geom_point(aes(x=mean(i_qaly), y=mean(i_cost)), colour="red",size=1,shape = 7,stroke=2) + 
  geom_hline(yintercept = 0, size = 0.9, linetype = 1, colour="grey25") + 
  geom_vline(xintercept = 0, size = 0.9, linetype = 1, colour ="grey25") + 
  xlab("Incremental QALYs") + 
  ylab("Incenmental Cost") + 
  theme_classic() + 
  theme(axis.line.x = element_line(color="blue", size = 0, linetype=2),
               axis.line.y = element_line(color="blue", size = 0, linetype=2),
        axis.title = element_text(size = 14),axis.text = element_text(size = 12)) + 
  geom_abline(slope = 160000,intercept = 0,linetype = 2,linewidth=0.7,colour="blue")

psa

psa_mean<-mean(results_psa$i_cost)/mean(results_psa$i_qaly)

threh <- 160000
psa_prob_1 <- (length(results_psa$icer[results_psa$i_qaly>=0&results_psa$icer<threh])+length(results_psa$icer[results_psa$i_qaly<0&results_psa$icer>threh]))/5000
psa_prob_2 <- length(results_psa$i_qaly[results_psa$i_qaly>=0])/5000


## CEAC plot ##
num <- length(seq(from = 0, to = 740000, by = 2000))
CEAC <- as.data.frame(matrix(data=NA, nrow=num, ncol=4))
colnames(CEAC) <- c("Threshold","ProbCE","ProbCEbase")
CEAC$Threshold <- seq(from = 0, to = 740000, by = 2000)
for (i in 1:num) {       
  
  INMB <- results_psa$i_qaly*CEAC$Threshold[i] - results_psa$i_cost
  
  CEAC$ProbCE[i] <- sum(INMB >=0)/length(results_psa$i_qaly)
  
  CEAC$ProbCEbase[i] <- 1 - CEAC$ProbCE[i]
  
}

new_data <- as.data.frame(matrix(data=NA, nrow=2*num, ncol=3))
colnames(new_data) <- c("CEthreshold","ProbCE","Decision")
new_data$CEthreshold[1:num] <- CEAC$Threshold
new_data$ProbCE[1:num] <- CEAC$ProbCE
new_data$Decision[1:num] <- 'Cost-effective'

new_data$CEthreshold[(num+1):(num*2)] <- CEAC$Threshold
new_data$ProbCE[(num+1):(num*2)] <- CEAC$ProbCEbase
new_data$Decision[(num+1):(num*2)] <- 'Not Cost-effective'

# plot the CEAC 
p_CEAC <- ggplot(new_data,aes(CEthreshold,ProbCE)) + 
  geom_line(aes(x = CEthreshold,y = ProbCE,group=Decision,linetype= Decision, colour=Decision),size=0.8) + 
  scale_color_manual(values=c("red","black")) + 
  scale_size_manual(values=c(10,1)) + 
  scale_linetype_manual(values=c(1,2)) + 
  scale_x_continuous(breaks = c(0,100000,200000,300000,400000,500000,600000,686000),expand = c(0,0))+
  scale_y_continuous(limits = c(0,1),breaks = c(0,0.2,0.4,0.6,0.68,0.9,1.0),labels = scales::percent_format(accuracy = 1),expand = c(0,0))+
  annotate("segment",x=160000,xend=160000,y=0,yend=0.6778,colour = "blue",size=0.4,linetype = "dashed")+
  annotate("segment",x=0,xend=160000,y=0.6778,yend=0.6778,colour = "blue",size=0.4,linetype = "dashed")+
  annotate("segment",x=686000,xend=686000,y=0,yend=0.9,colour = "blue",size=0.4,linetype = "dashed")+
  annotate("segment",x=0,xend=686000,y=0.9,yend=0.9,colour = "blue",size=0.4,linetype = "dashed")+
  xlab("Willingness-to-pay Threshold (Thai Bhat per QALY)") + 
  ylab("Probability of PET/CT being Cost Effective") + 
  theme_classic()+
  theme(axis.title = element_text(size = 14),axis.text = element_text(size = 12),
        legend.text = element_text(size = 11))


p_CEAC






