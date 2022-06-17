##########################################################
###########################Model_PE_VT####################
##########################################################


rm(list=ls())
library(BayesianTools)
library(foreach)
library(doParallel)
library(emdbook)
library(doSNOW)
library(rootSolve)
library(deSolve)
library(gtools)
library(reshape)
library(ithir)
library(devtools)
library(Rmisc)

#==============================
## Model_PE_VT
#==============================
setwd("xxxx/")

##steady.1D
Model_DM_mcmc<- function(pars, initOM, ANPP, BNPP, temp, tempref, DM, PE, nn){
 
   Model_DM <- function(time, initOM, parms=NULL,DM, PE){
    
    N = length(depth)
     
    CA <- initOM[1:N]              
    CS <- initOM[(N+1):(2*N)]
    CP <- initOM[(2*N+1):(3*N)]
    
    Ad_cabl = fr_ANPP * INPUT_ablitter
    
    if(DM == TRUE){
      Tr_ca  = D_ca * diff(diff(c(CA[1], CA, 0.8*CA[N]))/ dx) / dx
      Tr_cs  = D_cs * diff(diff(c(CS[1], CS, 0.8*CS[N]))/ dx) / dx
      Tr_cp  = D_cp * diff(diff(c(CP[1], CP, 0.8*CP[N]))/ dx) / dx

      ## advection
      Ad_ca  = A_ca * diff(c(CA[1], CA)) / dx 
      Ad_cs  = A_cs * diff(c(CS[1], CS)) / dx
      Ad_cp  = A_cp * diff(c(CP[1], CP)) / dx
      
      DM_ca = Tr_ca - Ad_ca
      DM_cs = Tr_cs - Ad_cs
      DM_cp = Tr_cp - Ad_cp
      
      DM_ca[1] = Ad_cabl*fr_agl_ca
      DM_cs[1] = Ad_cabl*(1-fr_agl_ca)
      DM_cp[1] = 0
      
    }else{
      DM_ca = DM_cs = DM_cp = rep(0,N)
    }

    VT = DM_ca + DM_cs + DM_cp
    
    ## priming effect
    if (PE == TRUE) {
      Cf = INPUT_be + VT + CA
      
      PE_CA = 1+ V_max_CA * INPUT_be/(Km_CA + INPUT_be)
      PE_CS = 1+ V_max_CS * (DM_ca + INPUT_be + CA)/(Km_CS + DM_ca + INPUT_be + CA)
      PE_CP = 1+ V_max_CP * (DM_ca + INPUT_be + CA + CS)/(Km_CP + DM_ca + INPUT_be + CA + CS)
      
    }else{
      PE_CA = PE_CS = PE_CP = 1
      }
    
    ## dynamics of biom and hum
    dCA = INPUT_be* fr_agl_ca + k_cs * CS * PE_CS  * fr_cs_ca + k_cp * CP * PE_CP  * fr_cp_ca -
          k_ca * CA * PE_CA  + DM_ca
    
    dCS = INPUT_be * (1-fr_agl_ca)  + k_ca * CA * PE_CA  * fr_ca_cs -
          k_cs * CS * PE_CS + DM_cs
    
    dCP = k_ca * CA* PE_CA  * fr_ca_cp + k_cs * CS * PE_CS  * fr_cs_cp -
          k_cp * CP * PE_CP + DM_cp
    
    return(list(c(dCA, dCS, dCP), vt= VT, pe = PE_C))
   }
   
  ## parameters
  pars  = as.numeric(pars)
  D_ca0 = pars[1]
  D_cs0 = pars[2]
  D_cp0 = pars[3]
  A_ca0 = pars[4]
  A_cs0 = pars[5]
  A_cp0 = pars[6]
  ddm = pars[7]
  V_max0 = pars[8]
  Km_cf0 = pars[9]
  k_ca00 = pars[10]
  k_cs00 = pars[11]
  k_cp00 = pars[12]
  dzk = pars[13]
  fr_cs_ca = pars[14]
  fr_cp_ca = pars[15]
  fr_ca_cp = pars[16]
  fr_cs_cp = pars[17]
  fr_ca_cs = pars[18]
  fr_ANPP = 0.1
  fr_agl_ca = 0.9 ### The fraction of NPP to  CA pool 
  

  ## temperature function
  gas = 0.008314
  # Ea_ca = 30
  # Ea_cs = 70 #70 for a tannin compound from Davidson and Janssens, 2006; 135.5 from Leifeld and Lutzow,2014,#activation energy for slow pool decompostion
  # Ea_cp = 135 #Leifeld and Lutzow,2014
  # 
  Ea_ca = 47
  Ea_cs = 47 
  Ea_cp = 47 #Allison SD, Wallenstein MD, Bradford MA. Soil-carbon response to warming dependent on microbial physiology. Nature Geoscience 3, 336-340 (2010).
  
  ft_ca = exp((-Ea_ca/gas*(1/(temp+273)-1/(tempref+273)))) #
  ft_cs = exp((-Ea_cs/gas*(1/(temp+273)-1/(tempref+273)))) #¡¤
  ft_cp = exp((-Ea_cp/gas*(1/(temp+273)-1/(tempref+273)))) #¡¤
  
  k_ca0 = k_ca00 * ft_ca
  k_cs0 = k_cs00 * ft_cs
  k_cp0 = k_cp00 * ft_cp
  
  ## depth effect of parameters
  D_ca = D_ca0 * exp(-depth/ddm)
  D_cs = D_ca0 * exp(-depth/ddm)
  D_cp = D_ca0 * exp(-depth/ddm)
  
  A_ca = A_ca0 * exp(-depth/ddm)
  A_cs = A_cs0 * exp(-depth/ddm)
  A_cp = A_cp0 * exp(-depth/ddm)
  
  V_max_CA = V_max_CS = V_max_CP = V_max0
  Km_CA = Km_CS = Km_CP = Km_cf0 
  
  k_ca = k_ca0 * exp(-depth/dzk)
  k_cs = k_cs0 * exp(-depth/dzk)
  k_cp = k_cp0 * exp(-depth/dzk)
  
  ## input
  INPUT_ablitter = ANPP / dx #litter_below ##gC/m2/year
  INPUT_be = BNPP / dx
  
  ## solve model steady state
  sol <- steady.1D(y = initOM, func = Model_DM, parms = NULL, nspec = 3,
                   names = c("CA","CS","CP"), DM=DM, PE=PE) 
  
  layer.soil <- rep(paste0("L", c(1:nn)), times = c(4,4,4,4,4,10,10))
  
  if(!attr(sol,"steady")){
    SOC_dx005 <- rep(-9999,times= nn)
    VT <- rep(NA, times = nn)
    pe <- rep(NA, times = nn)
  }else if(any(sol$y<=0) | any(sol$y %in% NaN) | any(sol$y %in% NA)){
    SOC_dx005 <- rep(-9999,times=nn)
    VT <- rep(NA, times = nn)
    pe <- rep(NA, times = nn)
  }else{
    df11 <- data.frame(layer.soil, SOC00 = rowSums(sol$y *0.05 ))
    SOC_dx005 = tapply(df11$SOC00, list(df11$layer.soil), sum)

    df.VT <- data.frame(layer.soil, VT = sol$vt*dx)
    VT = tapply(df.VT$VT, list(df.VT$layer.soil), mean)
    #
    df.PE <- data.frame(layer.soil, PE = sol$pe)
    pe = tapply(df.PE$PE, list(df.PE$layer.soil), mean)

  }
  # return(data.frame(SOC_dx005, VT, pe))
  return(SOC_dx005)
  
 
}


### run MCMC
#========================================================           
cores=detectCores()
cl <- makeCluster(cores[1]-5) #not to overload your computer
registerDoSNOW(cl)

paramers = read.csv("data/Prior_paras.csv") ##prior parameters
# pars_selec <- c(1:3, 5:15) 
# paramers[pars_selec, c(2:4)] = log(paramers[pars_selec, c(2:4)]) ##log translation
pars = paramers$best

nProfile =2000

Input_data <- read.csv( "E:/Mingming/SOC-DM-model/Century_model/data/Input_anpp_bnpp_soc.csv")
envfactor <- read.csv("E:/Mingming/SOC-DM-model/result/envfactor.csv")
envtemp <- envfactor$AMT

dx = 0.05 ## discretization of soil depth
depth = c(seq(0,2, dx)[-1]-dx/2)
N = length(depth)

DMPE= c(T,T,F,F, T,F,T,F)
initOM0 = rep(10,3*N)## initial value

pb <- txtProgressBar(max = nProfile, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

settings <- list(iterations = 1000000,  burnin=100000,  adaptation = 20000)
prior <- createUniformPrior(lower = paramers$lower, best = paramers$best, upper = paramers$upper)


for (ii in 1:4) {
  
  DM = DMPE[ii]
  PE = DMPE[ii+4]
  
  system.time(foreach(iProfile = 1:nProfile, .packages=c("BayesianTools","emdbook","rootSolve"),
                      .combine = list,.options.snow = opts) %dopar%{
                        
                        data0 = Input_data[iProfile,]
                        location =as.numeric(data0[1])
                        BNPP = as.numeric(data0[c(15:54)])
                        ANPP = as.numeric(data0[7])
                        obs  = as.numeric(data0[8:14])
                        
                        likelihood <- function(parm, sum = T){
                          # set parameters that are not calibrated on default values
                          x = parm
                          # x[pars_selec] = exp(parm[pars_selec])
                          result <- Model_DM_mcmc(pars=x, initOM=initOM0, ANPP=ANPP, BNPP=BNPP,
                                                    temp=20,tempref=20, DM=DM, PE=PE)
                          diff_mod_obs <- result - obs
                          # run the model
                          llValues = (dnorm(diff_mod_obs, sd = 0.2*obs, log = TRUE))
                          if (sum == T) return(sum(llValues))
                          else return(llValues)
                        }
                        
                        bayesianSetup = createBayesianSetup(likelihood, prior, best = paramers$best)
                        out_Century <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DREAMzs", settings = settings)
                        saveRDS(out_Century, file = paste0("I:/out_MCMC_3CP_TT/MCMC_0.05dx_",
                                                        substr(as.character(DM),1,1),substr(as.character(PE),1,1),
                                                        "/","MCMC","_3CP_",location,"_",
                                                        substr(as.character(DM),1,1),substr(as.character(PE),1,1),
                                                        ".Rdata"))
                        stopParallel(bayesianSetup)
                      })
  print(paste0( substr(as.character(DM),1,1),substr(as.character(PE),1,1)))
}
stopCluster(cl)

#======================================================== 


