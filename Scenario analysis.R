##########################################################
########Warming and carbon input increase scenarios#######
##########################################################

library(ggplot2)
library(gtools)
library(caret)
library(foreach)
library(doSNOW)
library(doParallel)
library(readr)
library(metR)
library(cowplot)
library(RColorBrewer)

rm(list = ls())

####response of SOC on warming to steady sate
FUN_socwarming <- function(iProfile, ins_npp, ins_temp, pars_selec=c(1:3, 5:14) ){
  
  out_Century_mcmc<- readRDS(paste(path, "MCMC_0.05dx_",
                                   substr(as.character(DM),1,1),substr(as.character(PE),1,1),
                                   "/", filenames[ip[iProfile]],sep = ""))
  gelmanDiagnostics(out_Century_mcmc)
  pars_1000 = getSample(out_Century_mcmc, parametersOnly = T,numSamples = npars,thin = 1,
                        start = 200000, coda = T)
  pars_1000 = do.call(rbind.data.frame, pars_1000)
  # pars_1000[,pars_selec] = exp(pars_1000[,pars_selec])
  
  data0 = Input_data[ip[iProfile],]
  location =as.numeric(data0[1])
  BNPP = as.numeric(data0[c(15:54)])
  ANPP = as.numeric(data0[7])
  obs  = as.numeric(data0[8:14])
  
  ggja<- function(x){
    initOM0 = rep(10,3*length(depth))
    SOCpool0<- Century_DM_mcmc(pars_1000[x,], initOM=initOM0,
                               ANPP, BNPP, temp=envtemp[ip[iProfile]], 
                               tempref=envtemp[ip[iProfile]],
                               DM=DM, PE=PE)
    
    if(all(is.na(SOCpool0))){
      Re_Cpool <- matrix(NA, N, 3)
      Re_TC <- matrix(NA, N, 1)
      deltaSOC <-  matrix(NA, N, 1)
      re = NA
    }else{
      ##warming 2¡æ
      anpp <- ANPP*ins_npp
      bnpp <- BNPP*ins_npp
      initOM2 <- melt(SOCpool0)$value
      TEMP = envtemp[ip[iProfile]] + ins_temp
      SOCpool2w<- Century_DM_mcmc(pars_1000[x,], initOM=initOM2,
                                  ANPP=anpp, BNPP=bnpp, 
                                  temp= TEMP, 
                                  tempref=envtemp[ip[iProfile]] ,
                                  DM=DM, PE=PE)
      
      if(all(is.na(SOCpool2w))){
        Re_Cpool <- matrix(NA, N, 3)
        Re_TC <- matrix(NA, N, 1)
        deltaSOC <-  matrix(NA, N, 1)
        re = NA
      }else{
        deltaSOC <- rowSums(SOCpool2w-SOCpool0)
        Re_Cpool <- (SOCpool2w-SOCpool0)/SOCpool0
        Re_TC <- (rowSums(SOCpool2w)-rowSums(SOCpool0))/rowSums(SOCpool0)
        re = (sum(rowSums(SOCpool2w)*dx) - sum(rowSums(SOCpool0)*dx))/ sum(rowSums(SOCpool0)*dx)
      }
    }
    return(cbind(deltaSOC, Re_Cpool, Re_TC, re))
    # return(re)
  }
  
  SOC0 <- sapply(1:npars, ggja)
  ReCpool1000means <- rowMeans(SOC0, na.rm = T)
  ReCpool1000sd <- apply(SOC0, 1, function(x){sd(x, na.rm = T)})
  # ReCpool1000means <- mean(SOC0, na.rm = T)
  # ReCpool1000sd <- sd(SOC0, na.rm = T)
  Cpools<- rep(c("delta.SOC","CA","CS","CP","TC", "allpro"), each=N)
  return(data.frame(Cpools, ReCpool1000means,ReCpool1000sd))
}

npars =1000
Input_data <- read.csv( "xxxx/Input_anpp_bnpp_soc.csv")

envfactor <- read.csv("xxxx/envfactor.csv")
envtemp <- envfactor$AMT

nProfile =2000
ip<-c(1:nProfile)
ip<- ip[-c(1171,1742)]

dx = 0.05 ## discretization of soil depth
depth = c(seq(0,2, dx)[-1]-dx/2)
N = length(depth)

DMPE= c(T,T,F,F, T,F,T,F)
initOM0 = rep(10,3*N)

nProfile = 2000
npars =1000
ip<-c(1:nProfile)
ip<- ip[-c(1171,1742)]

ins_npp_seq <- seq(1,1.4, 0.1)
ins_temp_seq = c(0, 1.5, 2.0, 3.0, 4.0, 5.0)  ##IPCC
ch.temp.npp = expand.grid(ins_npp_seq, ins_temp_seq)

cores=detectCores()
cl <- makeCluster(cores[1]-20) #not to overload your computer
registerDoSNOW(cl)

pb <- txtProgressBar(max = nProfile, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

for (jj in 9:dim(ch.temp.npp)[1]) {
  
  ins_npp = ch.temp.npp[jj,1]
  ins_temp = ch.temp.npp[jj,2]
  path = "xxxx/"
  
  
  for (ii in 1:4) {
    DM = DMPE[ii]
    PE = DMPE[ii+4]
    
    print(paste0(ins_temp, substr(as.character(DM),1,1),substr(as.character(PE),1,1)))
    
    a <- list.files(paste("MCMC_result/MCMC_0.05dx_",
                          substr(as.character(DM),1,1),substr(as.character(PE),1,1),sep = ""))
    filenames <- mixedsort(a)
    
    SOC100yr2000Profile <- foreach(iProfile = 1:length(ip) , .packages = c("BayesianTools","rootSolve","reshape","deSolve"),
                                   combine = list, .options.snow = opts) %dopar% {
                                     FUN_socwarming(iProfile, ins_npp, ins_temp)
                                   }
    
    saveRDS(SOC100yr2000Profile, file = paste("xxxx/re_SOC_3cp_",ins_temp,"_",ins_npp,
                                              "_",substr(as.character(DM),1,1),substr(as.character(PE),1,1),".Rdata", sep = ""))
  }
}
stopCluster(cl)
# FUN_socwarming(5)


#### response of SOC to warming and carbon input increase 
dir.sce = "xxxx/"
file.names = mixedsort(list.files(dir.sce))

re.sce = list()
for (ii in 1:length(file.names)) {
  
  dat = readRDS(paste0(dir.sce, file.names[ii]))
  str.scen = strsplit(file.names[ii], '_')
  
  VT = gsub(".Rdata", "", str.scen[[1]][6])
  inT = as.numeric(str.scen[[1]][4])
  inCinput =   as.numeric(str.scen[[1]][5])
  
  re.TC.pro = data.frame(re.TC =sapply(1:length(dat), function(xx){
    unique(dat[[xx]]$ReCpool1000means[161:200])
  }))
  
  re.TC.pro$VT = VT  
  re.TC.pro$inT = inT
  re.TC.pro$inCinput = inCinput
  
  re.sce[[ii]] = re.TC.pro
}

re.sce = do.call(rbind, re.sce)
# re.sce$inCinput = (re.sce$inCinput-1)*100

re.sce.mean = aggregate(re.sce$re.TC, list(re.sce$VT, re.sce$inT, re.sce$inCinput), function(x)mean(x, na.rm = T))
re.sce.mean$Group.3 = factor((re.sce.mean$Group.3 - 1)*100)
re.sce.mean$Group.2 = factor(re.sce.mean$Group.2)
re.sce.mean$Group.1 = factor(re.sce.mean$Group.1, labels = c("M","M+PE","M+VT","M+PE+VT"))

re.sce.mean.TT = re.sce.mean[which(re.sce.mean$Group.1 == "M+PE+VT"), ]
re.sce.mean.TT.3 = data.frame(apply(re.sce.mean.TT, 2, function(x)rep(x, each = 3)))
re.sce.mean.TT.3$x = as.numeric(re.sce.mean.TT.3$x)
re.sce.mean.fTT = re.sce.mean[-which(re.sce.mean$Group.1 == "M+PE+VT"), ]

bias =( re.sce.mean.fTT$x - re.sce.mean.TT.3$x)*100
re.sce.mean.fTT$bias = bias
re.sce.mean.fTT$Vars = paste0(re.sce.mean.fTT$Group.1, "-", re.sce.mean.TT.3$Group.1)

p_scenario_all <-
  ggplot(re.sce.mean.fTT)+
  geom_tile(aes(Group.2, Group.3, fill = bias))+
  geom_text(aes(Group.2, Group.3, label = round(bias,0)), size=3)+
  facet_wrap(~Vars, ncol=2)+
  xlab("Warming level (¡ãC)")+
  ylab(("Carbon input increase (%)"))+
  scale_fill_gradient2(nam = "¦¤SOC (%)",low = "#D62828", mid = "white", high = "#1B98E0")+
  theme_bw()
p_scenario_all

p_scenario_all <-
  ggplot(re.sce.mean)+
  geom_tile(aes(Group.2, Group.3, fill = 100*x))+
  geom_text(aes(Group.2, Group.3, label = round(100*x,0)), size=3)+
  facet_wrap(~Group.1, ncol=2)+
  xlab("Warming level (¡ãC)")+
  ylab(("Carbon input increase (%)"))+
  scale_fill_gradient2(nam = "¦¤SOC (%)",low = "#D62828", mid = "white", high = "#1B98E0")+
  theme_bw()
p_scenario_all
ggsave(p_scenario_all, filename="I:/Fig_model_VT/Fig_3cp/Fig.4.Scenario_all_response.pdf", 
       height = 4.13, width = 6.25, unit="in")

re.sce.mean$Group.2 = as.numeric(as.character(re.sce.mean$Group.2))
re.sce.mean$Group.3 = factor(re.sce.mean$Group.3)

p_Tin_deltaSOC <- 
ggplot(re.sce.mean)+
  geom_hline(yintercept = 0, color = "grey", linetype=2)+
  geom_vline(xintercept = 2, color = "grey", linetype=2)+
  geom_line(aes(Group.2, x*100, linetype = Group.1, color = Group.3), size=0.5)+
  xlab("Warming level (¡ãC)")+
  ylab(("¦¤SOC (%)"))+
  scale_color_manual(name = "Carbon input (%)",
                     values = brewer.pal(9, "Set1")[1:5])+
  scale_linetype_manual(name = "Model variates", values = c("twodash","dashed","dotted","solid"))+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color= "black"),
        )

fun.slope = function(x, y){
  res = summary(lm(y~x))
  res$coefficients[2,1]
}

re.sce.mean1 = re.sce.mean
re.sce.mean1$Group.2 = as.numeric(as.character(re.sce.mean1$Group.2))
re.sce.mean1$Group.3 = as.numeric(as.character(re.sce.mean1$Group.3))
mod = unique(re.sce.mean1$Group.1)
cinput = unique(re.sce.mean1$Group.2)
var1 = expand.grid(mod, cinput)
var1$slope = NA
for (ii in 1:nrow(var1)) {
    dat0=re.sce.mean1[which(re.sce.mean1$Group.1 == var1[ii,1]&re.sce.mean1$Group.2==var1[ii,2]),]
    x=dat0$Group.3
    y=dat0$x*100
    var1$slope[ii] = fun.slope(x, y)
}

P1 = ggplot(var1[1:4,], aes(Var1, slope))+
  geom_bar(stat="identity", fill = "#00A8E8", width = 0.4)+
  scale_y_continuous(expand = c(0,0.01))+
  ylab(expression(paste("%¦¤SOC %C+"^-1)))+
  xlab("")+
  theme_half_open()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(color= "black", size=10, angle=45, hjust = 1),
        axis.title = element_text(size=10))

mod = unique(re.sce.mean1$Group.1)
warming = unique(re.sce.mean1$Group.3)
var2 = expand.grid(mod, warming)
var2$slope = NA
for (ii in 1:nrow(var)) {
  dat0=re.sce.mean1[which(re.sce.mean1$Group.1 == var2[ii,1]&re.sce.mean1$Group.3==var2[jj,2]),]
  x=dat0$Group.2
  y=dat0$x*100
  var2$slope[ii] = fun.slope(x, y)
}

P2 = ggplot(var2[1:4,], aes(Var1, slope))+
  geom_bar(stat="identity", fill = "#D62828", width = 0.4)+
  scale_y_reverse(expand = c(0,0.1))+
  ylab(expression(paste("%¦¤SOC ¡ãC"^-1)))+
  xlab("")+
  theme_half_open()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(color= "black", size=10, angle=45, hjust = 1),
        axis.title = element_text(size=10))
  
P_sensitivity = plot_grid(P2, P1, ncol = 2)

df.Re.MC.3CP.2T12C = read.csv("I:/df.Re.MC.3CP.2T12C.csv")
df.Re.MC.3CP.2T12C$Mod[which(df.Re.MC.3CP.2T12C$Mod%in%c("F","FF"))] = "M"
df.Re.MC.3CP.2T12C$Mod[which(df.Re.MC.3CP.2T12C$Mod%in%c("T","TF"))] = "M+VT"
df.Re.MC.3CP.2T12C$Mod[which(df.Re.MC.3CP.2T12C$Mod%in%c("FT"))] = "M+PE"
df.Re.MC.3CP.2T12C$Mod[which(df.Re.MC.3CP.2T12C$Mod%in%c("TT"))] = "M+VT+PE"

df.Re.3CP.2T12C = df.Re.MC.3CP.2T12C[which(df.Re.MC.3CP.2T12C$bMod=="Carbon pool model"),]
P_MC3CP_Re_2T12C_7L<-
  ggplot(df.Re.3CP.2T12C, aes(Depth, Re*100))+
  geom_ribbon(aes(ymin = (Re-sd)*100, ymax = (Re+sd)*100, fill = Mod), alpha = 0.4)+
  geom_line(aes(color = Mod, linetype=Mod), size=1)+
  # facet_wrap(~bMod)+
  geom_hline(yintercept = 0, linetype=2, color = "grey")+
  scale_x_reverse()+
  # scale_x_continuous(expand = c(0,0))+
  coord_flip()+
  xlab("Depth (m)")+
  ylab("¦¤SOC (%)")+
  scale_fill_manual(name="Model variates", values = c( "#E76F51","#FDCA40","#00A8E8","#2BA84A"))+
  scale_color_manual(name="Model variates", values = c("#E76F51","#FDCA40","#00A8E8", "#2BA84A"))+
  scale_linetype_manual(name = "Model variates", values = c("twodash","dashed","dotted","solid"))+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none",
        legend.background = element_blank())

p_Re <-
  ggdraw()+
  draw_plot(p_Tin_deltaSOC, 0,0, 0.65,1)+
  draw_plot(P_MC3CP_Re_2T12C_7L, 0.65,0,0.35,1)

ggdraw(p_Re)+
  draw_plot(P_sensitivity, 0.12,0.66,0.4,0.32)

ggsave(p_Re, filename="I:/Fig_model_VT/Fig_3cp/P_Re.pdf", 
       height =4.26 , width = 8.74, unit="in")

  
class(re.sce.mean$Group.3)

re.sce.mean$Group.2= as.numeric(as.character(re.sce.mean$Group.2))
re.sce.mean$Group.3= as.numeric(as.character(re.sce.mean$Group.3))

value = c(30,20,20,-10,-20,-30)

p_scenario_faceline <-
ggplot(re.sce.mean,
       aes(x=Group.2, y=Group.3, z=x*100), order = 5)+
  geom_contour_fill(binwidth = 1, show.legend = T)+
  stat_contour(binwidth=10, color = "grey", size = 0.1)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  facet_wrap(~Group.1, ncol=2)+
  xlab("Warming level (¡ãC)")+
  ylab(("Carbon input increase (%)"))+
  scale_fill_gradient2(name = "¦¤SOC (%)",low = "#D62828", mid = "white", high = "#1B98E0")+
  theme_bw()
p_scenario_faceline
ggsave(p_scenario_faceline, filename="I:/Fig_model_VT/Fig_3cp/p_scenario_faceline.pdf", 
       height = 4.13, width = 6.25, unit="in")

cols.re = ((colorRampPalette(c("#D62828", "grey", "#1B98E0"))(7)))

TinCin0 = read.csv("I:/Tin-Cin-zero.csv")
TinCin0$deltaSOC = factor(TinCin0$deltaSOC)
p_Tin_Cin0 <- 
ggplot(TinCin0, aes(Tin, Cin))+
  geom_line(aes(linetype = Mod, color = deltaSOC), size=0.5)+
  # scale_fill_gradient2(name = "¦¤SOC (%)",low = "#D62828", mid = "white", high = "#1B98E0")+
  scale_color_manual(name = "¦¤SOC (%)", values = cols.re)+
  scale_linetype_manual(values = c("twodash","dashed","dotted","solid"))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  xlab("Warming level (¡ãC)")+
  ylab(("Carbon input increase (%)"))+
  theme_bw()+
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(),
        # legend.position = c(x=0.8, y=0.3),
        legend.background = element_blank())
p_Tin_Cin0
ggsave(p_Tin_Cin0, filename="xxx/Fig_p_Tin_Cin0_line.pdf", 
       height = 4.70, width = 6.17, unit="in")
 

p_scenario_faceline_0 <- 
  ggdraw()+
  draw_plot(p_scenario_faceline, 0, 1/3, 1, 2/3)+
  draw_plot(p_Tin_Cin0, 0, 0, 0.455, 1/3)+
  draw_label(label = c("a"), x=c(0.009), y=c(0.99))+
  draw_label(label = c("b"), x=c(0.009), y=c(1/3))
ggsave(p_scenario_faceline_0, filename="xxx/Fig_3cp/p_scenario_faceline_0.pdf", 
       height = 7.39, width = 7.05, unit="in")





### bias between M and M+PE+VT among biomes
scen.2t2c.TT = re.sce[which(re.sce$VT == "TT" & re.sce$inT == 2.0 & re.sce$inCinput == 1.2),]
scen.2t2c.FF = re.sce[which(re.sce$VT == "FF" & re.sce$inT == 2.0 & re.sce$inCinput == 1.2),]
scen.2t2c.TT$re.bias = (scen.2t2c.TT$re.TC -scen.2t2c.FF$re.TC)/scen.2t2c.FF$re.TC
scen.2t2c.TT$re.bias[which(abs(scen.2t2c.TT$re.bias )>10)] = NA
scen.2t2c.TT$biomes = factor(biomes)

eco.type = c("TS forests"  ,           "Temperate forests",        
             "Boreal forests",          "Temperatre grasslands",
             "TS grasslands/savannas", "Med shrublands", 
             "Tundra"  ,               "Deserts",            
             "Croplands")

scen.2t2c.TT$biomes = factor(scen.2t2c.TT$biomes, labels = eco.type)
scen.2t2c.TT.maen = aggregate(scen.2t2c.TT$re.bias, list(scen.2t2c.TT$biomes), function(x)mean(x, na.rm  = T))

p_re_bias_biomes <- 
ggplot(scen.2t2c.TT, aes(biomes, re.bias))+
  geom_boxplot(fill = "#DCDCDD",outlier.colour = "grey",show.legend = F, lwd=0.1, outlier.size = 1)+
  geom_hline(yintercept = 0, linetype=2, color = "black")+
  geom_point(data = scen.2t2c.TT.maen, aes(Group.1, x),  color = "#DD1C1A", size=2.5)+
  xlab("Biomes")+
  ylab("Relative bias of response")+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.background = element_blank(), 
        plot.background = element_blank())
p_re_bias_biomes
ggsave(p_re_bias_biomes, filename="I:/Fig_model_VT/Fig.5.p_re_bias_biomes.pdf", 
       height = 4.3, width = 4.5, unit="in")

dir.sce = "I:/Scensior analysis/"
file.names = mixedsort(list.files(dir.sce))



########################################################################
dir.sce = "I:/Scensior analysis/"
file.names = mixedsort(list.files(dir.sce))

re.sce = list()
dx = 0.05 ## discretization of soil depth
depth = c(seq(0,2, dx)[-1]-dx/2)
N = length(depth)
layer.soil <- rep(paste0("L", c(1:7)), times = c(4,4,4,4,4,10,10))

for (ii in 1:length(file.names)) {
  
  dat = readRDS(paste0(dir.sce, file.names[ii]))
  str.scen = strsplit(file.names[ii], '_')
  
  VT = gsub(".Rdata", "", str.scen[[1]][6])
  inT = as.numeric(str.scen[[1]][4])
  inCinput =   as.numeric(str.scen[[1]][5])
  
  re.TC =lapply(1:length(dat), function(xx){
    bjlad = dat[[xx]][121:160, ]
    bjlad$depth = depth
    bjlad$layer = layer.soil
    return(bjlad)
  })
  re.TC.pro = do.call(rbind.data.frame, re.TC)
  
  re.TC.pro$VT = VT  
  re.TC.pro$inT = inT
  re.TC.pro$inCinput = inCinput
  
  re.sce[[ii]] = re.TC.pro
}

re.sce = do.call(rbind.data.frame, re.sce)
re.sce$ReCpool1000means[which(abs(re.sce$ReCpool1000means) >= 1, )] = NA
re.sce$ReCpool1000sd[which(abs(re.sce$ReCpool1000means) >= 1, )] = NA

# fun_uncer <- function(XX){mean(sample(XX, 1998, replace = T), na.rm = T)}
# fun.bootstrap <- function(yy){aggregate(re.sce$ReCpool1000means, 
#                         list(re.sce$VT, re.sce$inT, re.sce$inCinput, re.sce$depth), 
#                         function(XX)fun_uncer(XX))
#   }
# 
# bootstrap200 <- lapply(1:200, fun.bootstrap)
# bootstrap200 = do.call(rbind.data.frame,bootstrap200)
# 
# re.sce.mean = aggregate(bootstrap200$x, list(bootstrap200$Group.1, bootstrap200$Group.2, 
#                                              bootstrap200$Group.3, bootstrap200$Group.4), function(x)mean(x, na.rm=T))
# re.sce.sd = aggregate(bootstrap200$x, list(bootstrap200$Group.1, bootstrap200$Group.2, 
#                                              bootstrap200$Group.3, bootstrap200$Group.4), sd)
# names(re.sce.mean) = c("Mod","Tin","Cin","Depth","Re")
# re.sce.mean$sd = re.sce.sd$x
# head(re.sce.mean)
# re.sce.mean$Cin = (re.sce.mean$Cin - 1) *100

MCP.re.sce.mean.2T12C = re.sce.mean[which(re.sce.mean$Tin==2& re.sce.mean$Cin==1.2), ]

re.sce.mean$Tin = factor(re.sce.mean$Tin, labels = paste("+", unique(re.sce.mean$Tin), "¡ãC"))
re.sce.mean$Cin = factor(re.sce.mean$Cin, labels = paste("+", unique(re.sce.mean$Cin), "(%)"))

P_Re_Tin_Cin_7L<-
ggplot(re.sce.mean, aes(Depth, Re*100, color = Mod))+
  geom_hline(yintercept = 0, color = "grey", linetype = 2)+
  geom_line(size=0.5)+
  scale_x_reverse()+
  ylab("¦¤SOC (%)")+
  xlab("Soil depth (m)")+
  coord_flip()+
  facet_grid(Tin~Cin)+
  theme_bw()+
  scale_color_manual(name = "Model variates",
                     values = c("#2BA84A","#FDCA40","#00A8E8", "#E76F51"), 
                     label = c("M","M+PE","M+VT","M+VT+PE"))+
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        legend.position = "bottom",
        legend.direction =  "horizontal",
        axis.text = element_text(color = "black"))
ggsave(P_Re_Tin_Cin_7L, filename="I:/Fig_model_VT/Fig_3cp/Fig.S_Re_Tin_Cin_7L.pdf", 
       height = 8.24, width = 8.24, unit="in")

df.re = aggregate(re.sce[,2:3], list(re.sce$depth,re.sce$layer, re.sce$VT, re.sce$inT, re.sce$inCinput), function(x)mean(x, na.rm = T))

names(df.re) <- c("Depth","Layer","Mod","Tin", "Cin", "Re", "sd")

re.sce.mean = df.re
re.sce.mean$Tin = factor(re.sce.mean$Tin, labels = paste("+", unique(re.sce.mean$Tin), "¡ãC"))
re.sce.mean$Cin = factor(re.sce.mean$Cin, labels = paste("+", unique(re.sce.mean$Cin), "(%)"))

huafb = re.sce.mean[which(re.sce.mean$Tin == "+ 2 ¡ãC" & re.sce.mean$Cin == "+ 1.2 (%)"),]
huafb2w = re.sce.mean[which(re.sce.mean$Tin == "+ 2 ¡ãC"& re.sce.mean$Cin == "+ 1 (%)"),]
huafb12C = re.sce.mean[which(re.sce.mean$Tin == "+ 0 ¡ãC"&re.sce.mean$Cin == "+ 1.2 (%)"),]


huafb$Re[which(huafb$Mod == "FF")][8] = mean(huafb$Re[which(huafb$Mod == "FF")][1:7])
huafb$Re[which(huafb$Mod == "FF")][38:40] = mean(huafb$Re[which(huafb$Mod == "FF")][1:7])
huafb$sd[which(huafb$Mod == "FF")][8] = mean(huafb$sd[which(huafb$Mod == "FF")][1:7])
huafb$sd[which(huafb$Mod == "FF")][38:40] = mean(huafb$sd[which(huafb$Mod == "FF")][1:7])

huafb2w$Re[which(huafb2w$Mod == "FF")][8] = mean(huafb2w$Re[which(huafb2w$Mod == "FF")][1:7])
huafb2w$Re[which(huafb2w$Mod == "FF")][38:40] = mean(huafb2w$Re[which(huafb2w$Mod == "FF")][1:7])
huafb2w$sd[which(huafb2w$Mod == "FF")][8] = mean(huafb2w$sd[which(huafb2w$Mod == "FF")][1:7])
huafb2w$sd[which(huafb2w$Mod == "FF")][38:40] = mean(huafb2w$sd[which(huafb2w$Mod == "FF")][1:7])

huafb12C$Re[which(huafb12C$Mod == "FF")][8] = mean(huafb12C$Re[which(huafb12C$Mod == "FF")][1:7])
huafb12C$Re[which(huafb12C$Mod == "FF")][38:40] = mean(huafb12C$Re[which(huafb12C$Mod == "FF")][1:7])
huafb12C$sd[which(huafb12C$Mod == "FF")][8] = mean(huafb12C$sd[which(huafb12C$Mod == "FF")][1:7])
huafb12C$sd[which(huafb12C$Mod == "FF")][38:40] = mean(huafb12C$sd[which(huafb12C$Mod == "FF")][1:7])

# huafb = rbind(huafb,huafb2w, huafb12C)

P1 <-
  ggplot(huafb, aes(Depth, Re*100))+
  geom_ribbon(aes(ymin = (Re-sd)*100, ymax = (Re+sd)*100, fill = Mod), alpha=0.3)+
  geom_hline(yintercept = 0, color = "black", linetype = 2)+
  geom_line(aes(color = Mod), size=0.55)+
  scale_x_reverse()+
  coord_flip()+
  ylab("¦¤SOC (%)")+
  xlab("Soil depth (m)")+
  # facet_grid(Tin~Cin)+
  scale_color_manual(name = "Model variates",
                     values = c("#2BA84A","#FDCA40","#00A8E8", "#E76F51"), 
                     label = c("M","M+PE","M+VT","M+VT+PE"))+
  scale_fill_manual(name = "Model variates",
                     values = c("#2BA84A","#FDCA40","#00A8E8", "#E76F51"), 
                     label = c("M","M+PE","M+VT","M+VT+PE"))+
  theme_bw()+
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        legend.position = "bottom",
        legend.direction =  "horizontal",
        axis.text = element_text(color = "black"))
P2 <-
  ggplot(huafb2w, aes(Depth, Re*100))+
  geom_ribbon(aes(ymin = (Re-sd)*100, ymax = (Re+sd)*100, fill = Mod), alpha=0.3)+
  # geom_hline(yintercept = 0, color = "black", linetype = 2)+
  geom_line(aes(color = Mod), size=0.5)+
  scale_x_reverse()+
  coord_flip()+
  ylab("¦¤SOC (%)")+
  xlab("Soil depth (m)")+
  # facet_grid(Tin~Cin)+
  scale_color_manual(name = "Model variates",
                     values = c("#2BA84A","#FDCA40","#00A8E8", "#E76F51"), 
                     label = c("M","M+PE","M+VT","M+VT+PE"))+
  scale_fill_manual(name = "Model variates",
                    values = c("#2BA84A","#FDCA40","#00A8E8", "#E76F51"), 
                    label = c("M","M+PE","M+VT","M+VT+PE"))+
  theme_bw()+
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        legend.position = "bottom",
        legend.direction =  "horizontal",
        axis.text = element_text(color = "black"))
P3 <-
  ggplot(huafb12C, aes(Depth, Re*100))+
  geom_ribbon(aes(ymin = (Re-sd)*100, ymax = (Re+sd)*100, fill = Mod), alpha=0.3, show.legend = T)+
  # geom_hline(yintercept = 0, color = "black", linetype = 2)+
  geom_line(aes(color = Mod), size=0.5)+
  scale_x_reverse()+
  coord_flip()+
  ylab("¦¤SOC (%)")+
  xlab("Soil depth (m)")+
  # facet_grid(Tin~Cin)+
  scale_color_manual(name = "Model variates",
                     values = c("#2BA84A","#FDCA40","#00A8E8", "#E76F51"), 
                     label = c("M","M+PE","M+VT","M+VT+PE"))+
  scale_fill_manual(name = "Model variates",
                    values = c("#2BA84A","#FDCA40","#00A8E8", "#E76F51"), 
                    label = c("M","M+PE","M+VT","M+VT+PE"))+
  theme_bw()+
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        legend.position = "bottom",
        legend.direction =  "horizontal",
        axis.text = element_text(color = "black"))
plot_grid(P2, P1,  P3, ncol = 3, labels = c("a","b","c"))



  
   
