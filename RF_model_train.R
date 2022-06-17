##########################################################
########RF model, importance, Mapping for VT and PE#######
##########################################################


rm(list = ls())
library(lattice)
library(ggplot2)
library(caret)
library(ranger)
library(e1071)
library(tidyverse)
library(rpart)
library(C50)
library(reshape)
library(car)
library(gtools)
library(foreach)
library(doParallel)
library(doSNOW)
library(dplyr)
library(cowplot)

dir = "xxx/"
filenames = list.files(dir)

###Vertical carbon tronsport
dat.3cp = readRDS("xxxx/Result_3CP_DM_PE_PARS_TT.Rdata")
VT.3cp = data.frame(dat.3cp[[3]])
VT.3cp = (t(t(VT.3cp)*c(0.2,0.2,0.2,0.2,0.2,0.5,0.5)))
VT.3cp = (t(t(VT.3cp)/c(4,4,4,4,4,10,10)))
names(VT.3cp) = paste0("L", 1:7)

####Priming effect and saturation index
PE.3cp = data.frame(dat.3cp[[2]])
PEmax = data.frame(dat.3cp[[7]])$X8+1


## environmental factors
out.input = c(1171, 1742)
dir.env = "xxxx/Envfactor/"
file.env = mixedsort(list.files(dir.env))
envfactor = lapply(1:length(file.env), function(x)read.csv(paste0(dir.env, file.env[x]))[-out.input, ])
env.VT = lapply(1:length(envfactor), function(x){
  data.frame(VT3cp = VT.3cp[,x], PE.3cp =PE.3cp[,x] ,PE_satindex = PE_satindex[,x],envfactor[[x]][,-1])
})

### select variables
df =env.VT[[1]]
df$biomess = factor(df$biomess)
impVars = df[,!names(df)%in%c("X","VT3cp","PE.3cp","npp","biomess","PE_satindex",
                              "long","lat", "CFRAG", "ORGC", "BULK")]

form_all <- as.formula(VT3cp~SDTO+STPC+CLPC+
                         TAWC+CECS+CECc+SDTO+
                         TOTN+CNrt+PHAQ+ELCO+
                         bio1+bio2+
                         bio3+bio4+bio5+
                         bio6+bio7+bio8+
                         bio9+bio10+bio11+
                         bio12+bio13+bio14+
                         bio15+bio16+bio17+
                         bio18+bio19)

all_vifs_test <- try(car::vif(lm(form_all, data = df)), silent=TRUE)
if (class(all_vifs_test) == "try-error"){
  lm_alias <- alias(lm(form_all, data=df))
  broken <- data.frame(lm_alias$Complete)
  broken_var <- row.names(broken)
  nam_var <- names(impVars) %in% broken_var
  impVars <- impVars[!nam_var]
  form_all_new <- as.formula(paste("VT3cp ~ ",
                                   paste(names(impVars), collapse=" + "),
                                   sep=""))
  all_vifs <- car::vif(lm(form_all_new, data=df))
} else {
  all_vifs <- all_vifs_test
}

if(any(all_vifs>10)){
  all_vifs <- as.data.frame(all_vifs)
  while((nrow(all_vifs) > 2)& (max(all_vifs[, 1]) > 5)  &
        (class(all_vifs) != "try-error")) {
    remove_var <- rownames(all_vifs)[which(all_vifs[, 1] == max(all_vifs[, 1]))]
    impVars <- impVars[!names(impVars) %in% remove_var]
    fullForm <- paste ("VT3cp ~ ", paste (names(impVars), collapse=" + "), sep="")
    fullMod <- lm(as.formula(fullForm), data=df)
    all_vifs <- try(as.data.frame(car::vif(fullMod)), silent=TRUE)
  }
  vif_filtered_variables <- names(fullMod$model)[!names(fullMod$model) %in%
                                                   "VT3cp"]
} else {
  all_vifs <- as.data.frame(all_vifs)
  vif_filtered_variables <- rownames(all_vifs)[!names(all_vifs) %in% "VT3cp"]
}


#### train randforest model
# my tuneGrid object:
tgrid <- expand.grid(
  mtry = 2:4,
  splitrule = c("variance", "extratrees"),
  min.node.size = 2:6
)

#10 fold cross
control <- trainControl(method = "cv",number = 10,repeats = 3)
dir.mod = "I:/SOC-VT-PE-Model_revised/Mod_RF_VT_PE_PESI/"

cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoSNOW(cl)

r2 = array(NA, dim = c(7,3))
for(h in 1:length(env.VT)){
  
  ###VT for three_pool model
  dat.3cp = env.VT[[h]][, names(env.VT[[h]]) %in% c("VT3cp", Vars_VT)]
  names(dat.3cp) = c("VT","SDTO", "CLPC","TAWC", "CECc", "ELCO", 
                        "Bio1",  "Bio2",  "Bio4", 
                      "Bio5",  "Bio8",  "Bio12", "Bio14", "Bio15", "Bio18",
                     "Bio19")
  
  corcoef = apply(dat.3cp[,-1], 2, function(x)cor(x, dat.3cp$VT))
  corcoef.sign = corcoef
  corcoef.sign[corcoef.sign < 0 ] = "-"
  corcoef.sign[corcoef.sign > 0 ] = "+"

  model.3cp <- train(y = dat.3cp$VT,
                     x = dat.3cp[, !names(dat.3cp) %in%  "VT"],
                     method = "ranger",
                     trControl = control,
                     tuneGrid = tgrid,
                     importance = 'impurity',
                     keep.inbag=TRUE)
  
  importance.3cp = varImp(model.3cp)$importance
  importance.3cp$R2 = model.3cp$finalModel$r.squared
  importance.3cp$Soillayer = paste0("L",h)
  importance.3cp$corcoef = corcoef
  importance.3cp$corcoef.sign = corcoef.sign
  importance.3cp$soilall = sum(importance.3cp$Overall[1:8])/sum(importance.3cp$Overall)
  
  write.csv(importance.3cp, paste0(dir.mod, "Importance/Importance_3cp.VT_L", h, ".csv"))
  saveRDS(model.3cp,paste0(dir.mod,"Mod_VT/RF_3cp.VT_L",h,".Rdata"))
  r2[h,1] = model.3cp$finalModel$r.squared
  
  ####PE model
  dat.3cp.PE = env.VT[[h]][, names(env.VT[[h]]) %in% c("PE.3cp", Vars.PE)]
  names(dat.3cp.PE) = c("PE","CNrt",  "SDTO", "CLPC",
                        "TAWC", "CECc", "ELCO",   "PHAQ", "Bio1",  "Bio2",  "Bio4", 
                        "Bio5",  "Bio8",  "Bio12", "Bio14", "Bio15", "Bio18",
                        "Bio19")

  corcoef = apply(dat.3cp.PE[,-1], 2, function(x)cor(x, dat.3cp.PE$PE))
  corcoef.sign = corcoef
  corcoef.sign[corcoef.sign < 0 ] = "-"
  corcoef.sign[corcoef.sign > 0 ] = "+"

  model.3cp.PE <- train(y = dat.3cp.PE$PE,
                        x = dat.3cp.PE[, !names(dat.3cp.PE) %in%  "PE"],
                        method = "ranger",
                        trControl = control,
                        tuneGrid = tgrid,
                        importance = 'impurity',
                        keep.inbag=TRUE)

  importance.3cp.pe = varImp(model.3cp.PE)$importance
  importance.3cp.pe$R2 = model.3cp.PE$finalModel$r.squared
  importance.3cp.pe$Soillayer = paste0("L",h)
  importance.3cp.pe$corcoef = corcoef
  importance.3cp.pe$corcoef.sign = corcoef.sign
  importance.3cp.pe$soilall = sum(importance.3cp.pe$Overall[1:8])/sum(importance.3cp.pe$Overall)
  
  write.csv(importance.3cp.pe, paste0(dir.mod, "Importance/Importance_3cp.PE_L", h, ".csv"))
  write_rds(model.3cp.PE,paste0(dir.mod,"Mod_PE/RF_3cp.PE_L",h,".rds"))
  r2[h,2] = model.3cp.PE$finalModel$r.squared
  
  # 
  ####PESI
  dat.mc = env.VT[[h]][, names(env.VT[[h]]) %in% c("PE_satindex", Vars.PE)]
  names(dat.mc) = c("PESI","CNrt",  "SDTO", "CLPC",
                    "TAWC", "CECc", "ELCO",   "PHAQ", "Bio1",  "Bio2",  "Bio4", 
                    "Bio5",  "Bio8",  "Bio12", "Bio14", "Bio15", "Bio18",
                    "Bio19")
  
  corcoef = apply(dat.mc[,-1], 2, function(x)cor(x, dat.mc$PESI))
  corcoef.sign = corcoef
  corcoef.sign[corcoef.sign < 0 ] = "-"
  corcoef.sign[corcoef.sign > 0 ] = "+"

  model.mc <- train(y = dat.mc$PESI,
                    x = dat.mc[, !names(dat.mc) %in%  "PESI"],
                    method = "ranger",
                    trControl = control,
                    tuneGrid = tgrid,
                    importance = 'impurity',
                    keep.inbag=TRUE)

  importance.mc = varImp(model.mc)$importance
  importance.mc$R2 = model.mc$finalModel$r.squared
  importance.mc$Soillayer = paste0("L",h)
  importance.mc$corcoef = corcoef
  importance.mc$corcoef.sign = corcoef.sign
  importance.mc$soilall = sum(importance.mc$Overall[1:8])/sum(importance.mc$Overall)
  
  
  write.csv(importance.mc, paste0(dir.mod, "Importance/Importance_PESI_L", h, ".csv"))
  write_rds(model.mc,paste0(dir.mod,"Mod_PESI/RF_PESI_L",h,".rds"))
  r2[h,3] = model.mc$finalModel$r.squared

  print(h)
}
stopCluster(cl)

##### plot the importance of environment factor
library(circlize)
library(cowplot)
library(RColorBrewer)

dir.mod = "I:/SOC-VT-PE-Model_revised/Mod_RF_VT_PE_PESI/"
file.imp = list.files(paste0(dir.mod, "Importance/"))

dat = lapply(1:length(file.imp),function(x){
  dat = read.csv(paste0(dir.mod, "Importance/", file.imp[x]))
  return(dat)
})


# soil_layer<-c("0-0.2", "0.2-0.4", "0.4-0.6","0.6-0.8", "0.8-1", "1-1.5","1.5-2") 
soil_layer = paste0("L", 1:7)
env.VT_names <-  c("SDTO", "CLPC",    "TAWC",  "CECc",  "ELCO",  "Bio1",  "Bio2",  "Bio4", 
                "Bio5",  "Bio8",  "Bio12", "Bio14", "Bio15", "Bio18",
                "Bio19")
env.PE_names <-  c("CNrt",  "SDTO","CLPC", "TAWC",  "CECc",  "ELCO",  "PHAQ",  "Bio1",  "Bio2",  "Bio4", 
                   "Bio5",  "Bio8",  "Bio12", "Bio14", "Bio15", "Bio18",
                   "Bio19")
factors.VT = c("R2", "Rsoil",env.VT_names, "Rclim","Layers")
factors.PE = c("R2", "Rsoil",env.PE_names, "Rclim","Layers")



##Importance of VT
# nf <- layout(matrix(c(1, 2, 3),1), widths = c(1,1,1), heights = c(1,1,1))
# layout.show(nf)
###3cp model

par(mar = c(1.5, 1.5, 1.5, 1.5), pin = c(4.5,4.5))
circos.clear()
circos.initialize(factors.VT, xlim = c(0, 1))
for(i in 1:9) {
  circos.track(ylim = c(0, 1), track.height=0.05)
}


for (j in 1:length(factors.VT)) {
  circos.update(sector.index = factors.VT[j], track.index = 2, bg.col = NA, bg.border = NA) 
  circos.update(sector.index = factors.VT[j], track.index = 1, bg.col = NA, bg.border = NA)
  circos.text(0.5, 0.7,factors.VT[j], sector.index =factors.VT[j] , track.index = 2)
}

# circos.update(sector.index = factors[1], track.index = 2, bg.col = NA, bg.border = NA) 
# circos.update(sector.index = factors[1], track.index = 1, bg.col = NA, bg.border = NA)
# circos.text(0.5, 0.7,expression(italic(R)^"2"), sector.index =factors[1] , track.index = 2)

for (i in 3:9) {
  circos.update( sector.index =factors.VT[1], track.index = i, bg.col = NA, bg.border = NA) 
  circos.update(sector.index = factors.VT[length(factors.VT)], track.index = i, bg.col = NA, bg.border = NA)
  circos.text(0.5,0.5,round(unique(dat[[i+7-2]]$R2),2), sector.index =factors.VT[1] , track.index = i, cex=0.8)
  circos.text(0.5,0.5,soil_layer[i-2], sector.index =factors.VT[length(factors.VT)] , track.index = i, cex=0.8)
  
  circos.update( sector.index =factors.VT[2], track.index = i, bg.col = NA, bg.border = NA) 
  circos.update(sector.index = factors.VT[length(factors.VT)-1], track.index = i, bg.col = NA, bg.border = NA)
  circos.text(0.5,0.5,round(unique(dat[[i+7-2]]$soilall),1), sector.index =factors.VT[2] , track.index = i, cex=0.8)
  circos.text(0.5,0.5,1-round(unique(dat[[i+7-2]]$soilall),1), sector.index =factors.VT[length(factors.VT)-1] , track.index = i, cex=0.8)
}


# highlight.sector(factors.VT[2:3], track.index = 1, text = "Vegetation",col = "forestgreen",
#                  facing = "bending.inside", niceFacing = TRUE, text.vjust = "4mm")

highlight.sector(factors.VT[2:7], track.index = 1, text = "Soil", col = "#B2967D",
                 facing = "bending.inside", niceFacing = TRUE, text.vjust = "4mm")

highlight.sector(factors.VT[8:18], track.index = 1, text = "Climate", col = "#A4DEF9",
                 facing = "bending.inside", niceFacing = TRUE, text.vjust = "4mm")

min.imp = min(sapply(1:7+7, function(x)min(dat[[x]]$Overall)))
max.imp = max(sapply(1:7+7, function(x)max(dat[[x]]$Overall)))
col_fun1 = colorRamp2(seq(min.imp, max.imp, (max.imp-min.imp)/8), brewer.pal(9, "OrRd"))


for (i in 3:length(factors.VT)) {
  for (j in 1:7) {
    highlight.sector(factors.VT[i], col = col_fun1(dat[[j+7]]$Overall[i-2]*unique(dat[[j+7]]$R2)),
                     track.index = j+2, border = NA, lwd = 0.1, 
                     text = dat[[j+7]]$corcoef.sign[i-2], text.vjust = 0.5,cex=0.8)
  }
}
text(-0.9,1, "VT", font = 2)

lgd = Legend(title = "Importance", col_fun = col_fun1, title_position = "topcenter")
draw(lgd, x = unit(0.5, "npc"), just="center")

###################################################################PE
par(mar = c(1.5, 1.5, 1.5, 1.5), pin = c(4.5,4.5))
circos.clear()
circos.initialize(factors.PE, xlim = c(0, 1))
for(i in 1:9) {
  circos.track(ylim = c(0, 1), track.height=0.05)
}

LL=0
for (j in 1:length(factors.PE)) {
  circos.update(sector.index = factors.PE[j], track.index = 2, bg.col = NA, bg.border = NA) 
  circos.update(sector.index = factors.PE[j], track.index = 1, bg.col = NA, bg.border = NA)
  circos.text(0.5, 0.7,factors.PE[j], sector.index =factors.PE[j] , track.index = 2)
}

# circos.update(sector.index = factors.PE[1], track.index = 2, bg.col = NA, bg.border = NA) 
# circos.update(sector.index = factors.PE[1], track.index = 1, bg.col = NA, bg.border = NA)
# circos.text(0.5, 0.7,expression(italic(R)^"2"), sector.index =factors.PE[1] , track.index = 2)

for (i in 3:9) {
  circos.update( sector.index =factors.PE[1], track.index = i, bg.col = NA, bg.border = NA) 
  circos.update(sector.index = factors.PE[length(factors.PE)], track.index = i, bg.col = NA, bg.border = NA)
  circos.text(0.5,0.5,round(unique(dat[[i+LL-2]]$R2),2), sector.index =factors.PE[1] , track.index = i, cex=0.8)
  circos.text(0.5,0.5,soil_layer[i-2], sector.index =factors.PE[length(factors.PE)] , track.index = i, cex=0.8)
  
  circos.update( sector.index =factors.PE[2], track.index = i, bg.col = NA, bg.border = NA) 
  circos.update(sector.index = factors.PE[length(factors.PE)-1], track.index = i, bg.col = NA, bg.border = NA)
  circos.text(0.5,0.5,round(unique(dat[[i+LL-2]]$soilall),1), sector.index =factors.PE[2] , track.index = i, cex=0.8)
  circos.text(0.5,0.5,1-round(unique(dat[[i+LL-2]]$soilall),1), sector.index =factors.PE[length(factors.PE)-1] , track.index = i, cex=0.8)
}


# highlight.sector(factors.PE[2:3], track.index = 1, text = "Vegetation",col = "forestgreen",
#                  facing = "bending.inside", niceFacing = TRUE, text.vjust = "4mm")

highlight.sector(factors.PE[2:9], track.index = 1, text = "Soil", col = "#B2967D",
                 facing = "bending.inside", niceFacing = TRUE, text.vjust = "4mm")

highlight.sector(factors.PE[10:20], track.index = 1, text = "Climate", col = "#A4DEF9",
                 facing = "bending.inside", niceFacing = TRUE, text.vjust = "4mm")

min.imp = min(sapply(1:7+LL, function(x)min(dat[[x]]$Overall)))
max.imp = max(sapply(1:7+LL, function(x)max(dat[[x]]$Overall)))
col_fun1 = colorRamp2(seq(min.imp, max.imp, (max.imp-min.imp)/8), brewer.pal(9, "OrRd"))


for (i in 3:length(factors.PE)) {
  for (j in 1:7) {
    highlight.sector(factors.PE[i], col = col_fun1(dat[[j+LL]]$Overall[i-2]*unique(dat[[j+LL]]$R2)),
                     track.index = j+2, border = NA, lwd = 0.1, 
                     text = dat[[j+LL]]$corcoef.sign[i-2], text.vjust = 0.5,cex=0.8)
  }
}
text(-0.9,1, "PE", font = 2)

lgd = Legend(title = "Importance", col_fun = col_fun1, title_position = "topcenter")
draw(lgd, x = unit(0.5, "npc"), just="center")

#######################################################PE saturation index
###################################################################PE
par(mar = c(1.5, 1.5, 1.5, 1.5), pin = c(4.5,4.5))
circos.clear()
circos.initialize(factors.PE, xlim = c(0, 1))
for(i in 1:9) {
  circos.track(ylim = c(0, 1), track.height=0.05)
}

LL=14
for (j in 1:length(factors.PE)) {
  circos.update(sector.index = factors.PE[j], track.index = 2, bg.col = NA, bg.border = NA) 
  circos.update(sector.index = factors.PE[j], track.index = 1, bg.col = NA, bg.border = NA)
  circos.text(0.5, 0.7,factors.PE[j], sector.index =factors.PE[j] , track.index = 2)
}

# circos.update(sector.index = factors.PE[1], track.index = 2, bg.col = NA, bg.border = NA) 
# circos.update(sector.index = factors.PE[1], track.index = 1, bg.col = NA, bg.border = NA)
# circos.text(0.5, 0.7,expression(italic(R)^"2"), sector.index =factors.PE[1] , track.index = 2)

for (i in 3:9) {
  circos.update( sector.index =factors.PE[1], track.index = i, bg.col = NA, bg.border = NA) 
  circos.update(sector.index = factors.PE[length(factors.PE)], track.index = i, bg.col = NA, bg.border = NA)
  circos.text(0.5,0.5,round(unique(dat[[i+LL-2]]$R2),2), sector.index =factors.PE[1] , track.index = i, cex=0.8)
  circos.text(0.5,0.5,soil_layer[i-2], sector.index =factors.PE[length(factors.PE)] , track.index = i, cex=0.8)
  
  circos.update( sector.index =factors.PE[2], track.index = i, bg.col = NA, bg.border = NA) 
  circos.update(sector.index = factors.PE[length(factors.PE)-1], track.index = i, bg.col = NA, bg.border = NA)
  circos.text(0.5,0.5,round(unique(dat[[i+LL-2]]$soilall),1), sector.index =factors.PE[2] , track.index = i, cex=0.8)
  circos.text(0.5,0.5,1-round(unique(dat[[i+LL-2]]$soilall),1), sector.index =factors.PE[length(factors.PE)-1] , track.index = i, cex=0.8)
}


# highlight.sector(factors.PE[2:3], track.index = 1, text = "Vegetation",col = "forestgreen",
#                  facing = "bending.inside", niceFacing = TRUE, text.vjust = "4mm")

highlight.sector(factors.PE[2:9], track.index = 1, text = "Soil", col = "#B2967D",
                 facing = "bending.inside", niceFacing = TRUE, text.vjust = "4mm")

highlight.sector(factors.PE[10:20], track.index = 1, text = "Climate", col = "#A4DEF9",
                 facing = "bending.inside", niceFacing = TRUE, text.vjust = "4mm")

min.imp = min(sapply(1:7+LL, function(x)min(dat[[x]]$Overall)))
max.imp = max(sapply(1:7+LL, function(x)max(dat[[x]]$Overall)))
col_fun1 = colorRamp2(seq(min.imp, max.imp, (max.imp-min.imp)/8), brewer.pal(9, "OrRd"))


for (i in 3:length(factors.PE)) {
  for (j in 1:7) {
    highlight.sector(factors.PE[i], col = col_fun1(dat[[j+LL]]$Overall[i-2]*unique(dat[[j+LL]]$R2)),
                     track.index = j+2, border = NA, lwd = 0.1, 
                     text = dat[[j+LL]]$corcoef.sign[i-2], text.vjust = 0.5,cex=0.8)
  }
}
text(-0.9,1, "PESI", font = 2)

lgd = Legend(title = "Importance", col_fun = col_fun1, title_position = "topcenter")
draw(lgd, x = unit(0.5, "npc"), just="center")
dev.off()

