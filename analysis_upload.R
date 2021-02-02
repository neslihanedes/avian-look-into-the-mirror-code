### analysis of Mirror Experiment
# created by Hiroshi Matsui, PhD 

# set configration
rm(list=ls())
library(tidyverse)
library(rstan)
library(lme4)
library(data.table)
library(car)
library(multcomp)
library(umap)

# import data
setwd(paste(dirname(rstudioapi::getActiveDocumentContext()$path),"/rawdata", sep=""))

d = fread("d.csv", data.table = FALSE)
d_control = fread("d_control.csv")
d_gencam = fread("d_gencam.csv")

# the first two seconds were removed
d = d[d$time > 2,] 
d_control = d_control[d_control$time > 2,] 
d_gencam = d_gencam[d_gencam$time > 2,] 


TM = theme(legend.position = "none", 
           axis.text = element_text(size=rel(3), colour="black", family="Times New Roman"),
           panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), strip.text = element_text(size=rel(2), family = "Times New Roman"),
           panel.border = element_rect(fill = NA, size=1)
)




#---------------------------------------------------#
### Result 1 : latency to the feeder ###
#---------------------------------------------------#
setwd(paste(dirname(rstudioapi::getActiveDocumentContext()$path), "/rawdata", sep=""))
colPallet = c("darkgreen", "blue", "red")



latency = fread("latency.csv", data.table = FALSE)
latency = latency[!is.na(latency$mirror),] 
latency[is.na(latency$exp),]$exp = 2
latency = latency[latency$exp == 1,]

latency = tidyr::gather(latency, "cond", "latency", -id, -exp, na.rm=TRUE)
latency$cond_exp = interaction(latency$exp, latency$cond, drop = TRUE)
latency$id_cond_exp = interaction(latency$id, latency$cond, latency$exp, drop = TRUE)
latency$session = rep(1:6, 24)





latency.mean = split(latency, latency$cond) %>% lapply(., function(d){
  data.frame(latency = tapply(d$latency, d$session, mean),
             se = tapply(d$latency, d$session, sd)/sqrt( length(unique(d$id))  ),
             session = 1:6, cond = d$cond[1])
}) %>% do.call(rbind,.)


### Figrue 2a ###
p = ggplot(latency.mean) + 
  geom_errorbar(aes(x=session, ymin = latency - se, ymax = latency + se, colour = cond), lwd=2, width = .3) + 
  geom_line(aes(x=session, y=latency, colour = cond), lwd=1.5) + 
  geom_point(aes(x=session, y=latency), size=7) + 
  geom_point(aes(x=session, y=latency, colour = cond), size=5) + 
  scale_colour_manual(values = colPallet) + xlab(NULL) + ylab(NULL) + 
  scale_x_continuous(limits=c(0,6.2), breaks=1:6) +
  TM
p

# statistical analysis
latency$fsession = as.factor(latency$session)
mod = lmer(log(latency) ~ cond * fsession + (1|id), data=latency)
Anova(mod)
effectsize::eta_squared(mod)

# post-hoc analysis for interaction
mod1 = lmer(log(latency) ~ cond  + (1|id), data=latency[latency$session==1,])
mod2 = lmer(log(latency) ~ cond  + (1|id), data=latency[latency$session==2,])
mod3 = lmer(log(latency) ~ cond  + (1|id), data=latency[latency$session==3,])
mod4 = lmer(log(latency) ~ cond  + (1|id), data=latency[latency$session==4,])
mod5 = lmer(log(latency) ~ cond  + (1|id), data=latency[latency$session==5,])
mod6 = lmer(log(latency) ~ cond  + (1|id), data=latency[latency$session==6,])

# and multiple comparison 
library(multcomp)
summary(glht(mod1, linfct=mcp(cond= "Tukey")))
summary(glht(mod2, linfct=mcp(cond= "Tukey")))
summary(glht(mod3, linfct=mcp(cond= "Tukey")))
summary(glht(mod4, linfct=mcp(cond= "Tukey")))
summary(glht(mod5, linfct=mcp(cond= "Tukey")))
summary(glht(mod6, linfct=mcp(cond= "Tukey")))

# difference was observed only in the first session 


#----------------------------------------------------------------------------------------------------#
#---------------------------------------------------#
### Result 2 : total number of pecks
#---------------------------------------------------#

setwd(paste(dirname(rstudioapi::getActiveDocumentContext()$path), "/rawdata", sep=""))
peck = fread("totalpeck.csv", data.table = FALSE)
peck = split(peck, peck$id) %>% lapply(., function(d){
  d$day = 1:nrow(d)
  return(d)
}) %>% do.call(rbind,.)
peck = peck[peck$exp==1,]

peck = tidyr::gather(peck, "cond", "total", -id, -exp, -day, na.rm=TRUE)

d.peck = split(peck, peck$cond) %>% lapply(., function(d){
  split(d, d$cond) %>% lapply(., function(d){
    data.frame(mean = tapply(d$total, d$day, mean), 
               se = tapply(d$total, d$day, sd)/ sqrt(length(unique(d$id))),
               
               cond = d$cond[1], session=1:6)
  } ) %>% do.call(rbind,.)
}) %>% do.call(rbind,.)


peck$fsession = as.factor(peck$day)
mod = lmer(total ~ cond * fsession + (1|id), data=peck)
Anova(mod)
mod = lmer(total ~ cond + fsession + (1|id), data=peck)

summary(glht(mod, linfct=mcp(cond= "Tukey")))

effectsize::eta_squared(mod)


#---------------------------------------------------#
### Result 3 : time-series analysis
#---------------------------------------------------#

d.orientation = split(d, d$cond) %>% 
  lapply(., function(d) {
    df = data.frame(cos = tapply(d$cos, d$time, mean, na.rm=TRUE), 
                    time = sort(unique(d$time)),
                    sd = tapply(d$cos, d$time, sd, na.rm=TRUE),
                    se = tapply(d$cos, d$time, sd, na.rm=TRUE)/sqrt(length(unique(d$id))),
                    cond = d$cond[1])
  }) %>% do.call(rbind,.)

d.orientation$ymin = d.orientation$cos - d.orientation$se
d.orientation$ymax = d.orientation$cos + d.orientation$se
d.orientation = d.orientation[d.orientation$time <= 600,]

d.orientation_control = d.orientation[d.orientation$cond == "control",]
d.orientation_withoutC = d.orientation[d.orientation$cond != "control",]
d.orientation_control = d.orientation_control[,-5]

lab = c("Mirror", "Unknown")
names(lab) = c("mirror", "stranger")



d.orientation.id = 
  split(d, d$id) %>% 
  lapply(., function(d) {
    split(d, d$cond) %>% lapply(., function(d){
      df = data.frame(cos = tapply(d$cos, d$time, mean, na.rm=TRUE), 
                      sd = tapply(d$cos, d$time, sd, na.rm=TRUE),
                      se = tapply(d$cos, d$time, sd, na.rm=TRUE)/sqrt(length(unique(d$session))),
                      time = sort(unique(d$time)),
                      id = d$id[1],
                      cond = d$cond[1])
      return(df)
    }) %>% do.call(rbind,.) }) %>% do.call(rbind,.)



# Figure 1c: head orientation
p = ggplot() + 
  geom_hline(yintercept = 0, lty=2) + 
  geom_ribbon(data=d.orientation_control, aes(x=time, ymin=ymin, ymax=ymax), alpha=0.3, fill=colPallet[1]) + 
  geom_line(data=d.orientation_control, aes(x=time, y=cos), lwd=1.3, alpha=0.5, colour=colPallet[1])  + 
  geom_ribbon(data=d.orientation_withoutC, aes(x=time, ymin=ymin, ymax=ymax, fill=cond), alpha=0.5) + 
  geom_line(data=d.orientation_withoutC, aes(x=time, y=cos, colour=cond), lwd=1.5)  + 
  scale_fill_manual(values=colPallet[-1]) + 
  scale_colour_manual(values=colPallet[-1]) + 
  scale_y_continuous(limit=c(-.7, 1.3), breaks = c(-.5, 0, .5, 1)) + 
  facet_grid(cond~., labeller = labeller(cond = lab) ) + xlab(NULL) + ylab(NULL) + TM 
p




get_d_orientation = function(LEN = 10){
  SEQ = seq(min(d$time), max(d$time), length=LEN)
  SEQ = SEQ[-c(1, length(SEQ))]
  
  d.orientation = split(d, d$ser) %>% 
    lapply(.,function(d){
      mcos = NULL
      for(i in 1:length(SEQ)){
        if(i == 1){
          mcos = c(mcos, mean(d[d$time <= SEQ[i], ]$cos, na.rm=TRUE) )
        }else if(i == length(SEQ)){
          mcos = c(mcos, mean(d[d$time > SEQ[i], ]$cos, na.rm=TRUE) )
        }else{
          mcos = c(mcos, mean(d[d$time > SEQ[i-1] & d$time <= SEQ[i], ]$cos, na.rm=TRUE) )
        }
      }
      df = data.frame(time = SEQ, cos = mcos, ser = d$ser[1], cond = d$cond[1], id = d$id[1])
      return(df)
    }) %>% do.call(rbind,.)
  
  d.orientation = d.orientation[!is.nan(d.orientation$cos),]
  return(d.orientation)
}


# fitting
d.orientation = get_d_orientation(12)
rstan_options(auto_write=TRUE)
options(mc.cores = parallel::detectCores())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

fit = stan("stan_exp1.stan", data=list(N=nrow(d.orientation),
                                       N_time=length(unique(d.orientation$time)),
                                       N_id = length(unique(d$id)),
                                       id = as.integer(as.factor(d.orientation$id)),
                                       y=d.orientation$cos,
                                       cond=as.integer(d.orientation$cond)-1,
                                       times=as.integer(as.factor(d.orientation$time)),
                                       upperBound = max(as.integer(as.factor(d.orientation$time)))),
           warmup = 350000, iter=400000, chain=4, thin=100)
summary(fit)$summary[,"Rhat"][which(summary(fit)$summary[,"Rhat"] >= 1.1)]
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
saveRDS(fit, file = "stan_exp1.obj")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
fit = readRDS("stan_exp1.obj") 

# analysis on fitting
res = rstan::extract(fit)
TM2 =  theme(plot.margin= unit(c(.1, 1, 0, 0), "lines"))
d.stan = data.frame(time = sort(unique(d.orientation$time)),
                    lower = apply(res$controlVSmirror, 2, function(p) quantile(p, .025)),
                    upper = apply(res$controlVSmirror, 2, function(p) quantile(p, .975)),
                    mean = apply(res$controlVSmirror, 2, mean))

p1 = ggplot(d.stan) + 
  geom_ribbon(aes(x=time, ymin=lower, ymax=upper), alpha=.3, fill=colPallet[2]) + 
  geom_line(aes(x=time, y=mean), alpha=1, colour=colPallet[2])+ TM + TM2 + 
  xlab(NULL) + ylab(NULL) + 
  scale_y_continuous(limits = c(0, .45), breaks = c(0, .2, .4))

d.stan = data.frame(time = sort(unique(d.orientation$time)),
                    lower = apply(res$controlVSstranger, 2, function(p) quantile(p, .025)),
                    upper = apply(res$controlVSstranger, 2, function(p) quantile(p, .975)),
                    mean = apply(res$controlVSstranger, 2, mean)) 

p2 = ggplot(d.stan) + 
  geom_ribbon(aes(x=time, ymin=lower, ymax=upper), alpha=.3, fill=colPallet[3]) + 
  geom_line(aes(x=time, y=mean), alpha=1, colour=colPallet[3])  + TM + TM2 + 
  xlab(NULL) + ylab(NULL) +
  scale_y_continuous(limits = c(0, .45), breaks = c(0, .2, .4))



d.stan = data.frame(time = sort(unique(d.orientation$time)),
                    lower = -apply(res$mirrorVSstranger, 2, function(p) quantile(p, .025)),
                    upper = -apply(res$mirrorVSstranger, 2, function(p) quantile(p, .975)),
                    mean = -apply(res$mirrorVSstranger, 2, mean))
p3 = ggplot(d.stan) + 
  geom_ribbon(aes(x=time, ymin=lower, ymax=upper), alpha=.3, fill=colPallet[3]) + 
  geom_line(aes(x=time, y=mean), alpha=1, colour=colPallet[3]) +
  scale_y_continuous(limits = c(0, .45), breaks = c(0, .2, .4)) + TM + TM2 + 
  xlab(NULL) + ylab(NULL) 



# figure 1d
p = gridExtra::grid.arrange(p1,p2,p3) 
p



quantile(res$diff[,1], c(0.025, .975))
quantile(res$diff[,2], c(0.025, .975))
quantile(res$mirrorVSstranger[,10], c(0.025, .975))


#---------------------------------------------------#
### Result 4 : activity rates
#---------------------------------------------------#
d$feeding = 0
d[d$distance_to_food < 100,]$feeding = 1

# d$velocity = d$velocity *60
# d$velocity = d$velocity %>% log
# get mean velocity
df.act = split(d, d$ser) %>% lapply(., function(d){
  d = d[d$feeding == 0,]
  d = d[d$head_likeli > .95 & d$body_likeli > .95,]
  df = data.frame(act = mean(d$velocity, na.rm=TRUE),
                  cond = d$cond[1],
                  id = d$id[1])
  return(df)
}) %>% do.call(rbind,.)
df.act = df.act[!is.na(df.act$act),]

df.act.vis = split(df.act, df.act$cond) %>% lapply(., function(d){
  split(d, d$id) %>% lapply(., function(d){
    data.frame(cond = d$cond[1], id = d$id[1], act = mean(d$act))
  }) %>% do.call(rbind,.)
}) %>% lapply(., function(d){
  data.frame(cond = d$cond[1], act = mean(d$act), se = sd(d$act)/sqrt(nrow(d)))
}) %>% do.call(rbind,.)

df.act.id.mean = split(df.act, df.act$id) %>% lapply(., function(d){
  data.frame(act = tapply(d$act, d$cond, mean), id = d$id[1], cond = c( "control"  ,  "mirror" , "stranger" ))
}) %>% do.call(rbind,.)
df.act.id.mean$cond = factor(as.character(df.act.id.mean$cond), levels = c("mirror", "stranger", "control"))

WTH = 0
p = ggplot() + 
  geom_bar(data= df.act.vis, aes(x=cond, y=act), stat="identity", colour="black", lwd=2, fill="white") + 
  geom_bar(data= df.act.vis, aes(x=cond, y=act, fill=cond), stat="identity", colour="black", alpha = .5) + 
  xlab(NULL) + ylab(NULL) + scale_fill_manual(values=c(colPallet[2:3], colPallet[1])) + 
  scale_colour_manual(values=c(colPallet[2:3], colPallet[1])) + 
  geom_errorbar(data = df.act.vis, aes(x=cond, ymin=act-se, ymax=act+se), width=.1, lwd = 1.5) + 
  geom_line(data = df.act.id.mean, aes(x = cond, y= act, group=id), lwd = 2, alpha = .5, position = position_dodge(width = WTH)) + 
  geom_point(data = df.act.id.mean, 
             aes(x = cond, y= act, group = id), size = 8, position = position_dodge(width = WTH)) + 
  geom_point(data = df.act.id.mean, 
             aes(x = cond, y= act, colour=cond, group = id), size = 5, position = position_dodge(width = WTH)) + 
  scale_x_discrete(labels = c("mirror", "unknown", "solitude")) + 
  scale_y_continuous(expand=c(0,0), limit = c(0, 2.4)) + TM + 
  theme(plot.margin= unit(c(.7, 1, .1, 1), "lines"))

p

ggplot(df.act) + 
  geom_boxplot(aes(x=cond, y=log(act), fill=cond)) + 
  scale_fill_manual(values = colPallet) + TM +
  xlab(NULL) + ylab(NULL) + ylim(c(-2, 3))

# applying linear mixed model
mod = lmer(log(act) ~ cond + (1|id) , data=df.act)
Anova(mod)
effectsize::eta_squared(mod)

mod = lmer(log(act) ~ cond + (1|id) - 1, data=df.act)
Anova(mod)
summary(glht(mod, linfct=mcp(cond= "Tukey")))

confint(mod, method = "Wald")