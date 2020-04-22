####################################
# Figure and some estimations (ergm, bergm, tergm)
####################################


#-------------
# load packages
#-------------

library(dplyr) # for laziness
library(countrycode) # for data prep (working with countries)
library(network) # to prepare networks for printing
library(ggplot2) # dependency
library(GGally) # ggnet is part of GGally
#devtools::install_github("briatte/ggnet")
library(ggnet) # Network visualization
library(ergm) # for simple network analysis
library(Bergm) #for baysian models
library(btergm) #for temporal models


#-----------
# load data
#-----------

dt <- read.csv("https://www.designoftradeagreements.org/media/filer_public/dd/0d/dd0d0c0e-6489-4120-9d09-242a9b9c52c6/list_of_treaties_01_05_dyads.csv") 

dt <- dt %>%
  mutate(iso1c=countrycode(iso1,"iso3n","iso3c"),
         iso2c=countrycode(iso2,"iso3n","iso3c")) %>%
  mutate(iso1c=ifelse(country1=="Kosovo","KSV",as.character(iso1c)),
         iso2c=ifelse(country1=="Kosovo","KSV",as.character(iso2c))) %>%
  select(iso1c,iso2c,year) %>% mutate(pta=1)

dt_2020 <- dt %>% dplyr::select(iso1c,iso2c) %>% unique

#-----------
# create network
#----------
units <- unique(c(dt_2020$iso1c,dt_2020$iso2c))
m <- matrix(0,nrow=length(units),ncol=length(units))
rownames(m) <- units
colnames(m) <- units

for(i in 1:nrow(dt_2020)){
  c <- which(colnames(m)==dt_2020$iso1c[i]) #columns are incoming nodes
  r <- which(row.names(m)==dt_2020$iso2c[i]) #rows are outgoing nodes
  m[r,c] <- 1
}

net <- network(m,directed=FALSE)

# alternative libraries for plotting networks are igraph, networkd3, visNetwork,...

ggnet2(net,mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.01),
       label=TRUE,
       label.color = "black",
       label.size = 1.5,
       
       edge.color = "gray",
       edge.size=0.2,
       
       node.color="cadetblue1",
       
       size = "degree",
       legend.position = "none")

#ggsave("./figs/pta_net.png")

#--------
# ERGM
#--------
unit <- letters
unit_info <- expand.grid('unit1' = unit, 'unit2' = unit, stringsAsFactors = FALSE) 
unit_info$e <- sample(c(1,0),nrow(unit_info),replace=TRUE)
unit_info$z <- rnorm(nrow(unit_info))
unit_info$x <- rnorm(nrow(unit_info))

head(unit_info)

m <- matrix(0,nrow=length(unit),ncol=length(unit))
rownames(m) <- unit
colnames(m) <- unit

for(i in 1:nrow(dt_2020)){
  c <- which(colnames(m)==unit_info$unit1[i]) #columns are incoming nodes
  r <- which(row.names(m)==unit_info$unit2[i]) #rows are outgoing nodes
  m[r,c] <- unit_info$e[i]
}

net <- network(m,directed=TRUE)
set.vertex.attribute(net, 'x', unit_info$x)
set.vertex.attribute(net, 'z', unit_info$z)

?ergm::`ergm-terms`


fit0 <- ergm(net ~ edges) #corrects for how tense the network is
summary(fit0)

summary(net ~ triangle ) #check number of triangles

fit1 <- ergm(net ~ edges + triangle + nodecov('x'))
summary(fit1)

mcmc.diagnostics(fit1) 
# All values from your sampling procedure should form a bell-curve. 
# If you have multiple peaks for a variable, this variable did not converge and the variable makes the model degenerate. 

gof.model1 <- btergm::gof(fit1, nsim = 50,parallel = "multicore",ncpus=5)
gof.model1
plot(gof.model1)
#--------
# BERGM
#--------
fit2 <- bergm(net ~ edges + triangle + nodecov('x')) #bergmM would be the option for missing data!
summary(fit2)
plot(fit2)
bgof(fit2)
#--------
# TERGM
#--------
#devtools::install_github('vincentarelbundock/SpatialHelper')
library(SpatialHelper)

unit <- letters
time <- 1:10
unit_time <- expand.grid('unit' = unit, 'time' = time, stringsAsFactors = FALSE) 
unit_time$x <- rnorm(nrow(unit_time))
unit_time$k <- rnorm(nrow(unit_time))
head(unit_time)


dyad_time <- expand.grid('unit1' = unit, 'unit2' = unit, 'time' = time, stringsAsFactors = FALSE) 
dyad_time$w <- rnorm(nrow(dyad_time))
dyad_time$e <- rnorm(nrow(dyad_time))
dyad_time$z <- as.numeric(ifelse(dyad_time$w + dyad_time$e > 0, 1, 0))

head(dyad_time)

# Convert the panel data to network data inside an environment
env <- panel_to_network(unit_time, dyad_time)

# Identify the dependent network
env <- dependent_network('z', env)

?btergm::`ergm-terms`

attach(env)
fit3 <- btergm(z  ~ edges + twopath + nodecov('x') + edgecov(w) + istar(2)+memory(type="stability",lab=1), R = 500)
detach(env)
summary(fit3)

# Model fit
attach(env)
gof.model3 <- btergm::gof(fit3, nsim = 50)
detach(env)
gof.model3
plot(gof.model3)


