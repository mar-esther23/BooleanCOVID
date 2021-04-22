setwd("~/Downloads/Boolean_COVID")

library(BoolNet)
library(BoolNetPerturb)

net <-  loadNetwork("data_raw/macrophage_net_210421.csv")
lab <- read.csv("data_raw/macrophage_label_210421.csv")
env <- read.csv("data_raw/macrophage_environment.csv")

print("WT")
attr <- getAttractors(net, method="sat.restricted", maxAttractorLength=2)
labels <- labelAttractors(attr, lab, net$genes)
print(table(labels))

run.net <- function(net, condition) {
  print(env[condition,])
  net_ <- fixGenes(net, colnames(env), env[condition,])
  attr_ <- getAttractors(net_)
  labels_ <- labelAttractors(attr_, lab, net_$genes)
  attr_ <- attractorToDataframe(attr_)
  attr_$involvedStates <- labels_
  print(aggregate(. ~ involvedStates, data=attr_, FUN=sum))
  cfm_ = cellFateMap(net_, label.rules=lab)
  print(table(cfm_[,c('initial','final')]))
}

for (condition in rownames(env)) {
  run.net(net, condition)
}

condition = "CoV-sev"
print(env[condition,])
net_ <- fixGenes(net, colnames(env), env[condition,])
attr_ <- getAttractors(net_)
labels_ <- labelAttractors(attr_, lab, net_$genes)
attr_ <- attractorToDataframe(attr_)
attr_$involvedStates <- labels_
print(aggregate(. ~ involvedStates, data=attr_, FUN=sum))
cfm_ = cellFateMap(net_, label.rules=lab)
print(cfm_[cfm_$final=="M1",])
