library(statnet)
setwd('~/Documents/noordin')

# Retrieve Data (data omitted due to liscensing)
network.data <- read.csv('noordin_terrorist_friendship_network.csv',as.is=T,header=T,row.names=1)
attributes <- read.csv('Noordin Attribute Data.csv',as.is=T,header=T)

net <- network(network.data, vertex.attr=attributes,
                    vertex.attrnames=colnames(attributes),directed=F)
# Stats Function
netstats <- function(network){
  show(paste('Density: ' , gden(network)))
  show(paste('Average Degree: ', mean(degree(network, gmode='graph'))))

  suppressMessages(library(sand))
  suppressMessages(library(intergraph))

  inetwork <- asIgraph(network)

  show(paste('Average Geodesic Distance ', mean_distance(inetwork, directed=F)))
  show(paste('Diameter: ', diameter(inetwork, directed=F)))
  show(paste('Global Clustering Coefficient: ',mean(transitivity(inetwork, type='localundirected'), na.rm=T)))

  show(paste('Degree Centralization: ',centr_degree(inetwork, mode='total')$centralization))
  show(paste('Closeness Centralization: ',centr_clo(inetwork, mode='total')$centralization))
  show(paste('Betweenness Centralization: ',centr_betw(inetwork, directed=F)$centralization))

  detach('package:sand', unload=TRUE)
  detach('package:igraph', unload=TRUE)

  invisible(capture.output(num_components <- components(network)))

  show(paste('Components: ',num_components))
  show(paste('Component Ratio: ',(num_components-1)/(network.size(network)-1)))
  show(paste('Connectedness: ',connectedness(network)))
  show(paste('Transitivity Index: ',gtrans(network)))
  show('Triad Census')
  show(triad.census(network))
  show('Clique Census')
  show(clique.census(net, mode='graph', tabulate.by.vertex=FALSE, enumerate=FALSE))
}

# Goodness-of-fit Metrics (Returns Data Frame)
GOFMetrics <- function(...) {
  models <- list(...)
  AIC <- round(sapply(list(...), AIC), 0)
  BIC <- round(sapply(list(...), BIC), 0)
  df <- sapply(models, function(x) {
    frame <- as.data.frame(anova(models[[1]],x,test='chi2'))
    ifelse(tracemem(x) == tracemem(models[[1]]), frame[2,1], frame[3,1])
    })
  deviance <- sapply(models, function(x) {
    frame <- as.data.frame(anova(models[[1]],x,test='chi2'))
    ifelse(tracemem(x) == tracemem(models[[1]]), frame[2,2], frame[3,2])
    })
  chi <- sapply(models, function(x) {
    frame <- as.data.frame(anova(models[[1]],x,test='chi2'))
    ifelse(tracemem(x) == tracemem(models[[1]]), frame[2,5], frame[3,5])
    })
  formula <- sapply(models, function(x) {
    return(paste(as.character(x$formula), collapse = ' '))
  })
  show(formula)
  return(data.frame(AIC,BIC,df,deviance,chi,formula))
}

# Basic Plot
gplot(net, gmode="graph", displaylabels=T, vertex.cex=.75,
      vertex.col="green", vertex.border="green", edge.col="blue",
      arrowhead.cex=.5, label.pos=5, label.cex=.5)

# Attributes
library(sand)
library(intergraph)

attributes$dc <- unname(degree(asIgraph(net), mode = "total"))
attributes$cc <- unname(closeness(asIgraph(net), mode = "total"))
attributes$bc <- unname(betweenness(asIgraph(net), directed = F))


net %v% "education" <- attributes$education
net %v% "contact" <- attributes$contact
net %v% "military" <- attributes$military
net %v% "nation" <- attributes$nation
net %v% "status" <- attributes$status
net %v% "role" <- attributes$role
net %v% "group" <- attributes$group
net %v% "noordin" <- attributes$noordin
net %v% 'degreecentrality' <- attributes$dc
net %v% 'closenesscentrality' <- attributes$cc
net %v% 'betweennesscentrality' <- attributes$bc

detach('package:sand', unload=TRUE)
detach('package:igraph', unload=TRUE)

# ERGM Models

null_model <- ergm(net ~ edges)
aa_model <- ergm(net ~ edges + nodecov('education') + nodefactor('role'))
aah_model <- ergm(net ~ edges + nodecov('education') + nodefactor('role') + absdiff('education') + nodematch('role'))
circuit_model <- ergm(net ~ edges + nodecov('education') + nodefactor('role') + absdiff('education') + nodematch('role') + gwesp(.7, fixed=F) + gwdsp(.7, fixed=F) + gwdegree(.7, fixed=F))

# GOF
df <- GOFMetrics(null_model,aa_model,aah_model,circuit_model)

gof_local_null_model <- gof(null_model ~ model)
gof_global_null_model <- gof(null_model ~ degree + espartners + dspartners + distance + triadcensus)

gof_local_aa_model <- gof(aa_model ~ model)
gof_global_aa_model <- gof(aa_model ~ degree + espartners + dspartners + distance + triadcensus)

gof_local_aah_model <- gof(aah_model ~ model)
gof_global_aah_model <- gof(aah_model ~ degree + espartners + dspartners + distance + triadcensus)

gof_local_circuit_model <- gof(circuit_model ~ model)
gof_global_circuit_model <- gof(circuit_model ~ degree + espartners + dspartners + distance + triadcensus)

# GOF Plots
par(mfrow=c(1,1))
plot(gof_local_null_model)
par(mfrow=c(2,3))
plot(gof_global_null_model)

par(mfrow=c(1,1))
plot(gof_local_aa_model)
par(mfrow=c(2,3))
plot(gof_global_aa_model)

par(mfrow=c(1,1))
plot(gof_local_aah_model)
par(mfrow=c(2,3))
plot(gof_global_aah_model)

par(mfrow=c(1,1))
plot(gof_local_circuit_model)
par(mfrow=c(2,3))
plot(gof_global_circuit_model)

# Coloring by role
role_colors <- get.vertex.attribute(net,"role")
education <- get.vertex.attribute(net, 'education')
colors <- c('gray60','yellow','red2','red1','red3','seagreen2','purple2','seagreen3','purple1','seagreen1','purple3','gray')
role_colors[role_colors == 1] = colors[1]
role_colors[role_colors == 2] = colors[2]
role_colors[role_colors == 3] = colors[3]
role_colors[role_colors == 4] = colors[4]
role_colors[role_colors == 5] = colors[5]
role_colors[role_colors == 6] = colors[6]
role_colors[role_colors == 7] = colors[7]
role_colors[role_colors == 8] = colors[8]
role_colors[role_colors == 9] = colors[9]
role_colors[role_colors == 10] = colors[10]
role_colors[role_colors == 11] = colors[11]
role_colors[role_colors == 12] = colors[12]

# Functions for display
scalevector <- function(x,min,max){
  y <- sqrt(x)
  y <- ((y)-min(y))/(max(y)-min(y))
  y[is.nan(y)] <- 0
  y <- y*(max-min)+min
  return(y)
}

colorvector <- function(x){
  return(rgb(1-(scalevector(x,0,1)), 1-(scalevector(x,0,1)), 1, 1))
}

# Plots
library(intergraph)
library(sand)

par(mfrow=c(1,1))

# Role and Betweenness Centralization
set.seed(50)
plot(asIgraph(net), vertex.label=NA, vertex.color = role_colors, vertex.size=scalevector(unname(betweenness(asIgraph(net), directed = F)),3,10))

# Education and Betweenness Centralization
set.seed(50)
plot(asIgraph(net), vertex.label=NA, vertex.color = colorvector(education), vertex.size=scalevector(unname(betweenness(asIgraph(net), directed = F)),3,10))

# Role and Education
set.seed(50)
plot(asIgraph(net), vertex.label=NA, vertex.color = role_colors, vertex.size=scalevector(education,3,10), edge.color = 'black')

#Noordin
set.seed(50)
plot(asIgraph(net), vertex.label=NA, vertex.color = attributes$noordin, vertex.size = 5)

detach('package:sand', unload=TRUE)
detach('package:igraph', unload=TRUE)

getdistribution <- function(attribute, min, max, base) {
  return(sapply(min:max, function(x) {
      mean(base[attribute == x])
    }))
}

plot(getdistribution(attributes$role,1,12,attributes$dc),pch = 16, main = 'Degree Centrality Plotted Against Role', xlab = 'Role', ylab = 'Degree Centrality')
plot(getdistribution(attributes$role,1,12,attributes$cc),pch =  16, main = 'Closeness Centrality Plotted Against Role', xlab = 'Role', ylab = 'Closeness Centrality')
plot(getdistribution(attributes$role,1,12,attributes$bc),pch =  16, main = 'Betweenness Centrality Plotted Against Role', xlab = 'Role', ylab = 'Betweenness Centrality')

par(xpd=FALSE)
plot(getdistribution(attributes$education,1,8,attributes$dc),pch = 16, main = 'Degree Centrality Plotted Against Education', xlab = 'Education', ylab = 'Degree Centrality')
abline(lm(attributes$dc ~ attributes$education))

plot(getdistribution(attributes$education,1,8,attributes$cc),pch = 16, main = 'Closeness Centrality Plotted Against Education', xlab = 'Education', ylab = 'Closeness Centrality')
abline(lm(attributes$cc ~ attributes$education))

plot(getdistribution(attributes$education,1,8,attributes$bc),pch = 16, main = 'Betweenness Centrality Plotted Against Education', xlab = 'Education', ylab = 'Betweenness Centrality')
abline(lm(attributes$bc ~ attributes$education))
