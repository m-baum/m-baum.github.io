#######################################################################
##                                                                   ##
##  Title  : R Replication script for                                ##
##           "Optimizing Consistency and Coverage in                 ##
##            Configurational Causal Modeling"                       ##
##  Authors: M. Baumgartner & M. Ambühl                              ##
##  Version: 10/07/2020                                              ##
##                                                                   ##
#######################################################################


library(QCA)
library(cna)
library(cnaOpt)
library(ggplot2)
library(dplyr)
library(microbenchmark)
library(reshape)


# Section 3
# ---------
# Data from Giugni and Yamasaki (2009, 476), Table 1a
data <- structure(list(P = c(1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 
   0L), O = c(0L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L), A = c(1L, 
   0L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 0L), C = c(0L, 1L, 0L, 0L, 
   0L, 0L, 0L, 1L, 1L, 1L, 0L)), class = "data.frame", row.names = c("E75", 
   "E87", "P81", "P90", "N77", "E80", "E92", "P75", 
   "P87", "P92", "N80"))

# Run ConCovOpt
ccoGY <- conCovOpt(data, outcome = "C", allConCov = T)

# Configuration table in Table 2 (result of step 1)
CT_GY <- attr(ccoGY, "configTable")
CT_GY

# Exo-groups (Step 2)
attr(ccoGY$C, "exoGroups")

# Rep-list (reproduction list) (result of step 3)
attr(ccoGY$C, "reprodList")

# All con-cov scores (result of steps 4-5)
attr(ccoGY$C, "allConCov")

# Con-cov optima (result of step 6)
ccoGY$C

# Con-cov maximum
(best <- selectMax(ccoGY))
# alternative con-cov maxima
selectMax(ccoGY, crit = quote(con^0.75* cov^0.25))
selectMax(ccoGY, crit = quote(con*0.25 +  cov*0.75))
selectMax(ccoGY, crit = quote(pmin(con,cov)))

# Rep-assignments realizing the con-cov optima 
CT_GY$PHI1 <- reprodAssign(best, "C", id=1)
CT_GY$PHI2 <- reprodAssign(best, "C", id=2) # con-cov maximum
CT_GY

# Building instances of the con-cov optima using CCMs
# QCA-PS (PHI1)
tt <- truthTable(data, incl.cut1 = .67, outcome = "C")
minimize(tt, include = "?", all.sol = T, details = T)
# QCA-CS (published model) (PHI1)
minimize(tt, all.sol = T, details = T)

# CCubes (PHI2)
minimize(tt, pi.cons = .75, sol.con = .8, details = T, method = "CCubes")

# CNA (PHI1 / PHI2)
cna(data, cov = .5, ordering = list("C"))
cna(data, con.msc = .75, con = .8, ordering = list("C"))



# Section 4
# ---------
# Data from Britt et al. (2000)
data <- structure(list(C = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 
    1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L), M = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 
    1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
    1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L), A = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L), G = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 
    1L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
    0L, 0L, 0L, 0L), T = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 
    0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 0L, 0L, 
    1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L)), class = "data.frame", row.names = c(NA, 142L
    ))

# Published model
# Run QCA at consistency threshold 0.875
tt <- truthTable(data, outcome = "T", incl.cut1 = .875)
minimize(tt, include="?", all.sol = T, details = T)

# Run ConCovOpt
ccoBritt <- conCovOpt(data, outcome = "T", allConCov = T)

# Configuration table in Table 2 (result of step 1)
CT_Britt <- attr(ccoBritt, "configTable")
CT_Britt

# Exo-groups (result of step 2)
attr(ccoBritt$T, "exoGroups")

# Rep-list (reproduction list) (result of step 3)
attr(ccoBritt$T, "reprodList")

# All con-cov scores (result of steps 4-5)
attr(ccoBritt$T, "allConCov")

# Con-cov optima (result of step 6)
ccoBritt$T

# Plot in Figure 1a
plot1a <- plot(ccoBritt) +
   geom_point(shape = 15, data = ccoBritt$T[3,], colour = "black", size = 3) +
   geom_point(shape = 17, data = ccoBritt$T[16,], color = "black", size = 3) +
   xlim(c(0.85,1)) + ylim(c(0.4,1)) + geom_point(colour = "black") + 
   geom_line(colour = "black") + theme(legend.position = "none") + 
   geom_text(aes(label = rownames(ccoBritt$T)), hjust = 1, 
   vjust = 2, colour = "black", size = 3)
plot1a

# Plot in Figure 1b
product <- ccoBritt$T
product$prod <-  product$con*product$cov
product$names <- rownames(product)
prodcut2 <- product[,4:5]
prodcut2$names <- as.numeric(prodcut2$names)
plot1b <- ggplot(data=prodcut2, aes(x = names, y = prod))  +
   geom_point() + geom_line() +
   geom_point(shape = 15, data = prodcut2[3,], colour = "black", size = 3) +
   geom_point(shape = 17, data = prodcut2[16,], color = "black", size = 3) +
   ylim(c(0.4,1)) + geom_text(aes(label = prodcut2$names), hjust = -.3, 
                               vjust = 2, colour = "black", size = 3)
plot1b

# Rep-assignment realized by the published model
pub <- selectMax(ccoBritt, cond = quote(con >= 0.97 & cov >= 0.7))
CT_Britt$published <- as.vector(reprodAssign(pub, outcome = "T"))
CT_Britt

# Rep-assignment reaching con-cov maximum
max <- selectMax(ccoBritt)
CT_Britt$max <- as.vector(reprodAssign(max, outcome = "T"))
CT_Britt

# CCM model reaching con-cov maximum
# found with CNA
cna(data, ordering = list("T"), strict = T, con = .89, cov = 1,
    maxstep=c(5,5,15), what = "mac", details = "i", inus.only = F)

# found with CCubes
minimize(tt, include = "?", pi.cons = .89, sol.cons = .89, sol.cov = 1,
         all.sol = T, details = T, method = "CCubes")


# Section 5
# ---------
# Data from Verweij and Gerrits (2015)
data <- structure(list(E = c(0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 0L, 0L, 0L, 0L, 1L, 1L, 0L), M = c(0L, 0L, 1L, 1L, 1L, 1L, 
    0L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L), I = c(2L, 1L, 
    0L, 2L, 0L, 0L, 0L, 0L, 1L, 2L, 1L, 1L, 0L, 1L, 1L, 0L, 2L, 1L
    ), S = c(1L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 1L, 
    0L, 0L, 0L, 1L, 0L)), class = "data.frame", row.names = c("TUN", 
    "BI1", "BI2", "WAT", "LAN", "ZON", "ENV", "LEE", "WES", "RIJ", 
    "TRA", "CR1", "CR2", "SOI", "BAD", "PRO", "CIV", "THE"))

# Recode the outcome
data$S <- ifelse(data$S==1, 3, 2)
data

# Build configuration table (Table 3)
CT_VG <- mvct(data)
CT_VG <- CT_VG[c(1,2,3,5,6,7,4,8,9,10),] # re-order the rows
CT_VG

# Run ConCovOpt (on the pre-built CT)
ccoVerweiy <- conCovOpt(CT_VG, outcome = "S=3", allConCov = T)

# Exo-groups (result of step 2)
attr(ccoVerweiy$"S=3", "exoGroups")

# Rep-list (reproduction list) (result of step 3)
attr(ccoVerweiy$"S=3", "reprodList")

# All con-cov scores (result of steps 4-5)
attr(ccoVerweiy$"S=3", "allConCov")

# Con-cov optima (result of step 6)
ccoVerweiy$"S=3"

# Rep-assignments 1-4
maxVerweiy <- selectMax(ccoVerweiy)
CT_VG$phi1 <- as.vector(reprodAssign(maxVerweiy, outcome = "S=3", id=1))
CT_VG$phi2 <- as.vector(reprodAssign(maxVerweiy, outcome = "S=3", id=2))
CT_VG$phi3 <- as.vector(reprodAssign(maxVerweiy, outcome = "S=3", id=3))
CT_VG$phi4 <- as.vector(reprodAssign(maxVerweiy, outcome = "S=3", id=4))
CT_VG

# Rep-assignment realizing the con-cov maximum
CT_VG$max <- as.vector(reprodAssign(maxVerweiy, outcome = "S=3"))
CT_VG

# Published model, built by QCA, conservative solution
tt <- truthTable(data, outcome = "S{3}", incl.cut1 = 0.75)
minimize(tt, all.sol = T, details = T) 

# Parsimonious QCA solutions
minimize(tt, all.sol = T, details=T, include = "?")

# Canonical DNF realizing PHI4
DNF_cano4 <- DNFbuild(maxVerweiy, outcome = "S=3", reduce = F)
DNF_cano4
# DNF_cano4 scores (0.75, 1) in accounting for S=3
mvcond(paste0(DNF_cano4," <-> S=3"), data)

# Canonical DNF of PHI1 to PHI3
DNF_cano1 <- DNFbuild(maxVerweiy, outcome = "S=3", reduce = F, id=1)
DNF_cano2 <- DNFbuild(maxVerweiy, outcome = "S=3", reduce = F, id=2)
DNF_cano3 <- DNFbuild(maxVerweiy, outcome = "S=3", reduce = F, id=3)

# Removing E=1 from DNF_cano4 does not affect the con-cov score
mvcond("E=0*M=0*I=2 + M=1*I=0 + M=0*I=0 + M=1*I=1 + M=1*I=2 + E=0*M=0*I=0 <-> S=3", data)

# Further removing I=2 from the first disjunct affects the con-cov score
mvcond("E=0*M=0 + M=1*I=0 + M=0*I=0 + M=1*I=1 + M=1*I=2 + E=0*M=0*I=0 <-> S=3", data)

# Reducing redundancies from DNF_cano4
DNF_free4 <- ereduce(DNF_cano4, CT_VG, full = F)
DNF_free4
# The same can be obtained using the DNFbuild() function 
DNFbuild(maxVerweiy, outcome = "S=3", reduce = "ereduce")
# DNF_free4 scores (0.75, 1) in accounting for S=3
mvcond(paste0(DNF_free4," <-> S=3"), data)

# Recovering DNF_free4 by means of CNA and CCubes
# CNA
mvcna(data, ordering=list("S"), strict=T, con=.66)
# CCubes
minimize(tt, pi.cons=.66, sol.cons=.75, sol.cov=1, details=T, method = "CCubes")

# Reducing redundancies from DNF_cano2
DNF_free2 <- ereduce(DNF_cano2, CT_VG, full = F)
DNF_free2
# DNF_free2 scores (0.86, 0.67) in accounting for S=3
mvcond(paste0(DNF_free2," <-> S=3"), data)

# Recovering DNF_free2 by means of CNA and CCubes
# CNA
asf(mvcna(data, ordering=list("S"), strict = T, con = .75, cov = .66))
# CCubes
minimize(tt, pi.cons = .75, sol.cons = .75, sol.cov = .66, details = T, method = "CCubes")



# Section 6
# ---------
# Data from Basurto (2013)
data <- d.autonomy[15:30, c("EM","SP","CO","AU")]
# Re-organize the row order to fit the order in Table 4 
data <- data[c(1,2,4,6,3,5,7,8,10,14,9,12,15,11,13,16),]
colnames(data) <- c("E", "S", "C", "A")
data

# Run ConCovOpt
ccoBas <- conCovOpt(fsct(data), outcome = "A", allConCov = T)

# Configuration table in Table 4 (result of step 1)
CT_Bas <- attr(ccoBas, "configTable")
CT_Bas

# Exo-groups (result of step 2)
attr(ccoBas$A, "exoGroups")

# Rep-list (reproduction list) (result of step 3)
attr(ccoBas$A, "reprodList")

# All con-cov scores (result of steps 4-5)
attr(ccoBas$A, "allConCov")

# Con-cov optima (result of step 6)
ccoBas$A

# Plot in Figure 2a
plot3 <- plot(ccoBas) +
   geom_point(shape = 15, data = data.frame(con = 0.786,
                                       cov = 0.936), colour = "black", size = 3) +
  geom_point(shape = 17, data = ccoBas$A[5,], color = "black", size = 3) +
  geom_point(shape = 18, data = ccoBas$A[3,], color = "black", size = 4) +
  xlim(c(0.77, 1.001)) + ylim(c(0.8, 1.001))+
  geom_point(colour = "black") + geom_line(colour = "black")+
  theme(legend.position = "none") + geom_text(aes(label = rownames(ccoBas$A)),
                                      hjust = 1, vjust = 2, colour = "black", size = 3)
plot3

# Plot in Figure 2b
product <- ccoBas$A
product <- rbind(product,data.frame(con = 0.786, cov = 0.936, id=90))
product$prod <-  product$con*product$cov
product$names <- rownames(product)
prodcut2 <- product[,4:5]
prodcut2$names <-as.numeric(prodcut2$names)
plot4 <- ggplot(data = prodcut2[1:7,], aes(x = names, y = prod)) +
   geom_point() + geom_line()+
   geom_point(shape = 18, data=prodcut2[3,], colour = "black", size = 4) +
   geom_point(shape = 17, data=prodcut2[5,], color = "black", size = 3) +
   geom_point(shape = 15, data=prodcut2[8,], color = "black", size = 3) +
   ylim(c(0.7,1.0))+ 
   geom_text(aes(label=1:7), hjust = -.3, vjust=2, colour = "black", size = 3)
plot4

# Published model, intermediate solution at a consistency threshold of 0.79
tt <- truthTable(data, outcome = "A", incl.cut1 = .79)
minimize(tt, all.sol = T, details = T, dir.exp = c(1,1,1), include = "?")
# Parsimonious QCA solutions at a consistency threshold  of 0.79
minimize(tt, all.sol = T, details = T,  include = "?")

# QCA issues varying models under different threshold settings
tt1 <- truthTable(data, outcome = "A", incl.cut1 = .85)
minimize(tt1, all.sol = T, details = T, include = "?")
tt2 <- truthTable(data, outcome = "A", incl.cut1 = .9)
minimize(tt2, all.sol = T, details = T, include = "?")
tt3 <- truthTable(data, outcome = "A", incl.cut1 = .93)
minimize(tt3, all.sol = T, details = T, include = "?")

# Finding a model realizing PHI1 (cf. Table 4) by CCMs
# CNA
fscna(data, ordering = list("A"), strict = TRUE, con = .93, cov = .91)
# CCubes
tt <- truthTable(data, outcome = "A", incl.cut1 = .93)
minimize(tt, pi.cons = .93, sol.cons = .93, sol.cov = .91, details = T, method = "CCubes") 

# Computational limits
# CS data
# Simulate ideal data with over 500 cases for a random target
fullData <- ct2df(full.ct(10))
target <- randomAsf(fullData, outcome = "D", compl = 3)
x <- ct2df(selectCases(target, fullData))
nrow(x)
# Generate 23 imperfect pairs
y <- setdiff(fullData,x)
data <- rbind(y[sample(1:nrow(y), 23, replace=F), ], x)
# Process data with conCovOpt 
conCovOpt(data, outcome = "D", approx = F)

# FS data without approximation (i.e. without search heuristic)
# Simulate ideal data for a random target
fullData <- ct2df(full.ct(7))
target <- randomAsf(fullData, outcome = "D", compl = 3)
x <- ct2df(selectCases(target, fullData))
# Generate 10 imperfect pairs to yield a total of 74 cases
y <- setdiff(fullData,x)
y <- rbind(y[sample(1:nrow(y), 10, replace=F), ], x)
nrow(y)
# Fuzzify the data
data <- makeFuzzy(y, fuzzvalues = seq(0, 0.4, 0.1))
# Process data with conCovOpt 
conCovOpt(data, type = "fs", outcome = "D", approx = F)

# FS data with approximation (i.e. with search heuristic)
# Simulate ideal data for a random target
fullData <- ct2df(full.ct(12))
target <- randomAsf(fullData, outcome = "D", compl = 3)
x <- ct2df(selectCases(target, fullData))
# Add 30 imperfect pairs to yield a total of 2074 cases
y <- setdiff(fullData,x)
y <- rbind(y[sample(1:nrow(y), 30, replace=F), ], x)
nrow(y)
# Fuzzify the data
data <- makeFuzzy(y, fuzzvalues = seq(0, 0.4, 0.1))
# Process data with conCovOpt 
conCovOpt(data, type = "fs", outcome = "D", approx = TRUE)



# Section 7
# ---------
# Illustrate the danger of overfitting
# Draw a random data generating structure
set.seed(59)
# Generate a random target.
fullData <- full.ct(5)
target <- randomAsf(fullData, outcome = "E", compl = 2:3)
target
# Generate the complete, ideal data.
x <- ct2df(selectCases(target, fullData))
# Introduce fragmentation by eliminating 30% of the cases.
x <- x[-sample(1:nrow(x), nrow(x)*0.3), ] 
# Introduce random noise by adding 20% random cases.
data <- rbind(ct2df(fullData[sample(1:nrow(fullData), nrow(x)*0.2), ]), x)  
# Run CNA on data at a conventional threshold setting
cna(data, ordering=list("E"), strict = T, con = .75, cov = .75)
# Calculate con-cov optima for outcome E.
(cco <- conCovOpt(data, "E"))
# Select the con-cov maximum for outcome E.
max <- selectMax(cco)
# Find all redundancy-free DNFs returning the con-cov maximum.
(DNF_max <- DNFbuild(max, outcome = "E", reduce = "ereduce"))
condTbl(paste0(DNF_max, "<->E"), data)
# Test whether the redundancy-free DNF yield at least one submodel of the target structure.
any(is.submodel(paste0(DNF_max, "<->E"), target))
# Whenever that test returns FALSE, overfitting has occurred. By re-running
# lines 424-442 repeatedly, the reader can get a sense of how 
# frequently overfitting occurs.


