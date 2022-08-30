#######################################################################
##                                                                   ##
##  Title  : R Replication script for the BJPS manuscript            ##
##           "Boolean difference-making: a modern regularity         ##
##           theory of causation"                                    ##
##  Authors: Michael Baumgartner & Christoph Falk                    ##
##  Version: 26/09/2019                                              ##
##                                                                   ##
#######################################################################

# R version used 3.6
# Required R package (package version used: 2.2.1)
library(cna)

# Fundamentals
# ------------
# Table 1b:
dat1 <- allCombs(c(2,2,2,2,2)) -1
(tab1b <- selectCases("(A*B + C <-> D)*(c + a <-> E)", dat1))
# Minimally sufficient conditions for D:
ana1 <- cna(tab1b, what = "mac")
subset(msc(ana1), outcome == "D")
# RDN-biconditional for D:
subset(asf(ana1), outcome == "D")

# Structural redundancies
# -----------------------
# RDNB-conjunction for Table 1b:
ana1
# Structurally minimal conjunction of RDN-biconditionals:
minimalizeCsf(ana1)

# Test loop estimating the frequency of structural redundancies
# (different runs will generate different results; redundancy
# ratios will vary accordingly):
n <- 100
score <- vector("list", n)
for(i in 1:n){
   cat(i, "\n")
   x <- randomCsf(full.tt(8), n.asf = 3, compl = 2)
   y <- selectCases(x)
   score[[i]] <- csf(cna(y, details = T, rm.dup.factors = F), 1)
}
eval <- Filter(function(x) dim(x)[1] > 0, 
               lapply(score, function(x) subset(x, x$redundant == TRUE)))
# Structural redundancy ratio:
length(eval)/n

# Permanence
# ----------
# Table 1c:
(tab1c <- tt2df(tab1b)[, -2])
ana2 <- cna(tab1c)
# Structurally minimal conjunction of RDN-biconditionals:
minimalizeCsf(csf(ana2)$condition, dat1)

# Ambiguities
# -----------
# Table 2a:
(tab2a <- selectCases("(A + B <-> C)*(C + D <-> E)", dat1))
cna(tab2a)
# Tables 3a/b:
dat2 <- allCombs(c(2,2,2,2)) -1 
(tab3a <- selectCases("A*b + a*B + A*c <-> D", dat2))
# Two structurally minimal RDN-biconditionals:
cna(tab3a)
# Ambiguity resolution through factor set expansion:
(tab3b <- selectCases("A*b*E + a*B + B*c <-> D", dat1))
cna(tab3b)

# A new regularity theory
# ----------------------
# Table 4a:
(tab4a <- selectCases("(a + b + c <-> D)*(A + C <-> E)*(A + B <-> E)", dat1))
ana1 <- cna(tab4a, details = T)
# All RDN-biconditionals entailed by Table 4a:
asf(ana1)$condition
# All minimal theories for Table 4a:
(mt <- as.vector(minimalizeCsf(subset(csf(ana1, Inf),
                               exhaustiveness == 1)$condition, dat1)$condition))
# Table 4c:
dat3 <- allCombs(rep(2,7)) -1
(tab4c <- selectCases("(A + B*F <-> D)*(C + B*f <-> E)*(D + E <-> G)",
                      dat3))
cna(tab4c)