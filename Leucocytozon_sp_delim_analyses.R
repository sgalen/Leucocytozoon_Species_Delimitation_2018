# This script contains the code needed to replicate the UniFrac and 
# GMYC analyses from Galen et al. "Integrating coalescent species delimitation with analysis of host specificity reveals extensive cryptic diversity despite minimal mitochondrial divergence in the malaria parasite genus Leucocytozoon"

# libraries needed to complete analyses 

library(dplyr)
library(ape)
library(vegan)
library(phyloseq)
library(bGMYC)
library(splits)

#############
#############
############# UniFrac Analysis
#############
#############

### data file can be accessed from dryad

data  <- read.delim("Leuco_sp_delim_for_analysis.txt", header=T, stringsAsFactors=F)

### prepare data for weighted UniFrac analysis 

# select only haplotypes that have been recorded at least twice
target <- c("CB1", "BT1", "TRPIP2", "PERCAN04", "CATMIN05", "CATUST28", "TUMIG15", "CATGUT02", "TUMIG11", "CATUST09", "ACAFLA03", "CATUST14", "PERCAN01", "ROFI06", "ZOLEU02", "CNEORN01", "CATMIN01", "LOXLEU05")

# prepare matrix
samples  <- rep(data$host_genus_species, 2)
lineages  <- c(data$leuco_lin_1, data$leuco_lin_2)
df_for_pres_abs_pseq  <- data.frame(lineages, samples)
df_for_pres_abs_pseq$lineages <- as.character(df_for_pres_abs_pseq$lineages)
df_for_pres_abs_pseq$samples <- as.character(df_for_pres_abs_pseq$samples)
df_for_pres_abs_pseq  <- filter(df_for_pres_abs_pseq, lineages != "", lineages != "cannot_parse")
df_for_pres_abs_pseq <- filter(df_for_pres_abs_pseq, lineages %in% target)
data_pres_ab_pseq  <- table(df_for_pres_abs_pseq$lineages, df_for_pres_abs_pseq$samples)

# make abundance matrix of leuco species X host species 
x <- matrix(data=data_pres_ab_pseq, nrow=nrow(data_pres_ab_pseq), ncol=ncol(data_pres_ab_pseq))
rownames(x) <- rownames(data_pres_ab_pseq)
colnames(x) <- colnames(data_pres_ab_pseq)

# transpose matrix so that hosts are rows ('taxa' in the unifrac) and parasites are columns (samples in the unifrac)
x_t <- t(x)
# check that it looks correct (hosts as rows and parasites as columns)
head(x_t)

# drop tips from host tree to match host-parasite matrix
bird_tree <- read.nexus("bird_tree.tre")
names_to_drop <- setdiff(bird_tree$tip.label, rownames(x_t))
reduced_bird_tree <- drop.tip(bird_tree, names_to_drop)
plot(reduced_bird_tree)

# build phyloseq object from matrix and phylogeny
my_otu_table <- otu_table(x_t, taxa_are_rows = TRUE)
my_phyloseq <- phyloseq(my_otu_table, reduced_bird_tree)

# calculate weighted unifrac
my_weighted_unifrac <- UniFrac(my_phyloseq, weighted=T)

# do PCOA from weighted unifrac
my_weighted_unifrac_pcoa <- pcoa(my_weighted_unifrac)
biplot(my_weighted_unifrac_pcoa)

### weighted UniFrac randomizations

# this function takes the observed host-parasite matrix (my_matrix),
# two putative species or haplotypes that you want to test for a difference between
# (lin_1 and lin_2), and the number of randomizations you want to perform (rand_num)
# it will output the true weighted unifrac distance between the two species,
# the mean weighted unifrac from the randomized distribution, and the 
# p value of where the true distance falls relative to the random distribution
unifrac_rand <- function(my_matrix, lin_1, lin_2, rand_num) {
  my_true_otu_table <- otu_table(my_matrix, taxa_are_rows = TRUE)
  my_true_phyloseq <- phyloseq(my_true_otu_table, reduced_bird_tree)
  my_true_weighted_unifrac <- UniFrac(my_true_phyloseq, weighted=T)
  my_true_weighted_unifrac_matrix <- as.matrix(my_true_weighted_unifrac)
  my_true_distance <- my_true_weighted_unifrac_matrix[lin_1, lin_2]
  
  unifrac_distances <- numeric(rand_num)
  rand_mat <- permatfull(my_matrix, fixedmar="both", times=rand_num)
  
  for (i in 1:length(rand_mat$perm)) {
    my_rand_mat <- rand_mat$perm[[i]]
    my_otu_table <- otu_table(my_rand_mat, taxa_are_rows = TRUE)
    my_phyloseq <- phyloseq(my_otu_table, reduced_bird_tree)
    my_weighted_unifrac <- UniFrac(my_phyloseq, weighted=T)
    my_weighted_unifrac_matrix <- as.matrix(my_weighted_unifrac)
    my_dist <- my_weighted_unifrac_matrix[lin_1, lin_2]
    unifrac_distances[i] <- my_dist
  }
  
  hist(unifrac_distances, breaks=200)
  print(c("Weighted UniFrac Distance:", my_true_distance))
  print(c("Mean from Distribution:", mean(unifrac_distances)))
  print(c("P-value:", length(unifrac_distances[unifrac_distances>my_true_distance])/rand_num))
  return(unifrac_distances)
}

# unifrac significance test for select pairs from this study
cb1_bt1 <- unifrac_rand(x_t, "CB1", "BT1", 1000)
cb1_trpip2 <- unifrac_rand(x_t, "CB1", "TRPIP2", 1000)
trpip2_bt1 <- unifrac_rand(x_t, "TRPIP2", "BT1", 1000)
zoleu02_cneorn01 <- unifrac_rand(x_t, "ZOLEU02", "CNEORN01", 1000)
acafla03_catust14 <- unifrac_rand(x_t, "ACAFLA03", "CATUST14", 1000)
acafla03_percan01 <- unifrac_rand(x_t, "ACAFLA03", "PERCAN01", 1000)
catust14_percan01 <- unifrac_rand(x_t, "CATUST14", "PERCAN01", 1000)
catust28_catmin05 <- unifrac_rand(x_t, "CATUST28", "CATMIN05", 1000)
rofi06_loxleu05 <- unifrac_rand(x_t, "ROFI06", "LOXLEU05", 1000)

#############
#############
############# GMYC Analysis
#############
#############

# read in concatenated BEAST tree
all_strict_coal_reduced <- read.tree("concatenated_strict_coal_gblocks_taxa_removed.tre")
# drop outgroup
all_strict_coal_reduced <- drop.tip(all_strict_coal_reduced, "SPIARB01_SPIARB01")

# single threshold model (concatenated tree)
all_sc_single_result <- gmyc(all_strict_coal_reduced)
summary(all_sc_single_result)
plot(all_sc_single_result)
spec.list(all_sc_single_result)

# multiple threshold model (concatenated tree)
all_sc_multi_result <- gmyc(all_strict_coal_reduced, method="multiple")
summary(all_sc_multi_result)
plot(all_sc_multi_result)
spec.list(all_sc_multi_result)

####
#### single gene analyses
####

# read in ATPASE2 tree
atp_strict_coal_reduced <- read.tree("atp_strict_coal_reduced_taxa_removed.tre")
# drop outgroup
atp_strict_coal_reduced <- drop.tip(atp_strict_coal_reduced, "SPIARB01_SPIARB01")

# single threshold model (ATPASE2)
atp_sc_single_result <- gmyc(atp_strict_coal_reduced)
summary(atp_sc_single_result)
plot(atp_sc_single_result)
spec.list(atp_sc_single_result)

# multiple threshold model (ATPASE2)
atp_sc_multi_result <- gmyc(atp_strict_coal_reduced, method="multiple")
summary(atp_sc_multi_result)
plot(atp_sc_multi_result)
spec.list(atp_sc_multi_result)


# read in RPB1 tree
ddrp_strict_coal_reduced <- read.tree("ddrp_strict_coal_reduced_taxa_removed.tre")
# drop outgroup
ddrp_strict_coal_reduced <- drop.tip(ddrp_strict_coal_reduced, "SPIARB01_SPIARB01")

# single threshold model (RPB1)
ddrp_sc_single_result <- gmyc(ddrp_strict_coal_reduced)
summary(ddrp_sc_single_result)
plot(ddrp_sc_single_result)
spec.list(ddrp_sc_single_result)

# multiple threshold model (RPB1)
ddrp_sc_multi_result <- gmyc(ddrp_strict_coal_reduced, method="multiple")
summary(ddrp_sc_multi_result)
plot(ddrp_sc_multi_result)
spec.list(ddrp_sc_multi_result)


# read in POLD1 tree
dpol_strict_coal_reduced <- read.tree("dpol_strict_coal_reduced_taxa_removed.tre")
# drop outgroup
dpol_strict_coal_reduced <- drop.tip(dpol_strict_coal_reduced, "SPIARB01_SPIARB01")

# multiple threshold model (POLD1)
dpol_sc_single_result <- gmyc(dpol_strict_coal_reduced)
summary(dpol_sc_single_result)
plot(dpol_sc_single_result)
spec.list(dpol_sc_single_result)

# multiple threshold model (POLD1)
dpol_sc_multi_result <- gmyc(dpol_strict_coal_reduced, method="multiple")
summary(dpol_sc_multi_result)
plot(dpol_sc_multi_result)
spec.list(dpol_sc_multi_result)


# read in KPNB1 tree
impo_strict_coal_reduced <- read.tree("impo_strict_coal_reduced_taxa_removed.tre")
# drop outgroup
impo_strict_coal_reduced <- drop.tip(impo_strict_coal_reduced, "SPIARB01_SPIARB01")

# single threshold model (KPNB1)
impo_sc_single_result <- gmyc(impo_strict_coal_reduced)
summary(impo_sc_single_result)
plot(impo_sc_single_result)
spec.list(impo_sc_single_result)

# multiple threshold model (KPNB1)
impo_sc_multi_result <- gmyc(impo_strict_coal_reduced, method="multiple")
summary(impo_sc_multi_result)
plot(impo_sc_multi_result)
spec.list(impo_sc_multi_result)


# read PRPF6 tree
prp_strict_coal_reduced <- read.tree("prp_strict_coal_reduced_taxa_removed.tre")
# drop outgroup
prp_strict_coal_reduced <- drop.tip(prp_strict_coal_reduced, "SPIARB01_SPIARB01")

# single threshold model (PRPF6)
prp_sc_single_result <- gmyc(prp_strict_coal_reduced)
summary(prp_sc_single_result)
plot(prp_sc_single_result)
spec.list(prp_sc_single_result)

# multiple threshold model (PRPF6)
prp_sc_multi_result <- gmyc(prp_strict_coal_reduced, method="multiple")
summary(prp_sc_multi_result)
plot(prp_sc_multi_result)
spec.list(prp_sc_multi_result)


# read in SEC24A tree
sec24_strict_coal_reduced <- read.tree("sec24_strict_coal_reduced_taxa_removed.tre")
# drop outgroup
sec24_strict_coal_reduced <- drop.tip(sec24_strict_coal_reduced, "SPIARB01_SPIARB01")

# single threshold model (SEC24A)
sec24_sc_single_result <- gmyc(sec24_strict_coal_reduced)
summary(sec24_sc_single_result)
plot(sec24_sc_single_result)
spec.list(sec24_sc_single_result)

# multiple threshold model (SEC24A)
sec24_sc_multi_result <- gmyc(sec24_strict_coal_reduced, method="multiple")
summary(sec24_sc_multi_result)
plot(sec24_sc_multi_result)
spec.list(sec24_sc_multi_result)


# read in SF3B1 tree
sf3_strict_coal_reduced <- read.tree("sf3_strict_coal_reduced_taxa_removed.tre")
# drop outgroup
sf3_strict_coal_reduced <- drop.tip(sf3_strict_coal_reduced, "SPIARB01_SPIARB01")

# single threshold model (SF3B1)
sf3_sc_single_result <- gmyc(sf3_strict_coal_reduced)
summary(sf3_sc_single_result)
plot(sf3_sc_single_result)
spec.list(sf3_sc_single_result)

# multiple threshold model (SF3B1)
sf3_sc_multi_result <- gmyc(sf3_strict_coal_reduced, method="multiple")
summary(sf3_sc_multi_result)
plot(sf3_sc_multi_result)
spec.list(sf3_sc_multi_result)

####
#### bGMYC analyses
####

# bGMYC uses a distribution of trees, so we must read in the posterior distribution from BEAST
# read trees from concatenated analysis
concatenated_trees <- read.nexus("concatenated_strict_coal_gblocks.trees")
concatenated_trees <- lapply(concatenated_trees, drop.tip, tip=c("SPIARB01_SPIARB01", "PRS4416_ACAFLA03", "PRS4431_PERCAN01"))


result.single <- bgmyc.singlephy(concatenated_trees[[10000]], mcmc=50000, burnin=1, thinning=10, t1=2, t2=100, start=c(1,1,25))
plot(result.single)

result.multi <- bgmyc.multiphylo(concatenated_trees[9900:10000], mcmc=50000, burnin=25000, thinning=100)
plot(result.multi)
result.spec <- bgmyc.spec(result.multi)
result.probmat <- spec.probmat(result.multi)
plot(result.probmat, concatenated_trees[[10001]])
bgmycrates <- checkrates(result.multi)
plot(bgmycrates)

# read in ATPASE2 trees
atp_trees <- read.nexus("atp_strict_coal.trees")
atp_trees <- lapply(atp_trees, drop.tip, tip=c("SPIARB01_SPIARB01", "PRS4416_ACAFLA03", "PRS4431_PERCAN01"))
result.multi <- bgmyc.multiphylo(atp_trees[9900:10000], mcmc=50000, burnin=25000, thinning=100, py1 = 0, py2 = 2, pc1 = 0, pc2 = 2, t1 = 2, t2 = 51, scale = c(20, 10, 5), start = c(1, 0.5, 25), sampler = bgmyc.gibbs, likelihood = bgmyc.lik, prior = bgmyc.prior)
result.spec.atp <- bgmyc.spec(result.multi)
#result.probmat.atp <- spec.probmat(result.multi)

# read in RPB1 trees
ddrp_trees <- read.nexus("ddrp_strict_coal.trees")
ddrp_trees <- lapply(ddrp_trees, drop.tip, tip=c("SPIARB01_SPIARB01", "PRS4416_ACAFLA03", "PRS4431_PERCAN01"))
result.multi.ddrp <- bgmyc.multiphylo(ddrp_trees[9900:10000], mcmc=50000, burnin=25000, thinning=100, py1 = 0, py2 = 2, pc1 = 0, pc2 = 2, t1 = 2, t2 = 51, scale = c(20, 10, 5), start = c(1, 0.5, 25), sampler = bgmyc.gibbs, likelihood = bgmyc.lik, prior = bgmyc.prior)
result.spec.ddrp <- bgmyc.spec(result.multi.ddrp)
#result.probmat.ddrp <- spec.probmat(result.multi)

# read in POLD1 trees
dpol_trees <- read.nexus("dpol_strict_coal.trees")
dpol_trees <- lapply(dpol_trees, drop.tip, tip=c("SPIARB01_SPIARB01", "PRS4416_ACAFLA03", "PRS4431_PERCAN01"))
result.multi.dpol <- bgmyc.multiphylo(dpol_trees[9900:10000], mcmc=50000, burnin=25000, thinning=100, py1 = 0, py2 = 2, pc1 = 0, pc2 = 2, t1 = 2, t2 = 51, scale = c(20, 10, 5), start = c(1, 0.5, 25), sampler = bgmyc.gibbs, likelihood = bgmyc.lik, prior = bgmyc.prior)
result.spec.dpol <- bgmyc.spec(result.multi.dpol)
#result.probmat.dpol <- spec.probmat(result.multi)

# read in KPBN1 trees
impo_trees <- read.nexus("impo_strict_coal.trees")
impo_trees <- lapply(impo_trees, drop.tip, tip=c("SPIARB01_SPIARB01", "PRS4416_ACAFLA03", "PRS4431_PERCAN01"))
result.multi.impo <- bgmyc.multiphylo(impo_trees[9900:10000], mcmc=50000, burnin=25000, thinning=100, py1 = 0, py2 = 2, pc1 = 0, pc2 = 2, t1 = 2, t2 = 51, scale = c(20, 10, 5), start = c(1, 0.5, 25), sampler = bgmyc.gibbs, likelihood = bgmyc.lik, prior = bgmyc.prior)
result.spec.impo <- bgmyc.spec(result.multi.impo)
#result.probmat.impo <- spec.probmat(result.multi)

# read in PRPF6 trees
prp_trees <- read.nexus("prp_strict_coal.trees")
prp_trees <- lapply(prp_trees, drop.tip, tip=c("SPIARB01_SPIARB01", "PRS4416_ACAFLA03", "PRS4431_PERCAN01"))
result.multi.prp <- bgmyc.multiphylo(prp_trees[9900:10000], mcmc=50000, burnin=25000, thinning=100, py1 = 0, py2 = 2, pc1 = 0, pc2 = 2, t1 = 2, t2 = 51, scale = c(20, 10, 5), start = c(1, 0.5, 25), sampler = bgmyc.gibbs, likelihood = bgmyc.lik, prior = bgmyc.prior)
result.spec.prp <- bgmyc.spec(result.multi.prp)
#result.probmat.prp <- spec.probmat(result.multi)

# read in SEC24A trees
sec24_trees <- read.nexus("sec24_strict_coal.trees")
sec24_trees <- lapply(sec24_trees, drop.tip, tip=c("SPIARB01_SPIARB01", "PRS4416_ACAFLA03", "PRS4431_PERCAN01"))
result.multi.sec24 <- bgmyc.multiphylo(sec24_trees[9900:10000], mcmc=50000, burnin=25000, thinning=100, py1 = 0, py2 = 2, pc1 = 0, pc2 = 2, t1 = 2, t2 = 51, scale = c(20, 10, 5), start = c(1, 0.5, 25), sampler = bgmyc.gibbs, likelihood = bgmyc.lik, prior = bgmyc.prior)
result.spec.sec24 <- bgmyc.spec(result.multi.sec24)
#result.probmat.sec24 <- spec.probmat(result.multi)

# read in SF3B1 trees
sf3_trees <- read.nexus("sf3_strict_coal.trees")
sf3_trees <- lapply(sf3_trees, drop.tip, tip=c("SPIARB01_SPIARB01", "PRS4416_ACAFLA03", "PRS4431_PERCAN01"))
result.multi.sf3 <- bgmyc.multiphylo(sf3_trees[9900:10000], mcmc=50000, burnin=25000, thinning=100, py1 = 0, py2 = 2, pc1 = 0, pc2 = 2, t1 = 2, t2 = 51, scale = c(20, 10, 5), start = c(1, 0.5, 20), sampler = bgmyc.gibbs, likelihood = bgmyc.lik, prior = bgmyc.prior)
result.spec.sf3 <- bgmyc.spec(result.multi.sf3)
#result.probmat.sf3 <- spec.probmat(result.multi)
