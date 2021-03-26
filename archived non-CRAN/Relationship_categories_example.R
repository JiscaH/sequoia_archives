# 1: pedigree statistics, similar to R package pedantics' pedigreeStats()
# 2: Find all pairs of individuals which are PO, FS, HS, GP, ...


library(sequoia)
source("E:/Sequoia/Rpackage/sequoia/R/RelationshipCategories.R")  # will be include in future version of sequoia

#===================================================
# generate example pedigree  with dummy parents
data(Ped_HSg5, LH_HSg5, package="sequoia")
GenoM <- SimGeno(Ped=Ped_HSg5, ParMis=0.4)
SEQ <- sequoia(GenoM = GenoM, LifeHistData = LH_HSg5, MaxSibIter=10)



#===================================================
# pedigree stats

PedStats(Ped = SEQ$Pedigree,
         SNPd = SEQ$PedigreePar$id)  # vector with genotyped individuals
# $PedStats
#                       GG  GD  GT DG DD DT   TT
# Individuals           NA  NA 920 NA NA 80 1000
# Maternities          572 328 900 44 16 60  960
# Paternities          482 418 900 34 26 60  960
# MaternalGrandmothers 460 265 725 28 15 43  768
# MaternalGrandfathers 384 341 725 24 19 43  768
# PaternalGrandmothers 500 225 725 28 15 43  768
# PaternalGrandfathers 488 237 725 28 15 43  768
#
# $MeanSibshipSize
#             G  D     T
# Full      4.1 NA  4.17
# Maternal  8.0  8  8.00
# Paternal 12.0 12 12.00

# Abbreviations in $Pedstats as in PedCompare()
# $MeanSibshipSize split into sibship with a genotyped parent (G) or with a dummy parent (D),
# and across both (T). Full - G is for both dam & sire genotyped.

PedStats(Ped = SEQ$Pedigree,
         SNPd = NA)    # no info on who is genotyped -> totals only.
# $PedStats
#                        TT
# Individuals          1000
# Maternities           960
# Paternities           960
# MaternalGrandmothers  768
# MaternalGrandfathers  768
# PaternalGrandmothers  768
# PaternalGrandfathers  768
#
# $MeanSibshipSize
#              T
# Full      4.17
# Maternal  8.00
# Paternal 12.00



#===================================================
# pairwise relationships

RelCat <- GetRelCat(Ped_HSg5)
RelCat2 <- GetRelCat(SEQ$Pedigree[,1:3])

# Note: fairly slow on very large pedigrees
# returns dataframe with ID1 - ID2 - RC

# abbreviations used in RC:
# 1st degree:
# MPO   maternal parent-offspring (i.e, mother-offspring)
# PPO   paternal parent-offspring
# FS    full siblings

# 2nd degree:
# MHS   maternal half-siblings
# PHS   patnernal half-siblings
# MGM   maternal grandmother - grandoffspring
# MGF   maternal grandfather - grandoffspring
# PGM   paternal grandmother - grandoffspring
# PGF   paternal grandfather - grandoffspring
# MFA   maternal full avuncular  (full sib of mother)
# PFA   paternal full avuncular  (full sib of father)
# DFC1  double full first cousin (two sets of parents are FS)
# DGGP  double great-grandparent

# 3rd degree:
# MMHA  maternal-maternal half avuncular (MHS of mother)
# MPHA  PHS of mother
# PMHA  MHS of father
# PPHA  PHS of father
# GGP   great-grandparent
# FC1   full first cousin (one set of parents are FS)
# DHC1  double half first cousin (two sets of parents are HS)

# if a pair is related via more than 1 way, only the closest relationship
#   is (currently) given



head(RelCat)
#      ID1    ID2  RC
# 1 a01016 a00011 MPO
# 2 a01008 a00013 MPO
# 3 b01015 a00011 MPO
# 4 b01007 a00013 MPO
# 5 b01041 a00018 MPO
# 6 a01042 a00018 MPO


table(RelCat$RC)
#  MPO  PPO   FS  MHS  PHS  MGM  MGF  PGM  PGF  MFA  PFA DFC1 DGGP MMHA MPHA
#  960  960 3200 3520 7360  768  768  736  728 2480 2480   64  376 2816 5856
# PMHA PPHA  GGP  FC1 DHC1
# 2624 5424 3820 9056  704



#================
# Compare this to pedigree relatedness

library(pedantics)
PedStats <- pedigreeStats(Ped_HSg5, graphicalReport=FALSE, includeA=TRUE)  # can take quite a while!
#save(PedStats, file="Pedstats.Rdata")
#load("Pedstats.Rdata")

Rel.ped <- as.data.frame.table(PedStats$Amatrix)
names(Rel.ped) <- c("ID1", "ID2", "R.PED")

# For large pedigrees, take a manageable portion of related individuals only,
# otherwise merging becomes extremely slow.
# Alternatively, use package data.table
Pairs <- with(Rel.ped, Rel.ped[ID1 != ID2 & R.PED >= 1/8, ])

Pairs <- merge(Pairs, RelCat, all.x=TRUE)
Pairs <- Merge(Pairs, setNames(RelCat[, c("ID1", "ID2", "RC")],
                               c("ID2", "ID1", "RC")),
               by=c("ID1", "ID2"), all.x=TRUE)
# 'Merge' is defined in RelationshipCategories.R, and is like normal merging
# but with filling in the NA's of the non-by columns

with(Pairs, boxplot(R.PED ~ RC))


CHK <- Pairs[is.na(Pairs$RC) & Pairs$R.PED>=0.25, ]  # non-considered double/triple relationships



# and similar can be done for genomic relatedness
# (see vignette part about comparing pedigree & genomic relatedness)


