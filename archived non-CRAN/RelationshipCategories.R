


#===============================================================
#' Pedigree statistics
#'
#' Basic summary statistics for a pedigree

PedStats <- function(Ped, SNPd=NULL) {

  PedIN <- AddParPed(Ped[,1:3])
  names(PedIN) <- c("id", "dam", "sire")

  for (x in c("id", "dam", "sire")) PedIN[,x] <- as.character(PedIN[,x])
  Ped <- merge(PedIN, setNames(PedIN, c("dam", "MGM", "MGF")), all.x=TRUE)
  Ped <- merge(Ped,
               setNames(PedIN, c("sire", "PGM", "PGF")), all.x=TRUE)

  SNPd <- na.exclude(SNPd)


  PedStats <- matrix(NA, 7, 7,
                     dimnames=list(c("Individuals", "Maternities", "Paternities",
                                "MaternalGrandmothers", "MaternalGrandfathers",
                                "PaternalGrandmothers", "PaternalGrandfathers"),
                                c("GG", "GD", "GT", "DG", "DD", "DT", "TT")))
  MeanSibshipSize <- matrix(NA, 3, 3,
                            dimnames= list(c("Full", "Maternal", "Paternal"),
                                           c("G", "D", "T")))

  PedStats["Individuals", c("GT", "DT", "TT")] <- c(table(Ped$id %in% SNPd)[2:1], nrow(Ped))
  PedStats["Maternities", ] <- GDT("dam", Ped, SNPd)
  PedStats["Paternities", ] <- GDT("sire", Ped, SNPd)
  PedStats["MaternalGrandmothers", ] <- GDT("MGM", Ped, SNPd)
  PedStats["MaternalGrandfathers", ] <- GDT("MGF", Ped, SNPd)
  PedStats["PaternalGrandmothers", ] <- GDT("PGM", Ped, SNPd)
  PedStats["PaternalGrandfathers", ] <- GDT("PGF", Ped, SNPd)

  MeanSibshipSize["Maternal", "T"] <- mean(table(Ped$dam))
  MeanSibshipSize["Maternal", "G"] <- mean(table(Ped$dam[Ped$dam %in% SNPd]))
  MeanSibshipSize["Maternal", "D"] <- mean(table(Ped$dam[!Ped$dam %in% SNPd]))

  MeanSibshipSize["Paternal", "T"] <- mean(table(Ped$sire))
  MeanSibshipSize["Paternal", "G"] <- mean(table(Ped$sire[Ped$sire %in% SNPd]))
  MeanSibshipSize["Paternal", "D"] <- mean(table(Ped$sire[!Ped$sire %in% SNPd]))

  Ped$dam_sire <- with(Ped, ifelse(!is.na(dam) & !is.na(sire), paste(dam, sire, sep="_"), NA))
  MeanSibshipSize["Full", "T"] <- mean(table(Ped$dam_sire))
  MeanSibshipSize["Full", "G"] <- mean(table(Ped$dam_sire[Ped$dam %in% SNPd &
                                                             Ped$sire %in% SNPd]))
  MeanSibshipSize <- round(MeanSibshipSize, 2)

  if (length(SNPd)<=1) {
    PedStats <- PedStats[,"TT", drop=FALSE]
    MeanSibshipSize <- MeanSibshipSize[, "T", drop=FALSE]
  }

  list(PedStats=PedStats, MeanSibshipSize=MeanSibshipSize)
}


GDT <- function(x, ped=Ped, snpd=SNPd) {
  c(GG = sum(ped[,"id"] %in% snpd & ped[,x] %in% snpd),
    GD = sum(ped[,"id"] %in% snpd & !ped[,x] %in% snpd & !is.na(ped[,x])),
    GT = sum(ped[,"id"] %in% snpd & !is.na(ped[,x])),
    DG = sum(!ped[,"id"] %in% snpd & ped[,x] %in% snpd),
    DD = sum(!ped[,"id"] %in% snpd & !ped[,x] %in% snpd & !is.na(ped[,x])),
    DT = sum(!ped[,"id"] %in% snpd & !is.na(ped[,x])),
    TT = sum(!is.na(ped[,x])))
}


#===============================================================
#' Relationship categories
#'
#' For each pair of individuals, infer the type of relationship
#'

GetRelCat <- function(Ped) {
  PedIN <- Ped
  names(PedIN) <- c("id", "dam", "sire")
  for (x in c("id", "dam", "sire")) PedIN[,x] <- as.character(PedIN[,x])
  Ped <- merge(PedIN, setNames(PedIN, c("dam", "MGM", "MGF")), all.x=TRUE)
  Ped <- merge(Ped,
               setNames(PedIN, c("sire", "PGM", "PGF")), all.x=TRUE)
  rm(PedIN)

  RelCat <- list()
  # 1st degree
  RelCat$MPO <- setNames(Ped[!is.na(Ped$dam), c("id", "dam")], c("ID1", "ID2"))
  RelCat$PPO <- setNames(Ped[!is.na(Ped$sire), c("id", "sire")], c("ID1", "ID2"))
  RelCat$FS <- getAllSibs(Ped, cat="F")
  # 2nd degree
  RelCat$MHS <- getAllSibs(Ped, cat="M")
  RelCat$PHS <- getAllSibs(Ped, cat="P")
  for (x in c("MGM", "MGF", "PGM", "PGF")) {
    RelCat[[x]] <- setNames(Ped[!is.na(Ped[,x]), c("id", x)], c("ID1", "ID2"))
  }
  RelCat$MFA <- getAllAU(Ped, cat="MF")
  RelCat$PFA <- getAllAU(Ped, cat="PF")
  RelCat$DFC1 <- getAllCC(Ped, cat="DF")
  RelCat$DGGP <- getAllGGP(Ped, cat="D")
  # 3rd degree
  for (x in c("MM", "MP", "PM", "PP")) {
    RelCat[[paste0(x,"HA")]] <- getAllAU(Ped, cat=x)
  }
  RelCat$GGP <- getAllGGP(Ped, cat="A")
  RelCat$FC1 <- getAllCC(Ped, cat="F")
  RelCat$DHC1 <- getAllCC(Ped, cat="DH")  # double half cousins


  RCDF <- ldply(RelCat, .id="RC")
  RCDF <- RCDF[!duplicated(RCDF[,c("ID1", "ID2")]), ]  # !! TODO?
  RCDF$RC <- factor(RCDF$RC, levels=names(RelCat), ordered=TRUE)
  RCDF[, c("ID1", "ID2", "RC")]
}



#==============
# add non-assigned likely relatives  [TODO]







#################################################################################

getSibs <- function(x, ped, cat="F") {
  y <- NULL
  if (cat=="F") {
    if (!is.na(ped$sire[x]) & !is.na(ped$dam[x])) {
    y <- with(ped, id[which(sire==sire[x] & dam==dam[x] & id != id[x])]) }
  } else if (cat=="P" & !is.na(ped$sire[x])) {
    y <- with(ped, id[which(sire==sire[x] & id != id[x] &
                                  (dam!=dam[x] | is.na(dam[x]) | is.na(dam)))])
  } else if (cat=="M" & !is.na(ped$dam[x])) {
    y <- with(ped, id[which(dam==dam[x] & id != id[x] &
                            (sire!=sire[x] | is.na(sire[x]) | is.na(sire)))])
  }
  if (length(y)>0) {
    out <- data.frame(ID1 = ped$id[x], ID2 = y)
  }
}

getAllSibs <- function(ped, cat="F") {
  tmp <- lapply(1:nrow(ped), getSibs, ped, cat)
  unique(ldply(tmp[!sapply(tmp, is.null)]))
}


#===============
# FA, HA

getAU <- function(x, ped, cat="PP") {
  y <- NULL
  if (cat=="PF" & !is.na(ped$PGF[x]) & !is.na(ped$PGM[x])) {
    y <- with(ped, id[which(sire==PGF[x] & dam==PGM[x] & !id %in% c(sire[x], dam[x], id[x]))])
  } else if (cat=="MF" & !is.na(ped$MGF[x]) & !is.na(ped$MGM[x])) {
    y <- with(ped, id[which(sire==MGF[x] & dam==MGM[x] & !id %in% c(sire[x], dam[x], id[x]))])
  } else if (cat=="PP" & !is.na(ped$PGF[x])) {
    y <- with(ped, id[which(sire==PGF[x] & id!=sire[x] & id!=id[x])])  # exclude inbreds
  } else if (cat=="MP" & !is.na(ped$MGF[x])) {
    y <- with(ped, id[which(sire==MGF[x] & id!=dam[x] & id!=id[x])])
  } else if (cat=="PM" & !is.na(ped$PGM[x])) {
    y <- with(ped, id[which(dam==PGM[x] & id!=sire[x] & id!=id[x])])
  } else if (cat=="MM" & !is.na(ped$MGM[x])) {
    y <- with(ped, id[which(dam==MGM[x] & id!=dam[x] & id!=id[x])])
  }
  if (length(y)>0) {
    data.frame(ID1 = ped$id[x], ID2 = y)
  }
}

getAllAU<- function(ped, cat="PF") {
  tmp <- lapply(1:nrow(ped), getAU, ped, cat)
  unique(ldply(tmp[!sapply(tmp, is.null)]))
}


getCC <- function(x, ped, cat="P") {  # first cousins (from FS/PHS/MHS)
  y <- NULL
  gp <- na.exclude(unlist(ped[x, c("MGF","PGF","MGM","PGM")]))
  if (length(gp)==0) return()
  if (cat=="DF") {  # double full first cousins (1/4)
    y <- with(ped, id[which((PGF %in% gp & PGM %in% gp &
                              MGF %in% gp & MGM %in% gp) &
                                sire!=sire[x] & dam!=dam[x])])
  } else if (cat=="F") { # full first cousins (1/8)
    y <- with(ped, id[which(((PGF %in% gp & PGM %in% gp) |
                              (MGF %in% gp & MGM %in% gp)) &
                                sire!=sire[x] & dam!=dam[x])])
  } else if (cat=="DH") { # double half first cousins (1/8)
    y <- with(ped, id[which(((PGF %in% gp & MGF %in% gp) |
                              (MGF %in% gp & PGM %in% gp)) &
                                sire!=sire[x] & dam!=dam[x])])
  } else if (cat=="H") { # half (1/16)
    y <- with(ped, id[which((PGF %in% gp | MGF %in% gp |
                               PGM %in% gp | MGM %in% gp) &
                                sire!=sire[x] & dam!=dam[x])])
  }
  if (length(y)>0) {
    data.frame(ID1 = ped$id[x], ID2=y)
  }
}

getAllCC <- function(ped, cat="F") {
  tmp <- lapply(1:nrow(ped), getCC, ped, cat)
  df <- unique(ldply(tmp[!sapply(tmp, is.null)]))
  if (nrow(df)>0) {
    df2 <- rbind(df, setNames(df[,c("ID2","ID1")], c("ID1","ID2")))
    df2[duplicated(df2),]  # exclude parents of either ID1 or ID2 are HS
  } else {
    df
  }
}


#==========
getGGP <- function(x, ped, cat="A") {  # great-grandparent  (1/8)
  y <- NULL
  if (!is.na(ped$sire[x])) {
    y <- with(ped, c(y, unlist(ped[which(id==sire[x]), c("MGM", "MGF", "PGM", "PGF")])))
  }
  if (!is.na(ped$dam[x])) {
    y <- with(ped, c(y, unlist(ped[which(id==dam[x]), c("MGM", "MGF", "PGM", "PGF")])))
  }
  if (length(y)>0) {
    if (cat=="A") y <- unique(na.exclude(y))
    if (cat=="D") y <- unique(y[duplicated(y, incomparables=NA)])
  }
  if (length(y)>0) {
    data.frame(ID1 = ped$id[x], ID2=y)
  }
}

getAllGGP <- function(ped, cat="A") {
  tmp <- lapply(1:nrow(ped), getGGP, ped, cat=cat)
  unique(ldply(tmp[!sapply(tmp, is.null)]))
}


#============================================================
# utils (in R package sequoia)
Merge <- function(df1, df2, by, ...) {
  commonNames <- names(df1)[which(colnames(df1) %in% colnames(df2))]
  commonNames <- commonNames[!commonNames %in% by]
  dfmerged <- merge(df1,df2,by=by,...)
  for(i in commonNames){
    left <- paste0(i, ".x")
    right <- paste0(i, ".y")
    dfmerged[is.na(dfmerged[left]),left] <- dfmerged[is.na(dfmerged[left]),right]
    dfmerged[right]<- NULL
    colnames(dfmerged)[colnames(dfmerged) == left] <- i
  }
  dfmerged
}



AddParPed <- function(PedIN) {
  Ped <- unique(PedIN)
  UID <- unique(c(as.character(Ped[,1]),
                  as.character(Ped[,2]),
                  as.character(Ped[,3])))
  UID <- stats::na.exclude(UID)
  if (length(UID) > nrow(Ped)) {
    AddPed <- data.frame(id=setdiff(UID, Ped[,1]),
                         dam=NA,
                         sire=NA,
                         stringsAsFactors=FALSE)
    Ped <- merge(AddPed, PedIN, all=TRUE)  # presume ancestors
  }
  Ped
}

