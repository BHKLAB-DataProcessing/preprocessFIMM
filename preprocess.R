library(PharmacoGx)
library(readxl)


options(stringsAsFactors=FALSE)

matchToIDTable <- function(ids,tbl, column, returnColumn="unique.cellid") {
	sapply(ids, function(x) {
                          myx <- grep(paste0("((///)|^)",Hmisc::escapeRegex(x),"((///)|$)"), tbl[,column])
                          if(length(myx) > 1){
                            stop("Something went wrong in curating ids, we have multiple matches")
                          }
			  if(length(myx) == 0){return(NA_character_)}
                          return(tbl[myx, returnColumn])
                        })
}
symSetDiff <- function(x,y) return(setdiff(union(x,y), intersect(x,y)))

## Load in annotation files

cell.all <- read.csv("/pfs/downAnnotations/cell_annotation_all.csv")
drug.all <- read.csv("/pfs/downAnnotations/drugs_with_ids.csv")

# message("pfs files:")
# message(list.files("/pfs/"))

# message("input files:")

# message(list.files("/pfs/FIMMdata/"))


# file.copy("/pfs/FIMMdata/nature20171-s1.xls", "/pfs/out/nature20171-s1.xls")
# file.copy("/pfs/FIMMdata/nature20171-s2.xlsx", "/pfs/out/nature20171-s2.xlsx")

sens.processed <- as.data.frame(read_excel("/pfs/downFIMM/nature20171-s1.xls"))
sens.raw <- as.data.frame(read_excel("/pfs/downFIMM/nature20171-s2.xlsx"))


## Creating curation cell and drug table to resolve mismatches between processed and raw data

curationCell <- cell.all[!is.na(cell.all[,"FIMM.cellid"]),c("unique.cellid", "FIMM.cellid")]

curationDrug <- drug.all[!is.na(drug.all[,"FIMM.drugid"]),c("unique.drugid", "FIMM.drugid")]

curationTissue <- cell.all[!is.na(cell.all[,"FIMM.cellid"]),c("unique.tissueid", "FIMM.tissueid")]

rownames(curationCell) <- curationCell[,"unique.cellid"]

rownames(curationTissue) <- curationCell[,"unique.cellid"]

rownames(curationDrug) <- curationDrug[,"unique.drugid"]


# cell.ids <- colnames(sens.processed)[-c(1:4)]
# tissue.ids <- sens.processed[2,-(1:4)]

# drug.ids <- sens.processed[-c(1,2,3),3]

sens.data <- sens.processed[-1:-3,]

rownames(sens.data) <- sens.data[,3]

sens.data <- sens.data[,-1:-4]

library(reshape2)

sens.data <- apply(sens.data, c(1,2), as.numeric)

sens.data.m <- melt(sens.data)

colnames(sens.data.m) <- c("drugid", "cellid", "aac_published")

sens.data.m[,"cellid"] <- matchToIDTable(sens.data.m[,"cellid"], curationCell, "FIMM.cellid", "unique.cellid")

sens.data.m[,"drugid"] <- matchToIDTable(sens.data.m[,"drugid"], curationDrug, "FIMM.drugid", "unique.drugid")

## Remove NAs from the melted table 
sens.data.m <- sens.data.m[complete.cases(sens.data.m),]


## No replicates here?

exp.ids.profiles <- paste0(sens.data.m[,1], "_", sens.data.m[,2])

rownames(sens.data.m) <- exp.ids.profiles



sens.raw.data <- sens.raw[-1,]

colnames(sens.raw.data) <- sens.raw[1,]

colnames(sens.raw.data)[2:3] <- c("drugid", "cellid")


sens.raw.data[,"cellid"] <- matchToIDTable(sens.raw.data[,"cellid"], curationCell, "FIMM.cellid", "unique.cellid")

sens.raw.data[,"drugid"] <- matchToIDTable(sens.raw.data[,"drugid"], curationDrug, "FIMM.drugid", "unique.drugid")



exp.ids.raw <- paste0(sens.raw.data[,"drugid"], "_", sens.raw.data[,"cellid"])
rownames(sens.raw.data) <- exp.ids.raw


setequal(exp.ids.profiles, exp.ids.raw)

setdiff(exp.ids.profiles, exp.ids.raw)
setdiff(exp.ids.raw, exp.ids.profiles)

## There is more raw data than processed, so we add the two missing experiments 

new.profs <- sens.data.m[setdiff(exp.ids.raw, exp.ids.profiles),]

rownames(new.profs) <- setdiff(exp.ids.raw, exp.ids.profiles)
new.profs[,c("drugid", "cellid")] <- sens.raw.data[setdiff(exp.ids.raw, exp.ids.profiles),c("drugid", "cellid")]

sens.data.m <- rbind(sens.data.m, new.profs)

## reorder together 

sens.data.m <- sens.data.m[rownames(sens.raw.data),]

ndose <- 5

sens.raw.array <- array(dim=c(length(exp.ids.raw), ndose, 2), dimnames = list(exp.ids.raw, paste0("dose", seq_len(ndose)), c("Dose","Viability")))

for(exp in exp.ids.raw){

	sens.raw.array[exp,,"Dose"] <- as.numeric(sens.raw.data[exp, paste0("C", seq_len(ndose))])/1000
	sens.raw.array[exp,,"Viability"] <- 100 - as.numeric(sens.raw.data[exp, paste0("D", seq_len(ndose))])

}

sens.info <- sens.raw.data[,1:3]

sens.prof <- sens.data.m


## Extract provided cell info
cell.info <- data.frame(colnames(sens.processed)[-1:-4],t(sens.processed[1:3,-1:-4]))

colnames(cell.info) <- c("Cell line", sens.processed[1:3,4])

cell.info$cellid <- matchToIDTable(cell.info[,"Cell line"], curationCell, "FIMM.cellid", "unique.cellid")

rownames(cell.info) <- cell.info$cellid

cell.info$tissueid <- curationTissue[rownames(cell.info), "unique.tissueid"]

## extract provided drug info

drug.info <- sens.processed[-1:-3,1:3][,3:1]

drug.info$drugid <- matchToIDTable(drug.info[,"Drug name"], curationDrug, "FIMM.drugid", "unique.drugid")

rownames(drug.info) <- drug.info$drugid

rownames(drug.all) <- drug.all$unique.drugid

drug.info <- cbind(drug.info, drug.all[rownames(drug.info), c("smiles", "inchikey", "cid")])

saveRDS(drug.info, file="/pfs/out/drug.info.rds")
saveRDS(cell.info, file="/pfs/out/cell.info.rds")

saveRDS(curationCell, file="/pfs/out/curationCell.rds")
saveRDS(curationDrug, file="/pfs/out/curationDrug.rds")
saveRDS(curationTissue, file="/pfs/out/curationTissue.rds")

saveRDS(sens.info, file="/pfs/out/sens.info.rds")
saveRDS(sens.prof, file="/pfs/out/sens.prof.rds")

saveRDS(sens.raw.array, file="/pfs/out/sens.raw.rds")

sens.raw <- sens.raw.array

dir.create("/pfs/out/slices/")

sens.raw.x <- parallel::splitIndices(nrow(sens.raw), floor(nrow(sens.raw)/100))

for(i in seq_along(sens.raw.x)){

  slce <- sens.raw[sens.raw.x[[i]],,]
  saveRDS(slce, file=paste0("/pfs/out/slices/raw_sens_", i, ".rds"))

}





