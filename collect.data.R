#' set some needed directories KORA data directory
KORA.DIR <- "/storage/groups/groups_epigenereg/analyses/PV_K14115g_Heinig/"
F.METH <- paste0(KORA.DIR, "data/20160204/KORAF4_illuminamethylation450k_qn_bmiq_n1727/KF4_beta_qn_bmiq.RData");
F.EXPR <- paste0(KORA.DIR, "data/20160204/Expression/kora_f4_normalized.Rdata");
F.COVA <- paste0(KORA.DIR, "data/20160204/Expression/technical_covariables_kora_f4.Rdata");
F.IDMAP <- paste0(KORA.DIR, "data/20160204/individuals_covariates.csv")

# path to indexed genotype file (tabix)
F.SNPS <- paste0(KORA.DIR, "results/20160204/genoF4/dosage_combined/MAF001/full_sorted.bgz")

#' set some needed directories KORA data directory
KORA.DIR <- "/media/data/Masterarbeit/data/F4/"
F.METH <- paste0(KORA.DIR, "KORAF4_illuminamethylation450k_qn_bmiq_n1727/KF4_beta_qn_bmiq.RData");
F.EXPR <- paste0(KORA.DIR, "Expression/kora_f4_normalized.Rdata");
F.COVA <- paste0(KORA.DIR, "Expression/technical_covariables_kora_f4.Rdata");
F.IDMAP <- paste0(KORA.DIR, "individuals_covariates.csv")

# path to indexed genotype file (tabix)
F.SNPS <- paste0(KORA.DIR, "results/20160204/genoF4/dosage_combined/MAF001/full_sorted.bgz")


# sanity check when sourcing the script
if(!dir.exists(KORA.DIR)){
    message("Directory ", KORA.DIR, " does not exist. Be careful using the methods!");
} else {
    message("Using ", KORA.DIR, " as data directory.");
}


#'
#' Collects methylation, expression and genotype data for a KORA cohort
#' 
#' @param snp.ranges SNP ranges for which to get the genotype data
#' @param meth.probes Optional list of methylation array probe ids to be retrieved
#' @param expr.probes Optional list of expression array probe ids to be retrieved
#' @param cache.global Flag whether to cache the expression and methylation data
#' in the global environment and reuse them when calling the method again.
#' 
#' @return A data matrix containing methylation, expression and genotype data 
#' for a specific sentinel
#' 
#' @autho Johann Hawe
#' 
#' @date 02/06/2017
#' 
#' @export
#' 
collect.data <- function(snp.ranges=NULL, meth.probes=NULL, expr.probes=NULL, 
                         cache.global=F) {

	if(is.null(snp.ranges)){
	 stop("SNP ranges need to be supplied when loading KORA data.")
	}
	
	cat("Loading KORA data.\n")

	# load id mapping to intersect expr/meth/geno data
	ID.MAP <- read.table(F.IDMAP, 
					header=T, sep=";", stringsAsFactors = F)
	# convert ids to character (instead of int) for easier subsetting of data matrices
	ID.MAP$axio_s4f4 <- as.character(ID.MAP$axio_s4f4)
	ID.MAP$expr_s4f4ogtt <- as.character(ID.MAP$expr_s4f4ogtt)
	ID.MAP$meth_f4 <- as.character(ID.MAP$meth_f4)
	# drop individuals with NAs (e.g. due to missing BMI)
	ID.MAP <- na.omit(ID.MAP)
	# load dummy snp to get the individuals which have geno data available
	dummy <- GRanges("chr4", ranges=IRanges(156902056,width=1))
	dummy <- t(scanSNPs(dummy)$snps);

	if(!exists("f4.norm") | !exists("beta")) {
	 cat("Preparing raw data.\n");
	 # defined in sourced kora common script
	 load(F.EXPR); load(F.METH);
	 
	 if(cache.global){
	   assign("f4.norm", f4.norm, .GlobalEnv)
	   assign("beta", beta, .GlobalEnv)
	 }
	} 

	# gets as 687 individuals, having all data available (some ids of the id 
	# map are not contained within the data frame...)
	toUse <- which(ID.MAP$expr_s4f4ogtt %in% colnames(f4.norm) & 
				 ID.MAP$meth_f4 %in% colnames(beta) & 
				 ID.MAP$axio_s4f4 %in% rownames(dummy));

	ID.MAP <- ID.MAP[toUse,];
	ID.MAP$utbmi <- NULL
	ID.MAP$ul_wbc <- NULL
	ID.MAP$utalteru <- NULL

	cat("Using ", nrow(ID.MAP), " samples.\n");

	# sort our input data s.t. each row corresponds to the same individual 
	# using the created ID mapping table
	if(!is.na(meth.probes)){
	 meth <- t(beta[,ID.MAP$meth_f4]);
	 if(!is.null(meth.probes)){
	   meth <- meth[,meth.probes,drop=F];
	 }
	 # get rid of zero-variance probes
	 meth <- meth[,apply(meth,2,var, na.rm=T)!=0, drop=F]
	}

	if(!is.na(expr.probes)){
	 # use only those individuals for which we have all data available
	 expr <- t(f4.norm[,ID.MAP$expr_s4f4ogtt]);
	 if(!is.null(expr.probes)){
	   expr <- expr[,expr.probes,drop=F];
	 }
	 # remove zeor-variance probes...
	 expr <- expr[,apply(expr,2,var)!=0, drop=F]
	}
	if(!is.na(snp.ranges)){
	 geno <- get.genotypes(snp.ranges, snp.ids, cohort)
	 geno <- geno[ID.MAP$axio_s4f4,,drop=F]
	}
	# get the methylation PCA results
	load(paste0(DATA.DIR, "/kora/methylation/control_probe_pcs_n1727.RData"))
	pcs <- pcs[ID.MAP$meth_f4,]
	colnames(pcs) <- paste(colnames(pcs), "cp", sep="_")
	meth.pcs <- pcs
	rm(pcs);

	# load technical covariates for expression data
	load(paste0(DATA.DIR, "/kora/expression/technical_covariables_kora_f4.Rdata"));
	rownames(covars.f4) <- covars.f4[,1];
	covars.f4 <- covars.f4[ID.MAP$expr_s4f4ogtt,c(2:6)];
	covars.f4$sex <- as.factor(covars.f4$sex);

	# load houseman blood count data for methylation
	houseman <- read.table(paste0(DATA.DIR, 
						"/kora/methylation/Houseman/KF4_QN_BMIQ_estimated_cell_distribution_meanimpute473_lessThanOneFALSE.csv"), 
					  sep=";", header=T,row.names=1);
	houseman <- houseman[ID.MAP$meth_f4,];

	# create initial data frame with covariates
	covars.f4 <- cbind(as.data.frame(covars.f4), houseman, meth.pcs);

	# replace some colnames to easily match with the lolipop data
	cc <- colnames(covars.f4)
	cc[grepl("storage.time",cc)] <- "batch2"
	cc[grepl("plate",cc)] <- "batch1"
	colnames(covars.f4) <- cc
	covars.f4[,"batch1"] <- factor(covars.f4[,"batch1"])

	data <- cbind.data.frame(covars.f4, stringsAsFactors=F)
	if(!is.na(meth.probes)) {
	 data <- cbind.data.frame(data, meth, stringsAsFactors=F)
	}
	if(!is.na(expr.probes)){
	 data <- cbind.data.frame(data, expr, stringsAsFactors=F)
	} 
	if(!is.na(snp.ranges)) {
	 data <- cbind.data.frame(data, geno, geno.ids=rownames(geno),
						 stringsAsFactors=F)
	}  
	 
	return(data);    
}


#'
#' Scans prepared genotype files for SNPs within the provided genomic ranges.
#' TODO: can be improved by using differenti method instead of scanTabix, then
#' the output would not have been to be parsed
#'
#' @param ranges GRanges object containing the ranges for which to get the SNPs
#'
scanSNPs <- function(ranges) {
	require(Rsamtools);
	require(data.table);

	# rather naive implementation for now	
	snpTabix <- TabixFile(file=F.SNPS);
	result <- scanTabix(snpTabix, param=ranges);
	counts <- countTabix(snpTabix, param=ranges);
	result <- result[names(counts[counts>0])];
	result <- unlist(result);

	if(!is.null(result)) {
		# get data frame
		if(length(result) > 4000){
			cuts <- seq(0,length(result), by=floor(length(result)/50));
			data <- c();
			for(i in 2:length(cuts)){
				temp <- fread(paste(result[(cuts[i-1]+1):cuts[i]],collapse="\n"), 
				              sep="\t", data.table=F);
				data <- rbind(data,temp);
				message(paste0("read ", i-1, "/", length(cuts)-1, " chunks."));
			}
			temp <- read.table(text=result[(cuts[length(cuts)]+1):length(result)], 
			                   sep="\t", colClasses="character");
			data <- rbind(data,temp);
			rm(temp);
		}
		else {
			data <- read.table(text=result, sep="\t", colClasses="character");
		}
		if(nrow(data) != length(result)){
			stop("sanity check failed.");
		}
		data <- unique(data);
		message(paste("Processed", nrow(data), "SNPs." ));

		rm(result);
		rm(counts);

		# process the genotype information to get integers/factors
		for(i in 6:ncol(data)){
			data[,i] <- factor(round(as.numeric(data[,i])));
		}

		#create colnames using individual codes
		ids <- read.table(paste0(KORA.DIR,"/results/current/genoF4/individuals.txt"), stringsAsFactors=F, colClasses="character");
		colnames(data)<- c("chr", "name", "pos", "orig", "alt", ids[,1])
		rownames(data) <- data$name;

		return(list(snpInfo=data[,c(1,3,4,5)], snps=data[,6:ncol(data)]));
	} else {
		cat("No SNPs found in specified regions.\n")
		return(list());
	}
}
