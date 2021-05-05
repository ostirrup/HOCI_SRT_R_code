
library("PostcodesioR")
library("bbmle")
library("boot")
library("gmt")
library("ape")
library("parallel")

setwd("//your_working_directory_containing_data_and_R_scripts")

##Read in patient meta-data, convert dates to new variables in R format 

dataset_meta<-read.csv("hociSiteMetadata_dummy.csv", stringsAsFactors=FALSE, header=TRUE)
dataset_meta$date_sample_var<-as.Date(dataset_meta$sampleDate, format="%d-%b-%y")
dataset_meta$date_adm_var<-as.Date(dataset_meta$admissionDate, format="%d-%b-%y")
dataset_meta$date_symptom_var<-as.Date(dataset_meta$symptomOnsetDate, format="%d-%b-%y")

dataset_meta$date_diag_var<- dataset_meta$date_sample_var
dataset_meta[!is.na(dataset_meta$date_symptom_var),"date_diag_var"]<- pmin(dataset_meta[!is.na(dataset_meta$date_symptom_var),"date_sample_var"], dataset_meta[!is.na(dataset_meta$date_symptom_var),"date_symptom_var"])

dataset_meta[is.na(dataset_meta$unitHistoryPost),"unitHistoryPost"]<-""

### Process postcode data

dataset_meta$longitude<-0
dataset_meta$latitude<-0

for (i in 1:nrow(dataset_meta)) {
	outcode_result_tmp<-outward_code_lookup(dataset_meta[i,"residentialOuterPostcode"])
	if(!is.null(outcode_result_tmp)){
		dataset_meta[i,"longitude"]<-outcode_result_tmp$longitude
		dataset_meta[i,"latitude"]<-outcode_result_tmp$latitude
		rm(outcode_result_tmp)} else {
		dataset_meta[i,"longitude"]<-NA
		dataset_meta[i,"latitude"]<-NA
		print(dataset_meta[i,"residentialOuterPostcode"])
		}
	}


#ggplot(data=dataset_meta,  aes(x=longitude, y=latitude)) +	geom_jitter(aes(x=longitude, y=latitude, color=admissionStatus), width = 0.05, height = 0.05)

## Process sequence data

seq_data<-read.FASTA("hociSiteSequences_MSA_dummy.fasta", type="DNA")
	
#Create subset of cases with both available meta-data and sequence data
dataset_meta$seq_ind<- (dataset_meta$sequenceID %in% names(seq_data))
table(dataset_meta$seq_ind)
dataset_meta<-dataset_meta[dataset_meta$seq_ind==1,]
seq_data<-seq_data[dataset_meta$sequenceID]

seq_data_matrix<-as.character(as.matrix(seq_data))
table(seq_data_matrix[1,], useNA="always")
table(as.character(seq_data[[1]]))
for(i in 1:nrow(seq_data_matrix)) {seq_data_matrix[i,seq_data_matrix[i,]=="-"]<-"n"}	#Replace any dashes with 'n'
seq_data<-as.DNAbin(as.list(data.frame(t(seq_data_matrix),stringsAsFactors=FALSE))) 
table(seq_data_matrix[1,], useNA="always")
table(as.character(seq_data[[1]]))
names(seq_data)<-gsub("\\.","-",names(seq_data))
head(names(seq_data))

## Drop any sequences with >10% 'n' codes from meta-data and sequence data objects
N_count_vec<-numeric()
for (i in 1:length(seq_data)) {
	N_count_vec[i]<- sum(as.character(seq_data[[i]])=="n")
	}
summary(N_count_vec)
sum(N_count_vec>(0.1*29907))

coverage_exclusions<-names(seq_data)[N_count_vec>(0.1*29907)]
dataset_meta<-dataset_meta[!dataset_meta$sequenceID %in% coverage_exclusions, ]
seq_data<-seq_data[dataset_meta$sequenceID]


## Create full pairwise distance matrix (can take a few minutes for hundreds or thousands of sequences)
date()
full_distance_matrix <- dist.dna(seq_data, model="N", pairwise.deletion = TRUE, as.matrix=TRUE)
date()

#SNP_difs<-ggplot(data=data.frame(dist=as.numeric(full_distance_matrix[lower.tri(full_distance_matrix, diag = FALSE)])), aes(x=dist)) +
#	geom_histogram(binwidth=1) + xlab("Pairwise SNP distance") + ylab("Count") +
#	theme_bw()

## Summarise proportion of pairwise SNP differences at different thresholds
median(as.numeric(full_distance_matrix[lower.tri(full_distance_matrix, diag = FALSE)]))
100*mean(as.numeric(full_distance_matrix[lower.tri(full_distance_matrix, diag = FALSE)])==0)
100*mean(as.numeric(full_distance_matrix[lower.tri(full_distance_matrix, diag = FALSE)])<=1)
100*mean(as.numeric(full_distance_matrix[lower.tri(full_distance_matrix, diag = FALSE)])<=2)
100*mean(as.numeric(full_distance_matrix[lower.tri(full_distance_matrix, diag = FALSE)])<=3)


### Geo calibration model

source("geo_cluster_fun_SNP_para.R")

community_inclusion_ind <- (dataset_meta$institutionID=="" | (dataset_meta$admissionStatus %in% c("not_a_patient","a_and_e_patient"))) & (dataset_meta$isHealthcareWorker!=TRUE | is.na(dataset_meta$isHealthcareWorker))

dataset_com <- dataset_meta[community_inclusion_ind ,]
nrow(dataset_com)
SNP_com_dist_mat<- full_distance_matrix[dimnames(full_distance_matrix)[[1]] %in% dataset_com$sequenceID, dimnames(full_distance_matrix)[[2]] %in% dataset_com$sequenceID]
str(SNP_com_dist_mat)
sum(dataset_com$sequenceID != dimnames(SNP_com_dist_mat)[[1]])
summary(dataset_com$date_sample_var)

SNP_threshold<-2
options(mc.cores=1)

date()
m0 <- mle2(geo_cluster_fun_SNP_para, start=list(log_tau=0, logit_beta=0), data=list(dataset_com=dataset_com, SNP_com_dist_mat=SNP_com_dist_mat, SNP_threshold=SNP_threshold))
date()

beta_estimate<-unname(inv.logit(m0@coef["logit_beta"]))
tau_estimate<-unname(exp(m0@coef["log_tau"]))

### Process data and submit cases to algorithm

dataset_meta$focus_ind <- (dataset_meta$admissionStatus=="inpatient" & (as.numeric(dataset_meta$date_diag_var-dataset_meta$date_adm_var)>2) & (dataset_meta$unitID!="None") & (dataset_meta$unitID!=""))
dataset_focus<-dataset_meta[dataset_meta$focus_ind==1,]
dataset_focus<-dataset_focus[order(dataset_focus$date_diag_var),]

dataset_focus$t_adm_to_diag <- as.numeric(dataset_focus$date_diag_var - dataset_focus$date_adm_var)
summary(dataset_focus$t_adm_to_diag)
quantile(dataset_focus$t_adm_to_diag, na.rm=TRUE)
table(dataset_focus$residentialOuterPostcode)

dataset_focus$posterior_w<- (-1)
dataset_focus$posterior_h<- (-1)
dataset_focus$posterior_v<- (-1)
dataset_focus$prior_prob<- (-1)	
dataset_focus$ward_cond_prop<- (-1)
dataset_focus$hos_cond_prop<- (-1)
dataset_focus$com_cond_prop<- (-1)
dataset_focus$com_match_list<- ""
dataset_focus$ward_match_list<- ""
dataset_focus$hos_match_list<- ""
dataset_focus$hos_locmatch_list<- ""

dataset_focus$ward_ref_n<- 0
dataset_focus$hos_ref_n<- 0
dataset_focus$com_ref_n<- 0

dataset_focus$HCW_match_n<- 0

dataset_focus$posterior_notes<- ""

# Model parameters
SNP_threshold<-2		#Threshold for a close sequence match
Pw<-0.7			#Prior probability of infection on ward/unit given post-admission infection if no visitors allowed
Pv<-0.2			#Prior probability of infection from visitor given post-admission infection if visitors allowed
tCutHos<-21			#Look back period for sequences within hospital (days)
tCutCom<-42			#Look back period for sequences within community (days)


source("gen_posterior_fun_SNP_M2_locmatch_2.R")	

nrow(dataset_focus)
for(i in 1:nrow(dataset_focus)) {
	dataset_focus<-gen_posterior_fun_SNP_M2_locmatch_2(i=i, beta=beta_estimate, tau=tau_estimate, Pa=0.4, c_days=3, dataset_full=dataset_meta, dataset_focus=dataset_focus, seq_data=seq_data, SNP_threshold=SNP_threshold, Pw=Pw, Pv=Pv,  tCutHos, tCutCom)
	print(i)
	}


dataset_focus$posterior_total <- dataset_focus$posterior_w + dataset_focus$posterior_h + dataset_focus$posterior_v

write.csv(dataset_focus, file="dataset_focus_results.csv")

# Output variables summary
#posterior_w: posterior probability of infection on ward
#posterior_h: posterior probability of infection elsewhere in hospital
#posterior_v: posterior probability of infection from visitor to ward
#prior_prob: prior probability of post-admission infection (based on time from admission to sample/symptom date)
#ward_cond_prop: proportion close match sequences on ward in last 3 weeks (matching sequences in ward_match_list, denominator ward_ref_n)
#hos_cond_prop: proportion close match sequences in rest of hospital in last 3 weeks (matching sequences in hos_match_list, denominator hos_ref_n)
#com_cond_prop: weighted proportion close match sequences in local community over last 6 weeks (centred on postcode) (matching sequences in com_match_list)
#posterior_notes : at the moment, just flags if no matching sequence anywhere in reference datasets (prior probabilities are then returned)
#posterior_total: posterior_w + posterior_h + posterior_v
#ward_match_list: Sequence matches with the same unitID.
#hos_locmatch_list: Any matches identified by cross-referencing pre- and post-sampling locations (not counting anything already in the ward_match_list).
#ward_ref_n: n sequences in unit reference set
#hos_ref_n: n sequences in hospital reference set
#com_ref_n: n sequences in community reference set
#HCW_match_n: n close sequence matches to HCWs



