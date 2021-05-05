
gen_posterior_fun_SNP_M2_locmatch_2<- function(i, beta, tau, Pa, c_days, dataset_full, dataset_focus, seq_data, SNP_threshold, Pw, Pv, tCutHos, tCutCom) {

	indexID<-dataset_focus[i, "sequenceID"]
	indexSampleDate<-dataset_focus[i, "date_sample_var"] 
	indexDiagnosisDate<-dataset_focus[i, "date_diag_var"] 
	indexLocVec<- unique(unlist(c(dataset_focus[i, "unitID"], strsplit(dataset_focus[i, "unitHistoryPre"],"\\|"))))

	community_ref_set_ind<- (dataset_full$date_sample_var<=indexSampleDate) &
							(dataset_full$sequenceID!=indexID)	& 
							(as.numeric(indexSampleDate - dataset_full$date_sample_var)<=tCutCom)	& (
									(dataset_full$institutionID=="") |
									(dataset_full$institutionID!="" & (dataset_full$admissionStatus %in% c("not_a_patient","outpatient","a_and_e_patient"))) |
									(dataset_full$institutionID!="" & dataset_full$admissionStatus=="inpatient" & as.numeric(dataset_full$date_diag_var - dataset_full$date_adm_var)<=2)
									)
	
	dataset_com<-dataset_full[community_ref_set_ind,]
	
	posterior_w<- (-1)
	posterior_h<- (-1)
	posterior_v<- (-1)

	posterior_notes<-""

	ward_match_list<-""
	hos_match_list<-""
	hos_locmatch_list<-""
	com_match_list<-""
	
	ward_ref_n<-0
	hos_ref_n<-0
	com_ref_n<-0
	
	HCW_match_n<-0
	
	ward_cond_prop<- (-1)
	hos_cond_prop<- (-1)
	com_cond_prop<- (-1)
	

	if (dataset_focus[i, "visitorsOnWard"]==TRUE) {
													Pw<-Pw*(1-Pv)
													Pv<-Pv
													}

	if (dataset_focus[i, "visitorsOnWard"]==FALSE) {
													Pw<-Pw
													Pv<-0
													}
													
	if (dataset_focus[i, "admissionStatus"]=="inpatient") {
			if(dataset_focus[i, "symptomatic"]==FALSE){
				prior_prob<- ((1-Pa)*plnorm((as.numeric(indexSampleDate - dataset_focus[i, "date_adm_var"])+c_days) , meanlog = 1.621, sdlog = 0.418, lower.tail = TRUE, log.p = FALSE) + Pa*plnorm(as.numeric(indexSampleDate - dataset_focus[i, "date_adm_var"]) , meanlog = 1.621, sdlog = 0.418, lower.tail = TRUE, log.p = FALSE))
				} else {
				prior_prob<-plnorm(as.numeric(indexDiagnosisDate - dataset_focus[i, "date_adm_var"]) , meanlog = 1.621, sdlog = 0.418, lower.tail = TRUE, log.p = FALSE)
				}
			} else if (dataset_focus[i, "admissionStatus"]=="outpatient") {
			prior_prob<- 0.5
			} else (stop("incorrect patient type"))
	
	##ward level calculations
	ward_ref_set_ind<-(dataset_full$date_sample_var<=indexSampleDate) &
							(dataset_full$sequenceID!=indexID)	& 
							(as.numeric(indexSampleDate - dataset_full$date_sample_var)<=tCutHos)	& 
							(dataset_full$institutionID==dataset_focus[i, "institutionID"])	& 
							(dataset_full$unitID==dataset_focus[i, "unitID"])
	ward_ref_data<- dataset_full[ward_ref_set_ind,]
	
	if(nrow(ward_ref_data)>0){
		ward_ref_data$seq_dif<- (-1)
		for (j in 1:nrow(ward_ref_data)){
				ward_ref_data[j,"seq_dif"]<- as.numeric(ape::dist.dna(seq_data[which(names(seq_data) %in% c(ward_ref_data[j,"sequenceID"], indexID))], model="N", pairwise.deletion = TRUE, as.matrix=TRUE)[2])
				}
		ward_ref_data$match_ind<-ward_ref_data$seq_dif<=SNP_threshold
		ward_cond_prop<-mean(ward_ref_data$match_ind)
		ward_match_list<-paste(ward_ref_data[ward_ref_data$match_ind,"sequenceID"],sep="", collapse=", ")
		} else {ward_cond_prop<-0}
	ward_ref_n<-nrow(ward_ref_data)

	
	##hos level calculations
	hos_ref_set_ind<-(dataset_full$date_sample_var<=indexSampleDate) &
							(dataset_full$sequenceID!=indexID)	& 
							(as.numeric(indexSampleDate - dataset_full$date_sample_var)<=tCutHos)	& 
							(dataset_full$institutionID==dataset_focus[i, "institutionID"])	&
							(dataset_full$unitID!=dataset_focus[i, "unitID"])	& (
								(dataset_full$admissionStatus=="inpatient" & as.numeric(dataset_full$date_diag_var - dataset_full$date_adm_var)>2) |
								(as.logical(dataset_full$isHealthcareWorker)==TRUE & dataset_full$hcwWorkplaceID==dataset_focus[i, "institutionID"]) |
								(dataset_full$admissionStatus=="outpatient") |
								(dataset_full$admissionStatus=="a_and_e_patient" & as.numeric(indexSampleDate - dataset_full$date_sample_var)>2) |
								(dataset_full$admissionStatus=="inpatient" & as.numeric(dataset_full$date_diag_var - dataset_full$date_adm_var)<=2 & as.numeric(indexSampleDate - dataset_full$date_sample_var)>2)
								)
	
	hos_ref_data<- dataset_full[hos_ref_set_ind,]
	if(nrow(hos_ref_data)>0){
		hos_ref_data$seq_dif<- (-1)
		for (j in 1:nrow(hos_ref_data)){
				hos_ref_data[j,"seq_dif"]<- as.numeric(ape::dist.dna(seq_data[which(names(seq_data) %in% c(hos_ref_data[j,"sequenceID"], indexID))], model="N", pairwise.deletion = TRUE, as.matrix=TRUE)[2])
				}
		hos_ref_data$match_ind<-hos_ref_data$seq_dif<=SNP_threshold
		hos_cond_prop<-mean(hos_ref_data$match_ind)
		HCW_match_n<-sum(as.logical(hos_ref_data[hos_ref_data$match_ind==1,"isHealthcareWorker"]))
		hos_match_list<-paste(hos_ref_data[hos_ref_data$match_ind,"sequenceID"],sep="", collapse=", ")
		hos_locmatch_list<-""
		counter<-1	
		for (seq_match in hos_ref_data[hos_ref_data$match_ind,"sequenceID"]) {
				seqLocVec<- unique(unlist(c(hos_ref_data[hos_ref_data$sequenceID==seq_match, "unitID"], strsplit(hos_ref_data[hos_ref_data$sequenceID==seq_match, "unitHistoryPre"],"\\|"), strsplit(hos_ref_data[hos_ref_data$sequenceID==seq_match, "unitHistoryPost"],"\\|"))))
				matching_locations<- indexLocVec[indexLocVec %in% seqLocVec]
				if(length(matching_locations)>0) {hos_locmatch_list[counter]<- paste(seq_match, " in ", matching_locations, sep="", collapse=". ")
													counter<-counter+1
													rm(matching_locations)}
				}
		hos_locmatch_list<-paste(hos_locmatch_list, sep="", collapse=". ")		
		} else {hos_cond_prop<-0}
	hos_ref_n<-nrow(hos_ref_data)
	

	##com level calculations
	com_ref_data<- dataset_com
	if(nrow(com_ref_data)>0){
		com_ref_data$distance_to_index<- (100)
		com_ref_data$seq_dif<- (-1)
		for (j in 1:nrow(com_ref_data)){
				com_ref_data[j,"seq_dif"]<- as.numeric(ape::dist.dna(seq_data[which(names(seq_data) %in% c(com_ref_data[j,"sequenceID"], indexID))], model="N", pairwise.deletion = TRUE, as.matrix=TRUE)[2])
				if(com_ref_data[j,"residentialOuterPostcode"]==dataset_focus[i,"residentialOuterPostcode"]) {
					com_ref_data[j,"distance_to_index"]<- 0
					} else if (is.na(com_ref_data[j,"latitude"])) {
					com_ref_data[j,"distance_to_index"]<- 100
					} else {com_ref_data[j,"distance_to_index"]<-geodist(Nfrom=dataset_focus[i,"latitude"], Efrom=dataset_focus[i,"longitude"], Nto=com_ref_data[j,"latitude"], Eto=com_ref_data[j,"longitude"], units="km")}
				}
		com_ref_data$weighting<-(1-beta)*exp(-tau*com_ref_data$distance_to_index) + beta
		com_ref_data$match_ind<-com_ref_data$seq_dif<=SNP_threshold
		com_cond_prop<- sum(com_ref_data$weighting*com_ref_data$match_ind)/sum(com_ref_data$weighting)
		com_match_list<-paste(com_ref_data[com_ref_data$match_ind,"sequenceID"],sep="", collapse=", ")
		} else {com_cond_prop<-0}
	com_ref_n<-nrow(com_ref_data)


	##prob calcs
	if (ward_cond_prop>0 | hos_cond_prop>0 | com_cond_prop>0) {
	post_prob_denom <- (1-prior_prob)*com_cond_prop + prior_prob*Pw*ward_cond_prop + prior_prob*Pv*com_cond_prop + prior_prob*(1-Pw-Pv)*hos_cond_prop
	posterior_w <- prior_prob*Pw*ward_cond_prop / post_prob_denom
	posterior_h <- prior_prob*(1-Pw-Pv)*hos_cond_prop / post_prob_denom
	posterior_v <- prior_prob*Pv*com_cond_prop / post_prob_denom
	if(post_prob_denom==0) {posterior_notes<-"Prior indicates post-admission infection, but no sequence matches in institution"}
	} else {
	posterior_w <- prior_prob*Pw
	posterior_h <- prior_prob*(1-Pw-Pv)
	posterior_v <- prior_prob*Pv
	posterior_notes<-"No sequence match in any dataset"
	}
	
######################
	
	dataset_focus[i,"prior_prob"]<-prior_prob
	dataset_focus[i,"posterior_w"]<-posterior_w
	dataset_focus[i,"posterior_h"]<-posterior_h
	dataset_focus[i,"posterior_v"]<-posterior_v
	dataset_focus[i,"posterior_notes"]<-posterior_notes
	dataset_focus[i,"ward_cond_prop"]<-ward_cond_prop
	dataset_focus[i,"hos_cond_prop"]<-hos_cond_prop
	dataset_focus[i,"com_cond_prop"]<-com_cond_prop
	

	dataset_focus[i,"ward_match_list"]<-ward_match_list
	dataset_focus[i,"hos_match_list"]<-hos_match_list
	dataset_focus[i,"hos_locmatch_list"]<-hos_locmatch_list
	dataset_focus[i,"com_match_list"]<-com_match_list
	
	dataset_focus[i,"ward_ref_n"]<-ward_ref_n
	dataset_focus[i,"hos_ref_n"]<-hos_ref_n
	dataset_focus[i,"com_ref_n"]<-com_ref_n
	
	dataset_focus[i,"HCW_match_n"]<-HCW_match_n


	return(dataset_focus)
}

