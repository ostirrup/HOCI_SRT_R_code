
geo_cluster_fun_SNP_para<- function(log_tau, logit_beta) {

	tau<-exp(log_tau)
	beta<-inv.logit(logit_beta)
	p_vec<-numeric()


	para_fun<-function(i){
		distance_to_index_vec<-numeric()
		weighting_vec<-numeric()
		prior_dataset_indvec<- (dataset_com$date_sample_var<=dataset_com[i,"date_sample_var"] & dataset_com$sequenceID!=dataset_com[i,"sequenceID"] & (as.numeric((dataset_com[i,"date_sample_var"]-dataset_com$date_sample_var)<=42) ) ) #NOTE '<=' HERE
		prior_dataset_indvec_which<-which(prior_dataset_indvec==1)
		if(sum(prior_dataset_indvec)>1){
			SNP_match_indvec_tmp<- SNP_com_dist_mat[i,]<=SNP_threshold
			SNP_match_indvec<- SNP_match_indvec_tmp[prior_dataset_indvec]
			if(sum(SNP_match_indvec>0)) {
				for(j in 1:length(prior_dataset_indvec_which)) {
					ind_tmp<- prior_dataset_indvec_which[j]
					if(dataset_com[i,"residentialOuterPostcode"]==dataset_com[ind_tmp,"residentialOuterPostcode"]) {
						distance_to_index_vec[j]<- 0
						} else {
						distance_to_index_vec[j]<- geodist(Nfrom=dataset_com[i,"latitude"], Efrom=dataset_com[i,"longitude"], Nto=dataset_com[ind_tmp,"latitude"], Eto=dataset_com[ind_tmp,"longitude"], units="km")}
					}
				weighting_vec<- (1-beta)*exp(-tau*distance_to_index_vec) + beta
				p_calc <- (sum(SNP_match_indvec*weighting_vec)) / (sum(weighting_vec))
				return(p_calc)
				} else {return(NA)}
			} else {return(NA)}
		}

	p_results<-mclapply(1:nrow(dataset_com), para_fun)
	p_results<-as.numeric(p_results)
	p_results<-p_results[!is.na(p_results)]

	log_lik<- sum(log(p_results))
	return(-log_lik)
	}

