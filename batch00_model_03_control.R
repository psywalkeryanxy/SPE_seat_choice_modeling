#xinyuanyan
#value-based DDM for seat choice

library(rstan)



setwd("D:\\博士论文写作202107\\1大论文中\\RESULTS\\6_Seat_modeling\\data_F")


#Step 1: Write Model in Stan 


hierRLDDM <-"

data {
	int<lower=1> N;									// number of data items
	int<lower=1> L;									// number of levels
	int<lower=1, upper=L> participant[N];			// level (participant)

	int<lower=-1,upper=1> accuracy[N];				// accuracy (-1, 1)
	real<lower=0> rt[N];							// rt

	vector[4] drift_priors;							// mean and sd of the mu_ prior and sd_ prior
	vector[4] threshold_priors;						// mean and sd of the mu_ prior and sd_ prior
	vector[4] ndt_priors;							// mean and sd of the mu_ prior and sd_ prior
	real<lower=0, upper=1> starting_point;			// starting point diffusion model not to estimate
}
parameters {
	real mu_drift;
	real mu_threshold;
	real mu_ndt;

	real<lower=0> sd_drift;
	real<lower=0> sd_threshold;
	real<lower=0> sd_ndt;

	real z_drift[L];
	real z_threshold[L];
	real z_ndt[L];
}
transformed parameters {
	real drift_ll[N];								// trial-by-trial drift rate for likelihood (incorporates accuracy)
	real drift_t[N];								// trial-by-trial drift rate for predictions
	real<lower=0> threshold_t[N];					// trial-by-trial threshold
	real<lower=0> ndt_t[N];							// trial-by-trial ndt

	real drift_sbj[L];
	real<lower=0> threshold_sbj[L];
	real<lower=0> ndt_sbj[L];

	real transf_mu_drift;
	real transf_mu_threshold;
	real transf_mu_ndt;

	transf_mu_drift = mu_drift;						// for the output
	transf_mu_threshold = log(1+ exp(mu_threshold));
	transf_mu_ndt = log(1 + exp(mu_ndt));

	for (l in 1:L) {
		drift_sbj[l] = mu_drift + z_drift[l]*sd_drift;
		threshold_sbj[l] = log(1 + exp(mu_threshold + z_threshold[l]*sd_threshold));
		ndt_sbj[l] = log(1 + exp(mu_ndt + z_ndt[l]*sd_ndt));
	}

	for (n in 1:N) {
		drift_t[n] = drift_sbj[participant[n]];
		drift_ll[n] = drift_t[n]*accuracy[n];
		threshold_t[n] = threshold_sbj[participant[n]];
		ndt_t[n] = ndt_sbj[participant[n]];
	}
}
model {
	mu_drift ~ normal(drift_priors[1], drift_priors[2]);
	mu_threshold ~ normal(threshold_priors[1], threshold_priors[2]);
	mu_ndt ~ normal(ndt_priors[1], ndt_priors[2]);

	sd_drift ~ normal(drift_priors[3], drift_priors[4]);
	sd_threshold ~ normal(threshold_priors[3], threshold_priors[4]);
	sd_ndt ~ normal(ndt_priors[3], ndt_priors[4]);

	z_drift ~ normal(0, 1);
	z_threshold ~ normal(0, 1);
	z_ndt ~ normal(0, 1);

	rt ~ wiener(threshold_t, ndt_t, starting_point, drift_ll);
}
generated quantities {
	vector[N] log_lik;

	{for (n in 1:N) {
		log_lik[n] = wiener_lpdf(rt[n] | threshold_t[n], ndt_t[n], starting_point, drift_ll[n]);
	}
}
}

"





#preparing data as list in R

#read data from csv files
#this data has been transformed
importdata <- read.csv(file = 'F_DDM_control_44.csv', stringsAsFactors = TRUE)




#RLdata#show the data part
RLdata_dict <-list(
  N = length(importdata$subid),

  L = length(unique(importdata$subid)),
  participant = importdata$subid,
  

  rt= importdata$rt,
  accuracy = importdata$accuracy,


  drift_priors=c(1,30,1,30),
  threshold_priors=c(1,3,0,3),
  ndt_priors =  c(1,3,0,3),
  starting_point=0.5)


for (reprep in 1:100){
  fit <- stan(model_code = hierRLDDM, data = RLdata_dict, iter = 2000, chains = 2)
  
  #m <- as.matrix(fit)
  #View(m)
  whichtime = Sys.Date()
  filename = sprintf("%s_seat_control_model_03_44_%s.rds",reprep,whichtime)
  saveRDS(fit,filename)
  
  #check
  
  
}

