# SASKMI
A SAS Macro to Perform Kaplan-Meier Multiple Imputation for Survival Analyses with Competing Events
# ABSTRACT
Analysis of cumulative incidence functions require special attention when there are competing events. A competing event is an event that either precludes or changes the probability of the event of interest from occurring. There are different approaches to handling competing events in survival analysis. One approach, developed by Fine and Gray (1999), cannot be applied directly, for example, to stratified data or to data with time-varying covariates. For these situations, Ruan and Gray (2008) proposed an alternative approach (Kaplan-Meier multiple imputation (KMI)) which recovers the missing censoring times for those who experienced a competing event. The missing censoring times are imputed from a non-parametric multiple-imputation approach based on the Kaplan-Meier estimator. In this paper, we introduce a user-friendly SAS Macro (%SASKMI) to implement this approach. %SASKMI generates a dataset with a new event status and imputed times that can be analyzed using a regular Cox model (PROC PHREG).  The output dataset follows PROC MI data standards so PROC MIANALYSIS can be used to summarize results. To demonstrate the effectiveness of the new macro, using a real data example we compare the effect estimates, standard errors, and run times obtained after applying %SASKMI to that of a standard Fine and Gray model as implemented in PROC PHREG. We find that %SASKMI performs similarly to the Fine and Gray option in PROC PHREG but with significant run time reduction.
# INSTRUCTION of USING %SASKMI
The %SASKMI macro implements the KMI approach to impute potential censoring times for those individuals with a competing event. The macro is called with the parameters shown in Table 1, some of which are required. It is recommended to include covariates, such as age, sex and race, that potentially correlate with the censoring distribution. When covariates are specified, then a proportional hazard model for the censoring distribution is fit and the censoring time imputed from the estimated model. Otherwise, the time is imputed from the KM estimator of the observed censoring distribution.
# DATA=	(Required) SAS input data set name
# EVENT= 	(Required) The variable name of the outcome of interest. This variable usually coded as 0, 1, or 2 where 0 denotes a censoring event, 1 denotes the outcome of interest and 2 denotes the competing event. 
# EVENTCODE=	(Required) The outcome of interest code in EVENT variable, default is 1. 
# CENSCODE=	(Required) The censoring code in EVENT variable, default is 0.
# TIME=	(Required) The time to event variable. 
# ADJVAR=	List of covariates to be used to model the censoring distribution. If no variables included, censoring times are then imputed from unadjusted Kaplan-Meier survivor estimator.
# CLASS=	(Required if category variable exist in ADJVAR) All categorical covariates from ADJVAR must be included here. 
# NIMP=	(Required) The number of imputations.
# SEED=	Random control seed. The default is 123. Since the KMI approach is based on randomly drawing censoring times from the estimated conditional censoring distribution, it’s extremely important to specify SEED in order to replicate results. 
# OUT=	The name of the output dataset. The default is “OUT”.

# EXAMPLE
The following code runs the %SASKMI macro with 10 imputations. The KMI model includes the baseline characteristic: AGE, RACE, and REGION. 
%SASKMI(
	DATA		= Gnstudy, 
	ADJVAR	= AGE RACE REGION,
	CLASS 	= RACE REGION,
	EVENT 	= CV_DEATH, 
	FAILCODE	= 1, 
  CENSCODE	= 0, 
	TIME 		= TIME2CV_DEATH, 
	NIMP 		= 10, 
	SEED 		= 123, 
	OUT 		= OUT
);

