/************************************************************************************************************
       NAME: SASKMI.SAS                                         				  	 					
      TITLE: A SAS Macro to Perform Kaplan-Meier Multiple Imputation for Survival Analyses with Competing Events  				
     AUTHOR: Jialin Han, Stanford University                                          	            
   Software: SAS Enterprise Guide 7.1									    						    
  CREATE ON: Dec 2020                                      					 	    
  UPDATE ON: Apr 2021	                                   					 	    
  UPDATE ON: June 2024	- Fix Macro reach character limit problem										    
 DESCRIPTION: This program perform KMI on survival analysis with competing events 			    
													    
    Copyright (C) <2021>  <Jialin Han>				                  		    
													    
    This program is free software: you can redistribute it and/or modify				    
    it under the terms of the GNU General Public License as published by			     	    
    the Free Software Foundation, either version 3 of the License, or					    
    (at your option) any later version.								    
													    
    This program is distributed in the hope that it will be useful,				            
    but WITHOUT ANY WARRANTY; without even the implied warranty of					    
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the					    
    GNU General Public License for more details.							    
													    
    You should have received a copy of the GNU General Public License					    
    along with this program.  If not, see <http://www.gnu.org/licenses/>.				    
***********************************************************************************************************/

%macro SASKMI(	data =  , time = , event = , 
			class = , 
			adjvar =  , 
			nimp = , eventcode = 1, censcode = 0, seed = 123, out = out);
/* Start timer */
%let _timer_start = %sysfunc(datetime());
%if &event =  or &time =  %then %do;
%put ERROR: Argument 'time' or 'event' is missing with no default;
%abort cancel;
%end;
/*censor time*/
proc sql noprint;
select count(distinct &time) into: n_cens from &data where &event = &censcode;
quit;
%if &n_cens <= 1 %then %do;
%put ERROR: Argument kmi can not make imputation based on less than one censoring time;
%abort cancel;
%end;
proc sql noprint;
select distinct &time into: cens_time SEPARATED BY ' ' from &data where &event = &censcode
order by &time;
create table cens as
select distinct &time from &data where &event = &censcode order by &time;
quit;
/*Keep an overall _id to main dataset*/
data x_x;
set &data;
_id = _n_;
run;
/*Data need to impute*/
data xx;
set x_x(where= (&event not in (&eventcode., &censcode.)));
drop &event;
run;
data xxx;
set x_x(where= (&event in (&eventcode., &censcode.)));
run;
data x;
set x_x;
if &event ^= &censcode then &event = &eventcode; /*Convert three indicator into two*/
run; 
/*
Cox model, treat censored as outcome
In order to get same estimate between R and SAS, 
need to specific efron in ties and FH in baseline
Model time * censor code
*/
proc phreg data=x noprint;
class &class/ PARAM=reference;
model &time*&event(&eventcode) = &adjvar/ TIES=EFRON;
baseline out=g(keep=surv _id) covariates=xx timelist=&cens_time survival=surv/method=FH;
run;
data g1;
set xx;
surv = 1;
keep surv _id;
run;
data gg;
set g g1;
run;
proc sort data=gg;
by _id descending surv;
run;
/*Time Need to Impute*/
proc sql noprint;
select count(*) into: n_it from xx;
/*Update June 17, 2024:
I decidede to not output &time as a macro as it will truncate at ~60000 charater,
Instead let's use it as a dataset*/
proc sql;
create table itimes as select &time from xx order by  _id;
/*select &time into: itimes SEPARATED BY ' ' from xx order by _id;*/
quit;
/*Find interval function*/
%let a = 0;
%let b = 0;
proc iml;
use x_x;
read all var{"&time", "&event"};
close x_x;
/*When the largest time is an event, add 1 day to the largest time into cenoring distribution*/
i = loc(&event = &censcode);
if max(&time) ^= max(&time[i]) then do;
ct = 0||{&cens_time}||(max(&time)+1);
call symputx("a", 1); 
call symputx("b", max(&time));
end;
else do;
ct = 0||{&cens_time};
end;
/*it = {&itimes};*/
use itimes;
read all var _all_ into it;
close itimes;
it=t(it) ;*make sure it = vecotr;
tmp = repeat(0, ncol(it));
do i = 1 to (ncol(ct)-1);
idx = loc(it > ct[i] & it <= ct[i+1]);
if ncol(idx) > 0 then tmp[idx] = (i+1);
end;
create tmp var{tmp};
append;
close;
run;
proc sql noprint;
select tmp into: tmp1 - : %sysfunc(compress(tmp&n_it)) from tmp;
quit;
/*%put "check time = &b.";*/
/*Start IML*/
proc iml;
call randseed(&seed);
use tmp;
read all var {"tmp"};
close tmp;
x = repeat({0 0 0},(&n_it*&nimp.));
use cens;
read all var{"&time"};
close cens;
use gg;
read all var {"surv","_id"};
close gg;
do i = 1 to &n_it;
	t = tmp[i];
	s1 = surv[((i-1)*(nrow(&time)+1)+1):(i*(nrow(&time)+1))];
	if t >= nrow(s1) then s2 = s1;
	else s2 = s1[1:t,]//repeat(s1[t],(nrow(s1)-t)); 
	spr = s1/s2;
	wp = -dif(spr,1);
	wp = wp[2:nrow(wp),];
	if &a. = 1 then do;
			wp = wp//spr[nrow(s1)];
/*If largest time is event then add one day*/
			time = {&cens_time}||(%sysevalf(&b+1));
	end;
	else time = &time;
	x[((i-1)*&nimp.+1):(i*&nimp.),2] = colvec(sample(time, &nimp, "Replace", wp));
	x[((i-1)*&nimp.+1):(i*&nimp.),1] = _id[i*(nrow(&time)+1)];
	x[((i-1)*&nimp.+1):(i*&nimp.),3] = T(1:&nimp.);
/*Random select from cens.time based on weights*/
/*proc surveyselect data=gg3 out=tt n=1 seed=&seed method=pps_wr noprint;*/
/*size wp;*/
/*run;*/
end;
create tt from x[colname = {"_id" "newtime" "nimp"}];
append from x;
close tt;
/*create t var{t};*/
/*append;*/
/*close;*/
/*create s1 var{s1};*/
/*append;*/
/*close;*/
/*create s2 var{s2};*/
/*append;*/
/*close;*/
/*create wp var{wp};*/
/*append;*/
/*close;*/
run;
/*End of IML*/
data xxx_n;
set xxx;
do i = 1 to &nimp;
output;
end;
newevent = &event;
newtime = &time;
rename i = nimp ;
run;
proc sort data=tt;
by _id;
run;
/*Fixed March 17, merge problem*/
data xx_n;
merge xx(in=a) tt(in=b);
by _id;
/*drop &time ;*/
newevent = 0;
&event = 2;
if a and b then output;
run;
data x_o;
set xx_n xxx_n;
run;
proc sort data=x_o tagsort;
by nimp _id;
run;
data &out;
set x_o;
if &event ^= 2 then do;
newevent = &event;
newtime = &time;
end;
rename nimp = _imputation_;
drop _id;
run;
/*In the end: design the output*/
title 'Kaplan-Meier Multiple Imputation for Competing Risks';
proc iml;
o = {"Data Set" "&data.",
	"Method" "KMI",
	"Number of Imputations" "&nimp.",
	"Time variable" "&time",
	"Event variable" "&event",
	"Censoring code" "&censcode",
	"Event code" "&eventcode",
	"Model Info" "&event*&time = &adjvar",
	"Seed for random number generator" "&seed"
};
use tt;
read all var {"newtime" "nimp"};
close tt;
use xxx;
read all var {"&time" "&event"};
close xxx;
o1 = repeat(0,(2+&nimp.),5);
idx1 = &time[loc(&event = &eventcode.)];
o1[1,1] = nrow(idx1);
o1[1,2] = mean(idx1);
o1[1,3] = min(idx1);
o1[1,4] = median(idx1);
o1[1,5] = max(idx1);
idx = &time[loc(&event = &censcode.)];
o1[2,1] = nrow(idx);
o1[2,2] = mean(idx);
o1[2,3] = min(idx);
o1[2,4] = median(idx);
o1[2,5] = max(idx);
Time = colvec(("Event time")||("Censoring time")||("Imputed time 1":"Imputed time &nimp."));
do i = 1 to &nimp;
idx = newtime[loc(nimp = i)];
o1[(i+2),1] = nrow(idx);
o1[(i+2),2] = mean(idx);
o1[(i+2),3] =  min(idx);
o1[(i+2),4] = median(idx);
o1[(i+2),5] = max(idx);
end;
col = {"Nobs" "Mean" "Min" "Median" "Max"};
create o1 from o1[colname = col  rowname = Time];
append from o1[ rowname = Time];
close;
print o[label="Model Information"];
quit;
title "Summary of Imputed Time";
proc print data=o1 noobs ;
run;

/*ods exclude all;*/
/*proc datasets lib= work;*/
/*delete o: x: g: tt tmp;*/
/*run;*/
/*ods exclude none;*/
/*ods _all_ close;*/
/* Stop timer */
data _null_;
file log;
  dur = datetime() - &_timer_start;
  put 30*'-' / ' TOTAL DURATION:' dur time13.2 / 30*'-';
run;
%mend;

