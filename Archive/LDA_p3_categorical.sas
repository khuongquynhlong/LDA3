%let path F:\Dropbox\Long\1. Master of Statistics\0. Program\3. Semester 3\6. Longditudial Data Analysis\5. Data and code\LDA3;

libname lda "&path";


/*Multiple Imputation*/ 
/*=============================================================================*/

/*Read long form from R*/
proc import datafile = "&path\df_long_im.csv" out=hearing_long_im
		dbms=csv replace;
		guessingrows=max;
	getnames=yes;
run;


data hearing_long_im;
	set hearing_long_im;
	sideclss = side2;
	timeclss = time;
	time2 = time**2;
	if y <=15 then ycat = 1;
	else if y <= 25 then ycat = 2;
	else ycat = 3;
run;

proc sort data=hearing_long_im;
	by _imputation_ id side time;
run;


proc print data=hearing_long_im(obs=100);
run;

/*Q5*/ 
/*=========================================================*/
/*Combine additional intercept for gmpinfo*/
/*
data gmpinfoint; 
	length _Imputation_ 8 Parameter $12 Effect $20;
	input _Imputation_ Parameter Effect;
	datalines;
	1 Intercept1 Intercept1
	1 Intercept2 Intercept2
	2 Intercept1 Intercept1
	2 Intercept2 Intercept2
	3 Intercept1 Intercept1
	3 Intercept2 Intercept2
	4 Intercept1 Intercept1
	4 Intercept2 Intercept2
	5 Intercept1 Intercept1
	5 Intercept2 Intercept2
	6 Intercept1 Intercept1
	6 Intercept2 Intercept2
	7 Intercept1 Intercept1
	7 Intercept2 Intercept2
	8 Intercept1 Intercept1
	8 Intercept2 Intercept2
	9 Intercept1 Intercept1
	9 Intercept2 Intercept2
	10 Intercept1 Intercept1
	10 Intercept2 Intercept2
;
run;

proc print data=gmpinfoint(obs=50); run;
*/

/*Q5 GEE with miltiple imputation*/ 
/*=========================================================*/
/*Time squared*/ 
proc gee data=hearing_long_im;
	class id timeclss side ycat;
	by _Imputation_;
	model ycat = time age side2 time*age time2*age / 
		dist =  multinomial link=cumlogit type3;
	repeated subject=id/ type=ind covb corrw 
		within=side*timeclss modelse ecovb;
	ods output GEEEmpPEst=lda.gmparms parminfo=lda.gmpinfo 
				modelinfo=modelinfo GEERCov=lda.gmcovb;
run;

proc import datafile = "&path\gmparms.csv" out=gmparms
		dbms=csv replace; guessingrows=max; getnames=yes;
run;
proc import datafile = "&path\gmcovb.csv" out=gmcovb
		dbms=csv replace; guessingrows=max; getnames=yes;
run;
proc import datafile = "&path\gmpinfo.csv" out=gmpinfo
		dbms=csv replace; guessingrows=max; getnames=yes;
run;

proc mianalyze parms=gmparms parminfo=gmpinfo covb=gmcovb;
	modeleffects prm6 prm7 time age side2 time*age time2*age;
run;




/*Time linear*/ 

proc gee data=hearing_long_im;
	class id timeclss side ycat;
	by _Imputation_;
	model ycat = time age side2 time*age time*side2 / 
		dist =  multinomial link=cumlogit type3;
	repeated subject=id/ type=ind corrw 
		within=side*timeclss modelse ecovb;
	ods output GEEEmpPEst=lda.gmparms parminfo=lda.gmpinfo 
				modelinfo=modelinfo GEERCov=lda.gmcovb;
run;


proc import datafile = "&path\gmparms.csv" out=gmparms
		dbms=csv replace; guessingrows=max; getnames=yes;
run;
proc import datafile = "&path\gmcovb.csv" out=gmcovb
		dbms=csv replace; guessingrows=max; getnames=yes;
run;
proc import datafile = "&path\gmpinfo.csv" out=gmpinfo
		dbms=csv replace; guessingrows=max; getnames=yes;
run;

proc mianalyze parms=gmparms parminfo=gmpinfo covb=gmcovb;
	modeleffects prm6 prm7 time age side2 time*age time*side2;
run;



proc gee data=hearing_long_im;
	class id timeclss side ycat;
	by _Imputation_;
	model ycat = time age side2/ 
		dist =  multinomial link=cumlogit type3;
	repeated subject=id/ type=ind corrw 
		within=side*timeclss modelse ecovb;
	ods output GEEEmpPEst=lda.gmparms parminfo=lda.gmpinfo 
				modelinfo=modelinfo GEERCov=lda.gmcovb;
run;

proc import datafile = "&path\gmparms.csv" out=gmparms
		dbms=csv replace; guessingrows=max; getnames=yes;
run;
proc import datafile = "&path\gmcovb.csv" out=gmcovb
		dbms=csv replace; guessingrows=max; getnames=yes;
run;
proc import datafile = "&path\gmpinfo.csv" out=gmpinfo
		dbms=csv replace; guessingrows=max; getnames=yes;
run;

proc mianalyze parms=gmparms parminfo=gmpinfo covb=gmcovb;
	modeleffects prm6 prm7 time age side2;
run;



/*Q5 WGEE*/ 
/*=========================================================*/

proc import datafile = "&path\hearing_timeround.csv" out=hearing_long
		dbms=csv replace;
		guessingrows=max;
	getnames=yes;
run;

proc sort data=hearing_long;
	by id side timeround;
run;

data hearing_long;
	set hearing_long;
	if side = "Left" then side2 = 0; 
	else side2 = 1;
	sideclss = side2;
	timeclss = timeround;
	time = timeround;
	time2 = timeround**2;
run;

proc print data=hearing_long(obs=20);
run;


/*Compute weights with marcro*/ 

/* Dropout marcro*/;


%macro dropout(data=,id=,time=,response=,out=);
%if %bquote(&data)= %then %let data=&syslast;
proc freq data=&data noprint;
tables &id /out=freqid;
tables &time / out=freqtime;
run;
proc iml;
reset noprint;
use freqid;
read all var {&id};
nsub = nrow(&id);
use freqtime;
read all var {&time};
ntime = nrow(&time);
time = &time;
use &data; 
read all var {&id &time &response};
n = nrow(&response);
dropout = j(n,1,0);
ind = 1;
do while (ind <= nsub);
  j=1;
  if (&response[(ind-1)*ntime+j]=.) then print "First Measurement is Missing";
  if (&response[(ind-1)*ntime+j]^=.) then
    do;
      j = ntime;
      do until (j=1);
        if (&response[(ind-1)*ntime+j]=.) then 
          do;
            dropout[(ind-1)*ntime+j]=1;
            j = j-1;
          end;
          else j = 1;
      end; 
    end;
  ind = ind+1;
end;
prev = j(n,1,1);
prev[2:n] = &response[1:n-1];
i=1;
do while (i<=n);
  if &time[i]=time[1] then prev[i]=.;
  i = i+1;
end;
create help var {&id &time &response dropout prev};
append;
quit;
data &out;
merge &data help;
run;
%mend;

%dropout(data=hearing_long,id=id,time=timeround,response=y,out=hearing_miss);


/*Psi model*/ 
proc genmod data=hearing_miss DESCENDING;
	model dropout = prev age side2 time / dist=binomial pred;
	ods output obstats=pred;
run;

proc print data=pred (obs=100);
run;

data pred;
	set pred;
	keep observation pred;
run;

data hearing_miss;
	merge pred hearing_miss;
run;

proc print data=hearing_miss(obs=100); run;



%macro dropwgt(data=,id=,time=,pred=,dropout=,out=);
%if %bquote(&data)= %then %let data=&syslast;
proc freq data=&data noprint;
tables &id /out=freqid;
tables &time / out=freqtime;
run;
proc iml;
reset noprint;
use freqid;
read all var {&id};
nsub = nrow(&id);
use freqtime;
read all var {&time};
ntime = nrow(&time);
time = &time;
use &data; 
read all var {&id &time &pred &dropout};
n = nrow(&pred);
wi = j(n,1,1);
ind = 1;
do while (ind <= nsub);
    wihlp=1;
    stay=1;
    /* first measurement */
    if (&dropout[(ind-1)*ntime+2]=1)
      then do;
          wihlp = pred[(ind-1)*ntime+2];
        stay=0;
      end;
    else if (&dropout[(ind-1)*ntime+2]=0)
      then wihlp = 1-pred[(ind-1)*ntime+2];
    /* second to penultimate measurement */
    j=2;
    do while ((j <= ntime-1) & stay);
      if (&dropout[(ind-1)*ntime+j+1]=1)
        then do;
          wihlp = wihlp*pred[(ind-1)*ntime+j+1];
          stay=0;
        end;
      else if (&dropout[(ind-1)*ntime+j+1]=0)
        then wihlp = wihlp*(1-pred[(ind-1)*ntime+j+1]);
      j = j+1;
    end;
    j=1;
    do while (j <= ntime);
      wi[(ind-1)*ntime+j] = wihlp;
      j = j+1;
    end;
    ind = ind+1;
end;
create help var {&id &time &pred &dropout wi};
append;
quit;
data &out;
merge &data help;
data &out;
set &out;
wi = 1/wi;
run;
%mend;


%dropwgt(data=hearing_miss,id=id,time=time,pred=pred,dropout=DROPOUT,out=hearing_weight)

proc print data=hearing_weight(obs=100); run;


proc print data=hearing_weight;
var id ylag age side2 time dropout prev pred wi;
run;









