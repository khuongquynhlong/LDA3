%let path F:\Dropbox\Long\1. Master of Statistics\0. Program\3. Semester 3\6. Longditudial Data Analysis\5. Data and code\LDA3;
%let path C:\Users\khuon\Dropbox\Long\1. Master of Statistics\0. Program\3. Semester 3\6. Longditudial Data Analysis\5. Data and code\LDA3;

libname lda "&path";


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



/*Multiple Imputation*/ 
/*=============================================================================*/

proc import datafile = "&path\hearing_timeround_wide.csv" out=hearing_wide
		dbms=csv replace;
		guessingrows=max;
	getnames=yes;
run;


proc sort data=hearing_wide;
	by side;
run;

data hearing_wide;
	set hearing_wide;
	if side = "Left" then side2 = 0; 
	else side2 = 1;
run;

proc mi data=hearing_wide seed=2023 simple nimpute=10 round=0.1 out=lda.hearing_wide_mono;
	mcmc impute = monotone;
	var side2 age t0 t2 t4 t6 t8 t10 t12 t14 t16 t18 t20 t22;
run;


proc print data=lda.hearing_wide_mono(obs=100);
run;


proc print data=hearing_wide;
run;

/*Imputation*/

proc mi data=hearing_wide seed=2023 out=lda.hearing_wide_im simple nimpute=10 round=0.1;
	var age t0 t2 t4 t6 t8 t10 t12 t14 t16 t18 t20 t22;
	by side;
run;


proc print data=lda.hearing_wide_im (obs=100);
run;



/*Q3: Factors influence missingness*/ 
/*=============================================================================*/

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
	time = timeround;
	sideclss = side2;
	timeclss = timeround;
	time2 = timeround**2;
	if y <=15 then ycat = 1;
	else if y <= 25 then ycat = 2;
	else ycat = 3;
run;


proc print data=hearing_long(obs=20);
run;

%dropout(data=hearing_long,id=id,time=timeround,response=y,out=hearing_miss);


proc print data=hearing_miss(obs=200);
run;


/*Psi model*/ 
proc genmod data=hearing_miss DESCENDING;
	model dropout = ylag age side2 / dist=binomial;
run;



/*Q4: Compare a direct likelihood analysis with multiple imputation*/ 
/*==================================================================*/


/*------------------Q4a:Direct Likelihood */ 
proc print data=hearing_long(obs=20);
run;

proc mixed data=hearing_long method=REML;
	class id timeclss sideclss;
	model y = age side2 timeround timeround*side2 timeround*age / covb s;
	repeated timeclss / type = sp(exp)(time) local subject=id r rcorr;
	random intercept / subject=id g s;
	random intercept timeround / type=un subject=sideclss(id) g s vcorr;
run;


/*------------------Q4b: LMM with MI*/ 
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
	by _Imputation_ id side time;
run;

proc mixed data=hearing_long_im method=REML;
	class id timeclss sideclss;
	by _Imputation_;
	model y = age side2 time time*side2 time*age/ covb s;
	repeated timeclss / type = sp(exp)(time) local subject=id;
	random intercept /  subject=id;
	random intercept time / type=un subject=sideclss(id);
	ods output solutionF = nlparms CovB=nlcovb;	
run;


proc mianalyze parms=nlparms covb(effectvar=rowcol)=nlcovb wcov bcov tcov;
	modeleffects intercept age side2 time time*side2 time*age;
run;


/*Q5: WGEE vs. GEE with miltiple imputation*/ 
/*==================================================================*/

/*------------------Q5a:WGEE */ 
/*Created weights*/
proc import datafile = "&path\df_long_im_mono.csv" out=hearing_long_im_mono
		dbms=csv replace;
		guessingrows=max;
	getnames=yes;
run;

data hearing_long_im_mono;
	set hearing_long_im_mono;
	sideclss = side2;
	timeclss = time;
	time2 = time**2;
	if y <=15 then ycat = 1;
	else if y <= 25 then ycat = 2;
	else ycat = 3;
run;


proc sort data=hearing_long_im_mono;
	by _Imputation_ id side time;
run;

proc print data=hearing_long_im_mono(obs=100); run;

data hearing_long_im_mono1;
	set hearing_long_im_mono;
	where _Imputation_ = 1;
run;
data hearing_long_im_mono2;
	set hearing_long_im_mono;
	where _Imputation_ = 2;
run;
data hearing_long_im_mono3;
	set hearing_long_im_mono;
	where _Imputation_ = 3;
run;
data hearing_long_im_mono4;
	set hearing_long_im_mono;
	where _Imputation_ = 4;
run;
data hearing_long_im_mono5;
	set hearing_long_im_mono;
	where _Imputation_ = 5;
run;
data hearing_long_im_mono6;
	set hearing_long_im_mono;
	where _Imputation_ = 6;
run;
data hearing_long_im_mono7;
	set hearing_long_im_mono;
	where _Imputation_ = 7;
run;
data hearing_long_im_mono8;
	set hearing_long_im_mono;
	where _Imputation_ = 8;
run;
data hearing_long_im_mono9;
	set hearing_long_im_mono;
	where _Imputation_ = 9;
run;
data hearing_long_im_mono10;
	set hearing_long_im_mono;
	where _Imputation_ = 10;
run;

%dropout(data=hearing_long_im_mono1,id=id,time=time,response=y,out=hearing_miss_mono1);
%dropout(data=hearing_long_im_mono2,id=id,time=time,response=y,out=hearing_miss_mono2);
%dropout(data=hearing_long_im_mono3,id=id,time=time,response=y,out=hearing_miss_mono3);
%dropout(data=hearing_long_im_mono4,id=id,time=time,response=y,out=hearing_miss_mono4);
%dropout(data=hearing_long_im_mono5,id=id,time=time,response=y,out=hearing_miss_mono5);
%dropout(data=hearing_long_im_mono6,id=id,time=time,response=y,out=hearing_miss_mono6);
%dropout(data=hearing_long_im_mono7,id=id,time=time,response=y,out=hearing_miss_mono7);
%dropout(data=hearing_long_im_mono8,id=id,time=time,response=y,out=hearing_miss_mono8);
%dropout(data=hearing_long_im_mono9,id=id,time=time,response=y,out=hearing_miss_mono9);
%dropout(data=hearing_long_im_mono10,id=id,time=time,response=y,out=hearing_miss_mono10);


data hearing_miss_mono;
	set hearing_miss_mono1 hearing_miss_mono2 hearing_miss_mono3 hearing_miss_mono4 
		hearing_miss_mono5 hearing_miss_mono6 hearing_miss_mono7 hearing_miss_mono8
		hearing_miss_mono9 hearing_miss_mono10;
run;

proc print data=hearing_miss_mono(obs=200); run;

/*Psi model*/ 
proc genmod data=hearing_miss_mono DESCENDING;
	model dropout = age side2 time/ dist=binomial pred;
	ods output obstats=pred;
run;

proc print data=pred (obs=100);
run;

data pred;
	set pred;
	keep observation pred;
run;

data hearing_miss_mono;
	merge pred hearing_miss_mono;
run;



data hearing_miss_mono1;
	set hearing_miss_mono;
	where _Imputation_ = 1;
run;
data hearing_miss_mono2;
	set hearing_miss_mono;
	where _Imputation_ = 2;
run;
data hearing_miss_mono3;
	set hearing_miss_mono;
	where _Imputation_ = 3;
run;
data hearing_miss_mono4;
	set hearing_miss_mono;
	where _Imputation_ = 4;
run;
data hearing_miss_mono5;
	set hearing_miss_mono;
	where _Imputation_ = 5;
run;
data hearing_miss_mono6;
	set hearing_miss_mono;
	where _Imputation_ = 6;
run;
data hearing_miss_mono7;
	set hearing_miss_mono;
	where _Imputation_ = 7;
run;
data hearing_miss_mono8;
	set hearing_miss_mono;
	where _Imputation_ = 8;
run;
data hearing_miss_mono9;
	set hearing_miss_mono;
	where _Imputation_ = 9;
run;
data hearing_miss_mono10;
	set hearing_miss_mono;
	where _Imputation_ = 10;
run;

%dropwgt(data=hearing_miss_mono1,id=id,time=time,pred=pred,dropout=DROPOUT,out=hearing_weight1)
%dropwgt(data=hearing_miss_mono2,id=id,time=time,pred=pred,dropout=DROPOUT,out=hearing_weight2)
%dropwgt(data=hearing_miss_mono3,id=id,time=time,pred=pred,dropout=DROPOUT,out=hearing_weight3)
%dropwgt(data=hearing_miss_mono4,id=id,time=time,pred=pred,dropout=DROPOUT,out=hearing_weight4)
%dropwgt(data=hearing_miss_mono5,id=id,time=time,pred=pred,dropout=DROPOUT,out=hearing_weight5)
%dropwgt(data=hearing_miss_mono6,id=id,time=time,pred=pred,dropout=DROPOUT,out=hearing_weight6)
%dropwgt(data=hearing_miss_mono7,id=id,time=time,pred=pred,dropout=DROPOUT,out=hearing_weight7)
%dropwgt(data=hearing_miss_mono8,id=id,time=time,pred=pred,dropout=DROPOUT,out=hearing_weight8)
%dropwgt(data=hearing_miss_mono9,id=id,time=time,pred=pred,dropout=DROPOUT,out=hearing_weight9)
%dropwgt(data=hearing_miss_mono10,id=id,time=time,pred=pred,dropout=DROPOUT,out=hearing_weight10)



data hearing_weight;
	set hearing_weight1 hearing_weight2 hearing_weight3 hearing_weight4 
		hearing_weight5 hearing_weight6 hearing_weight7 hearing_weight8
		hearing_weight9 hearing_weight10;
run;


proc print data=hearing_weight(obs=100); run;



/*WGEE*/ 
proc genmod data=hearing_weight;
	by _imputation_;
	scwgt wi;
	class id timeclss side ycat;
	model ycat = time age side2 time*side2 time*age / 
		dist =  multinomial link=cumlogit type3;
	repeated subject=id/ type=ind covb corrw 
		within=side*timeclss modelse ecovb;
	ods output GEEEmpPEst=lda.gmparms parminfo=lda.gmpinfo 
				modelinfo=modelinfo GEERCov=lda.gmcovb;
run;


proc import datafile = "&path\gmparms2.csv" out=gmparms
		dbms=csv replace; guessingrows=max; getnames=yes;
run;
proc import datafile = "&path\gmcovb2.csv" out=gmcovb
		dbms=csv replace; guessingrows=max; getnames=yes;
run;
proc import datafile = "&path\gmpinfo2.csv" out=gmpinfo
		dbms=csv replace; guessingrows=max; getnames=yes;
run;

proc mianalyze parms=gmparms parminfo=gmpinfo covb=gmcovb wcov bcov tcov;
	modeleffects prm6 prm7 time age side2 time*side2 time*age;
run;


proc mianalyze parms=gmparms parminfo=gmpinfo covb=gmcovb wcov bcov tcov;
	modeleffects time age side2 time*side2 time*age;
run;



/*------------------Q5b:GEE with MI*/ 

proc print data=hearing_long_im(obs=100); run;

/*Time*/ 
proc gee data=hearing_long_im;
	class id timeclss side ycat;
	by _Imputation_;
	model ycat = time age side2 time*side2 time*age / 
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

proc mianalyze parms=gmparms parminfo=gmpinfo covb=gmcovb wcov bcov tcov;
	modeleffects prm6 prm7 time age side2 time*side2 time*age;
run;


/*Q6: Sensitivity analysis*/ 
/*==================================================================*/
proc print data=hearing_wide(obs=100);
run;

/*Convert to Monotone missingness*/
proc mi data=hearing_wide seed=2023 simple nimpute=10 round=0.1 out=sen_wide_mono;
	mcmc impute = monotone;
	var side2 age t0 t2 t4 t6 t8 t10 t12 t14 t16 t18 t20 t22;
run;


proc print data=sen_wide_mono(obs=100);
run;


/*------------------ CCMV*/
proc mi data=sen_wide_mono seed=2023 simple out=sen_wide_ccnv nimpute=1;
	var side2 age t0 t2 t4 t6 t8 t10 t12 t14 t16 t18 t20 t22;
	monotone reg;
	mnar model (t0 t2 t4 t6 t8 t10 t12 t14 t16 t18 t20 t22 / modelobs=ccmv);
run;

proc print data=sen_wide_ccnv(obs=2000); run;

proc sort data=sen_wide_ccnv;
	by _imputation_ id side;
run;

/*Long transform*/
proc transpose data=sen_wide_ccnv out=sen_long_ccnv;
	by _imputation_ id side age side2;
run;

data sen_long_ccnv;
	set sen_long_ccnv; 
	time1=compress(_NAME_,'t');
	time=input(time1, 8.);
	y = col1;
	drop time1;
	drop col1;
	drop _NAME_;
	sideclss = side2;
	timeclss = time;
	time2 = time**2;
	if y <=15 then ycat = 1;
	else if y <= 25 then ycat = 2;
	else ycat = 3;
run;


proc mixed data=sen_long_ccnv method=REML;
	class id timeclss sideclss;
	by _Imputation_;
	model y = age side2 time time*side2 time*age/ covb s;
	repeated timeclss / type = sp(exp)(time) local subject=id;
	random intercept /  subject=id;
	random intercept time / type=un subject=sideclss(id);
	ods output solutionF = nlparms CovB=nlcovb;	
run;

proc mianalyze parms=nlparms covb(effectvar=rowcol)=nlcovb;
	modeleffects intercept age side2 time time*side2 time*age;
run;


/*------------------ NCMV*/
proc mi data=sen_wide_mono seed=2023 simple out=sen_wide_ncnv nimpute=1;
	var side2 age t0 t2 t4 t6 t8 t10 t12 t14 t16 t18 t20 t22;
	monotone reg;
	mnar model (t0 t2 t4 t6 t8 t10 t12 t14 t16 t18 t20 t22 / modelobs=ncmv);
run;

proc print data=sen_wide_ncnv(obs=2000); run;

proc sort data=sen_wide_ncnv;
	by _imputation_ id side;
run;

/*Long transform*/
proc transpose data=sen_wide_ncnv out=sen_long_ncnv;
	by _imputation_ id side age side2;
run;

data sen_long_ncnv;
	set sen_long_ncnv; 
	time1=compress(_NAME_,'t');
	time=input(time1, 8.);
	y = col1;
	drop time1;
	drop col1;
	drop _NAME_;
	sideclss = side2;
	timeclss = time;
	time2 = time**2;
	if y <=15 then ycat = 1;
	else if y <= 25 then ycat = 2;
	else ycat = 3;
run;


proc mixed data=sen_long_ncnv method=REML;
	class id timeclss sideclss;
	by _Imputation_;
	model y = age side2 time time*side2 time*age/ covb s;
	repeated timeclss / type = sp(exp)(time) local subject=id;
	random intercept /  subject=id;
	random intercept time / type=un subject=sideclss(id);
	ods output solutionF = nlparms CovB=nlcovb;	
run;

proc mianalyze parms=nlparms covb(effectvar=rowcol)=nlcovb;
	modeleffects intercept age side2 time time*side2 time*age;
run;
