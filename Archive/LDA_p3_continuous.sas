%let path F:\Dropbox\Long\1. Master of Statistics\0. Program\3. Semester 3\6. Longditudial Data Analysis\5. Data and code\LDA3;

libname lda "&path";


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
	time2 = timeround**2;
run;


proc print data=hearing_long(obs=20);
run;


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


proc print data=hearing_miss(obs=200);
run;


/*Psi model*/ 
proc genmod data=hearing_miss DESCENDING;
	model dropout = ylag age side2 / dist=binomial;
run;



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

proc mi data=hearing_wide seed=2023 simple nimpute=10 round=0.1 out=hearing_wide_mono;
	mcmc impute = monotone;
	var side2 age t0 t1 t2 t3 t4 t5 t6 t7 t8 t9 t10 t11 t12 t13 t14 t15 t16 t17 
	t18 t19 t20 t21 t22;
run;


proc print data=hearing_wide;
run;

/*Imputation*/

proc mi data=hearing_wide seed=2023 out=lda.hearing_wide_im simple nimpute=10 round=0.1;
	var age t0 t1 t2 t3 t4 t5 t6 t7 t8 t9 t10 t11 t12 t13 t14 t15 t16 t17 
	t18 t19 t20 t21 t22;
	by side;
run;


proc print data=lda.hearing_wide_im (obs=100);
run;


/*Read long form from R*/
proc import datafile = "&path\df_long_im.csv" out=hearing_long_im
		dbms=csv replace;
		guessingrows=max;
	getnames=yes;
run;




/*Q4*/ 
/*=============================================================================*/
/*Model with random intercept and slopes*/
proc print data=hearing_long(obs=100);
run;

proc mixed data = hearing_long method=reml;
	class side timeclss id;
	model y = age side2 timeround age*timeround side2*timeround age*time2 side2*time2 /s noint;
	random intercept timeround /type = cs subject=id g gcorr v vcorr;
	repeated timeclss/type = cs subject = id r rcorr;
run;

/*remove side2*time2*/
proc mixed data = hearing_long method=reml;
	class side timeclss id;
	model y = age side2 timeround age*timeround side2*timeround age*time2 /s noint;
	random intercept timeround /type = cs subject=id g gcorr v vcorr;
	repeated timeclss/type = cs subject = id r rcorr;
run;

/*remove side2*timeround*/
proc mixed data = hearing_long method=reml;
	class side timeclss id;
	model y = age side2 timeround age*timeround age*time2 /s noint;
	random intercept timeround /type = cs subject=id g gcorr v vcorr;
	repeated timeclss/type = cs subject = id r rcorr;
run;



proc mixed data = hearing_long method=reml;
	class side timeclss id;
	model y = age side2 timeround age*timeround age*time2 /s noint;
	random intercept/type = cs subject=id;
	random intercept/type = cs subject = side2(id);
	repeated timeclss/type = cs subject = side2(id) r rcorr;
run;





	












