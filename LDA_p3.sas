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


proc print data=hearing_miss(obs=20);
run;


/*Psi model*/ 
proc genmod data=hearing_miss DESCENDING;
	model dropout = ylag age side2 / dist=binomial;
run;



/*Q4*/ 
/*=============================================================================*/

proc import datafile = "&path\hearing_timeround_wide.csv" out=hearing_wide
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


