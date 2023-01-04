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
	time = timeround;
	time2 = timeround**2;
run;


proc print data=hearing_long(obs=20);
run;


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

/*Q4*/ 
/*=============================================================================*/
/*GLMM*/
/*GLMM Direct likelihood*/

proc glimmix data=hearing_long method=quad(qpoints=10 qcheck);
	class id timeclss ycat sideclss;
	model ycat = time age side2 time*age time2*age/ 
	dist =  multinomial link=cumlogit solution;
	random intercept /subject=id;
	random intercept time /subject=sideclss(id) type=chol;
run;



proc nlmixed data=hearing_long qpoints=5 maxiter=100 technique=newrap;
	parms int1=9.2100 int2=12.1362 beta1=0.2445 beta2=-0.1037 
			beta3=-0.1756 beta4=-0.00557 beta5=-0.00006
			sigmab1 = 2.6149 sigmab2 = 1.1680 sigmab3 = 0.1303 rho = 0.002596;
	eta = beta1*time + beta2*age + beta3*side2 + 
		  beta4*time*age + beta5*time2*age + b1 + b2 + b3*time;
	if ycat=1 then z = 1/(1+exp(-(int1+eta)));
	else if ycat=2 then z = 1/(1+exp(-(int2+eta))) - 1/(1+exp(-(int1+eta)));
	else z = 1 - 1/(1+exp(-(int2+eta)));
	if z > 1e-8 then ll = log(z);
	else ll = -1e100;
	model ycat ~ general(ll);
	random b1 ~ normal(0,sigmab1**2) subject=id; 
	random b2 b3 ~ normal([0, 0], 
					[sigmab2**2, rho*sigmab2*sigmab3, sigmab3**2]) 
	subject = side2(id);
	estimate 'sigmab1^2' sigmab1**2;
	estimate 'sigmab2^2' sigmab2**2;
	estimate 'sigmab3^2' sigmab3**2;
	estimate 'Cov23' rho*sigmab2*sigmab3;
run;



















