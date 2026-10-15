ods listing; %*OPEN;
run;

/*import data*/
proc import datafile="C:\Ags.txt" /*specify where the data is*/
out=Ags
dbms=TAB 
REPLACE; 
getnames=yes;
run;

/*checking data*/
proc print data=Ags;
run;

proc univariate data=Ags NORMAL PLOT;
var A;
qqplot A/ normal (mu=est sigma=est color=red l=1);
run;

/*MIXED MODEL FITTING*/ 
/*You have to make several runs to choose the best covariance structure*/
ods html style=journal;
ods graphics on;
proc mixed data=Ags plot=all;
class  Treatment Time Plant;
model A = Treatment|Time;
repeated Time/ type=cs sub=Plant;
lsmestimate Treatment*Time 
[1,1 1] [-1,2 1], [1,1 1][-1,3 1], [1,2 1] [-1,3 1],
[1,1 6] [-1,2 6], [1,1 6][-1,3 6], [1,2 6] [-1,3 6]/Adjust=BON ELSM; /*comparisons between specific levels of treatment*time for specific hyphotesis testing*/
run;
quit;
ods graphics off;
ods html close;


ods html style=journal;
ods graphics on;
proc mixed data=Ags plot=all;
class  Treatment Time Plant;
model A = Treatment|Time;
repeated Time/ type=csh sub=Plant;
lsmestimate Treatment*Time
[1,1 1] [-1,2 1], [1,1 1][-1,3 1], [1,2 1] [-1,3 1],
[1,1 6] [-1,2 6], [1,1 6][-1,3 6], [1,2 6] [-1,3 6]/Adjust=BON ELSM;
run;
quit;
ods graphics off;
ods html close;

ods html style=journal;
ods graphics on;
proc mixed data=Ags plot=all;
class  Treatment Time Plant;
model A = Treatment|Time;
repeated Time/ type=un sub=Plant;
lsmestimate Treatment*Time
[1,1 1] [-1,2 1], [1,1 1][-1,3 1], [1,2 1] [-1,3 1],
[1,1 6] [-1,2 6], [1,1 6][-1,3 6], [1,2 6] [-1,3 6]/Adjust=BON ELSM;
run;
quit;
ods graphics off;
ods html close;

ods html style=journal;
ods graphics on;
proc mixed data=Ags plot=all;
class  Treatment Time Plant;
model A = Treatment|Time;
repeated Time/ type=ar(1) sub=Plant;
lsmestimate Treatment*Time
[1,1 1] [-1,2 1], [1,1 1][-1,3 1], [1,2 1] [-1,3 1],
[1,1 6] [-1,2 6], [1,1 6][-1,3 6], [1,2 6] [-1,3 6]/Adjust=BON ELSM;
run;
quit;
ods graphics off;
ods html close;

/*BEST MODEL BY AIC*/
ods html style=journal;
ods graphics on;
proc mixed data=Ags plot=all;
class  Treatment Time Plant;
model A = Treatment|Time;
repeated Time/ type=toep sub=Plant;
lsmestimate Treatment*Time
[1,1 1] [-1,2 1], [1,1 1][-1,3 1], [1,2 1] [-1,3 1],
[1,1 8] [-1,2 8], [1,1 8][-1,3 8], [1,2 8] [-1,3 8]/Adjust=BON ELSM;
run;
quit;
ods graphics off;
ods html close;


