data gc;
input ikn:$ 3. servdate:yymmdd8. dayssupl quant_strength;
datalines;
1  20140710  4   30
1  20140929  10  30
1  20141224  5   40
1  20150210  5   100
1  20150222  15  60
1  20150226  25  80
1  20150603  30  100
1  20150811  30  90
1  20151221  20  50
1  20151230  30  50
1  20160201  30  50
1  20160310  30  50
1  20160410  30  70
1  20160520  50  80
1  20160715  60  200
2  20130801  7   300
2  20130802  7   300
2  20130809  7   100
2  20130812  7   100
2  20130816  7   100
2  20130930  30  50
2  20131115  14  60
2  20140106  20  10
2  20140309  40  20
2  20140608  40  100
2  20140718  40  200
3  20110615  14  460
3  20120106  20  80
4  20120630  50  500
4  20130101  10  120
5  20130204  10  200
5  20130212  10  200
6  20130519  10  300
6  20130601  10  200
7  20130819  10  300
7  20131001  30  210
7  20131101  10  200
;
run;

data gc;
set gc;
FORMAT servdate date9.;
daily = quant_strength/dayssupl;
run;

/************** Push dates back or forward first *************/
data gcpushback;
set gc;
by ikn;
retain cumshift;

prevdate=lag(servdate);
format prevdate date9.;
prevsupl=lag(dayssupl);
shift = (prevdate+prevsupl) - servdate;

if first.ikn then do;
	prevdate=.;
	prevsupl=.;
	shift=.;
	cumshift=0;
end;

servdate = servdate+cumshift;
if shift > 0 and shift <= 7 then do;
	servdate = servdate+shift;
	cumshift = sum(shift, cumshift);
end;

prevdate = lag(servdate);
run;


PROC SORT DATA=gcpushback(keep= ikn servdate dayssupl quant_strength daily) OUT=gcforward0;
  BY ikn DESCENDING servdate;
RUN ;


data gc_cleaned;
set gcforward0;
by ikn;

prevdate=lag(servdate);
format prevdate date9.;
prevsupl=lag(dayssupl);
shift = (servdate + dayssupl) - prevdate;

if first.ikn then do;
	prevdate=.;
	prevsupl=.;
	shift=.;
end;

if shift > 7 then do;
	dayssupl = prevdate - servdate;
	quant_strength = dayssupl*daily;
end;

run;

/* The gc_cleaned servdate is sorted in reverse order because of gc_forward, now sort it back */
PROC SORT DATA=gc_cleaned(keep= ikn servdate dayssupl quant_strength daily) OUT=gc_cleaned_order;
  BY ikn servdate;
RUN ;


data gc_cleaned_visual;
set gc_cleaned_order;
by ikn;
prevdate=lag(servdate);
format prevdate date9.;
prevsupl=lag(dayssupl);
shift = servdate - (prevdate+prevsupl);
enddate = servdate+dayssupl;
format enddate date9.;

if first.ikn then do;
	prevdate=.;
	prevsupl=.;
	shift=.;
end;
run;


/**************************************** Followdate determination *****************************************/
/************* If first >= 450, followdate is the second dispense date ****************/
/* First dispense */
data gc_copy1;
set gc_cleaned_order;
run;
proc sort nodupkey data=gc_copy1 out=gc_first;
by ikn;
run;

data gc_1st450;
set gc_first;
if quant_strength <=450 then delete;
run;

/* Second dispense */
data gc_copy2;
set gc_cleaned_order;
run;
proc sort nodupkey data=gc_copy2 dupout=gc_follow;
by ikn;
run;
proc sort nodupkey data=gc_follow out=gc_second;
by ikn;
run;

data gc_1st450_follow;
merge gc_1st450(in=a) gc_second(in=b);
by ikn;
if a;
keep ikn servdate;
rename servdate=followdate;
run;


data gc_final_1st450_0;
merge gc_cleaned_order(in=a) gc_1st450_follow(in=b);
by ikn;
if b;
run;

/* Get cumdose */
data gc_final_1st450_1;
set gc_final_1st450_0;
by ikn;
retain cumdose;
cumdose = sum(quant_strength, cumdose);
if first.ikn then cumdose = quant_strength; 
run;

/* To merge with gc_finalle450 */
data gc_final_1st450;
set gc_final_1st450_1;
by ikn;
cumdose = lag(cumdose);
if first.ikn then cumdose = .;
run;


/************** if first <=450, find the followdate in the rest dispense cumulatively ************/
data gc_cleaned_orderle450;
merge gc_1st450(in=a) gc_cleaned_order(in=b);
by ikn;
if a=0;
run;


data gc_cleaned_follow0;
set gc_cleaned_orderle450;
by ikn;
retain cumdose;
cumdose = sum(quant_strength, cumdose);
if first.ikn then cumdose = quant_strength; 
run;

data gc_cleaned_follow0;
set gc_cleaned_follow0;
by ikn;
cumdose = lag(cumdose);
if first.ikn then cumdose = .;
run;

data gc_cleaned_follow1;
set gc_cleaned_follow0;
by ikn;
if cumdose + quant_strength >= 450 then do;
    dayssupl2= dayssupl - ceil((450-cumdose)/daily);
	dayssupl = ceil((450-cumdose)/daily);
	quant_strength = dayssupl*daily;
	
	followdate = servdate + dayssupl;
	quant_strength2 = dayssupl2*daily;
	cumdose2 = cumdose+quant_strength;
end;
format followdate date9.;
run;


data gc_cleaned_follow11;
set gc_cleaned_follow1;
if followdate = . then delete;
run;

/* One servdate before follow date */
proc sort nodupkey data=gc_cleaned_follow11 (keep=ikn servdate) out = gc_cleaned_serv;
by ikn;
run;

data gc_cleaned_serv;
set gc_cleaned_serv;
rename servdate = preserv450;
run;

/* Follow Date */
proc sort nodupkey data=gc_cleaned_follow11 (keep=ikn followdate) out = gc_cleaned_follow;
by ikn;
run;

/* Original serv record modified before turning 450 */
proc sort nodupkey data=gc_cleaned_follow11 (keep=ikn servdate dayssupl quant_strength cumdose daily) out = gc_cleaned_split_modify;
by ikn;
run;

/* Merge in the follow date and prev servdate 450 to follow0, TODO: merge in preserv450 date and follow date info*/
data gc_followed0;
merge gc_cleaned_follow0 gc_cleaned_follow gc_cleaned_serv;
by ikn;
if servdate = preserv450 then delete;
run;


/* New follow date record at or after turned 450 */
proc sort nodupkey data=gc_cleaned_follow11 (keep=ikn followdate dayssupl2 quant_strength2 cumdose2 daily) out = gc_cleaned_split_record0;
by ikn;
run;

data gc_cleaned_split_record;
set gc_cleaned_split_record0;
rename followdate = servdate dayssupl2 = dayssupl quant_strength2 = quant_strength cumdose2 = cumdose;
run;


/* Append in preserv450 date, rest dates and follow date info*/
data gc_final0;
set gc_cleaned_split_modify gc_followed0 gc_cleaned_split_record;
by ikn;
run;

proc sort data=gc_final0 out=gc_final1;
by ikn servdate;
run;

data gc_finalle450;
merge gc_final1(drop=followdate preserv450) gc_cleaned_follow;
by ikn;
run;

/* To merge with gc_finalle450. p.s. seems like no need to proc sort by ikn. this by ikn already ordered it */
data gc_final_both450;
set gc_final_1st450 gc_finalle450;
by ikn;
run;


data gc_final;
set gc_final_both450;
enddate = servdate+dayssupl;
format enddate date9.;
admcensdate = followdate+365;
format admcensdate date9.;
if dayssupl=0 and quant_strength=0 then delete;
if servdate>admcensdate then delete;
if servdate<followdate then delete;
run;

/* Notice: patient 5 has followdate as . because cumulative dose did not reach 450, so automatically excluded */

/*** Now exclude records before followdate, and keep only ikn that has record after followdate using nodupkey ***/
data gc_final_copy;
set gc_final;
run;

/* This is the final ikn list, obtained by taking the first dispense info after followdate*/
proc sort nodupkey data=gc_final_copy out=gc_final_ikn_list(keep=ikn);
by ikn;
run;


/*************************** Fake data for death or fracture ***************************/
data gc_death;
input ikn:$ 3. deathdate:yymmdd8. fracdate:yymmdd8. ;
datalines;
1  20160723  20160623
3  .         20120306
4  20130201  .
5  20170809  .
6  .         20190802
7  20131105  20131110
;  
run;

data gc_death;
set gc_death;
format deathdate date9.;
format fracdate date9.;
run;

data gc_final_death0;
merge gc_final(in=a) gc_death(in=b);
by ikn;
if a;
run;
