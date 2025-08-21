/*HPTN 037*/
/*create analysis dataset*/
/*To replicate Soc Sci paper*/
/*Then estimate Individual and Disseminated effects*/
/*and individual and disseminated package component effects*/
/*add race variable - April 2017*/
options ps=60 ls=80 nodate pageno=1 nofmterr;
dm log "clear;" continue; dm out "clear;" continue;
libname hptn 'Z:\Desktop\hptn sasdata\sasdata';
/*load formats*/
*created in SAS 34 bit, cannot use in 64 bit;
*will manually enter later;
*options fmtsearch = (hptn.formats);

%INC '\\Mac\Google Drive\HPTN 037\macros\dupcheck.txt';
*add table 1 macro;
options nocenter ps=78 ls=125 replace nofmterr;
*mautosource sasautos= 'M://HPTN 037/macros';
/*read in sas datasets*/
/*Key variables*/
data keys;
set hptn.keys;
run;

proc freq data=keys;
tables SSstatus*index;
where site=222;
run;
%DupCheck(Dataset = keys, ByVar = uid , PrintVar = uid , Remove = );
/*demographics*/
data dm;
set hptn.dm;
run;

%DupCheck(Dataset = dm, ByVar = uid , PrintVar = uid , Remove = );
/*end of study inventory*/
data esi;
set hptn.esi;
run;

%DupCheck(Dataset = esi, ByVar = uid , PrintVar = uid , Remove = );
/*screening HIV status*/
data ss;
set hptn.ss;
run;
proc freq data=ss;
tables SSstatus;
run;
%DupCheck(Dataset = ss, ByVar = uid , PrintVar = uid SSstatus , Remove =);

/*risk assessment*/
data ra;
set hptn.ra;
run;

%DupCheck(Dataset = ra, ByVar = uid visnum, PrintVar = uid visnum visit RA1dt ra1drink ra1drunk, Remove =);

*keep first observation (ordered by calendar time);
proc sort data=ra;
by uid visnum ra1dt;
run;
proc sort data=ra nodupkey out=ra2;
by uid visnum;
run;
proc freq data=ra2;
tables visnum visit;
run;
%DupCheck(Dataset = ra2, ByVar = uid visnum, PrintVar = uid visnum visit RA1dt ra1drink ra1drunk, Remove =);

/*social impact assessment*/
data sia;
set hptn.sia;
run;

%DupCheck(Dataset = sia, ByVar = uid visnum, PrintVar = uid visnum visit SIAdt SIAgood, Remove =);

proc sort data=sia;
by  uid visnum siadt;
run;
*keep first observation (ordered by calendar time);
proc sort data=sia nodupkey out=sia2;
by uid visnum ;
run;
%DupCheck(Dataset = sia2, ByVar = uid visnum, PrintVar = uid visnum visit SIAdt SIAgood, Remove =);

/*lab information*/
data ll;
set hptn.ll;
if visnum=. then delete; *244 missing visit from lab data;
run;

proc sort data=ll;
by uid visnum studyday;
run;
*keep first observation (ordered by calendar time);
proc sort data=ll nodupkey out=ll2;
by uid visnum ;
run;
%DupCheck(Dataset = ll2, ByVar = uid visnum, PrintVar = uid visnum visit LLsmp1dt, Remove =);
proc print data=ll (obs=10);
run;

/*intervention uptake*/
/*variables on two forms*/
/*we do not analyze this*/
/*if we want to analyze, we will need to figure out the difference between ex and ext datasets*/
data ex;
set hptn.ex;
run;
proc sort data=ex;
by uid visnum EXdt;
run;
proc sort data=ex nodupkey out=ex2;
by uid visnum ;
run;
%DupCheck(Dataset = ex2, ByVar = uid visnum, PrintVar = uid visnum visit exnpeo exnconv, Remove =);

/*intervention uptake*/
data ext;
set hptn.ext;
run;

proc sort data=ext;
by uid visnum EXTdt;
run;
proc sort data=ext nodupkey out=ext2;
by uid visnum ;
run;
%DupCheck(Dataset = ext2, ByVar = uid visnum, PrintVar = uid visnum visit exTnpeo exTnconv, Remove =);
proc format ;
value yes1fmt 1='Yes'
              0='No';
			  run;

data exfin;
merge ex2(in=ini) ext2(rename=(exTnpeo=exnpeo exTnconv=exnconv exThiv=exhiv) in=inj);
by uid visnum;
/*Talked to five or more conversations above HIV risk reduction*/
if exnpeo>=5 then fivepeo=1;   else if 0<exnpeo<5 then fivepeo=0; else if exhiv = 2 then fivepeo=0;
label fivepeo ='Talked to five or more people about HIV in the last 6 months';
format fivepeo yes1fmt.;
/*had ten or more conversations above HIV risk reduction*/
if exnconv>=10 then tenconv=1; else if 0<exnconv<10 then tenconv=0; else if exhiv = 2 then tenconv=0;
label tenconv='Had ten or more convserations about HIV in the last 6 months';
format tenconv yes1fmt.;
/*had at least one conversation about HIV risk reduction*/
if exhiv=1 then myexhiv=1; else if exhiv=2 then myexhiv=0;
format myexhiv yes1fmt.;
label myexhiv='Had at least one convseration about HIV in the last 6 months';
if ini then inex=1; else inex=0; if inj then inext=1; else inext=0;
run;
proc freq data=exfin;
tables exhiv*exnpeo*fivepeo exhiv*exnconv*tenconv inex*inext/list missing;
run;
%DupCheck(Dataset = exfin, ByVar = uid visnum, PrintVar = uid visnum exnpeo exnconv, Remove =);
%DupCheck(Dataset = ra, ByVar = uid visit, PrintVar = uid visnum visit ra1drink ra1drunk, Remove = );
%DupCheck(Dataset = sia, ByVar = uid visit, PrintVar = uid visnum visit SIAdt SIAbad SIAgood, Remove = );
%DupCheck(Dataset = ll, ByVar = uid visit, PrintVar = uid visnum visit studyday lleia llwb1s llwb2s llwb3s, Remove = );
/*Combine all FUP data into one dataset*/
/*all fup data*/
data allfup;
merge ll2 sia2 ra2(in=ini) exfin;
by uid visnum;
/*indicator for being in Risk Assessment Data*/
if ini then inrisk=1; else inrisk=0;
format inrisk yes1fmt.;
label inrisk ='Has Risk Assessment Form?';
run;
proc freq data=allfup;
tables visit*visnum visnum inrisk*visnum/list missing;
run;
proc freq data=allfup;
tables visnum;
where inrisk=1;
run;
/*missed visit form*/
data mv;
set hptn.mv;
*create visnum variable;
     if visit=301 then visnum=6;
else if visit=401 then visnum=12;
else if visit=501 then visnum=18;
else if visit=601 then visnum=24;
else if visit=701 then visnum=30;
run;
proc freq data=mv;
tables visnum;
run;

%DupCheck(Dataset = mv, ByVar = uid visnum);

proc sort data=keys;
by uid;
run;
proc sort data=allfup;
by uid visnum;
run;
proc sort data=mv;
by uid visnum;
run;
proc format ;
value package 1='Peer Education Baseline'
              2='Peer Education Booster'
              3='HCT alone'
              4='None';
run;
/*index screening form*/
data is;
set hptn.is;
run;
%DupCheck(Dataset = is, ByVar = uid );

proc sort data=is;
by uid;
run;
proc sort data=dm;
by uid;
run;
proc print data=allfup;
var uid nkid;
run;
/*add single record per person datasets to analysis dataset*/
data pallfup;
*take each subjects network id at baseline for each followup visit;
merge allfup(drop=nkid  in=ini) keys(drop=visit) is(in=inj drop=visit) DM(drop=visit);
/*create unique network IDs*/
mynkid=cats(site,nkid);
if inj then myindex=1; else myindex=0;
by uid;
if ini;
run;

proc freq data=pallfup;
tables myindex visnum treat myindex*index;
run;
/*combine follow-up with missed visits form*/
proc sort data=pallfup;
by uid visnum;
run;
proc sort data=mv;
by uid visnum;
run;
data pallfupvis;
merge pallfup(in=ini) mv(in=inj) ;
by uid visnum;
*network ID does not change over time by design;
if ini then infup=1; else infup=0;
if inj then inmv=1; else inmv=0;
*create indicator for visit missing;
*There are 417 who have a risk assessment form and a MV form;
if inmv=1 and inrisk ne 1 then missvis=1; else missvis=0;
run;

proc freq data=pallfupvis;
tables missvis inmv*inrisk/list missing;
run;
proc sort data=pallfupvis;
by mynkid  visnum descending myindex ;
run;
/*Define a variable to determine if the index attended the visit*/
data pallfupvis2;
set pallfupvis;
by mynkid  visnum descending myindex ;
if first.visnum then totalindexvis=0;
retain totalindexvis;
totalindexvis=sum(totalindexvis,myindex);
if totalindexvis>0 then indexattend=1; else indexattend=0;
run;
/*For the package components the index had to attend the visit to be exposed*/
proc freq data=pallfupvis2;
tables indexattend*missvis/list missing;
where treat=1;
run;
/*clusters not randomized to treatment - index missed visits*/
/*did not change exposure to package classificaition*/
proc sort data=pallfupvis2;
by mynkid;
run;
proc print data=pallfupvis2;
var mynkid uid visnum myindex missvis totalindexvis indexattend treat;
*where missvis =0;
where mynkid='222272';
run;

/*Define package components received by the participant*/
/*Define package components received by the network members*/
proc sort data=pallfupvis2;
by mynkid uid visnum;
run;
/*variables in Table 3 over time*/
/*combine format libraries*/
proc format ;
value myvis 1='Baseline'
            2='6'
			3='12'
			4='18'
			5='24'
			6='30';
value PEboost 1='None'
              2='One'
			  3='Two';
			run;
data pallfupvis3;
set pallfupvis2;
/*define time-varying package based on protocol*/
*what participant received;
if treat=1 and missvis=0 and myindex=1 then do;
HCT=1;
if visnum=0 then do; package=1;  PEbase=1; PEboost=0; end;
else if visnum=6 then do; package=2;  PEbase=0; PEboost=1; end;
else if visnum=12 then do; package=2;  PEbase=0; PEboost=1; end;
else do; package=3;  PEbase=0; PEboost=0; end;
end;
else if treat=2 and missvis=0 and myindex=1 then do; HCT=1; package=3;  PEbase=0; PEboost=0; end;
else if myindex=0 and missvis=0 then do; HCT=1; package=3; PEbase=0; PEboost=0; end;
else if missvis=1 then do; HCT=.; PEbase=.; package=.; PEboost=.; end;
*what network received (both index and network members); 
if treat=1 and missvis=0 then do;
nHCT=1;
if visnum=0 and indexattend=1 then do; npackage=1;  nPEbase=1; nPEboost=0; end;
else if visnum=6 and indexattend=1 then do; npackage=2;  nPEbase=0; nPEboost=1; end;
else if visnum=12 and indexattend=1 then do; npackage=2;  nPEbase=0; nPEboost=1; end;
else do; npackage=3; nHCT=1; nPEbase=0; nPEboost=0; end;
end;
else if treat=2 and missvis=0 then do; npackage=3; nHCT=1; nPEbase=0; nPEboost=0; end;
else if myindex=0 and missvis=0 then do; npackage=3; nHCT=1; nPEbase=0; nPEboost=0; end;
else if missvis=1 then do; nHCT=0; nPEbase=0; npackage=4; nPEboost=0; end;
format package npackage package. HCT nHCT yes1fmt. PEbase nPEbase yes1fmt.;
label package='Component of package recieved by participant';
label HCT =' HCT recieved by participant';
label PEbase='Patient recieved education at baseline';
label PEboost ='Patient recieved booster at visit';
label npackage='Component of package recieved by network';
label nHCT =' HCT recieved by network';
label nPEbase='Network recieved education at baseline';
label nPEboost ='Network recieved booster at visit';
if treat=1 then mytreat=1; else mytreat=0;
visnumcat=visnum;
  if visnumcat=0 then visnumcat2=1;
else if visnumcat=6 then visnumcat2=2;
else if visnumcat=12 then visnumcat2=3;
else if visnumcat=18 then visnumcat2=4;
else if visnumcat=24 then visnumcat2=5;
else if visnumcat=30 then visnumcat2=6;
label visnumcat2='Visit';
format visnumcat2 myvis.;
run;
proc freq data=pallfupvis3;
tables mytreat npackage package nHCT HCT nPEbase PEbase nPEboost PEboost missvis visnumcat2;
run;
proc print data=pallfupvis3 (obs=200);
var uid mynkid visnum treat myindex indexattend missvis PEboost nPEboost;
where missvis=0;
run;
*Packages over time;
*%table1(data=pallfupvis3,
exposure=visnumcat2,
varlist=mytreat npackage package nHCT HCT nPEbase PEbase nPEboost PEboost,
cat=mytreat  nHCT HCT nPEbase PEbase nPEboost PEboost,
poly= npackage package,
rtftitle=Package components by visit,
landscape=F, ageadj=F,
nortf=F, file=table4_package);
/*count number of times participant/network recieved a booster*/
proc sort data=pallfupvis3;
by mynkid uid visnum;
run;
data allfupvis;
set pallfupvis3;
by mynkid uid visnum;
retain  PEboostTot nPEboostTot;
if first.uid then do;  PEboostTot=.;  nPEboostTot=.; end;
PEboostTot=sum(PEboostTot,PEboost);
nPEboostTot=sum(nPEboostTot, nPEboost);
run;
proc print data=allfupvis (obs=100);
var uid visnum missvis treat myindex nPEboost nPEboostTot PEboost PEboostTot;
where uid in (0000256,  0000288, 0000324, 0000631,0000777) and missvis=0;
run;

/*Some have visit data but also have missing visit form*/
/*they are missing treatment information, so not in analysis*/
proc freq data=allfupvis;
tables visnum*missvis infup*inmv*inrisk/list missing;
where treat ne .;
*where missvis=0;
run;

proc sort data=allfupvis;
by treat;
run;
proc freq data=allfupvis;
by treat;
tables package*visnum;
where inrisk=1;
run;
proc print data=allfupvis;
where package=. and treat ne .;
run;
%DupCheck(Dataset = allfupvis, ByVar = uid visnum, PrintVar = uid visnum visit , Remove = );

/*number of follow-up visits*/
/*number of subjects*/
/*Follow-up lower than what is reported in Soc Sci Paper Text*/
proc freq data=allfupvis;
tables visnum;
*where treat ne .;
where inrisk ne .;
run;

proc freq data=allfupvis ;
tables treat*missvis/list missing;
run;
/*create flowchart*/
proc freq data=allfupvis;
tables treat;
where visnum=0 and missvis=0;
run;
proc freq data=allfupvis;
tables visnum;
where site=222 and treat=1 and missvis=0;
run;
proc freq data=allfupvis;
tables visnum;
where site=570 and treat=1 and missvis=0;
run;
proc freq data=allfupvis;
tables visnum;
where site=222 and treat=2 and missvis=0 ;
run;
proc freq data=allfupvis;
tables visnum;
where site=570 and treat=2 and missvis=0;
run;
/*number of networks at each visit*/

proc sort data=allfupvis nodupkey out=checka;
by  visnum mynkid ;
where site=222 and treat=1 and missvis=0;
run;
proc freq data=checka;
tables visnum;
run;
proc sort data=allfupvis nodupkey out=checkb;
by  visnum mynkid ;
where site=570 and treat=1 and missvis=0;
run;
proc freq data=checkb;
tables visnum;
run;
proc sort data=allfupvis nodupkey out=checkc;
by  visnum mynkid ;
where site=222 and treat=2 and missvis=0;
run;
proc freq data=checkc;
tables visnum;
run;
proc sort data=allfupvis nodupkey out=checkd;
by  visnum mynkid ;
where site=570 and treat=2 and missvis=0;
run;

/*Create Variables in Table 2 of Social Sci manuscript*/
/*Define formats*/
proc format ;
value trt 1='Treatment'
          2='Control';
value site 222='Philadelphia'
           570='Chiang Mai';
value mysite 1='Philadelphia'
             2='Chiang Mai';
value yes1fmt 1='Yes'
              0='No';
value dc 1='Did not drink'
         2='Never got drunk'
		 3='Sometimes got drunk'
		 4='Always got drunk';
value injd 1='0-5'
           2='6-14'
		   3='15-29'
		   4='Everyday';
run;
proc sort data=allfupvis;
by uid;
run;
data allfupvis1;
set allfupvis;
by uid;
/*drinking cat var*/
   if RA1drunk=. then drinkcat=1;
else if RA1drunk=5 then drinkcat=2;
else if RA1drunk=1 then drinkcat=4;
else if RA1drunk  ne . then drinkcat=3;
format drinkcat dc.;
label drinkcat = 'Alcohol Use';
/*number of sex partners in the last month*/
  totalpart = sum(RA7nfpt,RA7nmpt);
  label totalpart='Total number of sexual partners in the last month';
  if totalpart>1 then moreone=1;
  else if .<totalpart<=1 then moreone=0;
  else if RA7vagsx=2 then moreone=0;
  else moreone=.;
  label moreone='More than one sex partner';
  format moreone yes1fmt.;
  if totalpart>=1 then atleastone=1; 
  else if moreone ne . then atleastone=0;
  format atleastone yes1fmt.;
  label atleastone='At least one sex partner';
 /*unprotected sex in the last week*/
  totalsex=sum(RA7nsxpp,RA7nsxop);
  label totalsex='Total number of sex acts in the last week';
  totalpwk=sum(RA7ncdpp,RA7ncdop);
  label totalpwk='Total number of protected sex acts in the last week';
  totalunpwk=totalsex-totalpwk;
  label totalunpwk='Total number of unprotected sex acts in the last week';
  if RA7vagsx=2 then unpwk=0;
  else if RA7nsxpp=0 and RA7nsxop=0 then unpwk=0;
  else if totalunpwk>0 then unpwk=1;
  else if totalunpwk=0 then unpwk=0; 
  else if totalsex=0 then unpwk=0; else unpwk=.;
  label unpwk='Unprotected sex in the last week';
  /*no condom use at all*/
  if totalpwk=0 then nocondom=1; else if unpwk ne . then nocondom=0;
  label nocondom='No condom use for any sex act last week';
  /*unprotected sex with non primary partner*/
     if RA7ncdop ne . and (RA7nsxop-RA7ncdop)>0 then unpnon=1; 
else if RA7ncdop ne . and (RA7nsxop-RA7ncdop)=0 then unpnon=0;
else if RA7nsxop=0 then unpnon=0;
else if RA7vagsx=2 then unpnon=0;
else if RA7nsxop<0 then unpnon=.;
else unpnon=.;
   label unpnon='Unprotected sex with non-primary partner last week';
/*no condom use at all with nonprimary partner*/
if RA7ncdop=0 then nocondomnop=1; else if unpnon ne . then nocondomnop=0;
 label nocondomnop='No condom use with nonprimary partner last week';
/*number of days injected*/
if RA3inj6m=1 and visnum=0 then do; *only among those who reported injection in the last 6 months;
     if .<RA3ndays<6    then injectdays=1;
else if 6<=RA3ndays<15  then injectdays=2;
else if 15<=RA3ndays<30 then injectdays=3;
else if   RA3ndays>=30  then injectdays=4;
else if  RA3injlm=2     then injectdays=1;
format injectdays injd.;
label injectdays='Number of days injected in last month';
if RA3ndays>14 then inject14=1; else if injectdays ne . then inject14=0;
format inject14 yes1fmt.;
label inject14='Injected more than 14 days in the last month';
end;
run;

/*add injection drug use at baseline variable*/
proc sort data=allfupvis1;
by uid visnum;
run;
data vis0;
set allfupvis1;
by uid visnum;
if visnum=0;
run;
proc sort data=allfupvis1;
by uid;
run;
proc sort data=vis0;
by uid;
run;
data allfupvis2;
merge allfupvis1(in=ini) vis0(keep=uid RA3inj6m RA3injlm rename=(RA3inj6m=injectbase RA3injlm=injectbaselm)) ;
by uid;
/*injection risk behaviors*/
if visnum>0 and injectbase=1 then do;
     if RA3injlm=1 then injectlmof=1;
else if RA3injlm=2 then injectlmof=0; 
else if inrisk=1 then injectlmof=0;
format injectlmof yes1fmt.;
label injectlmof ='Injected in last month (and injected at baseline)';
if RA3ndays>14 then inject14f=1; else if injectlmof ne . then if inrisk=1 then inject14f=0;
else if inrisk=1 then inject14f=0;
format inject14f yes1fmt.;
end;
label inject14f='Injected more than 14 days in the last month (and injected at baseline)';
if visnum>0 and injectbase=1 then do;
  if .<RA3ndays<6 then injectdays=1;
else if 6<=RA3ndays<15 then injectdays=2;
else if 15<=RA3ndays<30 then injectdays=3;
else if RA3ndays>=30 then injectdays=4;
else if  RA3injlm=2 then injectdays=1;
format injectdays injd.;
label injectdays = 'Number of days injected in last month (among baseline IDUs)';
if RA3ndays>14 then inject14=1; else if injectdays ne . then inject14=0;
format inject14 yes1fmt.;
label injectdays = 'Injected more than 14 days in the last month (among baseline IDUs)';
/*equipment sharing in the last month*/
if RA3nwatr>0 then indwater=1; else if RA3nwatr=0 or Ra3injlm=2 then indwater=0; else if inrisk=1 then indwater=0;
if RA3ncook>0 then indcook=1;  else if RA3ncook=0 or Ra3injlm=2 then indcook=0;  else if inrisk=1 then indcook=0;
if RA3ncotn>0 then indcotn=1;  else if RA3ncotn=0 or Ra3injlm=2 then indcotn=0;  else if inrisk=1 then indcotn=0;
if RA3nload>0 then indload=1;  else if RA3nload=0 or Ra3injlm=2 then indload=0;  else if inrisk=1 then indload=0;
if RA3nndl>0 then indnndl=1;   else if RA3nndl=0 or  Ra3injlm=2 then indnndl=0;  else if inrisk=1 then indnndl=0;
/*Those missing RA4dkwel either did not inject in the last month or last 6 months*/
/*Two did not respond and are classified as no*/
if RA4dkwel=1 or RA4dkwel=2 or RA4dkwel=3 or RA4dkwel=4 then inddkwel=1; else if inrisk=1 then inddkwel=0;
if RA5clean=2 then noclean=1; else if RA5clean=1 or RA3injlm=2 then noclean=0; else if inrisk=1 then noclean=0;
if RA4evpas=1 then indpass=1; else if inrisk=1 then indpass=0; 
if RA4evuse=1 then induse=1;  else if inrisk=1 then induse=0;
if RA4othrs=1 then indgall=1; else if inrisk=1 then indgall=0;
end;
label indwater='Risk Behavior: Shared rinse water';
label indcook='Risk Behavior: Shared cooker';
label indcotn='Risk Behavior: Shared cotton';
label indload='Risk Behavior: Front/back loading';
label indnndl='Risk Behavior: Injected with unclean syringe';
label inddkwel='Risk Behavior: Injected with people not known well';
label noclean='Risk Behavior: Cleaned needles'; *this outcome is not in Soc Sci Paper tables;
label indpass='Risk Behavior: Passed on a synringe';
label induse='Risk Behavior: Used a synringe after someone else';
label indgall='Risk Behavior: Injected in a shooting gallergy';
/*create binary indicator for randomized treatment*/
if treat=1 then mytreat=1; else mytreat=0;
label mytreat = 'Randomized Treatment';
format mytreat yes1fmt.;
/*There was no 3 month visit in the study, so delete*/
if visnum=3 then delete;
run;

proc freq data=allfupvis2;
tables RA4dkwel inddkwel RA5clean noclean 
RA4evpas indpass
 RA4evuse induse
 RA4othrs indgall;
 *where inrisk=1 and injectbase=1 and visnum>0;
 where visnum>0 and Ra3injlm=2 and inrisk=1 and injectbase=1;
 run;
 /*formats for other variables*/
proc format ;
value sex 1='Male'
          0='Female';
value ageg 1='18-20'
           2='21-30'
		   3='31-40'
		   4='40+';
value ed  1 = 'No schooling'                            
          2 = 'Primary schooling'                       
          3 = 'Secondary schooling'                
          4 = 'Completed secondary'           
          5 = 'Vocational or trade schooling'           
          6 = 'University or comm. college';         
value emp    1='Full time'                               
             2='Part-time/ocassional'                              
            3='Unemployed';
value marit 1='Single'
            2='Married'
			3='Living with partner'
			4='Separated/DIvorced/Widowed';
value myracegrp 1='White'
                2='Black'
				3='Natie Hawaiian/Pacific Islander'
				4='Asian'
				5='American Indian/Alaskan Native'
				6='Other'
				7='Multiple races'
				8='None';
			run;
proc freq data=allfupvis2;
tables DM1latin DM1white*DM1black*DM1nhopi*DM1asian*DM1aian*DM1photh/list;
run;
* Add Time-varying covariates;
data allfupvis3;
set allfupvis2;
/*age categories, marital satus, education and employment*/
     if 18<=age<=20 then agegrp=1;
else if 20<age<=30 then agegrp=2;
else if 30<age<=40 then agegrp=3;
else if age>40     then agegrp=4;
format agegrp ageg.;
label agegrp='Age (years)';

if DM1latin=1 then hispanic=1; else if DM1latin=2 then hispanic=0;
label hispanic='Hispanic';
format hispanic yes1fmt.;

     if DM1white=1 and DM1black=0 and DM1nhopi=0 and DM1asian=0 and DM1aian=0 and DM1photh=0 then myrace=1;
else if DM1black=1 and DM1white=0 and DM1nhopi=0 and DM1asian=0 and DM1aian=0 and DM1photh=0 then myrace=2;
else if DM1nhopi=1 and DM1white=0 and DM1black=0 and DM1asian=0 and DM1aian=0 and DM1photh=0 then myrace=3;
else if DM1asian=1 and DM1white=0 and DM1black=0 and DM1nhopi=0 and DM1aian=0 and DM1photh=0 then myrace=4;
else if DM1aian=1  and DM1white=0 and DM1black=0 and DM1nhopi=0 and DM1asian=0 and DM1photh=0 then myrace=5;
else if DM1photh=1 and DM1white=0 and DM1black=0 and DM1nhopi=0 and DM1asian=0 and DM1aian=0 then myrace=6;
else if DM1white=0 and DM1black=0 and DM1nhopi=0 and DM1asian=0 and DM1aian=0 and DM1photh=0 then myrace=8;
else if DM1photh ne . and DM1white ne . and DM1black ne .  and DM1nhopi ne .  and DM1asian ne .  and DM1aian ne . then myrace=7;
else myrace=.;
format myrace myracegrp.;
label myrace='Race';

if myrace =1 then nonwhite=0; else if myrace ne . then nonwhite=1; format nonwhite yes1fmt.;
label nonwhite='Nonwhite Race';

     if DM2educ=1 then myeduc=1;
else if DM2educ=2 then myeduc=2;
else if DM2educ=3 then myeduc=3;
else if DM2educ=4 then myeduc=4;
else if DM2educ=5 then myeduc=5;
else if DM2educ in (6,7,8) then myeduc=6;
format myeduc ed.;
label myeduc='Education';

     if DM2EMPL=1 then myemploy=1;
else if DM2EMPL in (2,3) then myemploy=2;
else if DM2EMPL=4        then myemploy=3;
format myemploy emp.;
label myemploy='Employment';

     if DM2marit=1 then mymarit=1;
else if DM2marit=2 then mymarit=2;
else if DM2marit=3 then mymarit=3;
else if DM2marit in (4,5,6) then mymarit=4;
format mymarit marit.;
label mymarit='Marital Status';
/*recode variables to 0,1*/
array allvars{20} DM1sex 
        RA2crack RA2cocai RA2tkamp RA2smamp RA2hero RA2benzo
        RA2prog RA1strt RA1jail
        
        RA3hero RA3spdbl RA3heram RA3cocai RA3amph 
        RA4evpas RA4evuse
        RA5clean RA3inj6m RA3injlm;
  do i=1 to 20;
          if allvars(i)=2 then allvars(i)=0;
	 else if allvars(i)=1 then allvars(i)=1;
	 end;
run;
proc freq data=allfupvis3;
tables nonwhite myrace hispanic;
where site=222;
run;
/*prepare variables for models*/
data panalydata;
set allfupvis3;
/*create indicator variables for covariates*/
if mymarit=1 then indmarit=1; else if mymarit ne . then indmarit=0;
if myeduc>4 then indeduc=1; else if myeduc ne . then indeduc=0;
if myemploy=3 then indunemp=1; else if myemploy ne . then indunemp=0;
*if heroin not used in last month or last 6 months then set variables to no;
if RA3inj6m=0 or RA3injlm=0 then do;
RA3hero=0; 
RA3spdbl=0; 
RA3heram=0; 
RA3cocai=0;  
RA3amph=0;  
end;
run;
proc sort data=panalydata;
by uid;
run;
data analydata;
set panalydata;
by uid;
*lag time-varying covariates by one visit;
*ensures correct temportal order (exposure and confounder occur before outcome);
array timevar(15) RA2crack_prev RA2cocai_prev RA2smamp_prev RA2hero_prev RA2benzo_prev
		RA2prog_prev RA1strt_prev RA1jail_prev
		drinkcat_prev injectdays_prev
		RA3hero_prev RA3spdbl_prev RA3heram_prev RA3cocai_prev RA3amph_prev;  
array original(15) RA2crack RA2cocai RA2smamp RA2hero RA2benzo
		RA2prog RA1strt RA1jail
		drinkcat injectdays
		RA3hero RA3spdbl RA3heram RA3cocai RA3amph;
array myretain(15) rRA2crack rRA2cocai rRA2smamp rRA2hero rRA2benzo
		rRA2prog rRA1strt rRA1jail
		rdrinkcat rinjectdays
		rRA3hero rRA3spdbl rRA3heram rRA3cocai rRA3amph;
do i=1 to 15;
retain rRA2crack rRA2cocai rRA2smamp rRA2hero rRA2benzo
		rRA2prog rRA1strt rRA1jail
		rdrinkcat rinjectdays
		rRA3hero rRA3spdbl rRA3heram rRA3cocai rRA3amph;
if first.uid then myretain(i)=.;
if original(i) >= 0 then myretain(i)=original(i);
timevar(i)=lag(myretain(i));
if first.uid then timevar(i)=original(i);
end;
run;

/*Define indicator variables for model*/
data analydata2a;
set analydata;
by uid visnum;
*create indicator variables for inject days previous;
if injectdays_prev=1 then injectdays1_prev=1; else if injectdays_prev ne . then injectdays1_prev=0;
if injectdays_prev=2 then injectdays2_prev=1; else if injectdays_prev ne . then injectdays2_prev=0;
if injectdays_prev=3 then injectdays3_prev=1; else if injectdays_prev ne . then injectdays3_prev=0;
if injectdays_prev=4 then injectdays4_prev=1; else if injectdays_prev ne . then injectdays4_prev=0;
*create indicator for inject days;
if injectdays=1 then injectdays1=1; else if injectdays ne . then injectdays1=0;
if injectdays=2 then injectdays2=1; else if injectdays ne . then injectdays2=0;
if injectdays=3 then injectdays3=1; else if injectdays ne . then injectdays3=0;
if injectdays=4 then injectdays4=1; else if injectdays ne . then injectdays4=0;

*define any injection risk behvaior;
if inddkwel=1 or indcook=1 or indwater=1 or indcotn=1 or indload=1 or indnndl=1 
or indpass=1 or induse=1 or indgall=1 then anyrisk=1;
else if inddkwel ne . or indcook ne . or indwater ne . or indcotn ne .  or indload ne .  
or indnndl ne .  or indpass ne . or induse ne . or indgall ne . then  anyrisk=0;
/*patient education at baseline*/
/*retain across follow-up visits*/
/*participant and network exposure*/
/*create variables that are any report of risk behavior over FUP*/
if first.uid then do; PEbaseret=PEbase; nPEbaseret=nPEbase; end;
retain PEbaseret nPEbaseret anyriskever anyinduse anyindwater anyindcook anyindcotn 
anyindload anyindpass anynoclean  anyinddkwel anyindgall;
if first.uid then do; anyriskever=0; anyinduse=0; anyindwater=0;
anyindcook=0; anyindcotn=0; anyindload=0; anyindpass=0;
anynoclean=0; anyinddkwel=0 ; anyindgall=0; end;
if anyrisk=1 and visnum>0 then anyriskever=1;
if induse=1 and visnum>0 then anyinduse=1;
if indwater=1 and visnum>0 then anyindwater=1;
if indcook=1 and visnum>0 then anyindcook=1;
if indcotn=1 and visnum>0 then anyindcotn=1;
if indload=1 and visnum>0 then anyindload=1;
if indpass=1 and visnum>0 then anyindpass=1;
if noclean=1 and visnum>0 then anynoclean=1;
if inddkwel=1 and visnum>0 then anyinddkwel=1;
if indgall=1 and visnum>0 then anyindgall=1;
run;
/*lag the exposure by one visit*/
/*ensures exposure happens before outcome*/
/*correct temportal order*/
data analydata2;
set analydata2a;
by uid visnum;
/*boosters at previous visit*/
/*participant and network exposure*/
if missvis=0 then do; prev_PEboost=lag(PEboost); prev_nPEboost=lag(nPEboost);end;
if first.uid then do; prev_PEboost=PEboost; prev_nPEboost=nPEboost; end;
/*Total number of boosters*/
/*participant and network exposure*/
if missvis=0 then do; prev_PEboosttot=lag(PEboosttot); prev_nPEboosttot=lag(nPEboosttot); end;
if first.uid then do; prev_PEboosttot=PEboosttot; prev_nPEboosttot=nPEboosttot; end;
run;
proc print data=analydata2;
	var uid visnum missvis indexattend PEboosttot prev_PEboosttot  nPEboosttot prev_nPEboosttot PEboost prev_PEboost;
	*where lagPEboosttot>prev_PEboosttot;
	where uid in (0000256,  0000288, 0000324, 0000631,0000777) and missvis=0;
	run;

proc sort data=analydata2;
by uid visnum;
run;
data analydata3;
set analydata2;
by uid visnum;
/*create baseline versions of time-varying covariates*/
if first.uid then do; basedkwel=inddkwel; basecook=indcook;  basewater=indwater;
    basecotn=indcotn; baseload=indload; basenndl=indnndl; basepass=indpass;
    baseuse=induse; basegall=indgall; baseatleastone=atleastone;
    baseunpnon=unpnon; basenocondomnop=nocondomnop; baseunpwk=unpwk;
    basenocondom=nocondom; baseanyrisk=anyrisk;  
    basecrack=RA2crack; basecocai=RA2cocai; basesmap=RA2smamp; basehero=RA2hero; basebenzo=RA2benzo;
    baseprog=RA2prog; basestrt=RA1strt; basejail=RA1jail;
	basedrink=drinkcat;    baseinjectdays =injectdays; 
	baseheroi=RA3hero; basespdl=RA3spdbl; baseheram=RA3heram; basecocaii=RA3cocai; baseamph=RA3amph; 
end;
retain basedkwel basecook basewater basecotn baseload basenndl basepass
    baseuse basegall baseatleastone baseunpnon basenocondomnop baseunpwk
    basenocondom baseanyrisk 
    basecrack basecocai basesmap basehero basebenzo
    baseprog basestrt basejail basedrink
    baseinjectdays baseheroi basespdl baseheram basecocaii baseamph;
	*create indicator for sometimes or always got drunk (baseline);
	if basedrink =3 or basedrink=4 then basedrunk=1; else if basedrink ne . then basedrunk=0;
	*create indicator for sometimes or always got drunk (time-varying);
	if drinkcat =3 or drinkcat=4 then drunk=1; else if drinkcat ne . then drunk=0;
	/*lag drunk varible*/
    drunk_prev =lag(drunk);
    if first.uid then drunk_prev=drunk;
	*create indicator for inject days;
if baseinjectdays=1 then baseinjectdays1=1; else if baseinjectdays ne . then baseinjectdays1=0;
if baseinjectdays=2 then baseinjectdays2=1; else if baseinjectdays ne . then baseinjectdays2=0;
if baseinjectdays=3 then baseinjectdays3=1; else if baseinjectdays ne . then baseinjectdays3=0;
if baseinjectdays=4 then baseinjectdays4=1; else if baseinjectdays ne . then baseinjectdays4=0;

	*indicator for all covariates not missing;
	if DM1sex ne .  and age ne . and indmarit ne . and indeduc ne . and 
      indunemp ne . and basecrack ne . and basecocai ne . and basebenzo ne . and
      baseprog ne . and basestrt in (0,1) and basejail in (0,1) and basedrunk ne . then nomiss=1; else nomiss=0;
run;
proc print data=analydata3;
var uid visnum drinkcat drunk drunk_prev;
run;
/*create indicator for one and two boosters (time-varying)*/
proc sort data=analydata3;
by uid visnum;
run;
data analydata4;
set analydata3;
by uid visnum;
if missvis=0 then do;
if  prev_PEboosttot=1 then PEboost1=1; else if  prev_PEboosttot ne . then PEboost1=0;
if  prev_PEboosttot=2 then PEboost2=1; else if  prev_PEboosttot ne . then PEboost2=0;
if  prev_nPEboosttot=1 then nPEboost1=1; else if  prev_nPEboosttot ne . then nPEboost1=0;
if  prev_nPEboosttot=2 then nPEboost2=1; else if  prev_nPEboosttot ne . then nPEboost2=0;
end;
/*lag the previous exposure variables*/
/*used to model P(exposure)*/
/*if missvis=0 then lagprev_PEboost=lag(prev_PEboost);
if indexattend=1 then lagprev_nPEboost=lag(prev_nPEboost); */
if missvis=0 then do;
lagprev_PEboost=lag(prev_PEboost);lagprev_nPEboost=lag(prev_nPEboost);
end;
if first.uid then do ; lagprev_PEboost=prev_PEboost; lagprev_nPEboost=prev_nPEboost;
end;
if missvis=0 then do ;
lagPEboosttot=lag(prev_PEboosttot); 
lagnPEboosttot=lag(prev_nPEboosttot); end;
if first.uid then do; lagPEboosttot=prev_PEboosttot; lagnPEboosttot=prev_nPEboosttot; end;
run;
data last;
set analydata4;
by uid visnum;
if visnum <24 then drop=1; else drop=0;
if last.uid and visnum>0;
run;
proc freq data=last;
tables visnum drop;
run;
proc lifetest data=last ;
      time visnum*drop(0);
     ;
   run;
   /*check against Latkin 2009*/
   data checklatkin;
   set analydata4;
   by uid visnum;
   if injectbase=1 then do;
   /*define these risk measures for the first visit*/
   if .<RA3ndays<6 then injectdays=1;
else if 6<=RA3ndays<15 then injectdays=2;
else if 15<=RA3ndays<30 then injectdays=3;
else if RA3ndays>=30 then injectdays=4;
else if  RA3injlm=2 then injectdays=1;
format injectdays injd.;
label injectdays = 'Number of days injected in last month (among baseline IDUs)';
if RA3ndays>14 then inject14=1; else if injectdays ne . then inject14=0;
format inject14 yes1fmt.;
label injectdays = 'Injected more than 14 days in the last month (among baseline IDUs)';
/*equipment sharing in the last month*/
if RA3nwatr>0 then indwater=1; else if RA3nwatr=0 or Ra3injlm=2 then indwater=0; else if inrisk=1 then indwater=0;
if RA3ncook>0 then indcook=1;  else if RA3ncook=0 or Ra3injlm=2 then indcook=0;  else if inrisk=1 then indcook=0;
if RA3ncotn>0 then indcotn=1;  else if RA3ncotn=0 or Ra3injlm=2 then indcotn=0;  else if inrisk=1 then indcotn=0;
if RA3nload>0 then indload=1;  else if RA3nload=0 or Ra3injlm=2 then indload=0;  else if inrisk=1 then indload=0;
if RA3nndl>0 then indnndl=1;   else if RA3nndl=0 or  Ra3injlm=2 then indnndl=0;  else if inrisk=1 then indnndl=0;
/*Those missing RA4dkwel either did not inject in the last month or last 6 months*/
/*Two did not respond and are classified as no*/
if RA4dkwel=1 or RA4dkwel=2 or RA4dkwel=3 or RA4dkwel=4 then inddkwel=1; else if inrisk=1 then inddkwel=0;
if RA5clean=0 then noclean=1; else if RA5clean=1 or RA3injlm=2 then noclean=0; else if inrisk=1 then noclean=0;
if RA4evpas=1 then indpass=1; else if inrisk=1 then indpass=0; 
if RA4evuse=1 then induse=1;  else if inrisk=1 then induse=0;
if RA4othrs=1 then indgall=1; else if inrisk=1 then indgall=0;
*define any injection risk behvaior;
if inddkwel=1 or indcook=1 or indwater=1 or indcotn=1 or indload=1 or indnndl=1 
or indpass=1 or induse=1 or indgall=1 then anyrisk=1;
else if inddkwel ne . or indcook ne . or indwater ne . or indcotn ne .  or indload ne .  
or indnndl ne .  or indpass ne . or induse ne . or indgall ne . then  anyrisk=0;
end;
   if missvis=1 then delete;
   if first.uid;
   run;
   proc freq data=checklatkin;
   tables nonwhite myrace hispanic;
   run;
proc means data= checklatkin;
var DM1sex age DM1age agegrp mymarit myeduc myemploy  nonwhite myrace hispanic
        RA2crack RA2cocai RA2smamp RA2hero RA2benzo
		RA2prog RA1strt RA1jail
		drinkcat moreone unpwk unpnon 
		RA3inj6m RA3injlm
		RA3hero RA3spdbl RA3heram RA3cocai RA3amph 
        injectdays  indwater indcook indcotn indload indpass induse
        noclean inddkwel;
		where mytreat=0;
		run;
proc freq data=checklatkin;
tables DM1sex agegrp mymarit myeduc myemploy nonwhite myrace hispanic
RA2crack RA2cocai RA2smamp RA2hero RA2benzo
RA2prog RA1strt RA1jail
		drinkcat moreone unpwk unpnon 
RA3inj6m RA3injlm;
where mytreat=1;
run;
proc freq data=checklatkin;
tables RA3hero RA3spdbl RA3heram RA3cocai RA3amph ;
where RA3injlm=1 and mytreat=1;
run;
proc freq data=checklatkin;
     tables   injectdays  indwater indcook indcotn indload indpass induse
        noclean RA5clean inddkwel ;
where injectbase=1 and mytreat=1;
run;
/*output final dataset*/
/*only keep Philadelphia site only*/
data hptn.analydataFINAL_Apr17;
set analydata4;
if site=222;
run;
/*risk behaviors only assessed among those who reported injection drug use 6 months prior to baseline*/
proc freq data=hptn.analydataFINAL_Apr17;
tables site indwater indcook indcotn indload indpass induse
        noclean inddkwel indgall injectbase;
		where injectbase=1;
run;
/*Check table 1*/
data test;
set hptn.analydataFINAL_Apr17;
by uid visnum;
if injectbase=1 then do;
 if .<RA3ndays<6 then injectdays=1;
else if 6<=RA3ndays<15 then injectdays=2;
else if 15<=RA3ndays<30 then injectdays=3;
else if RA3ndays>=30 then injectdays=4;
else if  RA3injlm=2 then injectdays=1;
format injectdays injd.;
label injectdays = 'Number of days injected in last month (among baseline IDUs)';
if RA3ndays>14 then inject14=1; else if injectdays ne . then inject14=0;
format inject14 yes1fmt.;
label injectdays = 'Injected more than 14 days in the last month (among baseline IDUs)';
/*equipment sharing in the last month*/
if RA3nwatr>0 then indwater=1; else if RA3nwatr=0 or Ra3injlm=2 then indwater=0; else if inrisk=1 then indwater=0;
if RA3ncook>0 then indcook=1;  else if RA3ncook=0 or Ra3injlm=2 then indcook=0;  else if inrisk=1 then indcook=0;
if RA3ncotn>0 then indcotn=1;  else if RA3ncotn=0 or Ra3injlm=2 then indcotn=0;  else if inrisk=1 then indcotn=0;
if RA3nload>0 then indload=1;  else if RA3nload=0 or Ra3injlm=2 then indload=0;  else if inrisk=1 then indload=0;
if RA3nndl>0 then indnndl=1;   else if RA3nndl=0 or  Ra3injlm=2 then indnndl=0;  else if inrisk=1 then indnndl=0;
/*Those missing RA4dkwel either did not inject in the last month or last 6 months*/
/*Two did not respond and are classified as no*/
if RA4dkwel=1 or RA4dkwel=2 or RA4dkwel=3 or RA4dkwel=4 then inddkwel=1; else if inrisk=1 then inddkwel=0;
if RA5clean=0 then noclean=1; else if RA5clean=1 or RA3injlm=2 then noclean=0; else if inrisk=1 then noclean=0;
if RA4evpas=1 then indpass=1; else if inrisk=1 then indpass=0; 
if RA4evuse=1 then induse=1;  else if inrisk=1 then induse=0;
if RA4othrs=1 then indgall=1; else if inrisk=1 then indgall=0;
*define any injection risk behvaior;
if inddkwel=1 or indcook=1 or indwater=1 or indcotn=1 or indload=1 or indnndl=1 
or indpass=1 or induse=1 or indgall=1 then anyrisk=1;
else if inddkwel ne . or indcook ne . or indwater ne . or indcotn ne .  or indload ne .  
or indnndl ne .  or indpass ne . or induse ne . or indgall ne . then  anyrisk=0;
end;
if missvis=1 then delete;
if first.uid;
run;
proc freq data=test;
tables site indwater indcook indcotn indload indpass induse
        noclean RA5clean inddkwel indgall ;
run;
proc print data=test;
var uid visnum anyrisk inddkwel indcook indwater indcotn indload indnndl indpass induse indgall;
run;
proc freq data=test;
tables anyrisk inddkwel*indcook*indwater*indcotn*indload*indnndl*indpass*induse*indgall/list missing;
run;
title 'Intevention';

proc means data= test;
var DM1sex age DM1age agegrp mymarit myeduc myemploy nonwhite myrace hispanic
        RA2crack RA2cocai RA2smamp RA2hero RA2benzo
		RA2prog RA1strt RA1jail
		drinkcat moreone unpwk unpnon 
		RA3inj6m RA3injlm
		RA3hero RA3spdbl RA3heram RA3cocai RA3amph 
        injectdays  indwater indcook indcotn indload indpass induse
        noclean inddkwel;
		where mytreat=1;
		run;
/*drugs injected last month among those reported injecting last month*/
proc freq data=test;
tables RA3hero RA3spdbl RA3heram RA3cocai RA3amph ;
where RA3injlm=1 and mytreat=1;
run;
proc freq data=test;
tables injectbase RA3injlm;
where mytreat=1;
run;
proc freq data=test;
     tables   injectdays  indwater indcook indcotn indload indpass induse indnndl
        noclean RA5clean inddkwel indgall;
where injectbase=1 and mytreat=1;
run;
title 'Control';
		proc means data= test;
var DM1sex age DM1age agegrp mymarit myeduc myemploy nonwhite myrace hispanic
        RA2crack RA2cocai RA2smamp RA2hero RA2benzo
		RA2prog RA1strt RA1jail
		drinkcat moreone unpwk unpnon 
		RA3inj6m RA3injlm
		RA3hero RA3spdbl RA3heram RA3cocai RA3amph 
        injectdays  indwater indcook indcotn indload indpass induse
        noclean inddkwel;
		where mytreat=0;
		run;

proc freq data=test;
tables RA3hero RA3spdbl RA3heram RA3cocai RA3amph ;
where RA3injlm=1 and mytreat=0;
run;

proc freq data=test;
tables injectbase RA3injlm;
where mytreat=0;
run;
/*injection drug use in last month among those who reported last six months at enrollment*/

proc freq data=test;
     tables   injectdays  indwater indcook indcotn indload indpass induse indnndl
        noclean RA5clean inddkwel indgall;
where injectbase=1 and mytreat=0;
run;

