libname gwas "\\dktgenna1.niddk.nih.gov\groupcollaboration\genstatgroup\Lauren_W\birthweight\bwt_GWAS\PECRB_EGG_bwtGWAS_noGAadj_compare";
libname pecrb "\\dktgenna1.niddk.nih.gov\groupcollaboration\genstatgroup\Lauren_W\birthweight\bwt_GWAS\bwt_noGAadj_GWAS";
libname ps "\\dktgenna1.niddk.nih.gov\groupcollaboration\genstatgroup\Lauren_W\birthweight\bwt_PS";
libname pscopy "\\dktgenna1.niddk.nih.gov\groupcollaboration\genstatgroup\Lauren_W\birthweight\bwt_PS\";
libname impgric '\\dktgenna1.niddk.nih.gov\GroupCollaboration\GenStatGroup\Pimachip_imputed_GRIC';

/* SECTION 1: IMPORT BWT GWAS SUMSTATS*/

	/*own birth weight from EGG EUR*/
	PROC IMPORT OUT= EGG_EUR_sumstats
	DATAFILE= "\\dktgenna1.niddk.nih.gov\groupcollaboration\genstatgroup\Lauren_W\birthweight\bwt_LD\EGG_EUR_sumstats.txt"
	DBMS=dlm 
	replace;
	GETNAMES=YES;
	RUN;/*13891969*/
	proc contents data=EGG_EUR_sumstats; run;
	/*extract genome-wide significant SNPs*/
	data ld.EGG_EUR_sumstats_gwsig; set EGG_EUR_sumstats;
	if p < 0.00000005; run; /*9408*/
	data pscopy.EGG_EUR_sumstats_gwsig; set ld.EGG_EUR_sumstats_gwsig;
	run;
	proc contents data=pscopy.EGG_EUR_sumstats_gwsig; run;

	/*upload genome-wide sig fetal GWAS SNPs from Warrington et al (2019)*/
	PROC IMPORT OUT= EGG_EUR_fetalbwt_gwsig
	DATAFILE= "\\dktgenna1.niddk.nih.gov\groupcollaboration\genstatgroup\Lauren_W\birthweight\bwt_GWAS\PECRB_EGG_bwtGWAS_noGAadj_compare\EGG_FetalGWAS_cleaned.txt"
	DBMS=tab 
	replace;
	GETNAMES=YES;
	RUN; /*200*/
	data gwas.EGG_EUR_fetalbwt_gwsig; set EGG_EUR_fetalbwt_gwsig;
	if p < 0.00000005; run; /*193*/
	data pscopy.EGG_EUR_fetalbwt_gwsig; set gwas.EGG_EUR_fetalbwt_gwsig;
	run;
	proc contents data=pscopy.EGG_EUR_fetalbwt_gwsig; run; /*has rsIDs --> can merge by that*/

	/*PECRB noGAadj own birth weight sumstats*/
	data PECRB_sumstats; set pecrb.bwt_imp_nogaadj_anno;
	drop chr; run; /*5109994*/
	data PECRB_fetalbwt_all; set PECRB_sumstats; 
	chr=input(chromosome, 16.); 
	if chromosome="X" then chr=23;
	pos=position37;
	drop snp; run;
	data chk_PECRB; set PECRB_fetalbwt_all; if chr=23; run; /*9235*/
	proc sort data=PECRB_fetalbwt_all; by snp138; run; /*5109994*/
	data pscopy.PECRB_fetalbwt_all; set PECRB_fetalbwt_all; run;
	proc contents data=pscopy.PECRB_fetalbwt_all; run; /*rsID = varname; no 'SNP' var --> can merge with Warrington data*/

	/*get impnums by merging with anno_filtered by chrom and pos
	data imp_anno; set impgric.anno; 
	keep impnum nchrom position37;
	run; /*10446452
	data PECRB_fetalbwt_mergeimp; set PECRB_fetalbwt_all;
	proc freq; tables chr; run;
	proc sort data=imp_anno; by nchrom position37; 
	data PECRB_fetalbwt_mergeimp; set PECRB_fetalbwt_mergeimp;
	nchrom=.; run;
	data PECRB_fetalbwt_mergeimp; set PECRB_fetalbwt_mergeimp;
	nchrom=chr; run;
	proc sort data=PECRB_fetalbwt_mergeimp; by nchrom position37;*/

	/*check number of SNPs missing impnums
	data missing_impnum; set PECRB_fetalbwt_mergeimp;
	if impnum=.; run; /*520092
	data PECRB_fetalbwt_all_impnums; merge PECRB_fetalbwt_mergeimp (in=bwt) imp_anno (in=anno);
	by nchrom position37; if bwt; run; /*5110119
	proc sort data=PECRB_fetalbwt_all_impnums; by impnum;
	proc sort data=PECRB_fetalbwt_all; by impnum;
	data PECRB_fetalbwt_no_impnums; merge PECRB_fetalbwt_all_impnums (in=impnum) PECRB_fetalbwt_all (in=all); 
	if all and not impnum; run;

/* SECTION 2: MERGE BWT SUMSTATS WITH PECRB GENO DATA*/

	/*for all EGG EUR variants from GWAS of own bwt that have p<5E-8*/
	proc sort data=gwas.EGG_EUR_sumstats_gwsig; by chr pos; /*9408*/
	proc sort data=PECRB_fetalbwt_all_impnums; by chr pos; 
	data gwsEGG_all_PECRB_overlap; merge ld.EGG_EUR_sumstats_gwsig (in=egg) PECRB_fetalbwt_all_impnums (in=pecrb);
	by chr pos; if pecrb; run; /*371*/
	data egg_all_chk; set gwsEGG_all_PECRB_overlap; 
	keep impnum; run;

	/*for all in the fetal GWAS identified by Warrington et al (2019) paper with gwsig (p<5E-8) in fetal GWAS*/
	data gwas.EGG_EUR_fetalbwt_gwsig; set gwas.EGG_EUR_fetalbwt_gwsig; 
	snp138=snp; proc contents; run;
	proc sort data=gwas.EGG_EUR_fetalbwt_gwsig; by snp138; run; /*193*/
	proc sort data=PECRB_fetalbwt_all; by snp138; run; /*5109994*/
	data gwsWarr_fetalbwt_PECRB_overlap; merge gwas.EGG_EUR_fetalbwt_gwsig (in=warr) PECRB_fetalbwt_all (in=pecrb); 
	by snp138; if warr and pecrb; run; /*126*/
	proc sort data=gwas.EGG_EUR_fetalbwt_gwsig; by chr pos; 
	proc sort data=PECRB_fetalbwt_all; by chr pos; run; /*5077595
	data gwsWarr_fetalbwt_addPECRBmiss; merge gwas.EGG_EUR_fetalbwt_gwsig (in=warr) pecrb_missID (in=pecrb); 
	by chr pos; if warr and pecrb; run; /*0 --> none of the SNPs without rsIDs for varname variable in PECRB data merge successfully 
	with EGG EUR data from Warrington, et al (2019) */

	/*copy all data to bwt_PS*/
	data ps.gwsEGG_all_PECRB_overlap; set gwsEGG_all_PECRB_overlap;
	run; /*371*/
	data ps.gwsWarr_fetalbwt_PECRB_overlap; set gwsWarr_fetalbwt_PECRB_overlap;
	run; /*126, when merged by rsIDs --> this will be used in the following section*/
	data ps.PECRB_fetalbwt_impnums; set PECRB_fetalbwt_all_impnums;
	run; /*5109994*/
/* OUTCOME OF PART2: USE 126-SNP POLYGENIC SCORE THAT IS THE OVERLAP BETWEEN PECRB & EGG EUR GW-SIG VARIANTS FROM WARRINGTON ET AL (2019)*/


/*SECTION 3: EXTRACT IMPUTED DATA FOR EACH OF 7701 GENOTYPED INDIVIDUALS FROM PECRB*/
	data snplist; set ps.gwsWarr_fetalbwt_PECRB_overlap;
	keep impnum; run; /*126*/
	proc sort data=snplist; by impnum; 	run;

	%include '\\dktgenna1.niddk.nih.gov\GroupCollaboration\GenStatGroup\Pimachip_imputed_GRIC\pickgeno_imputed7701.sas';
	/*%include calls the program listed below the pickgeno set statement below*/

	data '\\dktgenna1.niddk.nih.gov\userhome\wedekindle\birthweight\bwt_PS\bwt_PS_impgeno.sas7bdat'; set pickgeno; 
	run; /*970326*/
	data PS.bwt_PS_impgeno; set pickgeno; 
	proc contents; run; /*970326*/

		/*proc datasets nolist; delete pickgeno picklist badlist ;
		run;
		proc sort data=snplist; by impnum;
		run;

		data filteredlist; set impgric.anno_filtered; by nchrom;
		  if first.nchrom then n_filtered=0;
		  n_filtered=n_filtered+1;
		  retain n_filtered;
		  keep impnum nchrom n_filtered;
		  run;
		data picklist badlist; merge snplist(in=pick keep= impnum) filteredlist(in=ok); by impnum;
		if pick and ok then output picklist;
		if pick and not ok then output badlist;
		run;
		proc print data=badlist; var impnum; title 'SNPs not in filtered list';
		run;

		data nsnps; set picklist end=last1; if last1 then call symput('nsnps',left(_n_));
		run;
		%put &nsnps;
		run;
		*------------------------------------------------------------------------;

		%macro readgeno(ds,n,m);
		%let stop=%eval(&n * 7701);
		%let start=%eval(&stop-7700);
		proc datasets nolist; delete geno1;
		data geno1; set impgric.geno&ds(firstobs=&start obs=&stop);
		proc append base=pickgeno data=geno1;
		run;
		%mend;

		%macro getdat; options nonotes;
		%do i=1 %to &nsnps;
		  data snp1; set picklist; if _n_ eq &i; 
		    call symputx('c',nchrom);
		    call symputx('nf',n_filtered); call symputx('mnum',impnum);
		  run;
		  %readgeno(&c,&nf,&mnum)
		%end; options notes;
		%mend;

		%getdat
		run; 

		data chkimpnums; merge snplist(in=inlist) pickgeno(in=ingeno); by impnum;
		if (not inlist or not ingeno);
		if first.impnum then put impnum '*************MISMATCH ERROR**************';
		run;
		options notes;
		run;*/

	proc contents data=ps.gwsWarr_fetalbwt_PECRB_overlap; run;

/*change annotation prefixes to signify data source*/
	data gwsWarr_fetalbwt_PECRB_overlap; set ps.gwsWarr_fetalbwt_PECRB_overlap; 
	/*adding pecrb prefixes*/
	pecrb_gt=gt;
	pecrb_a1=allele1;
	pecrb_a2=allele2;
	pecrb_f1=est_maf;
	pecrb_ac=est_mac;
	pecrb_dosage_sd=dosage_sd;
	pecrb_chi=chi;
	pecrb_n=NAv;
	pecrb_b=bSNP;
	pecrb_b_se=bSNPse;
	pecrb_p=p_SNP_;
	pecrb_varexp=varexp;
	drop chromosome allele1 allele2 gt id est_maf est_mac dosage_sd chi NAv bSNP bSNPse p_SNP_ varexp snp138 a1 a2; 
	/*adding egg prefixes*/
	egg_a1=ea;
	egg_a2=nea;
	egg_f1=eaf;
	egg_b=beta;
	egg_b_se=se;
	egg_p=p;
	egg_n=n;
	drop ea nea eaf beta se p n;
	run; /*126*/
	proc contents data=gwsWarr_fetalbwt_PECRB_overlap; run;
/* compare PECRB and EGG (Warrington 2019) alleles and flip signs for EGG (Warrington) EAF, Beta
all are relative to (a1=PECRB a1) Phoenix data*/
	data compare_alleles; set gwsWarr_fetalbwt_PECRB_overlap;
	ps_a1="";
	ps_a2="";
	ps_f1=.;
	ps_b=.;
	ps_b_se=.;
	run;
	data compare_alleles; set compare_alleles; 
	ps_a1=pecrb_a1;
	ps_a2=pecrb_a2;
	ps_f1=pecrb_f1; /*because PS and PECRB both have same a1 and a2 assignments, then ps_f1 will also be pecrb_f1*/
	/*beta*/
	if pecrb_a1=egg_a1 and pecrb_a2=egg_a2 then ps_b=egg_b;
	if pecrb_a1=egg_a2 and pecrb_a2=egg_a1 then ps_b=-1*egg_b;
	/*beta standard error*/
	ps_b_se=pecrb_b_se;
	run;


/* START SECTION 4: DEVELOPING POLYGENIC SCORES PER PERSON*/
	data ps.gwsWarr_fetalbwt_PECRB_overlap; set compare_alleles; run; /*126*/

/*cut anno (merged AM-DIAGRAM-PECRB) data set and merge with genotyped data from PECRB (all imputed, 0.01 < AF < 0.99);
ps_beta the regression coefficient for each allele in each SNP*/
	data cut_anno; set ps.gwsWarr_fetalbwt_PECRB_overlap;
	keep impnum ps_a1 ps_a2 ps_b; 
	proc contents; run; /*126*/
	proc sort data=cut_anno; by impnum; 
	proc sort data=PS.bwt_PS_impgeno; by impnum;
	data scorecalc1; merge cut_anno (in=anno) PS.bwt_PS_impgeno (in=genos);
	by impnum; if anno; run; /*970326*/
	data scorecalc2; set scorecalc1; 
	SNPscore=.; 
	add_risk=.; run;
	data scorecalc3; set scorecalc2; 
	SNPscore=ps_b*add_a1; 
	if sign(ps_b)=1 then add_risk=add_a1; 
	if sign(ps_b)=-1 then add_risk=2-add_a1; 
	abs_wt=abs(ps_b); run; /*970326*/

/*merge PS.bwt_PS_impgeno (N=7701 all variants imputed) with allele scores by impnum,
preserving nih numbers to indicate geno and score components for each person*/
	proc sort data=scorecalc3; by NIH; run;
	/*sum to generate raw and weigted scores*/
	proc means noprint sum data=scorecalc3; by NIH;
	var add_risk SNPscore ps_b;
	output out=scores sum=rawriskscore wtdriskscore sum_wt; run;
	proc means data=scores; var rawriskscore wtdriskscore sum_wt; run;
	data PS.bwt_ps_by_NIH; set scores;
	keep NIH rawriskscore wtdriskscore sum_wt; run; /*7701; sum_abs_wt is constant across all NIH: -0.035776*/
/*corr*/
	proc corr data=PS.bwt_ps_by_NIH; 
	var rawriskscore; 
	with wtdriskscore;
	run; /*0.95, p <0.001*/

/*--> some people suggest standardizing polygenic scores by sum_abs_wt*/
/*simplify */
	proc contents data=ps.gwsWarr_fetalbwt_PECRB_overlap; run;
