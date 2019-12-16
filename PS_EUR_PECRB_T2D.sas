libname impres  "LOCATION_OF_PECRB_SUMMARY_STATISTICS";
libname PS "LOCATION_OF_PS_DATA";
libname jdata "LOCATION_OF_PECRB_ANNOTATION_FILES";
libname imp "LOCATION_OF_PECRB_IMPUTED_GENOTYPES";

/*Append 'chr' prefix to allow the full genotypes file to merge with the 3 groups of PECRB GWAS summary statistics:
7,996,323 variants in total*/
data fill_t2d; set impres.impres_diab;
length chrid $15;
length pos $12;
pos=left(put(position37,12.));
	if chr eq "1" then chrid=put("chr1:"||pos, 15.);
	if chr eq "2" then chrid=put("chr2:"||pos, 15.);
	if chr eq "3" then chrid=put("chr3:"||pos, 15.);
	if chr eq "4" then chrid=put("chr4:"||pos, 15.);
	if chr eq "5" then chrid=put("chr5:"||pos, 15.);
	if chr eq "6" then chrid=put("chr6:"||pos, 15.);
	if chr eq "7" then chrid=put("chr7:"||pos, 15.);
	if chr eq "8" then chrid=put("chr8:"||pos, 15.);
	if chr eq "9" then chrid=put("chr9:"||pos, 15.);
	if chr eq "10" then chrid=put("chr10:"||pos, 15.);
	if chr eq "11" then chrid=put("chr11:"||pos, 15.);
	if chr eq "12" then chrid=put("chr12:"||pos, 15.);
	if chr eq "13" then chrid=put("chr13:"||pos, 15.);
	if chr eq "14" then chrid=put("chr14:"||pos, 15.);
	if chr eq "15" then chrid=put("chr15:"||pos, 15.);
	if chr eq "16" then chrid=put("chr16:"||pos, 15.);
	if chr eq "17" then chrid=put("chr17:"||pos, 15.);
	if chr eq "18" then chrid=put("chr18:"||pos, 15.);
	if chr eq "19" then chrid=put("chr19:"||pos, 15.);
	if chr eq "20" then chrid=put("chr20:"||pos, 15.);
	if chr eq "21" then chrid=put("chr21:"||pos, 15.);
	if chr eq "22" then chrid=put("chr22:"||pos, 15.);
	proc contents;
run; /*23988969*/
/* extract 'All' group for merging with Mahajan PS, by chrid*/
data pecrb_t2d_snps; set fill_t2d;
	if group="FullPima";
	keep chrid chr position37 allele1 allele2 gt f1_nondb f1_db /*f1_all*/ b_dose1 se_dose1 p_dose1 impnum;
run; /*7996323*/
proc sort data=pecrb_t2d_snps; by impnum; run;
data imp_anno; set jdata.anno_filtered; 
proc contents; run; /*7996323*/
data pecrb_t2d_anno; merge imp_anno (in=anno) pecrb_t2d_snps (in=t2d);
by impnum; if anno and t2d; run; /*7996323*/


/*import 274 annotated SNPs that were extracted, with effect estimates, from Mahajan 2018
this will be used for annotation data (e.g. gene, rsID, EUR sumstats)*/
proc import 
	datafile="\\dktgenna1.niddk.nih.gov\userhome\wedekindle\T2D_PS\new_t2d_snps_n307.txt"
	out=PS.t2d_mahajan_2018_gwext
	dbms=dlm
	replace;
	delimiter=",";
	getnames=yes;
	guessingrows=5000000;
run; /*274*/
data t2d_mahajan_2018_gwext; set PS.t2d_mahajan_2018_gwext;
am_a1=".";
am_a2=".";
am_OR=.;
am_p=.;
am_ulor=.;
am_logor=.;
am_se_logor=.;
am_a1=riskal;
am_a2=nriskal;
am_OR=or;
am_p=pvalue;
am_ulor=ulor;
am_logor=logor;
am_se_logor=se_logor;
drop riskal nriskal or pvalue ulor logor se_logor; 
proc contents; run; /*274*/
/* add chrid column*/
data t2d_mahajan_2018_gwext_fill; set t2d_mahajan_2018_gwext;
length chrid $15;
length pos $12;
pos=left(put(position37,12.));
	if chrom eq 1 then chrid=put("chr1:"||pos, 15.);
	if chrom eq 2 then chrid=put("chr2:"||pos, 15.);
	if chrom eq 3 then chrid=put("chr3:"||pos, 15.);
	if chrom eq 4 then chrid=put("chr4:"||pos, 15.);
	if chrom eq 5 then chrid=put("chr5:"||pos, 15.);
	if chrom eq 6 then chrid=put("chr6:"||pos, 15.);
	if chrom eq 7 then chrid=put("chr7:"||pos, 15.);
	if chrom eq 8 then chrid=put("chr8:"||pos, 15.);
	if chrom eq 9 then chrid=put("chr9:"||pos, 15.);
	if chrom eq 10 then chrid=put("chr10:"||pos, 15.);
	if chrom eq 11 then chrid=put("chr11:"||pos, 15.);
	if chrom eq 12 then chrid=put("chr12:"||pos, 15.);
	if chrom eq 13 then chrid=put("chr13:"||pos, 15.);
	if chrom eq 14 then chrid=put("chr14:"||pos, 15.);
	if chrom eq 15 then chrid=put("chr15:"||pos, 15.);
	if chrom eq 16 then chrid=put("chr16:"||pos, 15.);
	if chrom eq 17 then chrid=put("chr17:"||pos, 15.);
	if chrom eq 18 then chrid=put("chr18:"||pos, 15.);
	if chrom eq 19 then chrid=put("chr19:"||pos, 15.);
	if chrom eq 20 then chrid=put("chr20:"||pos, 15.);
	if chrom eq 21 then chrid=put("chr21:"||pos, 15.);
	if chrom eq 22 then chrid=put("chr22:"||pos, 15.);
run;

/*file sent from Anubha with 307 variants included in the latest T2D PS*/
proc import 
	datafile="\\dktgenna1.niddk.nih.gov\userhome\wedekindle\T2D_PS\t2d_variants_am.txt"
	out=PS.t2d_ps_am
	dbms=dlm
	replace;
	getnames=yes;
	guessingrows=5000000;
run; /*307*/

proc sort data=PS.t2d_ps_am; by chrid;
proc sort data=pecrb_t2d_anno; by chrid;
data am_pecrb_t2d_snps; merge pecrb_t2d_anno (in=pecrb) PS.t2d_ps_am (in=am); 
by chrid; if pecrb and am; run; /*245 --> this is the final merge between Phoenix and Anubha data*/
/*add pecrb anno data*/
data am_pecrb_t2d_snps; set am_pecrb_t2d_snps;
pecrb_a1=".";
pecrb_a2=".";
pecrb_gt=".";
pecrb_f1=.;
pecrb_a1=allele1;
pecrb_a2=allele2;
pecrb_gt=gt;
pecrb_f1=f1_all;
rsid=snp138;
drop chr allele1 allele2 gt f1_all id; run;
proc contents data=am_pecrb_t2d_snps; run; /*245*/
proc sort data=t2d_mahajan_2018_gwext_fill; by chrid;
proc sort data=am_pecrb_t2d_snps; by chrid;
data am_pecrb_t2d_snps_anno; merge t2d_mahajan_2018_gwext_fill (in=anno) am_pecrb_t2d_snps (in=pick);
by chrid; if pick; run; /*245 obs, 32 var --> annotated*/
proc sort data=am_pecrb_t2d_snps_anno; by chrom nchrom position37; run;
proc contents data=am_pecrb_t2d_snps_anno; run;
data PS.am_pecrb_t2d_snps_anno; set am_pecrb_t2d_snps_anno;
pecrb_b=b_dose1;
pecrb_b_se=se_dose1;
pecrb_p=p_dose1;
drop chrom rsid f1_FullPima f1_db f1_nondb f1_other b_dose1 se_dose1 p_dose1 a1 a2 snp pos;
proc contents;
run;
proc contents data=PS.am_pecrb_t2d_snps_anno; run;
/*28 SNPs missing annotation data for EUR data*/

/*Bring in data from Mahajan et al, Nat Genet (2018b) paper on fine-mapping T2D loci --> from DIAGRAM Consortium website: 
'Full meta-analysis of all studies to generate associations for T2D (in Europeans only) unadjusted for BMI'*/
proc import 
	datafile="\\dktgenna1.niddk.nih.gov\userhome\wedekindle\T2D_PS\Mahajan.NatGenet2018b.T2D.European\Mahajan.NatGenet2018b.T2D.European.txt"
	out=PS.t2d_eur_mahajan_sumstats
	dbms=tab
	replace;
	getnames=yes;
	guessingrows=5000000;
run; /*23465132*/
proc contents data=PS.t2d_eur_mahajan_sumstats; run;
data PS.t2d_eur_mahajan_sumstats; set PS.t2d_eur_mahajan_sumstats;
position37=pos;
nchrom=chr;
drop pos chr;
run;

data am_pecrb_dataprocessing; set PS.am_pecrb_t2d_snps_anno; run; /*245*/
proc sort data=PS.am_pecrb_t2d_snps_anno; by nchrom position37;
proc sort data=PS.t2d_eur_mahajan_sumstats; by nchrom position37;
proc sort data=am_pecrb_dataprocessing; by nchrom position37;
data all_anno; merge am_pecrb_dataprocessing (in=overlap) PS.t2d_eur_mahajan_sumstats (in=diamante);
by nchrom position37; if overlap; run;
proc contents data=all_anno; run;
data compare; set all_anno;
keep snp138 gene position37 SNP pecrb_b pecrb_b_se am_OR am_se_logor beta SE pecrb_a1 pecrb_a2 am_a1 am_a2 EA NEA pecrb_p am_p Pvalue pecrb_f1 EAF impnum;
run;
proc contents data=all_anno; run;

/*1. fill missing genes using dbsnp*/
data missinggene; set compare;
if gene = ""; keep snp138 gene position37 SNP impnum; run; /*28*/
data all_anno; set all_anno; 
if snp138="rs67156297" then gene="";
if snp138="rs9429893" then gene="SRGAP2";
if snp138="rs1437466" then gene="";
if snp138="rs66815886" then gene="ADAMTS9-AS2";
if snp138="rs8192675" then gene="SLC2A2";
if snp138="rs73069940" then gene="CTBP1";
if snp138="rs13130484" then gene="";
if snp138="rs13133548" then gene="FAM13A";
if snp138="rs6813195" then gene="";
if snp138="rs3936510" then gene="C5orf67";
if snp138="rs9470794" then gene="ZFAND3";
if snp138="rs12663159" then gene="KCNK16/KCNK17/LOC105375047";
if snp138="rs6905288" then gene="";
if snp138="rs849133" then gene="JAZF1";
if snp138="rs2366214" then gene="UBE3C";
if snp138="rs1802295" then gene="VPS26A";
if snp138="rs703980" then gene="ZMIZ1";
if snp138="rs11819995" then gene="ETS1";
if snp138="rs2261181" then gene="RPSAP52";
if snp138="rs7955901" then gene="";
if snp138="rs1790116" then gene="PITPNM2";
if snp138="rs34165267" then gene="";
if snp138="rs8008910" then gene="NRXN3";
if snp138="rs12445430" then gene="NLRC3";
if snp138="rs8071043" then gene="ZZEF1";
if snp138="rs11651755" then gene="HNF1B";
if snp138="rs6567160" then gene="";
if snp138="rs3787497" then gene="GNAS/GNAS-AS1";
run;
data PS.am_pecrb_diagram_ps_rawdata; set all_anno; run;

/*2. compare PECRB and DIAGRAM alleles and flip signs for DIAGRAM EAF, Beta
all are relative to (A1=PECRB A1) Phoenix data*/
data compare_alleles; set all_anno;
ps_a1="";
ps_a2="";
ps_f1=.;
ps_b=.;
ps_b_se=.;
run;
data compare_alleles; set compare_alleles; 
ps_a1=pecrb_a1;
ps_a2=pecrb_a2;
ps_f1=pecrb_f1;
/*beta*/
if pecrb_a1=ea and pecrb_a2=nea then ps_b=Beta;
if pecrb_a1=nea and pecrb_a2=ea then ps_b=-1*Beta;
/*beta standard error*/
ps_b_se=pecrb_b_se;
run;

/*summarize data for extracting imputed genotypes*/
data PS.impnums_245PS; set PS.am_pecrb_diagram_ps_sumstats; 
keep impnum; run;
data impgeno; set PS.impgeno_245;
proc contents; run; /*245*7701=1,886,745*/

/* END SECTION 1: FOUND OVERLAP BETWEEN MAHAJAN 2019 T2D PS, PECRB T2D SUMSTATS FROM GRIC*/


/* START SECTION 2: DEVELOPING POLYGENIC SCORES PER PERSON*/

/* read in GRS SNPs proposed by LEW 2018 meta-analysis*/
data PS.am_pecrb_diagram_ps_sumstats; set compare_alleles;
run; /*245*/

/*cut anno (merged AM-DIAGRAM-PECRB) data set and merge with genotyped data from PECRB (all imputed, 0.01 < AF < 0.99);
ps_beta is corrected for allele comparison between Mahajan (2018) and PECRB already; beta=log(OR) according to:
https://www.diagram-consortium.org/downloads/Mahajan.et.al.2018b.European.GWAS.readme.pdf*/
data cut_anno; set PS.am_pecrb_diagram_ps_sumstats;
keep impnum ps_a1 ps_a2 ps_b; run;
proc sort data=cut_anno; by impnum; 
proc sort data=PS.impgeno_245; by impnum;
data scorecalc1; merge cut_anno (in=anno) PS.impgeno_245 (in=genos);
by impnum; if anno; run; /*1886745*/
data scorecalc2; set scorecalc1; 
SNPscore=.; 
add_risk=.; run;
data scorecalc3; set scorecalc2; 
SNPscore=ps_b*add_a1; 
if sign(ps_b)=1 then add_risk=add_a1; 
if sign(ps_b)=-1 then add_risk=2-add_a1; 
abs_wt=abs(ps_b); run;

/*merge gricgeno_imputed211 (N=7701 all variants imputed) with allele scores by impnum, 
preserving nih numbers to indicate geno for each person*/
proc sort data=scorecalc3; by NIH; 
/*sum to generate raw and weighted (by log(OR) = beta) scores*/
proc means noprint sum data=scorecalc3; 
by NIH; 
var add_risk SNPscore ps_b; 
output out=scores sum=rawriskscore wtdriskscore sum_wt;
run;
proc means data=scores; var rawriskscore wtdriskscore sum_wt;
run;
data PS.t2d_ps_by_NIH; set scores;
keep NIH rawriskscore wtdriskscore sum_wt; run; /*7701; sum_abs_wt is constant across all NIH, with sum_abs_wt=14.1761
--> some people suggest standardizing polygenic scores by sum_abs_wt*/
proc contents data=PS.am_pecrb_diagram_ps_sumstats; run;
data PS.t2d_ps_pecrb_am_diagram_anno; set PS.am_pecrb_diagram_ps_sumstats;
diag_b=beta;
diag_N=Neff;
diag_a1=ea;
diag_a2=nea;
diag_f1=eaf;
diag_p=Pvalue;
diag_b=beta;
diag_b_se=SE;
pecrb_Qinfo=Qinfo;
pecrb_q5ok=q5ok;
run;
data PS.t2d_ps_pecrb_am_diagram_anno; set PS.t2d_ps_pecrb_am_diagram_anno;
keep gene chrid diag_b diag_N diag_a1 diag_a2 diag_f1 diag_p diag_b diag_b_se pecrb_Qinfo pecrb_q5ok snp chrid nchrom pecrb_a1 pecrb_a2 pecrb_b pecrb_b_se pecrb_f1 pecrb_gt pecrb_p position37 ps_b ps_b_se ps_a1 ps_a2 ps_f1 snp138;
proc contents; run;
