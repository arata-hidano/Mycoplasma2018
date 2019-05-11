 /*===========================================================
Mycoplasma project
1. Demographic model - 2018/10
2. Starting date June 1
=============================================================*/


/*===Parameters to receive from R=========*/
//int herd_size;
int risk_date;
double Se and Sp;

int nrow_testschedule = 395; 
double bv_r1 = 1/21.0;
double bv_r2 = 1/15.0;
double bv_r3 = 1/45.0;
double bv_r4 = 1/100.0;
double Se_Bio = 0.4 ;
double Sp_Bio = 0.95 ;
double Se_ID = 0.8 ;
double Sp_ID = 0.95 ;
double PCR_live = 0.9;
double PCR_slaughter = 0.95;
#define herd_size 400
//#define mastitis_recovery_date 17
/*PARAMETER TO ESTIMATE*/
//S,E,I(clinical+shed),I(clinical+noshed),I(nonclinical+shed),I(nonclinical+noshed),R
// r1: rate from E to Isub
// r2: rate from Isub to Icli
// r3: rate from Icli to Isub 
// r4: rate from Isub to E
double bv_r1_min = 1/21.0 ;//minimum r1
double bv_r1_max = 1/1.0; //maximum r1
double bv_r2_min = 1/15.0; 
double bv_r2_max = 1/4.0 ;
double bv_r3_min = 1/60.0;
double bv_r3_max = 1/30.0;
double bv_r4_min = 1/15.0;
double bv_r4_max = 1/365.0 ;


double bv_f = 0.05 ; //scale factor of subclinical shedding compared to clinical
double beta_mm = 0.02  ; //from cow to cow
double beta_mc = 0.01; //from milkers to calves
double rate_sero_conversion = 1/21; //mode 21 days 

/*==================================================================
DEFINE PARAMETER & STRUCTURE
===================================================================*/

/* C LIBRARIES TO INCLUDE */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <malloc.h>
#include <ctype.h>

/* Define variables*/
/*Variables that can change according to the start date*/
#define PSC 45 //Planned start of calving, set as 15 July which is 45 days
#define PSM 127 //Planned start of mating, set as 07 October 12 weeks from PSC
int dry_day = 364; //May 31st
int r1_initial_age = 304; //As of June 1st, they are max 304 days old

/*Other parameters*/

#define num_total_herd 1
#define id_calf_group 0
#define id_calf_weaned 1
#define id_r1_group 2
#define id_r2_group 3
#define id_lact_group 4
#define id_dry_group 5
#define id_sick_group 6
int const id_bull_group = 7;
#define calf_female_prop 0.5
#define calf_keep_prop 0.3
#define replacement_prop 0.22
#define calv_3weeks 0.6 //around 60%
#define calv_6weeks 0.87 //can be lower between 73% to 88%
#define calv_9weeks 0.98
#define submission_prop 0.9
#define conception_AI 0.48
#define conception_bull 0.55
#define mating_week_AI 6 //can be 5 weeks accroding to Melvin
#define mating_week_bull 4 //can be 7 weeks according to Melvin
#define mating_period 10
#define scanning 272
#define test_date1 100
#define test_date2 200

/*MILK related parameters*/
#define withhold_colostrum 5 //days
#define calf_milk_consumption 4 //litres
#define av_milk_production 15 //litre per day
#define bobby_pickup 4 //every 4 day bobby gets picked up
#define time_first_heat_min 10  // days until the first oestrus minimum value
#define time_first_heat_max 49  // days until the first Oestrus maximum value
#define interval_heat_min 18  // interval between oestrus events minimum
#define interval_heat_max 24  // interval between oestrus events maximum
#define av_gestation 282
#define error_gestation 10
#define zero 0
#define heifer_puberty 365 // reaches puberty and starts heat
#define heifer_puberty_error 30 //margin for reaching puberty
//#define calf_mortality 0.041
//#define R1_mortality 0.017
//#define R2_mortality 0.017
//#define mixed_mortality 0.017
#define weaning_wks 10 
#define sim_years 10
#define num_extra_animal_per_year 120

#define column_cull_empty 1
#define column_sell_empty 2
/*Define tables imported*/
int num_cull_sell_steps = 62 ; //62 steps for culling and death
int num_mortality_steps = 19 ;
int n_column_cull_sell_rate = 5; // days, cull_rate_empty, sell_rate_empty,cull_rate,sell_rate
int n_column_mortality_rate = 2; // day and mortality

int num_column_NumGrpAnimal = 15000;
int temp_num_animal;
// and sum of rates in each mng, but later will include no.sus/exposed/infectious/latent
char CullSellRatesFile[] = "E:/ARATA/Documents/Research/Massey/Postdc_Massey/MBOVIS/Model/CullSellRatesFile2.csv" ;
char MortalityFile[] = "E:/ARATA/Documents/Research/Massey/Postdc_Massey/MBOVIS/Model/mortality.csv" ;
char NumberAnimalDataFile[] = "E:/ARATA/Documents/Research/Massey/Postdc_Massey/MBOVIS/Model/NumberAnimalDataFile.csv" ;
char TestScheduleFile[] = 'E:/ARATA/Documents/Research/Massey/Postdc_Massey/MBOVIS/Model/TestSchedule.csv';
char TestResultFile[] = 'E:/ARATA/Documents/Research/Massey/Postdc_Massey/MBOVIS/Model/TestResultFile.csv' ;
/*Define general parameters*/



/* Define calf related parameters*/


/* Define heifer related parameters*/ 
int first_heat_rand ;
int pregnant_status;
int target_calf_met = 0;
int target_milk_met = 0;
/* Define lactating and non-lactating related parameters*/
//int num_total_groups = num_total_herd * (id_bull_group + 1) ; // gives the total number of groups


int interval_first_heat = time_first_heat_max - time_first_heat_min + 1 ;
int interval_heat = interval_heat_max - interval_heat_min + 1 ;
int calving_duration = 70; //9 to 12 weeks, now set as 10 weeks
/*Define age distribution of the mixed-age cow*/
double prop_age_2 = 0.2 ;
double prop_age_4 = 0.335;
double prop_age_7 = 0.3 ;
double prop_age_8 = 0.165 ;

/* Define culling, mortality, sale parameters*/
//create a matrix storing age on row and time-period on column
int num_points = 26; //26 fortnights
int num_age_cat = 10 ;//0,1,2,3,4,5,6,7,8,9 and just repeat 9yrs rate

/* Define simulation related parameters*/
int sim_days = sim_years*365;
int test_freq, last_test_date;

//remtaining will give a birth between 63 (9weeks) and 70 (10weeks)

/* Define hospital group related parameters*/

/*=======DISEASE SPECIFIC PARAMETERS==========*/
int num_bovis_status = 5 ; //need to modify?
#define Sus 0
#define Exp 1
#define Sub_shed 2
#define Cli_shed 3



/*Specify column for List_mng_status*/
#define column_total 0
#define column_sum_demographic 1
#define column_mastitis_both 2
#define column_mastitis 3
#define column_treated 4
#define column_Sus 5
#define column_Exp 6
#define column_Sub_shed 7
#define column_Cli_shed 8
#define column_r_SE 9
#define column_r_EIs 10
#define column_r_IsIc 11
#define column_r_IcIs 12
#define column_r_IsE 13
#define column_sum_mbovis 14
#define column_sum_nonbovis 15
#define column_sum_detection 16

/*Specify column for Milk_numbers*/
#define column_nonshed_C 0 //N cows noshedding in colostrum, no detected
#define column_nonshed_NC 1 // N cows noshedding in non-colostrum, no detected
#define column_nonshed_D 2 // N cows nonshedding and detected CM
#define column_subshed_C 3 // N cows subshedding in colostrum, no detected
#define column_subshed_NC 4 // N cows subshedding in non-colostrum, no detected
#define column_subshed_D 5// N cows subshedding detected CM
#define column_cli_C 6// N cows clinically shedding in colostrum, non detected
#define column_cli_NC 7//N cows clinically shedding in non colostrum, non detected
#define column_cli_D 8 //N cows clinically shedding detected CM


/*Column for TestSchedule*/
int herd_id ;
#define test_date 1
#define slaughter_pcr_nomilk_calf 2
#define slaughter_elisa_calf 3
#define live_pcr_nomilk_calf 4
#define live_elisa_calf 5
#define slaughter_pcr_nomilk_adult 6
#define slaughter_elisa_adult 7
#define live_pcr_nomilk_adult 8
#define live_elisa_adult 9
#define live_pcr_milk 10
#define btank 11
#define bdiscard 12
#define elisa_type 13


#define test_pcr_live_nomilk 3
#define test_pcr_slaughter_nomilk 4
#define test_pcr_milk 5

/*Parameters associated with mastitis or other conditions being detected*/
/*This is nuisance parameter - don't want to estimate, shall I just predefine with no stochasticity?*/
#define rate_detect_clinical_mastitis 0.95
//#define rate_detect_subclinical_mastitis 0.96 //this is sensitivity of SCC test when herd test done
//#define rate_detect_lameness 0.6
//#define rate_subclinical_mastitis_5d 0.036
//#define rate_subclinical_mastitis_later 0.0015
//#define rate_clinical_mastitis 0.03 //incindece rate given already subclinical
//#define rate_lameness 0.001
#define length_CM_treatment 6 //treatment + withholding
#define rate_CM_2 0.031
#define rate_CM_3_4 0.016
#define rate_CM_5_7 0.031
#define rate_CM_8 0.047 //all rate /cow-day
///using these values, create a table of CM rate






/* STRUCTURE DECLARATIONS */  
struct animal_node {
      long long akey ; // animal id
      int age_day;   /* age in day*/
      int parity ;
      int group;         // production type, 0 = calf, 1 = R1, 2 = R2, 3 = Lact, 4 = Dary, 5 = Sick, 6 = bull 
      int sex;          /* sex, female = 1 male = 0*/
      int breed;
      int mbovis_status;    /* disease status check definition above*/
      int dried;
      //int p1_yes; //1 for animals moving to Isub, 0 for not
     // int p2_yes; //1 for animals moving to Icli, 0 for not
      int NB_mastitis_status ; //non-mastitis = 0, clinical mastitis = 1 DUE to non-bovis
      int mastitis_detected ; //0 not detected 1 detected
      double mastitis_rate ; //IR CM - 0 if not milking or when developed CM due to non-bovis
      int mbovis_sero_status ;//0 neg 1 pos (how about weak pos?)
      int mbovis_elisa_status; //0 neg 1 detected
      int mbovis_pcr_status; //0 neg 1 detected
    //  int lameness_status ; //no lame 0 lame 1 detected 2
      int pregnant_status; // non-pregnant = 0 preganant to AI = 1 pregnant to bull = 2
      int calving_date; //date of calving (Day xx)
      int non_colostrum ;//0 if it's in colostrum and 1 if it's not in colostrum
      int next_heat_date; //date of the next heat (Day xx)
      //int num_births; // record how many births it already gave to: deleted as now it has parity var
      double bovis_antibody;//MBOVIS antibody
      int present; //whether this animal exists on the farm or not
      long long current_pro_id ;
      int index_cull_sell ;
      int index_mortality ;
      double sum_markov_rate ;
      //double sub_mastitis_rate ; //incidence rate of having subclinical mastitis
      struct animal_node *previous_node ;/*pointer to the previous node*/
      struct animal_node *next_node; /* pointer to next animal*/
    
   }; 

/*==========EVENT NODE=====================================*/ 
struct event_node {
      int event_type;   
	  /*event type: 0 calf to R1 and R1 to R2 happens at the same time
	                1 calving
	                2 heat
	                3 culling/sale index change
	                4 mortality changes
					*/
      long long akey ; //animal id
      struct animal_node *animal;         
      struct event_node *next_node; /* pointer to next event*/
   }; 
   
/*=======================================================================
FUNCTION LIST
========================================================================*/
void add_animal_group() ;
void add_event_node();
void read_cull_sell_rate() ;
void read_mortality() ;
int write_number_animals() ;
double update_markov_date() ;
void remove_animal_group() ;
void visualize_list() ;
void read_testschedule();
int test() ;
int test_bulk();
int write_test_result() ;
/*Setup vector*/


/*======================================================================
Initialise calving date, age, pregnant status, disease status etc
=======================================================================*/

//int fortnight_num = floor(today_date/14) ; //0-364 days, Day 364 becomes week 26, need to make it wek 25

int mng_group ;
int age_cat ;
int current_age ;
int current_index_cull_sell ;
int current_index_mortality ;
int next_cull_change_date ;
int next_mortality_change_date ;
int current_grp ;
int* num_culled ;
int* num_sold ;
int* num_death ;
int var_cull = 0;
int var_sold = 0;
int var_death = 0;
int last_pick_up = 0;

int* num_bobby;
int var_bobby = 0;

double markov_rate ;
int index_column = 0;
int i, j, abc;
int calving_date, next_heat_date ;
int next_non_markov_date = 0;
double updated_date ;


/*=====Variable related to testing========================*/
int counter_slaughter_pcr_nomilk ;
int counter_slaughter_pcr_milk ;
int counter_slaughter_ELISA ;

int counter_live_pcr_nomilk;
int counter_live_pcr_milk ;
int counter_live_ELISA ;

int target_slaughter_calf;
int target_live_calf;
int target_slaughter_adult;
int target_live_adult;
int target_adult;
int slaughter, live, group ;

int N_PPS, N_H, counter_slaughter, counter_live; 
int counter_interval_sick, counter_interval_healthy;
int sampling_interval_sick, sampling_interval_healthy ;


int ncol_testschedule = 14;
int ncol_test_table = 10 ;//check 10 columns 
int shedder,total ;

/*===============================================================================*/
/*Main starts from here*/
int main(void){
	
int length_animal_pointer = herd_size*2 + num_extra_animal_per_year*sim_years ;	
srand((unsigned)time(NULL));
int current_akey = 0 ;
double today_date = 0 ;
int mbovis_status ;
int year = 0;
int calendar_day = 0;
int parity = 0;
int first_bobby_born = 0;
int id_lact_group2 = id_lact_group;
num_culled = &var_cull;
num_sold = &var_sold;
num_death = &var_death;
num_bobby = &var_bobby;
struct animal_node* fake_animal;
fake_animal = (struct animal_node*)malloc(sizeof( struct animal_node ));

printf("Starts");	
double** List_mng_status =  (double**)malloc( sizeof(double *) *id_bull_group); //modify if more than one herd exists	

double **cull_sell_rate = (double**)malloc( sizeof(double *) *num_cull_sell_steps );
double **mortality = (double**)malloc( sizeof(double *) *num_mortality_steps );
double **NumGrpAnimal = (double**)malloc( sizeof(double *) *(id_bull_group+2) );
double **table_CM_rate = (double**)malloc(sizeof(double*)*8);// 0-7 parity
int *milk_numbers = (double*)malloc(sizeof(double)*(column_cli_D+1)) ; //table for milk
int** GroupSampleSize = (int**)malloc(sizeof(int*)*(id_dry_group+1)) ;
//TASK_38
int* test_result_table = (int*)malloc(sizeof(int)*ncol_test_table) ;
int** TestSchedule =  (int**)malloc(sizeof(int*)*nrow_testschedule) ;
/*=======================================================================
Set up linked list for the management group
========================================================================*/
struct animal_node* FarmGroupList[id_bull_group+1];//Just a single pointer, if get confused, check F&Q why
struct animal_node* current_animal;
struct animal_node* next_animal;
struct animal_node* previous_animal;

struct event_node* event_day[sim_days];
struct event_node* new_event ;
struct event_node* next_event;
struct event_node* current_event ;
struct event_node* previous_event;

/*Set up dynamic memories*/

for(i = 0 ; i < 8; i++)
{
	table_CM_rate[i] = (double*)malloc( sizeof(double) * 2);//2 columns: [0] normal [1]under treatment or after drying off
	table_CM_rate[i][0] = 0 ;
	if(i<1)
	{
		table_CM_rate[i][1] = 0;
	}
	else if(i==1)
	{
		table_CM_rate[i][1] = rate_CM_2;
	}
	else if(i<=3)
	{
		table_CM_rate[i][1] = rate_CM_3_4;
	}
	else if(i<=6)
	{
		table_CM_rate[i][1] = rate_CM_5_7;
	}
	else
	{
		table_CM_rate[i][1] = rate_CM_8;
	}
}

for(i=0; i<=column_cli_D; i++)
{
	milk_numbers[i] = 0;
}



for(i = 0; i < num_cull_sell_steps; i++)
{
	cull_sell_rate[i] = (double*)malloc( sizeof(double) * n_column_cull_sell_rate);
	
}
for(i = 0; i< num_mortality_steps; i++)
{
	mortality[i] = (double*)malloc( sizeof(double) * n_column_mortality_rate);
}

for(i = 0; i < id_bull_group + 1 ; i++)
{
	List_mng_status[i] = (double*)malloc( sizeof(double) * (column_r_IsE+1));//num column i
	
	for(j = 0; j < column_r_IsE+1; j++)
	{
		List_mng_status[i][j] = 0;
	}
	
}
for(i = 0; i < id_bull_group + 2 ; i++)
{
	NumGrpAnimal[i] =  (double*)malloc( sizeof(double) * num_column_NumGrpAnimal);
	for(abc = 0 ;abc <num_column_NumGrpAnimal; abc++)
	{
		NumGrpAnimal[i][abc] = 0;
	}
}

for(i=0; i<=id_dry_group; i++)
{
	GroupSampleSize[i] = (int*)malloc(sizeof(int)*2) ; //number of column is 2: slaughter and alive
}
//now read the cull_sell_rate table	
printf("Start read");
read_cull_sell_rate(CullSellRatesFile,cull_sell_rate,num_cull_sell_steps) ;	
printf("A");
read_mortality(MortalityFile,mortality,num_mortality_steps) ;
//printf("File read");
double sum_age_prop_4 = prop_age_2 + prop_age_4 ;
double sum_age_prop_7 = prop_age_7 + sum_age_prop_4 ;

for(i = 0; i < sim_days; i++)
                {
                event_day[i] = NULL;
                }
//Initialise the linked list
for(i = 0; i < id_bull_group+1; i++)
      	  	{
     	      FarmGroupList[i] = NULL; // initialise the animal struct
			      	    }

for(i=0;i<nrow_testschedule;i++)
{
	TestSchedule = (int*)malloc( sizeof(int) * ncol_testschedule);
}
read_testschedule(TestScheduleFile,TestSchedule,nrow_testschedule) ;

/*==================================================================================================================*/






			      	    
/*==========FROM HERE FARM SPECIFIC========================================================================*/			      	    
herd_id = 0 ;
var_cull = 0;
var_death = 0;
var_sold = 0;
var_bobby = 0;			      	    
struct animal_node **animal_node_pointer = malloc( sizeof(struct animal_node*) * length_animal_pointer);
	   
/*======ADD TEST EVENT======================================================*/
//can I pass test table information from R?
test_freq = 0;
for(i=0;i<nrow_testschedule;i++)
{
	if(TestSchedule[i][0]==herd_id)
	{
	test_date = TestSchedule[i][1] ;//TASK_40
	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
    new_event->event_type = 12;//testing
    new_event->next_node = NULL ;
    add_event_node(event_day,test_date, new_event) ;
	test_freq++; // counter for total testing on this farm	
	
	}
	if(TestSchedule[i][0]>herd_id)
	{
		break;
	}
}
printf("Total test is %d",test_freq) ;
last_test_date = test_date ;
	   

/*======ADD RISK EVENT==========================================*/
//risk event date will be a parameter, random sample from some distribution
//TASK_40
	//risk_date =  ;
	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
    new_event->event_type = 16;//Risk events
    new_event->next_node = NULL ;
    add_event_node(event_day,risk_date, new_event) ;     	    
/*Day 0
Add R1/R2/Dry to the linked list*/
//For now, assume no calves and heifers die and farmers only keep 
//as many calves/heifers as needed (i.e. herd size * replacement rate)


/*============ R1 HEIFER========================================================*/
//adding R1 heifer
for(i=0; i< herd_size * calf_keep_prop; i++ )
{
	mng_group = id_r1_group ;
	 //R1 heifers should be between 284 and 304 days old
	current_age =  r1_initial_age - rand()%21 ;
	current_index_cull_sell = 6;
	current_index_mortality = 2 ;
	next_cull_change_date = cull_sell_rate[current_index_cull_sell+1][0] - current_age ;
	next_mortality_change_date = mortality[current_index_mortality+1][0] - current_age ;
		
	animal_node_pointer[current_akey] =malloc(sizeof(struct animal_node)) ;
	animal_node_pointer[current_akey]->akey = current_akey;
	animal_node_pointer[current_akey]->age_day = current_age ;
	animal_node_pointer[current_akey]->mastitis_rate = table_CM_rate[0][0] ;
	animal_node_pointer[current_akey]->index_cull_sell = current_index_cull_sell ;
	animal_node_pointer[current_akey]->index_mortality = current_index_mortality ;
	animal_node_pointer[current_akey]->group = mng_group ; //
	animal_node_pointer[current_akey]->pregnant_status = 0;
	animal_node_pointer[current_akey]->parity = 0 ;
	animal_node_pointer[current_akey]->dried = 1 ;
	animal_node_pointer[current_akey]->present = 1 ;
	animal_node_pointer[current_akey]->previous_node = NULL;
	animal_node_pointer[current_akey]->next_node = NULL;
	//animal_node_pointer[current_akey]->lameness_status = 0;
	animal_node_pointer[current_akey]->NB_mastitis_status = 0;
	animal_node_pointer[current_akey]->mbovis_status = Sus ;
	animal_node_pointer[current_akey]->mbovis_sero_status = 0;
	mbovis_status = animal_node_pointer[current_akey]->mbovis_status;
	
	markov_rate = 
	cull_sell_rate[current_index_cull_sell][column_cull_empty]+
	cull_sell_rate[current_index_cull_sell][column_sell_empty] + 
	mortality[current_index_mortality][1] ;
	//then mortality
	List_mng_status[mng_group][1] = List_mng_status[mng_group][1] + markov_rate ;
	
	animal_node_pointer[current_akey]->sum_markov_rate = markov_rate ;
	add_animal_group(FarmGroupList,mng_group, animal_node_pointer[current_akey]) ;
	/*============ADDING EVENTS===================================================*/
	/*Heat*/
	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
	new_event->event_type = 2 ;
	new_event->akey = current_akey;
	new_event->animal = animal_node_pointer[current_akey]; //let's see if this works
	new_event->next_node = NULL ;
	first_heat_rand = rand()%heifer_puberty_error+heifer_puberty;
	add_event_node(event_day,(first_heat_rand-current_age), new_event) ;
	/*============CHANGE CULLING iNDEX==========================*/
	//need to add an event that changes this animals culling rate
	//but not for mortality because mortality changes when they calve
	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
	new_event->event_type = 3;//change in culling index
	//  no need to specify the next mortality index change date because that happens after 
	// the animal becomes 2 yrs old
	new_event->akey = current_akey;
	new_event->animal = animal_node_pointer[current_akey]; //let's see if this works
	new_event->next_node = NULL ;
	add_event_node(event_day,next_cull_change_date, new_event) ;
	/*============CHANGE MORTALITY iNDEX==========================*/
	//need to add an event that changes this animals culling rate
	//but not for mortality because mortality changes when they calve
	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
	new_event->event_type = 4;//change in culling index
	//  no need to specify the next mortality index change date because that happens after 
	// the animal becomes 2 yrs old
	new_event->akey = current_akey;
	new_event->animal = animal_node_pointer[current_akey]; //let's see if this works
	new_event->next_node = NULL ;
	add_event_node(event_day,next_mortality_change_date, new_event) ;
	
	current_akey ++;
	List_mng_status[mng_group][0] ++;
	List_mng_status[mng_group][mbovis_status+5] ++; //sus, exp, Sub_shed, Cli_shed


}
printf("R1 added");

	
/*==================R2 HEIFER================================================================*/
for(i=0; i< herd_size * calf_keep_prop; i++ )
{
	mng_group = id_r2_group ;
	current_age = 365+r1_initial_age-rand()%21;
	//R2 heifers should be between 649 and 669 days old
//	if(current_age <= 700)
//	{
		current_index_cull_sell = 12;	
//	}
//	else
//	{
	//	current_index_cull_sell = 13;
//	}
	next_cull_change_date = cull_sell_rate[current_index_cull_sell+1][0] - current_age ;
	current_index_mortality = 2 ;
	next_mortality_change_date = mortality[current_index_mortality+1][0] - current_age ;
	
	animal_node_pointer[current_akey] =malloc(sizeof(struct animal_node)) ;
	animal_node_pointer[current_akey]->akey = current_akey;
	
	
	animal_node_pointer[current_akey]->age_day = current_age ;
	animal_node_pointer[current_akey]->mastitis_rate = table_CM_rate[0][0] ;
	animal_node_pointer[current_akey]->index_cull_sell = current_index_cull_sell ;
	animal_node_pointer[current_akey]->index_mortality = current_index_mortality ;
	animal_node_pointer[current_akey]->group = mng_group ; //
	pregnant_status = 1;
	animal_node_pointer[current_akey]->pregnant_status = pregnant_status;
	animal_node_pointer[current_akey]->parity = 0 ;
	animal_node_pointer[current_akey]->dried = 1 ;
	animal_node_pointer[current_akey]->present = 1 ;
	animal_node_pointer[current_akey]->previous_node = NULL;
	animal_node_pointer[current_akey]->next_node = NULL;
	//animal_node_pointer[current_akey]->lameness_status = 0;
	animal_node_pointer[current_akey]->NB_mastitis_status = 0;
	animal_node_pointer[current_akey]->mbovis_status = Sus ;
	animal_node_pointer[current_akey]->mbovis_sero_status = 0;
	
	mbovis_status = animal_node_pointer[current_akey]->mbovis_status;
	
	add_animal_group(FarmGroupList,mng_group, animal_node_pointer[current_akey]) ;
	
	
	markov_rate = 
	cull_sell_rate[current_index_cull_sell][column_cull_empty+2*pregnant_status]+
	cull_sell_rate[current_index_cull_sell][column_sell_empty+2*pregnant_status] + 
	mortality[current_index_mortality][1] ;
	
	
	List_mng_status[mng_group][1] = List_mng_status[mng_group][1] + markov_rate ;
	
	animal_node_pointer[current_akey]->sum_markov_rate = markov_rate ;
	
	if(i<herd_size * replacement_prop *calv_3weeks)
	{
		calving_date =  (int)today_date + rand()%21 + PSC ; //note today_date = 0
//		next_heat_date = calving_date + rand()%interval_first_heat + time_first_heat_min + rand()%interval_heat + interval_heat_min;
	}
	else if(i<herd_size * replacement_prop *calv_6weeks)
	{
		calving_date = (int)today_date + rand()%21 + 21 + PSC ;
//		next_heat_date = calving_date + rand()%interval_first_heat + time_first_heat_min + rand()%interval_heat + interval_heat_min;	
	}
	else if(i<herd_size * replacement_prop *calv_9weeks)
	{
		calving_date = (int)today_date + rand()%21 + 42 + PSC ;
//		next_heat_date = calving_date + rand()%interval_first_heat + time_first_heat_min + rand()%interval_heat + interval_heat_min;
	}
	else
	{
		calving_date = (int)today_date + rand()%7 + 63 + PSC ;
//		next_heat_date = calving_date + rand()%interval_first_heat + time_first_heat_min + rand()%interval_heat + interval_heat_min;
	}
	animal_node_pointer[current_akey]->calving_date = calving_date;
//	animal_node_pointer[current_akey]->next_heat_date = next_heat_date ;
	
	/*========== ADDING EVENTS========================================================*/
	//create calving event: calving, heat and mortality and culling index change
	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
	new_event->event_type = 1;//calving
	new_event->akey = current_akey;
	new_event->animal = animal_node_pointer[current_akey]; //let's see if this works
	new_event->next_node = NULL ;
	add_event_node(event_day,calving_date, new_event) ;
	//create heat event
//	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
//	new_event->event_type = 2;//heat
//	new_event->akey = current_akey;
//	new_event->animal = animal_node_pointer[current_akey]; //let's see if this works
//	new_event->next_node = NULL ;
//	add_event_node(event_day,next_heat_date, new_event) ;
	//create event that indicates when culling/mortality index should change
	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
	new_event->event_type = 3;//change in culling index
	new_event->akey = current_akey;
	new_event->animal = animal_node_pointer[current_akey]; //let's see if this works
	new_event->next_node = NULL ;
	add_event_node(event_day,next_cull_change_date, new_event) ;
	//mortality
	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
	new_event->event_type = 4;//change in mortality index
	new_event->akey = current_akey;
	new_event->animal = animal_node_pointer[current_akey]; //let's see if this works
	new_event->next_node = NULL ;
	add_event_node(event_day,next_mortality_change_date, new_event) ;
	
	//adding event done
	current_akey ++;
	List_mng_status[mng_group][0] ++;
	List_mng_status[mng_group][mbovis_status+5] ++;
}
printf("R2 added");
/*=====================MIXED AGE DRY COWS=============================================*/

int age_indicator ;
double temp_prop ;
int dry_calved = 0;
int dry_culled =0;
for(i=0; i< herd_size * (1-replacement_prop); i++ )
{
	mng_group = id_dry_group ;//dry group when started
	temp_prop = (double)rand()/(double)RAND_MAX ; //get a random value between 0 and 1
	
	animal_node_pointer[current_akey] =malloc(sizeof(struct animal_node)) ;
	animal_node_pointer[current_akey]->group = mng_group ;
	animal_node_pointer[current_akey]->akey = current_akey;
	pregnant_status = 1 ;
	animal_node_pointer[current_akey]->pregnant_status = pregnant_status;
	animal_node_pointer[current_akey]->dried = 1 ;
	animal_node_pointer[current_akey]->present = 1 ;
	animal_node_pointer[current_akey]->previous_node = NULL;
	animal_node_pointer[current_akey]->next_node = NULL;
	//animal_node_pointer[current_akey]->lameness_status = 0;
	animal_node_pointer[current_akey]->NB_mastitis_status = 0;
	animal_node_pointer[current_akey]->mbovis_status = Sus ;
	animal_node_pointer[current_akey]->mbovis_sero_status = 0;
	mbovis_status = animal_node_pointer[current_akey]->mbovis_status;
	
	add_animal_group(FarmGroupList,mng_group, animal_node_pointer[current_akey]) ;
//printf("Animal added to group") ;
		/*ASSIGN A RANDOM AGE*/
		if(temp_prop<prop_age_2)
		{
		current_age = r1_initial_age+365*2-rand()%21;
		animal_node_pointer[current_akey]->age_day =  current_age ;
		//age between 1014-1034
		if(current_age<=1021)
		{
			current_index_cull_sell = 17;
			current_index_mortality = 4;
		}
		else
		{
		current_index_cull_sell = 18 ;
		current_index_mortality = 4 ;	
		}
		
		age_indicator = 2 ;
		age_cat = 3 ;
		//printf("Age 2");
		}
		/*AGE 3/4*/
		else if(temp_prop<sum_age_prop_4)
		{
		age_indicator = rand()%2+3 ;//3 or 4
		current_age = r1_initial_age+365*age_indicator-rand()%21;//1379 and 1764
		animal_node_pointer[current_akey]->age_day = current_age; 
		
		age_cat = 4 ; 
		if(current_age<=1385)
		{
		current_index_cull_sell = 23 ;
		current_index_mortality = 6 ;
		}
		else if(current_age<=1427)
		{
		current_index_cull_sell = 24 ;
		current_index_mortality = 6 ;	
		}
		else if(current_age<=1773)
		{
		current_index_cull_sell = 25 ;
		current_index_mortality = 7 ;
		}
		else if(current_age<=1791)
		{
		current_index_cull_sell = 30 ;
		current_index_mortality = 8 ;
		}
		else
		{
		current_index_cull_sell = 31 ;
		current_index_mortality = 9 ;
		}
	//	printf("Age 4");
		}
		/*AGE 5/6/7*/
		else if(temp_prop<sum_age_prop_7)
		{
		age_indicator = rand()%3+5 ;//2109-2859
		age_cat = 5 ;
		current_age = r1_initial_age+365*age_indicator-rand()%21;
		animal_node_pointer[current_akey]->age_day =  current_age ;
		if(current_age<=2113)
		{
		current_index_cull_sell = 35 ;
		current_index_mortality = 10 ;	
		}
		else if(current_age<=2155)
		{
		current_index_cull_sell = 36 ;
		current_index_mortality = 10 ;	
		}
		else if(current_age<=2158)
		{
		current_index_cull_sell = 37 ;
		current_index_mortality = 11 ;
		}
		else if(current_age<=2519)
		{
		current_index_cull_sell = 42 ;
		current_index_mortality = 12 ;
		}
		else if(current_age<=2523)
		{
		current_index_cull_sell = 43 ;
		current_index_mortality = 13 ;
		}
		else if(current_age<=2883)
		{
		current_index_cull_sell = 48 ;
		current_index_mortality = 14 ;
		}
		else
		{
		current_index_cull_sell = 49 ;
		current_index_mortality = 15 ;
		}
	//	printf("Age 7");
		}
		/*AGE >=8 but now limited to 8*/
		else
		{
		age_indicator = 8 ;
		age_cat = 6 ;
		current_age = r1_initial_age+365*age_indicator-rand()%21; //3204-3224
		animal_node_pointer[current_akey]->age_day = current_age ;
		if(current_age<=3205)
		{
		current_index_cull_sell = 53 ;
		current_index_mortality = 16 ;
		}
		if(current_age<=3247)
		{
		current_index_cull_sell = 54 ;
		current_index_mortality = 16 ;
		}
		else
		{
		current_index_cull_sell = 56 ;
		current_index_mortality = 17 ;
		}
		printf("Age 8");	
		}
	//	printf("cull index start\n");
		next_cull_change_date = cull_sell_rate[current_index_cull_sell+1][0] - current_age ;
	//	printf("cull index set\n");
		if(current_index_mortality<18)
		{
		next_mortality_change_date = mortality[current_index_mortality+1][0] - current_age ;
	//	printf("mortality < 18") ;	
	//	printf("mortality date is %d",next_mortality_change_date) ;
	//	printf("current age is %d",current_age);
	//	printf("current index is %d",current_index_mortality);
		}
		else
		{
		next_mortality_change_date = mortality[16][0]+365 - current_age ;	
	//	printf("mortality over 18") ;
		}
		animal_node_pointer[current_akey]->index_cull_sell = current_index_cull_sell ;
		animal_node_pointer[current_akey]->index_mortality = current_index_mortality ;
	//0 calf, 1 R1, 2 R2, 3 2yr, 4 3-5yr, 5 6-7yr, 6 8yr
	markov_rate = 
	cull_sell_rate[current_index_cull_sell][column_cull_empty+pregnant_status*2]+
	cull_sell_rate[current_index_cull_sell][column_sell_empty+pregnant_status*2] + 
	mortality[current_index_mortality][1] ;
	animal_node_pointer[current_akey]->sum_markov_rate = markov_rate;
	
	List_mng_status[mng_group][1] = List_mng_status[mng_group][1] + markov_rate ;
	List_mng_status[mng_group][0] ++;
	List_mng_status[mng_group][mbovis_status+5] ++;
	
	//printf("Number is %d",(int)List_mng_status[mng_group][0]);
	/*REPRODUCTION PARAMETERS*/
	animal_node_pointer[current_akey]->parity = age_indicator - 1 ;// assuming cows have claved every year
	animal_node_pointer[current_akey]->mastitis_rate = table_CM_rate[age_indicator - 1][0] ;//max 7


	if(i<herd_size * (1-replacement_prop) *calv_3weeks)
	{
		calving_date =  (int)today_date + rand()%21 + PSC ;
	//	next_heat_date = calving_date + rand()%interval_first_heat + time_first_heat_min + rand()%interval_heat + interval_heat_min;
	}
	else if(i<herd_size * (1-replacement_prop) *calv_6weeks)
	{
		calving_date = (int)today_date + rand()%21 + 21 + PSC ;
	//	next_heat_date = calving_date + rand()%interval_first_heat + time_first_heat_min + rand()%interval_heat + interval_heat_min;	
	}
	else if(i<herd_size * (1-replacement_prop) *calv_9weeks)
	{
		calving_date = (int)today_date + rand()%21 + 42 + PSC ;
	//	next_heat_date = calving_date + rand()%interval_first_heat + time_first_heat_min + rand()%interval_heat + interval_heat_min;
	}
	else
	{
		calving_date = (int)today_date + rand()%7 + 63 + PSC ;
	//	next_heat_date = calving_date + rand()%interval_first_heat + time_first_heat_min + rand()%interval_heat + interval_heat_min;
	}
	animal_node_pointer[current_akey]->calving_date = calving_date;
	//animal_node_pointer[current_akey]->next_heat_date = next_heat_date ;
	//printf("Adding event");
	/*ADD EVENTS*/
	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
	
	new_event->event_type = 1;//calving
	new_event->akey = current_akey;
	//printf("address is %i akey is %lld\n",&new_event->akey,new_event->akey);
	new_event->animal = animal_node_pointer[current_akey]; //let's see if this works
	new_event->next_node = NULL ;
	add_event_node(event_day,calving_date, new_event) ;
	
//	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
//	new_event->event_type = 2;//heat
//	new_event->akey = current_akey;
	//printf("address is %i akey is %lld\n\n",&new_event->akey,new_event->akey);
	//system("pause");
//	new_event->animal = animal_node_pointer[current_akey]; //let's see if this works
//	new_event->next_node = NULL ;
//	add_event_node(event_day,next_heat_date, new_event) ;
	
		//create event that indicates when culling/mortality index should change
	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
	new_event->event_type = 3;//change in culling index
	new_event->akey = current_akey;
	new_event->animal = animal_node_pointer[current_akey]; //let's see if this works
	new_event->next_node = NULL ;
	add_event_node(event_day,next_cull_change_date, new_event) ;
	//mortality
	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
	new_event->event_type = 4;//change in mortality index
	new_event->akey = current_akey;
	new_event->animal = animal_node_pointer[current_akey]; //let's see if this works
	new_event->next_node = NULL ;
	add_event_node(event_day,next_mortality_change_date, new_event) ;	
	current_akey ++;
}
/* ADDING INITIAL ANIMALS DONE*/
//printf("Dry is %d",(int)List_mng_status[id_dry_group][0]);
//system("pause");

/*ADD Scanning and Dry-off event*/

for(i=0;i<sim_years;i++)
{
	//SCANNING
	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
	new_event->event_type = 7;//scanning
	new_event->next_node = NULL ;
	new_event->animal = fake_animal;
	add_event_node(event_day,scanning+365*i, new_event) ;
	//DRY
	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
	new_event->event_type = 9;//dry off day
	new_event->next_node = NULL ;
	new_event->animal = fake_animal;
	add_event_node(event_day,dry_day+365*i, new_event) ;
	
	
}

/*Create fake animal that moving event can point to*/

int num_pregnant = 0;
int num_submission = 0;
int num_heat = 0 ;
int num_add_heat = 0;
int dry_season = 1 ;
/*============================================================================================
Day procedes
=============================================================================================*/
while(today_date<sim_days || today_date<last_test_date+1)
{
	year = (int)floor(today_date/365);

	
	printf("today is %lf\n",today_date) ;
	printf("non-markov is %d",next_non_markov_date);
	//visualize_list(event_day,0) ;
/*====START OF DAY========================================================*/
if(next_non_markov_date<sim_days-1)
{
next_non_markov_date = ceil(today_date);	


//	if(next_non_markov_date==year*365)
//	{printf("R1 is %d\n",(int)List_mng_status[1][0]);
//	printf("R2 is %d\n",(int)List_mng_status[2][0]);
//	system("pause");
//	}
  while (event_day[next_non_markov_date] == NULL)
            {//loop2
		    next_non_markov_date++; // go to the next day
		    if(next_non_markov_date==sim_days-1)
		    {
		    	printf("non markov reached the last day\n");
		    	break;
			}
	 //   printf("Next non Markov is %d",next_non_markov_date) ;
	 //here I need to cancel events that are associated with animals that are removed
	 if(event_day[next_non_markov_date]!=NULL)
	 {
	 //	printf("NonM is %d",next_non_markov_date);
	 	current_event = event_day[next_non_markov_date];
	 	//trial
	 	current_animal = current_event->animal; 	
	 	previous_event = current_event;
	 	while(current_event->animal->present==0)
	 	{
	 	//	printf("This event does not exist anymore!");
	 	//	system("pause");
	 		current_event = current_event->next_node;
	 	//	if(current_event->event_type==2)
	 	//	{
	 	//		printf("Heat cancelled 1");
	 	//		system("pause");
		//	 }
	 		free(previous_event);
	 		previous_event = current_event;
	 		if(current_event==NULL)
	 		{
	 		//	printf("current event becomes null");
	 			break;
	 			
			 }
		 }
	 	event_day[next_non_markov_date] = current_event ;
	 //	printf("Head updated");
	while(current_event!=NULL)
	{
	//	printf("Type is %d",current_event->event_type) ;
	if(current_event->next_node == NULL)
	{
	//	printf("C");
		break;
	}
	while(current_event->next_node->animal->present==0)
		{
		//	if(current_event->next_node->event_type==2)
	 	//	{
	 			//printf("Heat cancelled 2");
	 			//system("pause");
			// }
			if(current_event->next_node->next_node!=NULL)
			{
			//	printf("A1\n");
				next_event = current_event->next_node->next_node;
			//	printf("Animal removed is %lld",current_event->next_node->animal->akey);
				free(current_event->next_node);
				//printf("This type is %d",current_event->next_node->event_type);
				current_event->next_node = next_event ;
			//	printf("A2\n");
			//	system("pause");
			}
			else
			{
			//	printf("B1\n");
				free(current_event->next_node);
				current_event->next_node = NULL;
			//	printf("B2\n");
				break;
			}
			
		}
		current_event = current_event->next_node ;
	}
	//printf("C2"); 	
	 	
	 }
            } // now get next non-markov date
    } //if(next_non_markov_date<sim_days-1) ends
  //  printf("NonM is %d",next_non_markov_date);
 // if(next_non_markov_date>=sim_days)
 // {
 // 	printf("Break");
 // 	break;
 // }
 // printf("updating markov date");
 printf("Before update");
  updated_date=update_markov_date(today_date,List_mng_status,cull_sell_rate,
  mortality, FarmGroupList,next_non_markov_date,
  num_culled,num_sold,num_death) ;
printf("After update");
/*==============MARKOV DID NOT HAPPEN=====================================================================*/
  if (updated_date==next_non_markov_date) // this means markov event did not happen
     {//LOOP NM1
     calendar_day = next_non_markov_date - year*365 ;
	//printf("Loop started\n");
     current_event = event_day[next_non_markov_date];
     while(current_event!=NULL)
     {
     	//printf("event is not null\n") ;
     //	printf("Event is %d", current_event->event_type) ;
     	/*
    EVENT 1: Calving 2:Heat 3:Change in culling index
        4: Change in mortality index
        5: Culling
        6: Sale
        7: Death
        8: Move to the next age grp
        9:
		when do we recalculate rate again?
		The end of the day! 
     	*/
/*======== CALVING==========================================*/
     	if(current_event->event_type==1  && current_event->animal->present==1)
     	{
     		
		    current_animal = current_event->animal ;
		    if(current_animal->pregnant_status==0)
		    {
		    	printf("Wait!She is not pregnant!\n");
		    	system("pause");
			}
     		current_grp = current_animal->group ;
     		
     		/*UPDATE PARITY*/
     		parity = current_animal->parity ;
     		if(parity<7)
     		{
     			parity++;
     			current_animal->parity = parity; //if less than 7 then update, otherwise remian the same
			 }
     		current_animal->mastitis_rate = table_CM_rate[parity][1] ;
     		current_animal->dried= 0;
     		
     	/*UPDATE colostrum status and add event to switch colostrum status*/
     	//UPDATE List_mng_status bovis numbers
     	//whether this animal remains in the herd or culled, need to minus from the original group
     	mbovis_status = current_animal->mbovis_status ;
     	List_mng_status[current_grp][mbovis_status+5]--;
     	List_mng_status[current_grp][0]--;
     	//whetehr adding to the new group depends on whether this animal remains so check below
     	//UPDATING milk_number table			
     		
     		
     //		printf("event is %d\n",current_event->event_type);
     		//if(current_grp==4)
     		//{
     		//	printf("Dry cow is %d\n",(int)List_mng_status[id_dry_group][0]);
     		//	dry_calved++;
			// }
			 if(current_animal->present == 0)
			 {
			 	printf("wait! Ghost cow calving!");
			 	printf("akey is %lld and grp is %d",current_animal->akey, current_animal->group);
			 		system("pause");
			 
			 }
			
     		//this animal anyway move into Lactation group[3]
     		
     		//if heifer -1 from heifer, if adult -1 from dry
     
     	
     	//		if(current_grp==4)
     	//	{
     	//		printf("Dry cow is %d\n",(int)List_mng_status[id_dry_group][0]);
     	//	}
     	//	printf("Dry calved cow is %d\n",dry_calved);
     		if(List_mng_status[current_grp][0]<0)
     		{
     			printf("Warning! Number gets below 0 in group %d",current_grp);
     			system("pause");
			 }
     		
     	//Now change this animal's group and related fields
		 //need to change total rate of demographic right?
		 //TASK_32	- no need this is done everytime when I update event rate
     	//remove R2 if there are 400 calved
		List_mng_status[id_lact_group][0]++;//whgether they are heifer or adult, goes to lact
     	List_mng_status[id_lact_group][mbovis_status+5]++;
     		if(mbovis_status<=1)
     		{
			 bovis_index=0 ;
			 }
			 else
			 {
			 	bovis_index = mbovis_status - 1 ;
			 }	
     		
     		if(current_animal->mastitis_detected==1)
     		{
     			printf("wait, this animal just calved so CM shouldn't exist");
     			system("pause");
			 }
			current_animal->non_colostrum = 0;//now in colostrum period
			milk_numbers[bovis_index*3+current_animal->non_colostrum]++; 
     	
     		
			if(next_non_markov_date+withhold_colostrum<sim_days)
			{
				new_event = (struct event_node*)malloc(sizeof( struct event_node ));
			    new_event->event_type = 13;//switching colostrum
			    new_event->akey = current_animal->akey;
			    new_event->animal = current_animal; //let's see if this works
			    new_event->next_node = NULL ;
				add_event_node(event_day,next_non_markov_date+withhold_colostrum, new_event) ;
			}
			
		 	current_animal->pregnant_status = 0;
     		current_animal->sum_markov_rate = 
					cull_sell_rate[current_animal->index_cull_sell][column_cull_empty+2*pregnant_status]+
					cull_sell_rate[current_animal->index_cull_sell][column_sell_empty+2*pregnant_status] + 
					mortality[current_animal->index_mortality][1] ;
     	
		 	remove_animal_group(FarmGroupList,current_grp,current_animal);
     		current_animal->group = id_lact_group ; //change the group
     		add_animal_group(FarmGroupList,id_lact_group,current_animal);
     		
     	if(target_milk_met ==0)
     		{
     			//heat will happne only if these animals are included for mating
     		next_heat_date = next_non_markov_date+
     		rand()%interval_first_heat + time_first_heat_min + rand()%interval_heat + interval_heat_min ;
     		new_event = (struct event_node*)malloc(sizeof( struct event_node ));
			new_event->event_type = 2;//heat
			new_event->akey = current_animal->akey;
			new_event->animal = current_animal; //let's see if this works
			new_event->next_node = NULL ;
			if(next_heat_date<sim_days)
			{
				num_add_heat++;
				add_event_node(event_day,next_heat_date, new_event) ;
			}	
			 }
     		/*ADd First heat event*/
     	
			
			
     	if(List_mng_status[0][0]>herd_size * calf_keep_prop)
     	{
     	(num_bobby*) ++; 
     	//late calvers still keep on calving (note no more R2 exists at this stage
     		
		 }
		 else if((double)rand()/(double)RAND_MAX<=calf_female_prop) //proportion of female and male 0.5
		 {
		 // Keep at farm new born calves
     	animal_node_pointer[current_akey] =malloc(sizeof(struct animal_node)) ;
		animal_node_pointer[current_akey]->akey = current_akey;	
		animal_node_pointer[current_akey]->age_day = 0; //just born
		animal_node_pointer[current_akey]->parity = 0;
		animal_node_pointer[current_akey]->mastitis_detected = 0;
		animal_node_pointer[current_akey]->dried = 1;
		animal_node_pointer[current_akey]->mastitis_rate = table_CM_rate[0][0];
		animal_node_pointer[current_akey]->NB_mastitis_status = 0;
		animal_node_pointer[current_akey]->mbovis_status = 0;
		animal_node_pointer[current_akey]->index_cull_sell = 0;
		animal_node_pointer[current_akey]->index_mortality = 0;
		animal_node_pointer[current_akey]->pregnant_status = 0;
		animal_node_pointer[current_akey]->present = 1 ;
		animal_node_pointer[current_akey]->group = id_calf_group;
		animal_node_pointer[current_akey]->next_node = NULL;
		animal_node_pointer[current_akey]->previous_node = NULL;
		
		next_cull_change_date = next_non_markov_date+cull_sell_rate[1][0] ;
		next_mortality_change_date = next_non_markov_date+mortality[1][0] ;	
		List_mng_status[id_calf_group][0]++;
		List_mng_status[id_calf_group][column_Sus]++;
		
		//Add weaning event
		new_event = (struct event_node*)malloc(sizeof( struct event_node ));
		new_event->event_type = 14; //moving to the next group
		new_event->next_node = NULL ;
		new_event->animal = animal_node_pointer[current_akey];
			if(next_non_markov_date+weaning_wks*7<sim_years*365)
		{
		add_event_node(event_day,next_non_markov_date+weaning_wks*7, new_event) ;//animals move to the next age group
		
		}
		
		if(List_mng_status[0][0]==1)
		{
			printf("The first calving date is %d in YEAR %d\n",calendar_day, year);
			dry_season = 0;
			system("pause");
		
		}
	//	printf("Calf size is %lf\n",List_mng_status[0][0]);
		//if this is the last animal to be kept
		if(List_mng_status[0][0]==herd_size * calf_keep_prop && target_calf_met==0)
		{
		target_calf_met = 1;
		printf("Target calf secured %d in YEAR %d\n",calendar_day, year);
	//	printf("Calendar is %d", calendar_day);
		system("pause")	;
		new_event = (struct event_node*)malloc(sizeof( struct event_node ));
		new_event->event_type = 8; //moving to the next group
		new_event->next_node = NULL ;
		new_event->animal = fake_animal;
	//	printf("Lact number is %d",(int)List_mng_status[id_lact_group][0]);
	//	system("pause");
	//	printf("setting up move date") ;
		if(next_non_markov_date+weaning_wks*7<sim_years*365)
		{
			add_event_node(event_day,next_non_markov_date+weaning_wks*7, new_event) ;//animals move to the next age group
		printf("Moving event is %d",next_non_markov_date+weaning_wks*7);
	
		}
	//	printf("last calf to add\n");
		}
		markov_rate = 
		cull_sell_rate[current_index_cull_sell][1]+
		cull_sell_rate[current_index_cull_sell][2] + 
		mortality[current_index_mortality][1] ;
		animal_node_pointer[current_akey]->sum_markov_rate = markov_rate ;
		//List_mng_status[0][1] = List_mng_status[0][1] + markov_rate ;
	
		add_animal_group(FarmGroupList,0, animal_node_pointer[current_akey]) ;
		
		/*ADD FIRST HEAT ASSUMING 12MO*/
		new_event = (struct event_node*)malloc(sizeof( struct event_node ));
		new_event->event_type = 2;
		first_heat_rand = rand()%heifer_puberty_error+heifer_puberty;
		new_event->akey = current_akey;
		new_event->animal = animal_node_pointer[current_akey]; //let's see if this works
		new_event->next_node = NULL ;
			if((first_heat_rand+next_non_markov_date)<sim_years*365)
		{
			add_event_node(event_day,first_heat_rand+next_non_markov_date, new_event) ;
		}
		
		new_event = (struct event_node*)malloc(sizeof( struct event_node ));
		new_event->event_type = 3;//change in culling index
		new_event->akey = current_akey;
		new_event->animal = animal_node_pointer[current_akey]; //let's see if this works
		new_event->next_node = NULL ;
			if(next_cull_change_date<sim_years*365)
		{
			add_event_node(event_day,next_cull_change_date, new_event) ;
		}
		
		//mortality
		new_event = (struct event_node*)malloc(sizeof( struct event_node ));
		new_event->event_type = 4;//change in mortality index
		new_event->akey = current_akey;
		new_event->animal = animal_node_pointer[current_akey]; //let's see if this works
		new_event->next_node = NULL ;
			if(next_mortality_change_date<sim_years*365)
			{
			add_event_node(event_day,next_mortality_change_date, new_event) ;
		}
		
	
		current_akey++;
		 }
		 else //this is male and becomes bobby, still not enough calves born
		 {
		 	(num_bobby*)++;
		 	if(first_bobby_born==0)
		 	{
		 	//arrange calf pick up
			//TASK_31
			new_event = (struct event_node*)malloc(sizeof( struct event_node ));
		    new_event->event_type = 15; //moving to the next group
		    new_event->next_node = NULL ;
		     new_event->animal = animal_node_pointer[current_akey];
		     last_pick_up = next_non_markov_date + 77 ; //11 weeks from teh first bobby born
			if(next_non_markov_date+bobby_pickup<sim_years*365)
		        {
		        add_event_node(event_day,next_non_markov_date+bobby_pickup, new_event) ;//animals move to the next age group
		
		        }
		        first_bobby_born = 1 ;
			 }
		 }
		 //if the number of lact reaches to the target
		 //remove all R2 when the calved cows reached 400
		 //this happens only once per lactation
		 //cows can calve after this but not heifers
		 //assume all bobbies from R born on this date
		 	if(List_mng_status[id_lact_group][0]==herd_size && target_milk_met==0)
     		{
     			printf("Day is %d",calendar_day);
     			//remove all R2 if remaining
     			if(List_mng_status[id_r2_group][0]<=0)
     			{
     				printf("No R2 animals remained\n");
				 }
				 else
				 {
				 	printf("R2 is %d",(int)List_mng_status[id_r2_group][0]);
				 	current_animal = FarmGroupList[id_r2_group];
				 	while(current_animal!=NULL)
				 	{
				 		(num_bobby*)++;
				 		current_animal->present = 0;
				 		mbovis_status = current_animal->mbovis_status;
				 		previous_animal = current_animal;
				 		current_animal = current_animal->next_node;
				 		remove_animal_group(FarmGroupList,id_r2_group,previous_animal);
				 		List_mng_status[id_r2_group][0]--;
				 		List_mng_status[id_r2_group][mbovis_status+5]--;
				 		//TASK_30 - change bovis status
				 		(*num_culled)++;
				 		
				 		
					 }
				 }
				
				if(List_mng_status[id_dry_group][0]<=0)
     			{
     				printf("No Dry animals remained");
				 }
				 //else
				// {
				// 	printf("Dry is %d",(int)List_mng_status[id_dry_group][0]);
				// 	current_animal = FarmGroupList[id_dry_group];
				// 	while(current_animal!=NULL)
				// 	{
			
				// 		current_animal->present = 0;
				// 		previous_animal = current_animal;
				// 		current_animal = current_animal->next_node;
				 		//printf("now akey is %lld",current_animal->akey);
				// 		List_mng_status[id_dry_group][0]--;
				// 		(*num_culled)++;
				 	
				// 		remove_animal_group(FarmGroupList,id_dry_group,previous_animal);
				 		
				 		
				//	 }
				//	 if(FarmGroupList[id_dry_group]==NULL)
				//	 {
				//	 	printf("Now all dry animals gone\n");
				//	 }
				// }
				 target_milk_met = 1 ;
				 printf("Target milking reached\n");
			//	 system("pause");
			 }
     	
		 } //calving done
		 	
	/*=============HEAT============================================*/
		 if(current_event->event_type==2 && calendar_day<=PSM+7*mating_period&& current_event->animal->present==1)
		 {
	//	 	printf("event is %d\n",current_event->event_type);
		 	//Heat
		//printf("Calendar is %d grp is %d",calendar_day,current_event->animal->group);
		//system("pause");
			num_heat++;
			current_animal = current_event->animal;
			if(current_animal->pregnant_status==1)
			{
			printf("Heat but pregnant grp is %d day is %d",current_animal->group, next_non_markov_date);
		 	printf("Calving date is %d",current_animal->calving_date);
		 	if(current_animal->present==0)
		 	{
		 		printf("Ghost heat");
			 }
			 system("pause");
			}
		 		if(((double)rand()/(double)RAND_MAX<=submission_prop ) && 
				 current_animal->pregnant_status==0 &&
				 calendar_day>=PSM)	//submission happen
		 		{
		 			num_submission++;
		 			double conception_rate;
		 			if(next_non_markov_date-365*year-PSM <=7*mating_week_AI)
		 			{
		 				conception_rate = conception_AI ;
					 }
					 else
					 {
					 	conception_rate = conception_bull ;
					 }
		 			if((double)rand()/(double)RAND_MAX <= conception_rate)
					 {
					 	current_animal->pregnant_status = 1 ;
					 	current_animal->sum_markov_rate = 
						cull_sell_rate[current_animal->index_cull_sell][column_cull_empty+2*pregnant_status]+
						cull_sell_rate[current_animal->index_cull_sell][column_sell_empty+2*pregnant_status] + 
						mortality[current_animal->index_mortality][1] ;
						
					 //	printf("animal pregnant");
					 //	system("pause");
					 	//add calving event
					 	calving_date = av_gestation+ rand()%error_gestation*2*zero - error_gestation*zero + next_non_markov_date ;
					 	current_animal->calving_date = calving_date;
						 new_event = (struct event_node*)malloc(sizeof( struct event_node ));
						new_event->event_type = 1;//calving
						num_pregnant++;
						new_event->akey = current_event->akey;
						new_event->animal = current_event->animal; //let's see if this works
						new_event->next_node = NULL ;
								if(calving_date<sim_years*365)
						{
						add_event_node(event_day,calving_date, new_event) ;
						}
						
					}
						
				 } // then otherwise set up next heat date
				 if(current_animal->pregnant_status==0)
				 {
				 		next_heat_date = rand()%interval_heat + interval_heat_min + next_non_markov_date;
				 		new_event = (struct event_node*)malloc(sizeof( struct event_node ));
						new_event->event_type = 2;//calving
						new_event->akey = current_event->akey;
						new_event->animal = current_event->animal; //let's see if this works
						new_event->next_node = NULL ;
						if(next_heat_date<sim_years*365)
			                  {
			            add_event_node(event_day,next_heat_date, new_event) ;
		                      }
					
				 }
			 }
	/*===============CHANGE IN CULLING INDEX===========================*/
	if(current_event->event_type==3&& current_event->animal->present==1)
		 {
	//	 	printf("event is %d\n",current_event->event_type);
		 	current_animal = current_event->animal;
		  current_index_cull_sell = current_animal->index_cull_sell;	
		  pregnant_status = current_animal->pregnant_status;
		// printf("index is %d\n",current_index_cull_sell);
		 //now if it does not exceed the limit
		 if(current_index_cull_sell<=60)
		 {
		 current_index_cull_sell++;
		 }
		  else
	     {
	     	current_index_cull_sell = 57;
		 }
		current_animal->index_cull_sell = current_index_cull_sell;
		current_index_mortality = current_animal->index_mortality;
		 //then re-calculate sum again		
		current_animal->sum_markov_rate = 
		cull_sell_rate[current_index_cull_sell][column_cull_empty+2*pregnant_status]+
		cull_sell_rate[current_index_cull_sell][column_sell_empty+2*pregnant_status] + 
		mortality[current_index_mortality][1] ;
	//	printf("index updated");
	//then have to add next change date
	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
	new_event->event_type = 3;//change in culling index
	new_event->akey = current_animal->akey;
	new_event->animal = current_animal; //let's see if this works
	new_event->next_node = NULL ;
	if(current_index_cull_sell<=60)
		 {
		 next_cull_change_date = next_non_markov_date+
		 cull_sell_rate[current_index_cull_sell+1][0] - cull_sell_rate[current_index_cull_sell][0] ;
		 }
		  else
	     {
	     	next_cull_change_date = next_non_markov_date+cull_sell_rate[57][0]+365-
			cull_sell_rate[current_index_cull_sell][0]  ;
		 }
	if(next_cull_change_date<sim_days)
			{
	
				add_event_node(event_day,next_cull_change_date, new_event) ;
			}
	     }
	/*===============CHANGE IN MORTALITY INDEX===========================*/
	if(current_event->event_type==4 && current_event->animal->present==1)
		 {
	//	 printf("event is %d\n",current_event->event_type);	
		 current_animal = current_event->animal;
		 current_index_mortality =current_animal->index_mortality;
		 pregnant_status = current_animal->pregnant_status;
		 //now if it does not exceed the limit
		 if(current_index_mortality<=17)
		 {
		 	current_index_mortality++;
		 }
		  else
	     {
	     	current_index_mortality = 16;
		 }
		current_animal->index_mortality = current_index_mortality;
		current_index_cull_sell = current_animal->index_cull_sell;
		 //then re-calculate sum again		
	current_animal->sum_markov_rate = 
		cull_sell_rate[current_index_cull_sell][column_cull_empty+2*pregnant_status]+
		cull_sell_rate[current_index_cull_sell][column_sell_empty+2*pregnant_status] + 
		mortality[current_index_mortality][1] ;	
		//mortality
	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
	new_event->event_type = 4;//change in mortality index
	new_event->akey = current_animal->akey;
	new_event->animal = current_animal; //let's see if this works
	new_event->next_node = NULL ;
		if(current_index_mortality<=17)
		 {
		 //	printf("index is %d",current_index_mortality);
		 next_mortality_change_date = next_non_markov_date+
		 mortality[current_index_mortality+1][0] - cull_sell_rate[current_index_mortality][0] ;
		 }
		  else
	     {
	     //	printf("index is %d",current_index_mortality);
	     	next_mortality_change_date = next_non_markov_date+cull_sell_rate[16][0]+365-
			cull_sell_rate[current_index_mortality][0]  ;
		 }
		 if(next_mortality_change_date<sim_days)
			{
	add_event_node(event_day,next_mortality_change_date, new_event) ;
	     }
	   //  printf("Event 4 done");
	 }
	/*=========Scanning==================================================*/     
	if(current_event->event_type==7)
	{
		//scan lact cows
		if(FarmGroupList[id_lact_group]!=NULL)
		{
			current_animal = FarmGroupList[id_lact_group];
			while(current_animal!=NULL)
			{
				next_animal = current_animal->next_node;
				mbovis_status = current_animal->mbovis_status ;
				if(current_animal->pregnant_status==0)
				{
				remove_animal_group(FarmGroupList,id_lact_group,current_animal);
				current_animal->present = 0;
				List_mng_status[id_lact_group][0]--;
				List_mng_status[id_lact_group][mbovis_status+5]--;
				//this changes the colostrum milk status
				//TASK_33	
     		
     				if(mbovis_status<=1)
     				{
			 		bovis_index=0 ;
					}
					else
					{
			 			bovis_index = mbovis_status - 1 ;
					}
     				///update milk_numbers table
     				if(current_animal->mastitis_detected==0)
     				{//milk_numbers table changes only when CM not detected
					milk_numbers[bovis_index*3+current_animal->non_colostrum]--; 
					}
					else
					{
					milk_numbers[bovis_index*3+2]--; 	
					}
			 
				(*num_culled)++;
				printf("empty culled at scanning\n");	
				}
				current_animal = next_animal ;
			}
		}
		//scen R2 too
		if(FarmGroupList[id_r2_group]!=NULL)
		{
			current_animal = FarmGroupList[id_r2_group];
			while(current_animal!=NULL)
			{
				next_animal = current_animal->next_node;
				if(current_animal->pregnant_status==0)
				{
				remove_animal_group(FarmGroupList,id_r2_group,current_animal);
				current_animal->present = 0;
				List_mng_status[id_r2_group][0]--;
				List_mng_status[id_r2_group][mbovis_status+5]--;
				//TASK_33
				(*num_culled)++;	
				}
				current_animal = next_animal ;
			}
		}
	}



	/*===========Move to next subgroup========================================*/
	if(current_event->event_type==8)
	{
		
	//	printf("Moving Now R2 is %d",(int)List_mng_status[id_r2_group][0]);
	//		printf("Now Lact is %d",(int)List_mng_status[id_lact_group][0]);
	//		system("pause");
	//	printf("event is %d\n",current_event->event_type);
		if(FarmGroupList[id_r2_group]!=NULL)
		{
			current_animal = FarmGroupList[id_r2_group];
		//	printf("There are still animals in R2! Day is %d\n",next_non_markov_date) ;
		//	system("pause");
			//these animals will be culled
			
			while(current_animal!=NULL)
			{
				next_animal = current_animal->next_node;
		//		if(current_animal->pregnant_status==0)
		//		{
		//			printf("This heifer empty\n");
		//		}
				remove_animal_group(FarmGroupList,id_r2_group,current_animal);
				current_animal->present = 0;
				List_mng_status[id_r2_group][0]--;
				(*num_culled)++;
				current_animal = next_animal ;
			}
		//	printf("Now R2 is %d",(int)List_mng_status[id_r2_group][0]);
	//		system("pause");
			FarmGroupList[id_r2_group] = NULL;
			
		}
	//	if(FarmGroupList[1]==NULL)
	//	{
	//		printf("There is no animals in R1!\n");
	//	}
		if(FarmGroupList[id_r2_group]==NULL) //now if R2 becomes 0
		{//
		
			current_animal = FarmGroupList[id_r1_group] ;
			while(current_animal!=NULL)
			{
				current_animal->group = id_r2_group;
				current_animal = current_animal->next_node;
			}
			current_animal = FarmGroupList[id_r1_group] ;//again get the first animal
			FarmGroupList[id_r2_group] = current_animal ;
			FarmGroupList[id_r1_group] = NULL;
			temp_num_animal = List_mng_status[id_r1_group][0]; //R1 heifers
			List_mng_status[id_r2_group][0] = temp_num_animal;
	//		printf("R2 heifer is %d",temp_num_animal);
			
			if(FarmGroupList[id_calf_weaned]==NULL)
			{
				printf("There is no animals in calf:YEAR %d\n",year) ;
				system("pause");
			}
			else
			{
				if(List_mng_status[id_calf_weaned][0]<=0)
				{
					printf("FarmGroupList has calves but not in List_mng");
					system("pause");
				}
				else
				{
					
					temp_num_animal = List_mng_status[id_calf_weaned][0];
					List_mng_status[id_r1_group][0] = temp_num_animal;
					printf("R1 heifer is %d",temp_num_animal);
					List_mng_status[id_calf_weaned][0]=0; //all calves gone
					current_animal = FarmGroupList[id_calf_weaned] ;
					while(current_animal!=NULL)
					{
					current_animal->group = id_r1_group;
					current_animal = current_animal->next_node;
					}
					current_animal = FarmGroupList[id_calf_weaned] ;//get the first animal again
					FarmGroupList[id_r1_group] = current_animal ;
					FarmGroupList[id_calf_weaned] = NULL ;
				}
				
			}
			//transfer List_mng_status information other than the total number of animals
			for(i=column_Sus; i<= column_Exp; i++)
			{
				List_mng_status[id_r2_group][i] = List_mng_status[id_r1_group][i]; //r1 to r2
				List_mng_status[id_r1_group][i] = List_mng_status[id_calf_weaned][i]; //r1 to r2
				List_mng_status[id_calf_weaned][i] = 0;
			}
		}
	}
	/*================Drying off=====================================*/
		if(current_event->event_type==9)
	{//cows go from 3 to 4 
	target_calf_met=0;
	target_milk_met=0;
	first_bobby_born = 0 ;
	dry_season = 1 ;
	int num_lact = 0;
	int num_dry = 0;
	int num_culled_nonpregnant = 0;
	int num_pregnant_at_dry = 0;
	for(i=0;i<=column_cli_D; i++)
	{
		milk_numbers[i] = 0 ; //reset to 0, if think about between dried animals transmission need to change this
	}
	
	//printf("pregnant is %d submission is %d heat is %d adding heat is %d",num_pregnant,num_submission,num_heat,num_add_heat);
	//system("pause");
	//	printf("event is %d\n",current_event->event_type);
		if(FarmGroupList[id_lact_group]==NULL)
		{
			printf("There are no milking animals!\n") ;
		}
		else
		{
			current_animal = FarmGroupList[id_lact_group] ;
			if(FarmGroupList[id_dry_group]!=NULL)
			{
				printf("There are already dry animals!\n") ;
			//	printf("num is %d",List_mng_status[id_dry_group][0]);
				//next_animal = FarmGroupList[id_dry_group];
				//while(next_animal!=NULL)
				//{
				//	printf("Akey and grp is %lld, %d status %d\n",next_animal->akey,next_animal->group,next_animal->pregnant_status);
				//	printf("Present is %d",next_animal->present);
				//	next_animal = next_animal->next_node;
				//}
			system("pause");
			}
			else
			{
				while(current_animal!=NULL)
				{
					current_animal->group = id_dry_group;
					//here cull animals if not-pregnant
					next_animal = current_animal->next_node;
					if(current_animal->pregnant_status==0)//if not pregnant
					{
						//Assumption_4
					//	printf("akey is %lld",current_animal->akey);
						remove_animal_group(FarmGroupList,id_lact_group,current_animal);
						current_animal->present = 0 ;
						(*num_culled)++;
						num_culled_nonpregnant++;
					//	printf("Non-pregnant culled");
					//	system("pause");
					//-1 this group because we need number of remaining lact animals
						List_mng_status[id_lact_group][0]--;
					}
					
					current_animal = next_animal ;
					num_lact++;
					parity = current_animal->parity ;
					//Assumption_5
					current_animal->mastitis_rate = table_CM_rate[parity][0] ;//no mastitis in dryoff
					current_animal->dried= 1 ;
				}
				current_animal = FarmGroupList[id_lact_group] ;//get the first animal again
				FarmGroupList[id_dry_group] = current_animal;
				while(current_animal!=NULL)
				{//what am I doin this? ah maybe I was checking if num_lact and num_dry being the same
					num_dry++;
					current_animal = current_animal->next_node ;
				}
			//	printf("Lact was %d, Dry is %d Cull is %d total cull is %d\n",num_lact, num_dry,num_culled_nonpregnant,*num_culled);
			//	printf("Death is %d Sale is %d\n",*num_death,*num_sold);
			//	system("pause");
				FarmGroupList[id_lact_group]=NULL;
				temp_num_animal = List_mng_status[id_lact_group][0];
				List_mng_status[id_lact_group][0]=0;
				List_mng_status[id_lact_group][1]=0;
				List_mng_status[id_dry_group][0]=temp_num_animal;
				
				for(i = column_Sus; i<= column_Cli_shed; i++)
				{
					List_mng_status[id_dry_group][i] = List_mng_status[id_lact_group][i];
					List_mng_status[id_lact_group][i] = 0;
				}
					if(temp_num_animal<0)
     					{
     						printf("Warning! Number in Dry gets below 0");
     						system("pause");
						 }
			}
		}	
	}
/*==========EVENT 10: HERD TEST===========================================*/
/*TASK_10 TOP
	if(current_event->event_type==10) //on herd test, (sub)clinical mstitis may be detected and treatment starts
	{
		//some subclinical and all clinical mstitis will be detected with a certain probability
		if(FarmGroupList[id_lact_group]==NULL)
		{
			printf("There are no milking animals!\n") ;
		}
		else
		{
			current_animal = FarmGroupList[id_lact_group] ;
			while(current_animal!=NULL)
			{
				if(current_animal->mastitis_detected==0)
				{
					if(current_animal->mastitis_status==1||current_animal->mbovis_status==Sub_shed)//subclinical
					{
						double random_value = (double)rand()/(double)RAND_MAX;
						if(random_value < rate_detect_subclinical_mastitis) //assuming that m bovis and others lift SCC similarly
						{
						current_animal->mastitis_detected	= 1; // now detected as (sub)mastitis
						}
					}
					else if(current_animal->mastitis_status==2||current_animal->mbovis_status==Cli_shed) //clinical will be detected
						{
						current_animal->mastitis_detected =1 ;
						}
					//once detected, they are isolated and milk goes to waste milk
				//but if it's m bovis they don't get cured by the normal treatment
					if(current_animal->mastitis_detected==1) //once mastitis detected, animals get separated or their milk go to waste milk
					{
					//add treatment event
						new_event = (struct event_node*)malloc(sizeof( struct event_node ));
						new_event->event_type = 11;//treatment
						new_event->akey = current_animal->akey;
						new_event->animal = current_animal; //let's see if this works
						new_event->next_node = NULL ;
						add_event_node(event_day,next_non_markov_date+mastitis_recovery_date, new_event) ;
					}	
				}
			
			current_animal = current_animal->next_node ;
			}
		
		}
	} TASK_10 DOWN*/
/*===============EVENT 11: Finishing treatment===========================*/

	if(current_event->event_type==11) //non-m-bovis mastitis will be treated
	{//TASK_14 completed
		//some subclinical and all clinical mstitis will be detected with a certain probability
		current_animal = current_event->animal ; //animal that got treated
		if(FarmGroupList[id_lact_group]==NULL)
		{
			printf("There are no milking animals!\n") ;
		}
		else if(current_animal->dried==1)
		{
		//if the animal is already dried off
		//set deetcted_status = 0? 
		//Assumption_7: CM due to non-bovis cured as soon as dried-off
		current_animal->mastitis_detected = 0;
		current_animal->NB_mastitis_status=0; //CM due to non-bovis cured
		current_animal->mastitis_rate = 0;
		
		//And no additional treatment for CM due to BOVIS
		}
		else
		{
		
			if(current_animal->NB_mastitis_status==0) //if CM was sue to bovis
			{
			//this means CM was due to BOVIS
			///then nothing changes, animals stay in this group, will receive treatment again
			//Assumption_3: assume animals under treatment doesn't get new non-bovis infection
			//milk_numbers table doesn't change because detection status remains the same
			new_event = (struct event_node*)malloc(sizeof( struct event_node ));
			new_event->event_type = 11;//treatment
			new_event->akey = current_animal->akey;
			new_event->animal = current_animal; //let's see if this works
			new_event->next_node = NULL ;
			add_event_node(event_day,next_non_markov_date+length_CM_treatment, new_event) ;
			
			
	    	}
	    	else 
			{
			//this means CM was either non-BOVIS or both non-BOVIS+BOVIS
				current_animal->NB_mastitis_status=0; //CM due to non-bovis cured	
			if(current_animal->mbovis_status!=3)
	    	{//if CM was only due to non-bovis
	    		mbovis_status = current_animal->mbovis_status;
				if(mbovis_status<=1)
				{
					bovis_index = 0;	
				}
				else
				{
				bovis_index = mbovis_status -1 ;
				}
	    		
	    		//TASK_19
	    		current_animal->mastitis_detected = 0; //not detected status anymore
	    		//back to milking
	    		parity = current_animal->parity;
	    		//Assumption_6: back to having non-zero CM rate
	    		current_animal->mastitis_rate = table_CM_rate[parity][1] ;
	    	
	    		milk_numbers[bovis_index*3+2]--;//used to be detected
				milk_numbers[bovis_index*3+current_animal->non_colostrum]++;//back to stat
			}//if still has CM due to bovis, then remained detected and more treatment
			else
				//CM both bovis and non-bovis
			{//milk_numbers table doesn't change because detection status remains the same
				new_event = (struct event_node*)malloc(sizeof( struct event_node ));
				new_event->event_type = 11;//treatment
				new_event->akey = current_animal->akey;
				new_event->animal = current_animal; //let's see if this works
				new_event->next_node = NULL ;
				add_event_node(event_day,next_non_markov_date+length_CM_treatment, new_event) ;	
			}
			
			}
		
			
			
		}
	}
	
	
/*===============EVENT 12: SURVEILLANCE TESTING: USE REAL DATA===============*/
//TASK_16

if(current_event->event_type==12) //testing for MBOVIS 
{


//test bulk or dicard
if(test_schedule[k][btank]==1)
{
	shedder = milk_numbers[4]+milk_numbers[7];
	total = milk_numbers[1]+shedder ;
	test_bulk(0,shedder,total,dilution,PCR_live, test_result_table);
}
if(test_schedule[k][bdisc]==1)
{
	shedder = milk_numbers[3]+milk_numbers[5]+milk_numbers[6]+milk_numbers[8];//colostrum or detected
	total = milk_numbers[0]+milk_numbers[2]+shedder ;
	test_bulk(1,shedder,total,dilution,PCR_live, test_result_table);
}
/*======Determine the observed sample size for slaughter and live samples.
Get the largest sample size for each sample type (PCR -nomilk or PCR - milk OR ELISA=====*/
//========CALF======================================================
group = 0;	
ELISA_type = test_schedule[k][elisa_type] ;
//target_live_calf and target_slaughter_calf is the largest sample size of all sample types
if(test_schedule[k][slaughter_pcr_nomilk_calf+4*group]>=test_schedule[k][slaughter_pcr_nomilk_calf+1+4*group])
{
	target_slaughter_calf = test_schedule[k][slaughter_pcr_nomilk_calf+4*group] ;
}
else
{
	target_slaughter_calf = test_schedule[k][slaughter_pcr_nomilk_calf+1+4*group];
}

//target_live_calf
if(test_schedule[k][live_pcr_nomilk_calf+4*group]>=test_schedule[k][live_pcr_nomilk_calf+1+4*group])
{
	target_live_calf = test_schedule[k][live_pcr_nomilk_calf+4*group];
}
else
{
	target_live_calf = test_schedule[k][live_pcr_nomilk_calf+1+4*group] ;
}

//========================Adult==================================================*/
group = 1 ;
if(test_schedule[k][slaughter_pcr_nomilk_calf+4*group]>=test_schedule[k][slaughter_pcr_nomilk_calf+1+4*group])
{
	target_slaughter_adult = test_schedule[k][slaughter_pcr_nomilk_calf+4*group] ;//never worry that it has calf in column - just specifying the first column
}
else
{
	target_slaughter_adult = test_schedule[k][slaughter_pcr_nomilk_calf+1+4*group];
}


//target_live_adult
if(test_schedule[k][live_pcr_nomilk_calf+4*group]>=test_schedule[k][live_pcr_nomilk_calf+1+4*group] && 
test_schedule[k][live_pcr_nomilk_calf+4*group]>=test_schedule[k][live_pcr_nomilk_calf+2+4*group])
{
	target_live_adult = test_schedule[k][live_pcr_nomilk_calf+4*group] ;
}
else if(test_schedule[k][live_pcr_nomilk_calf+1+4*group]>=test_schedule[k][live_pcr_nomilk_calf+4*group] && 
test_schedule[k][live_pcr_nomilk_calf+1+4*group]>=test_schedule[k][live_pcr_nomilk_calf+2+4*group])
{
	target_live_adult = test_schedule[k][live_pcr_nomilk_calf+1+4*group] ;
}
else
{
	target_live_adult = test_schedule[k][live_pcr_nomilk_calf+2+4*group] ;
}
//===================================ADULT ENDS

/*=========Determine how many samples to collect for slaughter and live in each management group======*/
target_adult = target_live_adult + target_slaughter_adult;
slaughter = 0;
live = 1 ;

//I prioritise slaughter sample
//so collect samples from calf, wean and R1, R2 as slaughter
//once slaughter sample size is met then finally start collecting live
//Decide how many samples to choose from each management group
for(i=0;i<id_bull_group;i++)
{
	GroupSampleSize[i][slaughter] = 0;
	GroupSampleSize[i][live] = 0;
	if(target_slaughter_calf +target_live_calf >0 && i <= id_r2_group) //if need calf sampling
	{
		if(target_slaughter_calf<=List_mng_status[i][0])
		{
			GroupSampleSize[i][slaughter] = target_slaughter_calf; //all calf samples are from id_calf_group
			target_slaughter_calf = 0;
			if(target_live_calf <= List_mng_status[i][0] - target_slaughter_calf )
			{
			GroupSampleSize[i][live] = target_live_calf;
			target_live_calf = 0;
			}
			else
			{
			GroupSampleSize[i][live] =	List_mng_status[i][0]- target_slaughter_calf;
			target_live_calf = target_live_calf - GroupSampleSize[i][live];
			}
		
		}
		else
		{
			GroupSampleSize[i][slaughter] = List_mng_status[i][0];
			target_slaughter_calf = target_slaughter_calf - List_mng_status[i][0] ;
			GroupSampleSize[i][live] = 0 ;
			//if require_calf_slaughter is <0 then all slaughter sample is collected but not live
			//if >0 then  List_mng_status[i][0] were all used for slaughter and not for live
		}
	}
	else if(i>=id_lact_group && target_adult>0)
	{
		if(FarmGroupList[i]!=NULL)
		{
			if(target_slaughter_adult+taget_live_adult<=List_mng_status[i][0])
			{
				GroupSampleSize[i][slaughter] = target_slaughter_adult ;
				GroupSampleSize[i][live] = taget_live_adult ;
			}
			else if(target_slaughter_adult <= List_mng_status[i][0])
			{
				GroupSampleSize[i][slaughter] =target_slaughter_adult;
				GroupSampleSize[i][live] = List_mng_status[i][0] - target_slaughter_adult;
			}
			else
			{
				GroupSampleSize[i][slaughter] = List_mng_status[i][0];
				GroupSampleSize[i][live] = 0 ;
			}
			target_adult = 0;
		}
	}

}
/*=========ENDS: Determine how many samples to collect for slaughter and live in each management group======*/

/*=======Do adult and calf sampling as a loop=====================================================*/
for(i=0;i<=Id_dry_group;i++)
{
	//TASK_37
N_PPS = 0;
N_H = 0;
counter_slaughter = 0;//this counts the number of sample collected for slaughter in a given mangement grp
counter_live = 0;//same for live samples
counter_interval_sick = 0;
counter_interval_healthy = 0;

if(GroupSampleSize[i][slaughter]+GroupSampleSize[i][live]>0)
{//if sampling from this group

	if(i==id_calf_group)
	{
	group = 0 ;
	counter_slaughter_pcr_nomilk = test_schedule[k][slaughter_pcr_nomilk_calf] ;
	counter_slaughter_ELISA = test_schedule[k][slaughter_elisa_calf] ;	
	
	counter_live_pcr_nomilk = test_schedule[k][live_pcr_nomilk_calf];
	counter_live_ELISA = test_schedule[k][live_elisa_calf] ;
	}
	else if(i>=id_lact_group)
	{
	group = 1 ;
	counter_slaughter_pcr_nomilk = test_schedule[k][slaughter_pcr_nomilk_adult] ;
	counter_slaughter_ELISA = test_schedule[k][slaughter_elisa_adult] ;
	
	counter_live_pcr_nomilk = test_schedule[k][live_pcr_nomilk_adult];
	counter_live_pcr_milk = test_schedule[k][live_pcr_milk];
	counter_live_ELISA = test_schedule[k][live_elisa_adult] ;
	}
//Count hoe many sick animals exist in this group
current_animal = FarmGroupList[i];
	while(current_animal!=NULL)
	{
		if(current_animal->mbovis_elisa_status==1||current_animal->mbovis_pcr_status==1||current_animal->mastitis_detected==1)
		{
		N_PPS++;	
		}
		current_animal=current_animal->next_node;
	}
	N_H = List_mng_status[i][0] - N_PPS ;

/*=======Two scenarios of sampling: =======*/


//Scenario 1: slaughter + live < sick: collecting all samples from sick animals
if(GroupSampleSize[i][slaughter]+ GroupSampleSize[i][live]<= N_PPS)//if sample size from slaughter is less than N_PPS, take all from this group
{
	
	//Determine sampling interval
	sampling_interval_sick = floor((GroupSampleSize[i][slaughter]+ GroupSampleSize[i][live])/N_PPS);
	if(sampling_interval_sick<1)
	{
	sampling_interval_sick = 1 ;
	}
current_animal = FarmGroupList[i];
while(current_animal!=NULL)
{
	//sampling sick animals for both slaughter and live
	if(current_animal->mbovis_elisa_status==1||current_animal->mbovis_pcr_status==1||current_animal->mastitis_detected==1)
	{
	counter_interval_sick++;
	if(counter_interval_sick<sampling_interval_sick)
		{
		current_animal = current_animal->next_node;
		}	
	else
		{
		counter_interval_sick = 0;
		if(counter_slaughter<GroupSampleSize[i][slaughter])
		{
			//still slaughter samples
				if(counter_slaughter_pcr_nomilk>0)
				{
				counter_slaughter_pcr_nomilk--;
				test(current_animal,test_pcr_slaughter_nomilk,slaughter,group,test_result_table) ;
				current_animal->present = 0;
				}
			
				if(counter_slaughter_ELISA>0 && ELISA_type<=2)
				{
				counter_slaughter_ELISA--;
				test(current_animal,ELISA_type,slaughter,group,test_result_table) ;
				current_animal->present = 0;
				}
		if(current_animal->present = 0) // means this animal was slaughtered
		{
			N_PPS--;//how many PPS remaining
			mbovis_status = current_animal->mbovis_status;
			List_mng_status[i][0] --;
			List_mng_status[i][mbovis_status+5] --;
					if(dry_season==0 && group ==1)
					{
						if(mbovis_status<=1)
						{
						mbovis_status = 0;	
						}
						else
						{
						mbovis_status = mbovis_status - 1 ;
						}
					milk_numbers[mbovis_status*3+current_animal->non_colostrum]--;
					}
		animal_to_remove = current_animal;
		current_animal = current_animal->next_node;
		remove_animal_group(FarmGroupList, i, animal_to_remove);	
		counter_slaughter++;		
		}		
				
		}//slaughter sampling done
		else if(counter_live<GroupSampleSize[i][live])
		{
				if(counter_live_pcr_nomilk>0)
				{
				counter_live_pcr_nomilk--;
				test(current_animal,test_pcr_live_nomilk,live,group,test_result_table) ;
				}
				if(counter_live_pcr_milk>0)
				{
				counter_live_pcr_milk--;
				test(current_animal,test_pcr_milk,live,group,test_result_table) ;
				}
				if(counter_live_ELISA>0 && ELISA_type<=2)
				{
				counter_live_ELISA--;
				test(current_animal,ELISA_type,live,group,test_result_table) ;
				}
				counter_live++;
				curent_animal = current_animal->next_node;
		}
		}
	}
	else
	{
		current_animal = current_animal->next_node;
	}
	if(counter_slaughter==GroupSampleSize[i][slaughter] && counter_live==GroupSampleSize[i][live])
	{
		break;
	}
	///while loop ends when enough sample collected or reached to the end
}
}//scneario 1 ends: when all samples can be taken from sick animals

//Scenario 2: slaughter sample is larger than N_PPS
//Need to sample from N_H for slaughter
else if(GroupSampleSize[i][slaughter] <= N_PPS + N_H)
{
//take all sick for slaughter
current_animal = FarmGroupList[i];
while(counter_slaughter < GroupSampleSize[i][slaughter]||counter_live<GroupSampleSize[i][live])
{
	if(current_animal->mbovis_elisa_status==1||current_animal->mbovis_pcr_status==1||current_animal->mastitis_detected==1)
	{
				if(counter_slaughter_pcr_nomilk>0)
				{
				counter_slaughter_pcr_nomilk--;
				test(current_animal,test_pcr_slaughter_nomilk,slaughter,group,test_result_table) ;
				current_animal->present = 0;
				}
				if(counter_slaughter_ELISA>0 && ELISA_type<=2)
				{
				counter_slaughter_ELISA--;
				test(current_animal,ELISA_type,slaughter,group,test_result_table) ;
				current_animal->present = 0;
				}
				//what do they do when animals found positive?
		if(current_animal->present = 0) // means this animal was slaughtered
		{
			N_PPS--;//how many PPS remaining
			mbovis_status = current_animal->mbovis_status;
			List_mng_status[i][0] --;
			List_mng_status[i][mbovis_status+5] --;
					if(dry_season==0 && group ==1)
					{
						if(mbovis_status<=1)
						{
						mbovis_status = 0;	
						}
						else
						{
						mbovis_status = mbovis_status - 1 ;
						}
					milk_numbers[mbovis_status*3+current_animal->non_colostrum]--;
					}
		animal_to_remove = current_animal;
		current_animal = current_animal->next_node;
		remove_animal_group(FarmGroupList, i, animal_to_remove);	
		counter_slaughter++;		
		}	
	}//all sick animals are slaughtered
	else
	{
		current_animal=current_animal->next_node;
	}
//when finishing collecting all sick animals
//then start sampling from healthy for slaughter
//can collect slaughter and live sample at the same time
if(N_PPS==0)
{
current_animal = FarmGroupList[i];//reset
counter_interval_healthy = 0;
sampling_interval_healthy = floor((GroupSampleSize[i][slaughter]+ GroupSampleSize[i][live]-counter_slaughter)/N_H);
if(sampling_interval_healthy<1)
{
sampling_interval_healthy = 1 ;
}

while(current_animal!=NULL)
{
//First collect remaining slaughter samples
//Then collect live samples later
counter_interval_healthy++;
if(counter_interval_healthy<sampling_interval_healthy)
{
	current_animal=current_animal->next_node;
}
else
{
	counter_interval_healthy = 0; //reset
	if(counter_slaughter<GroupSampleSize[i][slaughter])
	{
	//still slaughtering
				if(counter_slaughter_pcr_nomilk>0)
				{
				counter_slaughter_pcr_nomilk--;
				test(current_animal,test_pcr_slaughter_nomilk,slaughter,group,test_result_table) ;
				current_animal->present = 0;
				}
				if(counter_slaughter_ELISA>0 && ELISA_type<=2)
				{
				counter_slaughter_ELISA--;
				test(current_animal,ELISA_type,slaughter,group,test_result_table) ;
				current_animal->present = 0;
				}
		if(current_animal->present = 0) // means this animal was slaughtered
		{
			N_H--;//how many PPS remaining
			mbovis_status = current_animal->mbovis_status;
			List_mng_status[i][0] --;
			List_mng_status[i][mbovis_status+5] --;
					if(dry_season==0 && group ==1)
					{
						if(mbovis_status<=1)
						{
						mbovis_status = 0;	
						}
						else
						{
						mbovis_status = mbovis_status - 1 ;
						}
					milk_numbers[mbovis_status*3+current_animal->non_colostrum]--;
					}
		animal_to_remove = current_animal;
		current_animal = current_animal->next_node;
		remove_animal_group(FarmGroupList, i, animal_to_remove);	
		counter_slaughter++;		
		}
	}
	else if(counter_live<GroupSampleSize[i][live])
	{//live sample
				if(counter_live_pcr_nomilk>0)
				{
				counter_live_pcr_nomilk--;
				test(current_animal,test_pcr_live_nomilk,live,group,test_result_table) ;
				}
				if(counter_live_pcr_milk>0)
				{
				counter_live_pcr_milk--;
				test(current_animal,test_pcr_milk,live,group,test_result_table) ;
				}
				if(counter_live_ELISA>0 && ELISA_type<=2)
				{
				counter_live_ELISA--;
				test(current_animal,ELISA_type,live,group,test_result_table) ;
				}
				counter_live++;
				curent_animal = current_animal->next_node;
	}//live sample done
}

if(counter_slaughter==GroupSampleSize[i][slaughter] && counter_live==GroupSampleSize[i][live])
	{
		break;
	}		
}//while current_animal!=NULL Done
} //If N_PPS = 0 DONE
if(counter_slaughter==GroupSampleSize[i][slaughter] && counter_live==GroupSampleSize[i][live])
	{
		break;
	}	
}//All sampling done


} //else if Scenario 2 ends

} //sampling from this management group done

} //loop for going over different management group done
}
/*========EVENT 12: Sampling done=====================================*/

/*===============EVENT 13: Colostrum status change===============*/
	if(current_event->event_type==13)
	{
		current_animal = current_event->animal ;
		mbovis_status = current_animal->mbovis_status ;
     		
     		if(mbovis_status<=1)
     		{
			 bovis_index=0 ;
			}
			else
			{
			 	bovis_index = mbovis_status - 1 ;
			}
     		///update milk_numbers table
     		if(current_animal->mastitis_detected==0)
     		{//milk_numbers table changes only when CM not detected
			milk_numbers[bovis_index*3+current_animal->non_colostrum]--; 
			current_animal->non_colostrum = 1 ;
			milk_numbers[bovis_index*3+current_animal->non_colostrum]++;
			 }
				 
	}
/*===============EVENT 13: Colostrum status change ENDS===============*/

/*================EVENT 14: Weaning====================================*/
//TASK_30
	if(current_event->event_type==14)
	{
		current_animal = current_event->animal ;
		current_grp = current_animal->group ;
		mbovis_status = current_animal->mbovis_status;
		//remove this animal to the next group
		if(current_grp != id_calf_group)
		{
			printf("This is not calf!");
			system("pause");
		}
		else
		{
			
			List_mng_status[id_calf_group][0] --;
			List_mng_status[id_calf_weaned][0] ++;
			List_mng_status[id_calf_group][mbovis_status+column_Sus] --;
			List_mng_status[id_calf_weaned][mbovis_status+column_Sus] ++;
			//TASK_30
			///change the numebr of mbovis depending on this animal status
			remove_animal_group(FarmGroupList,current_grp,current_animal);
     		current_animal->group = id_calf_weaned ; //change the group
     		add_animal_group(FarmGroupList,id_calf_weaned,current_animal);
		}
     		
     		
				 
	}
/*===============EVENT 14: Weaning ends=================================*/

/*================EVENT 15: Bobby pick up==================================*/
//TASK_31
if(current_event->event_type==15)
	{
		if((num_bobby*)==0)
		{
			printf("No bobby exists");
		}
		else
		{
			(num_bobby*) = 0 ;
		}
		if(next_non_markov_date+bobby_pickup<sim_days && next_non_markov_date+bobby_pickup < last_pick_up)
		{
		new_event = (struct event_node*)malloc(sizeof( struct event_node ));
				new_event->event_type = 15;//treatment
				new_event->akey = current_animal->akey;
				new_event->animal = current_animal; //let's see if this works
				new_event->next_node = NULL ;
				add_event_node(event_day,next_non_markov_date+bobby_pickup, new_event) ;	
		}
		
	}
/*================EVENT 15: Bobby pick up ends==============================*/
	
/*=================EVENT 16: Risk Event========================*/
if(current_event->event_type==16)
   {//TASK_41
	//because now I don't know what's the actual risk animal
	//just assume one of milking animal is infected
	//if in the dry season, then assume one of dry adult
	if(dry_season==0)//milking season
	{
	group = id_lact_group ;
    }
    else
    {
    group = id_dry_group ;	
	}
	current_animal = FarmGroupList[group];
	random_number = rand%()List_mng_status[group][0] + 1;
	while(random_number>0)
	{
	random_number--;
	    if(random_number==0)
	    {
		current_animal->mbovis_status=1 ;
		break;
		}
	    else
	    {
	    current_animal= current_animal->next_node;	
	    }	
	}
	
   }

/*=================EVENT 16: Risk Event done====================*/
	
	/*============================================================================*/
     	//MOVE to the next event
     	struct event_node *previous_event;
		previous_event = current_event ; //rewire to the next event
	   current_event = current_event->next_node;
	   if (current_event!=NULL)
	   {
	 //  	printf("going to free event now") ;
	   
	   //printf("next event is %d, %lld", current_event->event_type, current_event->akey);
       }
       else
       {
      //	printf("next event is NULL\n") ;
	   }
	  // printf("now free event");
	   free(previous_event);///check if this works - previously this was in if(current_event!=NULL) blacket
	  // printf("event was %d",previous_event->event_type);
	  // system("pause");
	   event_day[next_non_markov_date] = current_event;
	   //printf("event pointing new event\n") ;
     }
    // printf("day is over\n");
    for(i=0;i<id_bull_group+1;i++)
	{
		//printf("number is %lf",List_mng_status[i][0]);
	//	printf("index column is %d",index_column) ;
	NumGrpAnimal[0][index_column] = today_date ;
	NumGrpAnimal[i+1][index_column] = (int)List_mng_status[i][0] ;
	}
	index_column++;	
	printf("index column is %d\n",index_column);
//	printf("index updated %d",index_column) ;
 	} //LOOP NM1
	today_date = updated_date;
	if((int)today_date==sim_days-1)
	{
		//when it reaches the final day of simulation
		printf("Reaches the end of the simulation") ;
		break;
	}
	//printf("Now today date is %lf\n",today_date) ;
	
	
}
printf("Loop done");
/*===========SIMULATION ENDS=============================================*/
printf("Total cull is %d\n",*num_culled);
printf("Total sale is %d\n",*num_sold);
printf("Total death is %d\n",*num_death);
system("pause");
	write_number_animals(NumberAnimalDataFile,NumGrpAnimal,id_bull_group,index_column) ;
	write_test_result(TestResultFile, test_result_table,nrow_testschedule,ncol_test_table) ;

} //main ends here


/*=================================================================================
FUNCTION
==================================================================================*/
 /* -------------------------------------------------------------------------- */
/* ADD animals TO production_pointer */
/* -------------------------------------------------------------------------- */
void add_animal_group(struct animal_node *List[], int manage_id, struct animal_node *node_to_add )
{     

 struct animal_node* current_node1;
 current_node1 = List[manage_id]; // thought X[a] = *(X+a) but seems X[a] is a pointer, so it's address
 if(current_node1 == NULL)
    {
        List[manage_id] = node_to_add;
    }
 else
    {
       current_node1->previous_node = node_to_add ;
       node_to_add -> next_node = current_node1;
       List[manage_id] = node_to_add;
       
    }

}

 /* -------------------------------------------------------------------------- */
/* add_stub: ADD event TO event_day POINTER */
/* -------------------------------------------------------------------------- */
void add_event_node(struct event_node *event_day[], int day, struct event_node *node_to_add )
{     

 struct event_node *current_node1;
 // printf("event day is %p", event_day[day]);
 
 current_node1 = event_day[day];
 if(current_node1 == NULL)
    {
    	// printf("current node is null and then connect to %lld \n", node_to_add -> akey);
        event_day[day] = node_to_add;
    //    printf("this was the first event of the day %d",day);
        // printf("event day is now pointing to %p", event_day[day]);
    }
 else
    {
    //	printf("current node is %lld", current_node1 -> akey);
       node_to_add -> next_node = current_node1;
    //   	printf("next node is now %lld", node_to_add -> akey);
       event_day[day] = node_to_add;
    //    printf("current node is now %lld",event_day[day]->akey);
    //    system("pause");
    }

}

/*=========Remove animal node from linked list================*/
void remove_animal_group(struct animal_node* FarmGroupList[],int current_grp,struct animal_node* current_animal)
{

if (current_animal -> previous_node != NULL)
	{ 
		struct animal_node *prev_animal;
        prev_animal = current_animal -> previous_node ; // get the node which moving one is connected
	    if (current_animal -> next_node != NULL) // if moving one's next node is conencted
	            {
	  	        struct animal_node *next_animal;
                next_animal = current_animal -> next_node ; // get the next animal
                prev_animal -> next_node = next_animal;    // reconenct previous_animal's next node
                next_animal -> previous_node = prev_animal;//similarly reconnect next_animal's previous node
		         //  printf("B ends");
				}	
		else // if next node is null
		        {
		   	    prev_animal -> next_node = NULL ;
		        }
		    //  printf("A ends");    
	}
else // if previous node is null 
	{
	    if (current_animal -> next_node != NULL) //and if next node is not null
	    {  
	  //  printf("D starts") ;
	  	struct animal_node *next_animal;
        next_animal = current_animal -> next_node ;
     //   printf("next animal akey is %lld",next_animal->akey) ;
        next_animal -> previous_node = NULL;
        FarmGroupList[current_grp] = next_animal;
      //  printf("D ends") ;
		}
		else 
		{
	//		printf("E starts") ;
		FarmGroupList[current_grp] = NULL ; 
		}
	//	       printf("C ends");  
	}
	current_animal->previous_node = NULL;
	current_animal->next_node = NULL;
}
/*---------------------------------------------------------------------------
ADD update_markov_date
--------------------------------------------------------------------------*/
double update_markov_date(double today_date, double **List_mng_status,
double** cull_sell_rate, double** mortality,
struct animal_node *FarmGroupList[], int next_non_markov_date,
int* num_culled, int* num_sold, int* num_death, double bv_f, double rate_sero_conversion, int id_lact_group2
int* num_bobby)
{
// Calculate the sum of rate of M events
	double day_to_markov;
	struct animal_node* current_animal;
	struct animal_node* next_animal;
	struct animal_node* previous_animal;
	double current_rate;
	double sum_rate = 0.0;
	double demographic_rate = 0.0 ;
	double total_sero_conversion = 0.0 ;
	double accumu_rate = 0.0 ;
	double r_SE,r_EIs, r_IsIc, r_IcIs, r_IsE;
	double av_milk_production = 15.0;
	double calf_milk_consumption = 4.0;
	//double sum_subclinical_rate = 0.0;
	//double sum_clinical_rate = 0.0;
	//double sum_lameness_rate = 0.0 ;
	//double sum_detection_subclinical = 0.0;
//	double sum_detection_clinical = 0.0 ;
	//double sum_detection_lameness = 0.0 ;

	int random_int;
	int i,j;
	int k = 0;
	int pregnant_status ;
	int change ;
	int num_CM_animal_undetected = 0;
	int counter ;
	
	//milk related
	double total_colostrum;
	double waste_milk;
	double healthy_milk;
	double required_milk;
	double shed_colostrum;
	double shed_waste;
	double shed_healthy;
/*================================================================================*/
/*===== Calculate the total rate for Markov events=================================*/
//No need to count the numebr of each status but calculate disease rate
	for(i=0; i<id_bull_group+1; i++)
	{//for each management group
	//printf("Management grp is %d",i);
	
	current_animal = FarmGroupList[i] ;
	List_mng_status[i][1] = 0; //reset
	List_mng_status[i][2] = 0; //reset: sero-conversion
	List_mng_status[i][14] = 0 ; //reset-MBOVIS disease transition rate
	List_mng_status[i][15] = 0 ; //reset - CM incidence rate (only milker)
	List_mng_status[i][16] = 0; //reset - CM detection rate (only milker)
/*================Calculating Sum of demographic event and CM==============*/
	while(current_animal!=NULL)
	{
		//printf("summing up rate");
		//visit each linked list and add up rate
		current_rate = current_animal->sum_markov_rate ;
	//	printf("Current rate is %lf\n",current_rate) ;
		List_mng_status[i][1] = List_mng_status[i][1]+current_rate;
		
	//***Sero-conversion event
		if(current_animal->mbovis_sero_status==0 && current_animal->mbovis_status!=0)
		{
		List_mng_status[i][2] = List_mng_status[i][2] + rate_sero_conversion ;	
		}
		
		
/*=======Calculating total rate of clinical mastitis=================*/		
	if(i==id_lact_group2 && dry_season==0)//meaning it's milking season
	{
		//calculating number of mastitis (do I need this because it should be stored in List_mng_status[3][3]
	//	if(current_animal->mastitis_status==1)
	//	{
		//	num_mastitis ++ ;///note num_mastitis is the number of non-bovis mastitis
		//}
	
		if(current_animal->NB_mastitis_status==0 && current_animal->mastitis_detected==0)
		{
			List_mng_status[i][15] = List_mng_status[i][15] + current_animal->mastitis_rate ;
		}
		//shall I count number of CM animals undetected one by one?
		if((current_animal->NB_mastitis_status==1||current_animal->mbovis_status==3)&&
		current_animal->mastitis_detected==0) //having CM but not detected
		{
			num_CM_animal_undetected++;
		}
		
	}
/*=======Calculating total rate of clinical mastitis done=================*/	
	//TASK_1 deleted TOP
	//	{

		
		//sum_subclinical_rate = sum_subclinical_rate + current_animal->sub_mastitis_rate;
		//if(current_animal->mastitis_status==1)//if subclinical
		//{
		//	sum_clinical_rate = sum_clinical_rate + @@@rate of clinical incindece
		//}
		
	////TASK_1 deleted DOWN	} 
				
	//	if(current_animal->next_node == NULL)
	//	{
	//		printf("Next node is NULL");
	//		}	
		current_animal = current_animal->next_node ;
	}
/*================Calculating Sum of demographic event done==========*/

/*================Calculating sum of Cow to calf transmission=======*/

/// TASK_2 deleted
		if(i==0 && dry_season==0)
		{
		//TASK_28 need to change milk supply - keep track of the total amount of colostrum and waste milk
		//keep them in tank and infection pressure changes every day
		//How do we do this? 
		//TASK_12 add milk to calf transmission	
		//transmission pressure depends on
		//1. number of calves and number of animals with colostrum and their bovis status
		//2. number of animals that are detected and their bovis status
	
		total_colostrum = (milk_numbers[0]+milk_numbers[3]+milk_numbers[6])*av_milk_production;
		waste_milk = (milk_numbers[2]+milk_numbers[5]+milk_numbers[8])*av_milk_production;
		healthy_milk = (milk_numbers[1]+milk_numbers[4]+milk_numbers[7])*av_milk_production;
		required_milk = (List_mng_status[0][0]+(num_bobby*))*calf_milk_consumption ;
		shed_colostrum = bv_f*milk_numbers[3] + milk_numbers[6];
		shed_waste = bv_f*milk_numbers[5] + milk_numbers[8];
		shed_healthy = bv_f*milk_numbers[4] + milk_numbers[7];
		
		//scenario1: when colostrum milk can afford all calves
		//if(total_colostrum>=required_milk)
		//if colostrum is sufficient
		//	{
		//	beta = shed_colostrum*beta_mc*calf_milk_consumption/total_colostrum;
		//	}
			//here beta is the amount of infection pressure individual receives
		//scenario2: when colostrum + waste milk can afford all calves
		//TASK_28
		if(total_colostrum+waste_milk>=required_milk)
		{
			beta = beta_mc*(shed_colostrum+shed_waste)*calf_milk_consumption/(total_colostrum+waste_milk) ;
		}
		//scenario3: when more milk is required
		else
		{
			beta = beta_mc*(shed_colostrum+shed_waste+	(required_milk-total_colostrum-waste_milk)*shed_healthy/healthy_milk
			)/(List_mng_status[0][0]+(num_bobby*));
		 //this comes from [1], [4], [7]
		}
		List_mng_status[i][9] = beta*List_mng_status[i][5];
		//for now calf there is no transition from exposed to infectious
		//Assumption_8
		List_mng_status[i][10] = 0;//exposed to Is
		List_mng_status[i][11] = 0;//Is to Ic
		List_mng_status[i][12] = 0];//Ic to Is
		List_mng_status[i][13] = 0;//Is to E 
		
		List_mng_status[i][14] = List_mng_status[i][9];//+
		//List_mng_status[i][10]+List_mng_status[i][11]+List_mng_status[i][12]+
	//	List_mng_status[i][13] ;	
	//TASK_12 DONE	
		}
/*==========TRANSMISSION FROM COW TO CALVES DONE===============================*/

/*==========TRANSMISSION BETWEEN COWS=====================================*/
		if(i==id_lact_group2 && dry_season==0)
		{
		//now have to add rate of subclinical mastitis happening if not under treatment (which is mastitis_detected==1)
		//num_no_mastitis = List_mng_status[i][0] - List_mng_status[i][3] ;
		//sum_clinical_rate = rate_subclinical_mastitis*num_no_mastitis;
		
		//number of non-mastitis animals (excluding animals that are already sub/clincial mas)
	//TASK_3 deleted
		//sum_clinical_rate = List_mng_status[i][2]*rate_clinical_mastitis;
	//TASK_4 deleted
		//sum_lameness_rate = (List_mng_status[i][0]-List_mng_status[i][4])*rate_lameness ;
		//List_mng_status[i][15] = sum_subclinical_rate + sum_clinical_rate + sum_lameness_rate ;
		
//TASK_5 deleted and modified
		//sum_detection_subclinical = (List_mng_status[i][2]+List_mng_status[i][7])*rate_detect_subclinical_mastitis ;
		//TASK_13
		//sum_detection_lameness = List_mng_status[i][4]*rate_detect_lameness ;
		///CM detection either non-bovis, bovis or both
		//this does not tell me how many detected so it's wrong
		
		List_mng_status[i][16] = num_CM_animal_undetected*rate_detect_clinical_mastitis ;
		//Assumption_9
		List_mng_status[i][9] = beta_mm*(milk_numbers[4]*bv_f+milk_numbers[7])
		*List_mng_status[i][5];
		
		List_mng_status[i][10] = bv_r1*List_mng_status[i][6];//exposed to Is
		List_mng_status[i][11] = bv_r2*List_mng_status[i][7];//Is to Ic
		List_mng_status[i][12] = bv_r3*List_mng_status[i][8];//Ic to Is
		List_mng_status[i][13] = bv_r4*List_mng_status[i][7];//Is to E 
		
		List_mng_status[i][14] = List_mng_status[i][9]+
		List_mng_status[i][10]+List_mng_status[i][11]+List_mng_status[i][12]+
		List_mng_status[i][13] ;
		} //calculating sum rate for milking group done
/*==========TRANSMISSION BETWEEN COWS DONE=====================================*/

/*================Calculating sum of m.bovis related Markov sum=======*/
		

		//sum altogether
		//TASK_6 no need to change
		sum_rate = 	sum_rate + List_mng_status[i][1] + List_mng_status[i][2]+ List_mng_status[i][14] +
		List_mng_status[i][15] + List_mng_status[i][16] ;
		demographic_rate = demographic_rate + List_mng_status[i][1] ;
		total_sero_conversion = total_sero_conversion + List_mng_status[i][2] ;
//	printf("sum is %lf\n",sum_rate);
	}
	
/*=============================================================================*/
//printf("sum rate is %lf",sum_rate);
//Now calculate a waiting time
	 double random_value = (double)(rand()+1)/((double)(RAND_MAX)+1);//avoids 0
//	 printf("random value is %f",random_value);
     day_to_markov =(-log(random_value))/sum_rate ; // Waiting time
     if(day_to_markov<0)
     {
     	printf("Error: markov date cannot be negative!") ;
     	system("pause") ;
	 }
	// printf("markov date is %lf",day_to_markov);
//	 printf("non markov is %d",next_non_markov_date);
/* NOW ASSESS WHETHER MARKOV EVENT OCCURS FIRST OR NOT*/
	if (next_non_markov_date>day_to_markov+today_date)
	{ // if markov comes first, choose markov events 
	   k++;
	   //printf("This is %d th markov",k);
	   today_date = day_to_markov+today_date ;
	  // printf("ceil days is %lf",today_date);
	   //system("pause");
	  
	random_value =  (double)(rand()+1)/(double)(RAND_MAX+1)*sum_rate;
	
	//Now determine this value is greater or smaller than sum of demographic rate
/*======================MARKOV is DISEASE EVENT===================*/
	if(random_value>demographic_rate)
	{
	//now do events other than demographic - CM occur + CM detection + MBOVIS transition
	i = 0;
	accumu_rate = 0.0 ;
	random_value = random_value - demographic_rate;
	int index_animal_to_chose = 0;
	struct animal_node* animal_to_chose ;
	
	if(random_value<List_mng_status[id_lact_group2][15]) //CM incidnece
	/*==============1. CM incidence===================================*/
	{
		//TASK_7 completed
	animal_to_chose = 	FarmGroupList[id_lact_group2];
	/* TASK_7 TOP
		if(random_value<sum_subclinical_rate)
		{
			//choose one random number from the number of animals without any mastitis
			index_animal_to_chose = rand()%num_no_mastitis + 1 ;
			counter = 0 ;
				while(counter<index_animal_to_chose)
				{
				if(animal_to_chose->mastitis_status==0)
				{
				counter++:	
				}
				if(counter==index_animal_to_chose)
				{
				animal_to_chose->mastitis_status = 1;//become subclinical mastitis
				List_mng_status[3][2]++; //increase the number of subclinical mastitis
				
				}
				else
				{
				animal_to_chose = animal_to_chose->next_node ; //go to next animal
				}	
				}
			
		} TASK_7 DOWN*/
		//TASK_8
		
		///how many animals don't have non-bovis CM?
		///Assumption_1: 
		//each cows have different rate so need to consider that
		random_value =  (double)(rand()+1)/(double)(RAND_MAX+1)*List_mng_status[id_lact_group2][15];
		
		while(random_value>accumu_rate)
		{
			//TASK_21
		 accumu_rate = accumu_rate+ animal_to_chose->mastitis_rate ;
		animal_to_chose = animal_to_chose->next_node;
		}
			if(animal_to_chose->NB_mastitis_status==1)
			{
				printf("This animal already has non Bovis CM!");
				system("pause");
			}
			if(animal_to_chose->mastitis_detected==1)
			{
				printf("This animal already is detected as CM!");
				system("pause");
			}
			
		
			
				animal_to_chose->NB_mastitis_status = 1;//become clinical mastitis
				animal_to_chose->mastitis_rate = 0; //now already CM so incidence rate becomes 0
				
				
				
				
		
		/* TASK_8: remove lameness
		else //else lameness
		{
			index_animal_to_chose = rand()%List_mng_status[3][4] + 1 ;
			///random number from the number of lameness
			counter = 0 ;
				while(counter<index_animal_to_chose)
				{
				if(animal_to_chose->lameness_status==0) //now assuming that m bovis not causing lameness
				{ //if subclinical but not due to M bovis
				counter++:	
				}
				if(counter==index_animal_to_chose)
				{
				animal_to_chose->lameness_status = 1;//become lameness
				List_mng_status[3][4]++; //increase the number of clinical mastitis
				}
				else
				{
				animal_to_chose = animal_to_chose->next_node ; //go to next animal
				}	
				}
		} TASK_8 remove lameness*/
	}
/*====================CM DETECTION========================================*/
	else if((random_value-List_mng_status[id_lact_group2][15])<List_mng_status[id_lact_group2][16])
	{
		//Assumption_2: 
	//random_value =  (double)(rand()+1)/(double)(RAND_MAX+1)*List_mng_status[3][16];	
	/*TASK_9 TOP
		if(random_value<sum_detection_subclinical)
		{//subclincial mastitis detection
		//
		//clinical mastitis will be detected - so pick up one CM animal 
		index_animal_to_chose = rand()%(List_mng_status[3][2]+List_mng_status[3][7]) + 1 ;	
		counter = 0 ;
				while(counter<index_animal_to_chose)
				{
				if(animal_to_chose->mastitis_status==1) //sub mastitis
				{ //if subclinical but not due to M bovis
				counter++:	
				}
				if(counter==index_animal_to_chose)
				{
				animal_to_chose->mastitis_detected = 1;//become detected
				///now the duration requried for this animal to recover depends on mbovis status
				if(animal_to_chose->mbovis_status==0)
				{
					//add event that this animal is going to recover in x days
					days_to_recover = @@@;
				}
				else
				{
					days_to_recover = @@@;
				}
				//List_mng_status[3][4]++; //increase the number of clinical mastitis
				}
				else
				{
				animal_to_chose = animal_to_chose->next_node ; //go to next animal
				}	
				}
		} TASK_9 */
		index_animal_to_chose = rand()%num_CM_animal_undetected + 1 ;
		counter = 0 ;
				while(counter<index_animal_to_chose)
				{
				if((animal_to_chose->NB_mastitis_status==1||animal_to_chose->mbovis_status==3)
				&& animal_to_chose->mastitis_detected==0) ///CM either non-bovis or bovis and not detected
				{ //if subclinical but not due to M bovis
				counter++:	
				}
				if(counter==index_animal_to_chose)
				{
				animal_to_chose->mastitis_detected = 1;//become detected
				animal_to_chose->mastitis_rate = 0; //CM rate was non-zero before if CM was only due to bovis
				//Update milk_numbers
				//animals have CM - either bovis, non-bovis or both
				mbovis_status = animal_to_chose->mbovis_status ;
				if(mbovis_status==0)
				{
				milk_numbers[2]++;
				milk_numbers[animal_to_chose->non_colostrum]--;
				}
				else
				{
				milk_numbers[(mbovis_status-1)*3+2]++;
				milk_numbers[(mbovis_status-1)*3+animal_to_chose->non_colostrum]--;	
				
				}
				
					//TASK_22: need to minus depending on its colostrum status
				
				//TASK_22: need to change milk_number[9], milk_number for non-shed and sub-shed too
				//depending on colostrum status
			
				
					//TASK_14 add treatment length and add event where treatment done and then starts again
				
					new_event = (struct event_node*)malloc(sizeof( struct event_node ));
					new_event->event_type = 11;//treatment
					new_event->akey = animal_to_chose->akey;
					new_event->animal = animal_to_chose; //let's see if this works
					new_event->next_node = NULL ;
					add_event_node(event_day,next_non_markov_date+length_CM_treatment, new_event) ;
					break; 
				}
				else
				{
				animal_to_chose = animal_to_chose->next_node ; //go to next animal
				}	
				}
	} //CM detection done
/*====================CM DETECTION DONE================================*/
		/* TASK_9 remove lameness detection
		else //lameness detection
		{
				index_animal_to_chose = rand()%(List_mng_status[3][4]) + 1 ;	
				counter = 0 ;
				while(counter<index_animal_to_chose)
				{
				if(animal_to_chose->lameness_status==1) //sub mastitis
				{ //if subclinical but not due to M bovis
				counter++:	
				}
				if(counter==index_animal_to_chose)
				{
				animal_to_chose->lameness_detected = 1;//become detected
				///now the duration requried for this animal to recover does not depend on mbovis status for lameness
				days_to_recover = @@;
				//add events
				//List_mng_status[3][4]++; //increase the number of clinical mastitis
				}
				else
				{
				animal_to_chose = animal_to_chose->next_node ; //go to next animal
				}	
				}
		} TASK_9 remove lameness detection*/
/*=================SERO-CONVERSION==================================*/
	else if((random_value-List_mng_status[id_lact_group2][15]-List_mng_status[id_lact_group2][16])<total_sero_conversion)
	{
	random_value = 	((double)(rand()+1)/(double)(RAND_MAX+1))*total_sero_conversion ;
	while(random_value>accumu_rate)
	{
		accumu_rate = accumu_rate+ List_mng_status[i][2] ;
		i++;
		//printf("i is %d\n", i) ;
	}	
	index_animal_to_chose = rand()%(int)(List_mng_status[i-1][2]/rate_sero_conversion)+1;
	animal_to_chose = FarmGroupList[i-1];
	counter = 0;
	while(counter<index_animal_to_chose)	
	{
		if(animal_to_chose->mbovis_sero_status==0 && animal_to_chose->mbovis_status!=0)
		{
			counter++;
		}
		if(counter==index_animal_to_chose)
		{
		animal_to_chose->mbovis_sero_status = 1 ;
		break;	
		}
		else
		{
		animal_to_chose = animal_to_chose->next_node ;
		}
	
	}
	}
	
	
/*=================Sero-conversion done===============================*/	
	else
/*====================MBOVIS STATUS TRANSITION=============================*/
	{
	random_value = 	random_value - List_mng_status[id_lact_group2][15] -List_mng_status[id_lact_group2][16] -total_sero_conversion;
	//i is already 0 in above
	///TASK_15: need to add cow to calf transmission
	while(random_value>accumu_rate)
	{
		accumu_rate = accumu_rate+ List_mng_status[i][14] ;//List_mng_status[i][14] is sum of disease markov in a given mng grp
		i++;
		//printf("i is %d\n", i) ;
	}
	if(i==4)//i-1 is 3 means it's milkers
	{
		
	random_value = ((double)(rand()+1)/(double)(RAND_MAX+1))* List_mng_status[i-1][14];
	//********Below obtain how many animals in each category, and get random number below pop size
	//********Then, select this nth animal having this status and change their status
	if(random_value<List_mng_status[i-1][column_r_SE]) //S to E
	{
	index_animal_to_chose = rand()%List_mng_status[i-1][column_Sus]+1 ;
	index = 5;
	change = 1;//change indicates how status changes e.g. from 5 to 6 means 1 plus
	//now jump to the animal linked list and change the status
	//no change in milk_numbers
	
	}
	else if(random_value<List_mng_status[i-1][column_r_SE+1])
	{
	index_animal_to_chose = rand()%List_mng_status[i-1][column_Sus+1]+1 ; //E to Is
	index = 6;
	change = 1;
	//now jump to the animal linked list and change the status
	
	}
	else if(random_value<List_mng_status[i-1][column_r_SE+2])
	{
	index_animal_to_chose = rand()%List_mng_status[i-1][column_Sus+2]+1 ; //Is to Ic
	index = 7;
	change = 1;
	//now jump to the animal linked list and change the status	
	
	}
	else if(random_value<List_mng_status[i-1][column_r_SE+3])
	{
	index_animal_to_chose = rand()%List_mng_status[i-1][column_Sus+3]+1 ; //Ic to Is
	index = 8;
	change = -1;
	//now jump to the animal linked list and change the status	
	}
	else if(random_value<=List_mng_status[i-1][column_r_SE+4])
	{
	index_animal_to_chose = rand()%List_mng_status[i-1][column_Sus+4]+1 ; //Is to E
	index = 7;
	change = -1;
	//now jump to the animal linked list and change the status	
	}
	else
	{
		printf("Out of range!");
		system("pause");
	}

	animal_to_chose = FarmGroupList[i-1];
	counter = 0; //if counter becomes equal to index_animal_to_chose, then that's the animal to choose
	//index_animal_to_chose indicates nth animal to choose
	//index indicates which disease status to choose (e.g. 6 means E staus)
	//index ranges 5 to 8, and status corresponds to 0 to 3
	mbovis_status = index-5 ;
	while(counter<index_animal_to_chose)
	{
		if(animal_to_chose->mbovis_status==mbovis_status) 
		{
		counter++:	
		}
		if(counter==index_animal_to_chose)
		{
		//TASK_19
	    //Update milk_numbers table
			if(mbovis_status<=1)
			{
				mbovis_status = 0; //E set to 0 for milk_numbers table changing purpose
			}
		  if(animal_to_chose->mastitis_detected==1)
		  {
		  	current_milk_numbers = mbovis_status*3 + 2 ;
		  }
		  else
		  {
		  	current_milk_numbers = mbovis_status*3 + animal_to_chose->non_colostrum ;
		  }
		
		if(index!=5)
		{///milk_numbers don't change if animal moving from S to E
		milk_numbers[current_milk_numbers]--;
		milk_numbers[current_milk_numbers+change*3]++;
		}
		
	    
		animal_to_chose->mbovis_status = animal_to_chose->mbovis_status + change;
		List_mng_status[i-1][index]--; //n of index status animal minus 1
		List_mng_status[i-1][index+change]++; // n of (index+change) status animal plus 1
		break;		
		}
		else
		{
		animal_to_chose = animal_to_chose->next_node ; //go to next animal
		}	
	 }//changing disease status done	
    }//if it is milkers done
    else if(i==1)
    {//calf infection
    index_animal_to_chose = rand()%List_mng_status[i-1][column_Sus]+1 ;
	animal_to_chose = FarmGroupList[i-1];
	counter = 0;
      while(counter<index_animal_to_chose)
      {
      	 if(animal_to_chose->mbovis_status==0)//if susceptible 
		{
		counter++:	
		}
		if(counter==index_animal_to_chose)
		{
			List_mng_status[0][column_Sus]--;
			List_mng_status[0][column_Sus+1]++;
			break;
	    }
	    else
	    {
	    animal_to_chose = animal_to_chose->next_node ;	
		}
	  }
	   
	} //if calf done
	else
	{
		printf("It's not either milker or calf!'");
		system("pause");
	}
	}
/*===========MBOVIS RELATED DISEASE EVENTS DONE============================*/	
	
	
	//event in k-1 happens
	} //end of disease markov
/*=====================MARKOV IS DEMOGRAPHIC EVENT========================================*/
	else
	{
		//now do demographic
		//TASK_19
	i = 0;
	accumu_rate = 0.0 ;
	//now choose one of management group for the event
	while(random_value>accumu_rate)
	{
		accumu_rate = accumu_rate+ List_mng_status[i][1] ;
		i++;
		//printf("i is %d\n", i) ;
	} //i is the group at which an event occurs	
	//printf("group decided");
	random_value = ((double)(rand()+1)/(double)(RAND_MAX+1))* List_mng_status[i-1][1];
	//now decides which animal is going to have an event
	accumu_rate = 0; //reset accumulate, this time use this for sum of rate over animals
	//j=0;
	if(FarmGroupList[i-1]==NULL)
	{
		printf("There is no animals in this group!") ;
	}
	else
	{
		current_animal = FarmGroupList[i-1] ;
		accumu_rate = current_animal->sum_markov_rate ;
		printf("choose animals");
	while(random_value>accumu_rate)
		{
		//	printf("curernt animal is %d",current_animal->akey) ;
		current_animal = current_animal->next_node ;
		if(current_animal==NULL)
		{
			printf("Current animal becomes null");
			system("pause");
		}
		accumu_rate = accumu_rate + current_animal->sum_markov_rate ;
		}
	
	//	printf("animal decided");
		//current_animal is the animal that will have an event
	//Choose event
	random_value = (double)(rand()+1)/(double)(RAND_MAX+1)* (current_animal->sum_markov_rate);
	//printf("random value is %lf\n",random_value);
	//printf("cull rate is %lf\n",cull_sell_rate[current_animal->index_cull_sell][1+2*(current_animal->pregnant_status)]);
	//printf("sale rate is %lf\n",cull_sell_rate[current_animal->index_cull_sell][2+2*(current_animal->pregnant_status)]);
	//printf("death is %lf\n",mortality[current_animal->index_mortality][1]);
	//system("pause");
	//obtain rates for culling, death and sale
	if(random_value<=cull_sell_rate[current_animal->index_cull_sell][1+2*(current_animal->pregnant_status)])
	{
		//cull
		(*num_culled)++;
		printf("this is cull");
		//system("pause");
	}
	else if(random_value<=(cull_sell_rate[current_animal->index_cull_sell][1+2*(current_animal->pregnant_status)]+cull_sell_rate[current_animal->index_cull_sell][1+2*(current_animal->pregnant_status)]))
	{
		//sale
		(*num_sold)++;
		printf("this is sold");
		//system("pause");
	}
	else
	{
		//death
		(*num_death)++;
		printf("this is death");
		//system("pause");
		
	}
	//printf("i is %d",i);
	List_mng_status[i-1][0]--;
	if(List_mng_status[i-1][0]<0)
	{
		printf("Warning! Animals below 0");
		system("pause");
	}
	//Need to update List_mng_status-bovis status and milk_numbers table
	mbovis_status = current_animal->mbovis_status ;
	List_mng_status[i-1][mbovis_status+5]-- ;//minus
	//milk_numbers changes only if it is milking animals
	if(i==4)
	{
		if(current_animal->mastitis_detected==1)
			{
			if(mbovis_status==0)
			{
			milk_numbers[2]--;
			}
			else
			{
			milk_numbers[(mbovis_status-1)*3+2]--;
			}
			}
		else
			{
			if(mbovis_status==0)
			{
			milk_numbers[current_animal->non_colostrum]--;
			}
			else
			{
			milk_numbers[(mbovis_status-1)*3+current_animal->non_colostrum]--;	
			}
			
			}	
	}

	printf("event decided");
	/*==========remove this animal=====================================*/
	if (current_animal -> previous_node != NULL)
	{ 
//	printf("A");
	// if moving animal's previous node is conencted to other animal
	    //    
		struct animal_node *prev_animal;
        prev_animal = current_animal -> previous_node ; // get the node which moving one is connected
	    if (current_animal -> next_node != NULL) // if moving one's next node is conencted
	            {//Removal situation D
	              //	printf("B starts");
	  	        struct animal_node *next_animal;
                next_animal = current_animal -> next_node ; // get the next animal
                prev_animal -> next_node = next_animal;    // reconenct previous_animal's next node
                next_animal -> previous_node = prev_animal;//similarly reconnect next_animal's previous node
		           //printf("B ends");
				}	
		else // if next node is null
		        {//Removal Situation B
		   	    prev_animal -> next_node = NULL ;
		   	    //printf("A ends");
		        }
		    //      
	}
else // if previous node is null 
	{//Situation C
	      //   printf("C starts");	
	    if (current_animal -> next_node != NULL) //and if next node is not null
	    {  
	  //  printf("D starts") ;
	  	struct animal_node *next_animal;
        next_animal = current_animal -> next_node ;
     //   printf("next animal akey is %lld",next_animal->akey) ;
        next_animal -> previous_node = NULL;
        FarmGroupList[i-1] = next_animal;
       // printf("D ends") ;
		}
		else 
		{//Situation A
	//		printf("E starts") ;
		FarmGroupList[i-1] = NULL ; 
	//	printf("C ends");
		}
		         
	}
// moving_animal->current_pro_id = -10 ; //means it's sent to slaughterhouse 
 current_animal->previous_node = NULL;
 current_animal->next_node = NULL ;
current_animal->present = 0;
	/*==Removing animal done=========================*/

//	printf("free now");
//	printf("akey is %lld",current_animal->akey) ;
//	free(current_animal) ; //free memory
//	printf("free done");
//	printf("akey is %lld",current_animal->akey) ;
//	system("pause");
	}
} //end of demographic markov event
	}  
	else
	{
		today_date = (double)next_non_markov_date ;
	 } 
return(today_date) ;	
}
// END OF update_markov_date

/*==============================================================================================*/

/* READ culling and selling data*/
void read_cull_sell_rate(char CullSellRatesFile[], double **cull_sell_rate, int num_cull_sell_steps)
{     
    /* OPEN INPUT FILE */
    printf("Opening file");
    FILE *temp = fopen(CullSellRatesFile,"r"); 
       
    /* DECLARE STORAGE VARIABLES */
    int line_num, day;
    double cull_rate_empty, sale_rate_empty,cull_rate,sale_rate;
    
    /* READ LINES OF FARM FILE */
    for(line_num = 0; line_num < num_cull_sell_steps; line_num++)
      { 
         fscanf(temp, "%d,%lf,%lf,%lf,%lf",&day, &cull_rate_empty, &sale_rate_empty, &cull_rate,&sale_rate);
         
             cull_sell_rate[line_num][0] = day;
             //printf("day is %d",day) ;
             cull_sell_rate[line_num][1] = cull_rate_empty/1000;
             cull_sell_rate[line_num][2] = sale_rate_empty/1000;
             cull_sell_rate[line_num][3] = cull_rate/1000;
             cull_sell_rate[line_num][4] = sale_rate/1000;
      }
      fclose(temp);
}

/*==========READ MORTALITY DATA====================================*/
void read_mortality(char MortalityFile[], double **mortality, int num_mortality_steps)
{     
    /* OPEN INPUT FILE */
    FILE *temp = fopen(MortalityFile,"r"); 
       
    /* DECLARE STORAGE VARIABLES */
    int line_num, day;
    double rate;
    
    /* READ LINES OF FARM FILE */
    for(line_num = 0; line_num < num_mortality_steps; line_num++)
      { 
         fscanf(temp, "%d,%lf",&day, &rate);
         
             mortality[line_num][0] = day;
             mortality[line_num][1] = rate/1000;
      }
      fclose(temp);
}
/*=======READ TEST SCHEDULE DATA===============================================================*/

void read_testschedule(char TestScheduleFile[], int **TestSchedule, int nrow_testschedule)
{     
    /* OPEN INPUT FILE */
    printf("Opening file");
    FILE *temp = fopen(TestScheduleFile,"r"); 
       
    /* DECLARE STORAGE VARIABLES */
    int line_num, id, date;
    int slaughter_PCR_nomilk_calf,slaughter_PCR_milk_calf,slaughter_ELISA_calf;
    int live_PCR_nomilk_calf,live_PCR_milk_calf,live_ELISA_calf;
    int slaughter_PCR_nomilk_adult,slaughter_PCR_milk_adult,slaughter_ELISA_adult;
    int live_PCR_nomilk_adult, live_PCR_milk_adult, live_ELISA_adult ;
    int btank, bdisc, ELISA_type ;
   
    
    /* READ LINES OF FARM FILE */
    for(line_num = 0; line_num < nrow_testschedule; line_num++)
      { 
         fscanf(temp, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d, %d,%d,%d,%d,%d,%d, %d"
		 ,&id, &date, 
		 &slaughter_PCR_nomilk_calf, &slaughter_PCR_milk_calf,&slaughter_ELISA_calf,
		 &live_PCR_nomilk_calf, &live_PCR_milk_calf,&live_ELISA_calf,
		 &slaughter_PCR_nomilk_adult, &slaughter_PCR_milk_adult,&slaughter_ELISA_adult,
		 &live_PCR_nomilk_adult, &live_PCR_milk_adult,&live_ELISA_adult,
		 &btank, &bdisc, &ELISA_type
		 );
         
             TestSchedule[line_num][0] = id;
             TestSchedule[line_num][1] = date;
             TestSchedule[line_num][2] = slaughter_PCR_nomilk_calf;
             TestSchedule[line_num][3] = slaughter_PCR_milk_calf;
             TestSchedule[line_num][4] = slaughter_ELISA_calf;
             
             TestSchedule[line_num][5] = live_PCR_nomilk_calf;
             TestSchedule[line_num][6] = live_PCR_milk_calf;
             TestSchedule[line_num][7] = live_ELISA_calf;
             
             TestSchedule[line_num][8] = slaughter_PCR_nomilk_adult;
             TestSchedule[line_num][9] = slaughter_PCR_milk_adult;
             TestSchedule[line_num][10] = slaughter_ELISA_adult;
             
             TestSchedule[line_num][11] = live_PCR_nomilk_adult;
             TestSchedule[line_num][12] = live_PCR_milk_adult;
             TestSchedule[line_num][13] = live_ELISA_adult;
             
             TestSchedule[line_num][14] = btank;
             TestSchedule[line_num][15] = bdisc;
             TestSchedule[line_num][16] = ELISA_type ;
             
      }
      fclose(temp);
}

/*========WRITE DOWN THE NUMBER OF ANIMALS IN EACH MANAGEMENT GROUP IN A GIVEN TIME===================*/
int write_number_animals(char* NumberAnimalDataFile,double** NumGrpAnimal,int id_bull_group,int n_column_output)
{
	FILE *Out = fopen(NumberAnimalDataFile, "w") ;
	int line_num, col_num;
	
	for (line_num = 0 ; line_num < id_bull_group+2; line_num ++)
	{
		for (col_num = 0;col_num <n_column_output+1 ; col_num++ )
		
		{
	    fprintf(Out,"%lf,",NumGrpAnimal[line_num][col_num]);
        }
        fprintf(Out,"\n");
    }
	fclose(Out);
	return 0;
}
/*============VISUALISE EVENTS=====================================================================*/
void visualize_list(struct event_node *event_day[], int day)  
{

 
 struct event_node *current_node1;
 current_node1 = event_day[day];
 //printf("current node is %lld",current_node1->akey);
 if(current_node1 != NULL)
    {
    //   printf("Day %d: ", day );
       while(current_node1 != NULL)
          {
          //	if (current_node1->des_pro_id==(1225*3+2)||current_node1->src_pro_id==(1225*3+2))
           //   {
			// printf("Event is %d", current_node1 -> event_type);
			 // system("pause") ;
		   // }
              
              current_node1 = current_node1 -> next_node; 
              
          }  
       printf("\n");
   }
   
 

}

/*=========TESTING ANIMALS=============================================================================*/
void test(struct animal_node* current_animal,int test_type,int slaughter_live, int group,int* test_result_table,
double Se_bio, double Sp_bio, double Se_ID, double Sp_ID, double PCR_live, double PCR_slaughter) 
{
	double random;
	double Se, Sp;
	random = (double)rand()/(double)RAND_MAX;
	int slaughter_PCR = 0;
	int ELISA_calf_bio = 1;
	int ELISA_calf_ID = 2;
	int ELISA_adult_bio = 3;
	int ELISA_adult_ID = 4;
	int PCR_calf = 5;
	int PCR_adult_nomilk = 6;
	int PCR_adult_milk = 7;

	//need to include Se/Sp into passing parameters
	if(test_type<=2)
	{
		if(test_type!=1) //ELISA bio
		{
			Se = Se_bio;
			Sp = Sp_bio;
			if((current_animal->mbovis_sero_status==1 && random<Se) ||
			(current_animal->mbovis_sero_status==0 && random>Sp))
			{
		
			current_animal->mbovis_elisa_status=1 ;
			test_result_table[group*2+ELISA_calf_bio] ++;
			} //detecting true pos
		}
		if(test_type!=0)
		{
			Se = Se_ID;
			Sp = Sp_ID;
			if((current_animal->mbovis_sero_status==1 && random<Se) ||
			(current_animal->mbovis_sero_status==0 && random>Sp))
			{
		
			current_animal->mbovis_elisa_status=1 ;
			test_result_table[group*2+ELISA_calf_ID] ++;
			}
		}
	
	
	}
	else //ANY PCR
	{
		if(test_type==3)//live PCR nomilk
		{
			Se = PCR_live;
			if(current_animal->mbovis_status>=1 && random < Se)
			{
			current_animal->mbovis_pcr_status=1 ;
			test_result_table[PCR_calf+group] ++;	
			}
			
		}
		else if(test_type==4)//slaughter  PCR
		{
			Se = PCR_slaughter;
			if(current_animal->mbovis_status>=1 && random < Se)
			{
			current_animal->mbovis_pcr_status=1 ;
			test_result_table[slaughter_PCR] ++;	
			}
			
		}
		else //milk PCR
		{
			Se = PCR_live;
			if(current_animal->mbovis_status>=2 && random < Se)
		  	{
		  	current_animal->mbovis_pcr_status=1 ;
			test_result_table[PCR_adult_milk] ++;
		  	}
		}
	

	
	}
}

/*=========Testing BULK MILK===================================================*/
void test_bulk(int tank_discard, int shedder, int total, double dilution, double PCR_live, int* test_result_table)
{
	double prop = double(shedder)/double(total); //how much diluted
	double random = (double)rand()/(double)RAND_MAX;
	if(prop > dilution) //more than threshold
	{
		if(random < PCR_live)
		{
			test_result_table[8+tank_discard]++;
		}
	}
}
/*==========Exporting test_result_table==============================================*/
int write_test_result(char* ResultFile,int* test_result_table,int nrow_testschedule,int n_column_output)
{
	FILE *Out = fopen(ResultFile, "w") ;
	int line_num, col_num;
	
	for (line_num = 0 ; line_num < nrow_testschedule; line_num ++)
	{
		for (col_num = 0;col_num <n_column_output ; col_num++ )
		
		{
	    fprintf(Out,"%d,",test_result_table[line_num][col_num]);
        }
        fprintf(Out,"\n");
    }
	fclose(Out);
	return 0;
}
