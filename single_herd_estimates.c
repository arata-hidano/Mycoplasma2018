 /*===========================================================
Mycoplasma project
Demographic model - Simulation Starting date June 1 2015. Infection comes in sampled risk date which is different for each farm



Major revision:

9 July 2019: One herd running so expland it to multiple herds at the same time
26 July 2019: TAG MTG
30 October 2019: MODIFY the code in order for the model to provide within-herd prevalence estimates with given parameters input




======A CHANGE MADE AFTER THE TAG MEETING===============================================================================
1. Don't feed waste milk to bobby calves

2. Sensitivty analysis may be necessary if farmers are using milk powder instead of bulk milk when waste milk insufficient


3. Sick cows can get infection from sick and Mbovis cows
Sick animals are currently not physically separated.
But milking order is different. Check line 3850.
Because sick cows are milked later than healthy cows, the contribution of sick cows to infection pressures on healthy cows is minimal.
And we stick to the stance that nose-to-nose transmission cannot be estimated. (but we set it to 0).
Transmission among sick cows is also minimal unless there are many sick cows and the same milk cup is used for more than one animals.
Milking cluster is disinfected after the morning use, making it unlikely that transmission occurs in PM milking after AM milking.

4. Can we think about respiratory infection? - currently not because
all calves receive the infection pressures both from respiratory and waste milk, the former depends on 
the number of infectious calves, the latter depends on the number of infectious adults contributing to milk.
Both have different beta, which cannot be determined unless there is data for weaned calves/heifer that shows nose-to-nose transmission.
Here we stick to a stance that we cannot really estimate nose-to-nose transmission parameters.
The approach will be to assume say, nose-to-nose transmission has 10% infection pressure of milk transmission and see how much it changes the results.
Also we can evaluate when the first calf starts shedding after being infected by milk - if this takes long calf-to-calf transmission is not
important before weaning. 

5. Waste milk tank - accumulate over time


=============================================================*/
/* C LIBRARIES TO INCLUDE */
#include <stdio.h>
//#include <R.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <malloc.h>
#include <ctype.h>

/*===Parameters to receive from R=========*/
//Parameters to estimate
//1. risk_date for 41 farms
//2. dilution
//3. bv_r1 to bv_r4; 4 parameters
//4. Se bio,ID,PCR live and slaughter; 4 parameters
//5. Sp for above 4
//6. scale factor
//7. beta from cow to cow and cow to calf; 2 parameters
//8. sero-conversion duration 
///58 parameters???

/*===Define whether we import pre-defined diagnostic testing dates====*/
int pre_defined_test = 0; //0 indicates we don't import pre-defined test schedule table; 1 otherwise
//int init_bovis_status = 2; //define the initial status of risk animal
int init_infected_number = 1 ; // maybe try 1,3,5,10,20
#define num_total_herd 1
//Specify the number of farms to be simulated
//the number is 41 when parametrised as of 9 July 2019, removed 2 farms that are clearly not spring calving
///*====Define variables related to recording the output=======*/
//int interval_disease_record = 5 ;//every 5 day, record the number of animals in each status
//int duration_recording = 365 ; //one year to record the disease status
//int column_disease_status = 4 ; // 4 disease status
//int row_disease_status = 73; //(365-5)/5+1


int herd_size;

/*PARAMETER TO ESTIMATE*/
int risk_date;
double dilution;
double bv_r1; // r1: rate from E to Isub
double bv_r2; // r2: rate from Isub to Icli
double bv_r3 ; // r3: rate from Icli to Isub 
double bv_r4 ; // r4: rate from Isub to E
double Se_bio ;
double Sp_bio;
double Se_ID ;
double Sp_ID ;
double PCR_live ;
double PCR_slaughter ;
double bv_f ;
double beta_mm ;
double beta_mc ;
double rate_sero_conversion;

/*=======Define waste milk management parameters==================*/
int waste_milk_accumulate = 1; //0 indicates waste milk is discarded everyday, 1 indicates they accumulate
int use_milk_powder = 0; //0 indicates insufficient calf milk comes from bulk milk 1 indicates insufficient comes from milk powder
int feed_waste_milk = 1 ; //1 indicates using waste milk, 0 otherwise
double* volume_waste_milk ; //the accumulated waste milk volume
double* density_mbovis_waste_milk ; //density of M bovis in waste milk
double* density_healthy ;

double av_milk_production = 15.0;
double calf_milk_consumption = 4.0;
double required_milk;
double topup_milk;
double healthy_milk;
double shed_healthy;

double intensity_mc = 0;
double effective_mbovis_waste_milk = 0;//the effective number of clinical shedders that shed into waste milk



/*==================================================================
DEFINE PARAMETER & STRUCTURE
===================================================================*/

int nrow_testschedule = 381; 
//To adjust this and export appropriate dataset see creating_asof_May_table.sql
int iterations = 1;
//int target_particle_size = 1;//4 pararell computing
/* Define variables*/
/*Variables that can change according to the start date*/
#define PSC 45 //Planned start of calving, set as 15 July which is 45 days
#define PSM 127 //Planned start of mating, set as 07 October 12 weeks from PSC
int dry_day = 364; //May 31st
int r1_initial_age = 304; //As of June 1st, they are max 304 days old


/*Other parameters*/
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

#define weaning_wks 10 
#define sim_years 10
#define prop_extra_animal_per_year 0.25
int num_extra_animal_per_year;
#define column_cull_empty 1
#define column_sell_empty 2
/*Define tables imported*/
int num_cull_sell_steps = 62 ; //62 steps for culling and death
int num_mortality_steps = 19 ;
int n_column_cull_sell_rate = 5; // days, cull_rate_empty, sell_rate_empty,cull_rate,sell_rate
int n_column_mortality_rate = 2; // day and mortality
int temp_num_animal;

//char CullSellRatesFile[] = "E:/Do/CullSellRatesFile2.csv" ;
//char MortalityFile[] = "E:/Do/mortality.csv" ;
//char TestScheduleFile[] = "E:/Do/TestSchedule.csv";

char CullSellRatesFile[] = "D:/ABCSMC/CullSellRatesFile2.csv" ;
char MortalityFile[] = "D:/ABCSMC/mortality.csv" ;
char TestScheduleFile[] = "D:/ABCSMC/TestSchedule.csv";


/*Define general parameters*/

/*Parameters only used if the number of animals to be recorded*/
//int num_column_NumGrpAnimal = 15000;
//char NumberAnimalDataFile[] = "C:/ABCSMC/NumberAnimalDataFile.csv" ;
//char TestResultFile[] = "C:/ABCSMC/TestResultFile.csv" ; //if to export the test results per farm

/* Define calf related parameters*/


/* Define heifer related parameters*/ 
int first_heat_rand ;
int pregnant_status;
int target_calf_met = 0;
int target_milk_met = 0;
/* Define lactating and non-lactating related parameters*/


int interval_first_heat = time_first_heat_max - time_first_heat_min + 1 ;
int interval_heat = interval_heat_max - interval_heat_min + 1 ;
int calving_duration = 70; //9 to 12 weeks, now set as 10 weeks

/*Define an initial age distribution of each age in milking group*/
double prop_age_2 = 0.2 ;
double prop_age_4 = 0.335;
double prop_age_7 = 0.3 ;
double prop_age_8 = 0.165 ;

/* Define simulation related parameters*/
int sim_days = sim_years*365;
int test_freq, last_test_date;
int ELISA_type ;


/*=======Define DISEASE SPECIFIC Column PARAMETERS==========*/
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
int test_date ;
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
#define bdisc 12
#define elisa_type 13


#define test_pcr_live_nomilk 3
#define test_pcr_slaughter_nomilk 4
#define test_pcr_milk 5

/*Parameters associated with mastitis or other conditions being detected*/
#define rate_detect_clinical_mastitis 0.95
#define length_CM_treatment 6 //treatment + withholding
#define rate_CM_2 0.031
#define rate_CM_3_4 0.016
#define rate_CM_5_7 0.031
#define rate_CM_8 0.047 //all rate /cow-day


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
   
      int NB_mastitis_status ; //non-mastitis = 0, clinical mastitis = 1 DUE to non-bovis
      int mastitis_detected ; //0 not detected 1 detected
      double mastitis_rate ; //IR CM due to non M bovis - 0 if not milking or if animals already developed CM due to non-bovis
      int mbovis_sero_status ;//0 neg 1 pos (weak pos considered to be positive)
      int mbovis_elisa_status; //0 neg 1 detected
      int mbovis_pcr_status; //0 neg 1 detected
    //  int lameness_status ; //no lame 0 lame 1 detected 2
      int pregnant_status; // non-pregnant = 0 preganant 1, not differentiating AI or bull
      int calving_date; //date of calving (Day xx)
      int non_colostrum ;//0 if it's in colostrum and 1 if it's not in colostrum
      int next_heat_date; //date of the next heat (Day xx)
      //int num_births; // record how many births it already gave to: deleted as now it has parity var
      
      int present; //whether this animal exists on the farm or already removed
      long long current_pro_id ;
      int index_cull_sell ;
      int index_mortality ;
      int tested ;
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
void test() ;
void test_bulk();
int write_test_result() ;
int write_disease_result() ;


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
int bovis_index;

int* num_bobby;
int var_bobby = 0;

double markov_rate ;
int random_number ;


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
int ncol_test_table = 10*2 ;//stores test positive and test conducted
double shedder,total ;
int current_particle = 0;
/*===============================================================================*/
/*Set up the input parameters*/
//int *risk_date_par[41];
//int *inf_par[5];
//double *scale_par[2];
//double *beta_par[2];
//double *Se_par[4];
//double *Sp_par[2];
//double *summary_statistics[10] ;

//Herd size is a deterministic number from ARDB
int herd_size_vector[] = {
892,
280,
400,
750,
368,
936,
730,
750,
600,
900,
157,
1303,
579,
585,
580,
700,
2020,
620,
783,
817,
443,
504,
300,
625,
1705,
550,
570,
787,
530,
370,
500,
230,
464,
1250,
150,
650,
750,
1040,
565,
800,
220		
};


/*===========FUNCTION TO SIMULAE A DISEASE STARTS==========================================================================================*/
/*Main starts from here*/
void herd_simulation(
int *risk_date_par,//disease entry date for each farm as a prior
int *inf_par, //inverse of bv_r1 to r4, and sero conversion as 5
double *scale_par,//bv_f(scale) as item 1, dilution as 2
double *beta_par, //1 is cow to cow beta_mm, 2 is cow to calf beta_mc
double *Se_par, //ELISA bio, ELISA ID, PCR live, PCR slaughter
double *Sp_par, //ELISA bio, ELISA ID
double *summary_statistics, //result of tests
int* seed,
int* duration,
int* in_herd_size,
int* init_bovis_status //status of animals that are introduced 1 Exposed and 2 Subclinical
){


//Rename some of input parameters for readability
dilution = scale_par[1];//2nd item
bv_f = scale_par[0];
rate_sero_conversion = 1/(double)inf_par[4];
Se_bio = Se_par[0] ;
Sp_bio = Sp_par[0];
Se_ID = Se_par[1] ;
Sp_ID = Sp_par[1] ;
PCR_live = Se_par[2] ;
PCR_slaughter = Se_par[3] ;

//Setting seed - seed is provided as an input for a reproducibility
//srand((unsigned)time(NULL));
srand((unsigned) *seed);

//Define some more parameters
int mbovis_status ;
int id_lact_group2 = id_lact_group;
int id_calf_group2 = id_calf_group;
num_culled = &var_cull;
num_sold = &var_sold;
num_death = &var_death;
num_bobby = &var_bobby;
struct animal_node* fake_animal;
struct animal_node* animal_to_remove;
fake_animal = (struct animal_node*)malloc(sizeof( struct animal_node ));
fake_animal->present=1 ;
animal_to_remove = (struct animal_node*)malloc(sizeof( struct animal_node ));

//printf("Starts");	

//Allocate dynamic memory	
double** List_mng_status =  (double**)malloc( sizeof(double *) *(id_bull_group+1)); //modify if more than one herd exists
double **cull_sell_rate = (double**)malloc( sizeof(double *) *num_cull_sell_steps );
double **mortality = (double**)malloc( sizeof(double *) *num_mortality_steps );
//double **NumGrpAnimal = (double**)malloc( sizeof(double *) *(id_bull_group+2) );
double **table_CM_rate = (double**)malloc(sizeof(double*)*8);// 0-7 parity
int *milk_numbers = (int*)malloc(sizeof(int)*(column_cli_D+1)) ; //table for milk
int** GroupSampleSize = (int**)malloc(sizeof(int*)*(id_dry_group+1)) ;
//TASK_38
int** test_result_table = (int**)malloc(sizeof(int*)*num_total_herd) ;
int** test_schedule =  (int**)malloc(sizeof(int*)*nrow_testschedule) ;
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
/*Those not require free*/

//Baseline mastitis rate table: age-dependent
for(i = 0 ; i < 8; i++)
{
	//table row represents age
	table_CM_rate[i] = (double*)malloc( sizeof(double) * 2);//2 columns: [0] Dry period or under-treatment [1]During lactation period
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

//Allocate memory for cull and mortality tables
for(i = 0; i < num_cull_sell_steps; i++)
{
	cull_sell_rate[i] = (double*)malloc( sizeof(double) * n_column_cull_sell_rate);
	
}
for(i = 0; i< num_mortality_steps; i++)
{
	mortality[i] = (double*)malloc( sizeof(double) * n_column_mortality_rate);
}

//now read the cull_sell_rate table	
read_cull_sell_rate(CullSellRatesFile,cull_sell_rate,num_cull_sell_steps) ;	
read_mortality(MortalityFile,mortality,num_mortality_steps) ;

//Set up List_mng_status table
for(i = 0; i < id_bull_group + 1 ; i++)
{
	List_mng_status[i] = (double*)malloc( sizeof(double) * (column_sum_detection+1));//num column i
}
//Set up NumGrpAnimal table if the number of animal in each management group to be recorded
//for(i = 0; i < id_bull_group + 2 ; i++)
//{
//	NumGrpAnimal[i] =  (double*)malloc( sizeof(double) * num_column_NumGrpAnimal);
//}

//Table to store the number of animals to test either as live or slaughter
for(i=0; i<=id_dry_group; i++)
{
	GroupSampleSize[i] = (int*)malloc(sizeof(int)*2) ; //number of column is 2: slaughter and alive
}



//Initialise test_result_table to store test results in simulation
for(i = 0;i<num_total_herd;i++ )
{
	test_result_table[i] = (int*)malloc( sizeof(int) * ncol_test_table);
	for(j=0;j<ncol_test_table;j++)
	{
	
	test_result_table[i][j]=0;	
	}
}
int ite,farm;

//these parameters necessary to change the proportion of calves to keep when they test bobby calves
double calf_keep_prop2;
double	calf_female_prop2;
/*===================================================================================================================*/
/*=======Starting particle===========================*/

/*=================Starting for each farm====================================*/

for(farm = 0; farm < num_total_herd; farm++)
{
//printf("farm is %d\n",farm);
target_calf_met = 0;
target_milk_met = 0;
//(*farm_done)++;
herd_id = farm ;
calf_keep_prop2	 = calf_keep_prop;
calf_female_prop2 = calf_female_prop;
///need to reset the value for each farm
//herd_size = herd_size_vector[farm] ;
herd_size = *in_herd_size;
//printf("HSIZE is %d\n",herd_size);
risk_date = risk_date_par[farm];


//Set up animal_node_pointer and associated parameters
num_extra_animal_per_year = (int)herd_size*prop_extra_animal_per_year;
int length_animal_pointer = herd_size*2 + num_extra_animal_per_year*sim_years ;	
struct animal_node **animal_node_pointer = malloc( sizeof(struct animal_node*) * length_animal_pointer);
int test_nth = 1;

/*Start of iteration: For ABC-SMC it does only 1 iteration per farm*/
for(ite=0;ite<iterations;ite++)
{
int current_akey = 0 ;
double today_date = 0 ;
int year = 0;
int calendar_day = 0;
int parity = 0;
int first_bobby_born = 0;
int index_column = 0;
int next_non_markov_date = 0;

/*===========RESET TABLES IN EACH ITERATION=========================*/

	for(i=0; i<=column_cli_D; i++)
	{
	milk_numbers[i] = 0;
	}


	for(i = 0; i < id_bull_group + 1 ; i++)
	{

	for(j = 0; j < column_sum_detection+1; j++)
	{
		List_mng_status[i][j] = 0;
	}
	
	}

//Reset NumGrpAnimal if number of animals in each management to be recorded
//for(i = 0; i < id_bull_group + 2 ; i++)
//{
//	for(abc = 0 ;abc <num_column_NumGrpAnimal; abc++)
//	{
//		NumGrpAnimal[i][abc] = 0;
//	}
//}

	
//	printf("Examine6\n");



//printf("File read");
double sum_age_prop_4 = prop_age_2 + prop_age_4 ;
double sum_age_prop_7 = prop_age_7 + sum_age_prop_4 ;


//Initialise the linked list
for(i = 0; i < sim_days; i++)
                {
                event_day[i] = NULL;
                }
//linked list for events
for(i = 0; i < id_bull_group+1; i++)
      	  	{
     	      FarmGroupList[i] = NULL; // linked list for animals
			      	    }
		      	    
/*==========FROM HERE FARM SPECIFIC========================================================================*/			      	    

var_cull = 0;
var_death = 0;
var_sold = 0;
var_bobby = 0;			      	    

	   
/*======ADD TEST EVENT======================================================*/
//can I pass test table information from R?
//currently not because .C function limits the number of parameters that can be passed
// Hence each time creating test schedule, which may be waste of computational time
//Ideally, use Rcpp 

//Read from test_schedule table and add testing events on this farm as an event
test_freq = 0;
if(pre_defined_test==1)
{
	

	
	//Set up test_schedule table which stores the actual testing information and read it from csv
	for(i=0;i<nrow_testschedule;i++)
		{
		test_schedule[i] = (int*)malloc( sizeof(int) * ncol_testschedule);
		}
		read_testschedule(TestScheduleFile,test_schedule,nrow_testschedule) ;
	

	for(i=0;i<nrow_testschedule;i++)
	{
		if(test_schedule[i][0]==herd_id)
		{
		test_date = test_schedule[i][1] ;//TASK_40
		//add some indicator to go back to this table
		new_event = (struct event_node*)malloc(sizeof( struct event_node ));
	    new_event->event_type = 12;//testing
	    new_event->akey = i ; //use akey as an index for test_schedule table
	    new_event->next_node = NULL ;
	    new_event->animal = fake_animal ;
	    add_event_node(event_day,test_date, new_event) ;
		test_freq++; // counter for total testing on this farm	
		
		}
		if(test_schedule[i][0]>herd_id)
		{
			break;
		}
	}	
	//printf("Total test is %d",test_freq) ;
	last_test_date = test_date ; //record last test date so that a simulation can end on this date
}
else
{
	last_test_date = PSC+ risk_date + *duration ;
}


	   

/*======ADD RISK EVENT==========================================*/
//risk event date will be a parameter, random sample from some distribution
//TASK_40
	
	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
    new_event->event_type = 16;//Risk events
    new_event->next_node = NULL ;
    new_event->animal = fake_animal;
    add_event_node(event_day,risk_date+PSC, new_event) ;     
	
		    


/*====ADDING INITIAL ANIMALS====================================================================*/

/*R1 HEIFER========================================================*/
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
	animal_node_pointer[current_akey]->dried = 1 ; //Simulation starts on 1 June, so it's a dry period
	animal_node_pointer[current_akey]->present = 1 ;
	animal_node_pointer[current_akey]->tested = 0;
	animal_node_pointer[current_akey]->previous_node = NULL;
	animal_node_pointer[current_akey]->next_node = NULL;
	//animal_node_pointer[current_akey]->lameness_status = 0;
	animal_node_pointer[current_akey]->NB_mastitis_status = 0;
	animal_node_pointer[current_akey]->mastitis_detected = 0;
	animal_node_pointer[current_akey]->mbovis_status = Sus ;
	animal_node_pointer[current_akey]->mbovis_sero_status = 0;
	animal_node_pointer[current_akey]->mbovis_elisa_status = 0;
	animal_node_pointer[current_akey]->mbovis_pcr_status = 0;
	
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
	new_event->animal = animal_node_pointer[current_akey]; 
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
	animal_node_pointer[current_akey]->tested = 0;
	animal_node_pointer[current_akey]->previous_node = NULL;
	animal_node_pointer[current_akey]->next_node = NULL;
	//animal_node_pointer[current_akey]->lameness_status = 0;
	animal_node_pointer[current_akey]->NB_mastitis_status = 0;
	animal_node_pointer[current_akey]->mbovis_status = Sus ;
	animal_node_pointer[current_akey]->mbovis_sero_status = 0;
	animal_node_pointer[current_akey]->mbovis_elisa_status = 0;
	animal_node_pointer[current_akey]->mbovis_pcr_status = 0;
	animal_node_pointer[current_akey]->mastitis_detected = 0;
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
	//printf("num is %lf", List_mng_status[mng_group][0]);
	List_mng_status[mng_group][0] ++;
	//printf("num now is %lf", List_mng_status[mng_group][0]);
	
	List_mng_status[mng_group][mbovis_status+5] ++;
}

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
		animal_node_pointer[current_akey]->tested = 0;
	animal_node_pointer[current_akey]->previous_node = NULL;
	animal_node_pointer[current_akey]->next_node = NULL;
	//animal_node_pointer[current_akey]->lameness_status = 0;
	animal_node_pointer[current_akey]->NB_mastitis_status = 0;
	animal_node_pointer[current_akey]->mbovis_status = Sus ;
	animal_node_pointer[current_akey]->mbovis_sero_status = 0;

	animal_node_pointer[current_akey]->mastitis_detected = 0;
	animal_node_pointer[current_akey]->mbovis_elisa_status = 0;
	animal_node_pointer[current_akey]->mbovis_pcr_status = 0;
	
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
	//	printf("Age 8");	
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
	//	printf("num is %lf", List_mng_status[mng_group][0]);
	List_mng_status[mng_group][0] ++;
	//	printf("num is %lf", List_mng_status[mng_group][0]);
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
//if(farm==36)
//{
//printf("Mixed added");	
//}

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
int* dry_season ;
int var_dry = 1;
int var_wet = 0;
dry_season = &var_dry ;
int k = 0;

/*Reset waste milk*/
double vwm = 0;
volume_waste_milk = &vwm; //they are also reset to 0 when the last calves get weaned
double dmv = 0;
density_mbovis_waste_milk = &dmv;
double dh = 0;
density_healthy = &dh;
intensity_mc = 0;

/*============================================================================================
Day procedes
=============================================================================================*/
while(today_date<sim_days && today_date<last_test_date+1)
{
	year = (int)floor(today_date/365);
	if(year==3)
	{
		if(farm==1||farm==34)
		{
		calf_keep_prop2	 = 1;
		calf_female_prop2 = 1;	
		}
		
	}
	/*Farm 1 and 34 are tested for a large number of calves in Year 3 - probably all calves including bobby calves are tested.
	So making these farms to keep all calves born in this particular year*/
//printf("Day is %f\n",today_date);
/*====START OF DAY========================================================*/
if(next_non_markov_date<sim_days-1)
{
next_non_markov_date = ceil(today_date);	



/*UPDATE NEXT NON MARKOV EVENT DATE*/
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
	 	//printf("NonM is %d",next_non_markov_date);
	 	current_event = event_day[next_non_markov_date];
	 //	printf("Type is %d",current_event->event_type) ;
	 	//trial
	 	current_animal = current_event->animal; 	
	 	previous_event = current_event;
	 	while(current_event->animal->present==0)//what did I do for fake animals?
	 	{
	 
	 		current_event = current_event->next_node;
	 

	 		free(previous_event);
	 		previous_event = current_event;
	 		if(current_event==NULL)
	 		{
	 		//	printf("current event becomes null");
	 			break;
	 			
			 }
		 }
	 	event_day[next_non_markov_date] = current_event ;
	 	//this is right - updating the root
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
	
			if(current_event->next_node->next_node!=NULL)
			{
			
				next_event = current_event->next_node->next_node;

				free(current_event->next_node);
			
				current_event->next_node = next_event ;
		
			}
			else
			{
		
					if(farm==39)

				free(current_event->next_node);
				current_event->next_node = NULL;
		
				break;
			}
			
		}
		current_event = current_event->next_node ;
	}
	
	 	
	 }
            } // now get next non-markov date
    } //if(next_non_markov_date<sim_days-1) ends
/*UPDATING NEXT NON MARKOV DATE DONE*******************************************************************/
 
 /*DETRMINE WHETHER MARKOV EVENT HAPPENS BEFORE THE NEXT MARKOV EVENT*/
  updated_date=update_markov_date(today_date,List_mng_status,cull_sell_rate,
  mortality, FarmGroupList,next_non_markov_date,
  num_culled,num_sold,num_death,bv_f,rate_sero_conversion,id_lact_group2,id_calf_group2,num_bobby,
  dry_season,milk_numbers, event_day,inf_par,beta_par,intensity_mc, farm) ;

//printf("After update");
/*==============MARKOV DID NOT HAPPEN=====================================================================*/
  if (updated_date==next_non_markov_date) // this means markov event did not happen
     {//LOOP NM1
     calendar_day = next_non_markov_date - year*365 ;
	//printf("Loop started\n");
     current_event = event_day[next_non_markov_date];
     while(current_event!=NULL)
     {

/*======== CALVING==========================================*/
if(current_event->event_type==1  && current_event->animal->present==1)
{
    	
		    current_animal = current_event->animal ;
//		    if(current_animal->pregnant_status==0)
//		    {
//		    	printf("Wait!She is not pregnant!\n");
//		    	system("pause");
//			}
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
     		
     	/*UPDATE colostrum status and add event to switch colostrum status
     	UPDATE List_mng_status bovis numbers
     	whether this animal remains in the herd or culled, need to minus from the original group*/
     	mbovis_status = current_animal->mbovis_status ;
     	List_mng_status[current_grp][mbovis_status+5]--;
     	List_mng_status[current_grp][0]--;
			
		//Update the number in the lactation group
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
     		

			current_animal->non_colostrum = 0;//now in colostrum period
			milk_numbers[bovis_index*3+current_animal->non_colostrum]++; 
     	
     	/*Add event when this animal comes back to milking after having colostrum period	*/
			if(next_non_markov_date+withhold_colostrum<sim_days)
			{
				new_event = (struct event_node*)malloc(sizeof( struct event_node ));
			    new_event->event_type = 13;//switching colostrum
			    new_event->akey = current_animal->akey;
			    new_event->animal = current_animal; //let's see if this works
			    new_event->next_node = NULL ;
				add_event_node(event_day,next_non_markov_date+withhold_colostrum, new_event) ;
			}
			
		 			current_animal->pregnant_status = 0;//not pregnant anymore
     				current_animal->sum_markov_rate = 
					cull_sell_rate[current_animal->index_cull_sell][column_cull_empty+2*pregnant_status]+
					cull_sell_rate[current_animal->index_cull_sell][column_sell_empty+2*pregnant_status] + 
					mortality[current_animal->index_mortality][1] ;
     	
		 			remove_animal_group(FarmGroupList,current_grp,current_animal);
     				current_animal->group = id_lact_group ; //change the group from dry to lactation
     				add_animal_group(FarmGroupList,id_lact_group,current_animal);
     	/*All animals calved before the milking number 
		then this animal will receive a submission to calve next year, so heat happens
		if the number of animals that calved is sufficient then this animal won't be inseminated 
		and sold or culled - here we assume they will be culled */	
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
     	
     	
		/*If the number of calves to keep is sufficient then new calves become bobby*/	
			
     	if(List_mng_status[0][0]>herd_size * calf_keep_prop2)
     	{
     	(*num_bobby) ++; 
     	/*late calvers still keep on calving but calves won't be kept.*/
     	
		}
		else if((double)rand()/(double)RAND_MAX<=calf_female_prop2) //proportion of female and male 0.5
		 {
		 /* This is female animal and going to be kept on farm*/
            	animal_node_pointer[current_akey] =malloc(sizeof(struct animal_node)) ;
		        animal_node_pointer[current_akey]->akey = current_akey;	
		        animal_node_pointer[current_akey]->age_day = 0; //just born
		        animal_node_pointer[current_akey]->parity = 0;
		        animal_node_pointer[current_akey]->mastitis_detected = 0;
		        animal_node_pointer[current_akey]->dried = 1;
			    animal_node_pointer[current_akey]->tested = 0;
		        animal_node_pointer[current_akey]->mastitis_rate = table_CM_rate[0][0];
		        animal_node_pointer[current_akey]->NB_mastitis_status = 0;
				animal_node_pointer[current_akey]->mbovis_status = 0;
				animal_node_pointer[current_akey]->mbovis_sero_status = 0;
	    		animal_node_pointer[current_akey]->mbovis_elisa_status = 0;
	    		animal_node_pointer[current_akey]->mbovis_pcr_status = 0;
		
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
				
		/*Add weaning event*/
				new_event = (struct event_node*)malloc(sizeof( struct event_node ));
				new_event->event_type = 14; //moving to the next group
				new_event->next_node = NULL ;
				new_event->animal = animal_node_pointer[current_akey];
				if(next_non_markov_date+weaning_wks*7<sim_years*365)
				{
				add_event_node(event_day,next_non_markov_date+weaning_wks*7, new_event) ;//animals move to the next age group
				}
		/*On the first day of calving, declare it's milking season	
		Set the volume_waste_milk and density_mbovis_waste_milk to 0
		Why are we not setting these volumes to 0 when the last calf gets weaned?
		Because there may be testing of waste milk after the weaning*/
			
				if((int)List_mng_status[0][0]==1 && *dry_season==1) /*If this is the first born in the season, decalre milking season starts*/
				{
				//	printf("The first calving date is %d in YEAR %d\n",calendar_day, year);
				
						dry_season = &var_wet;
						*volume_waste_milk = 0;
						*density_mbovis_waste_milk = 0;
						if(next_non_markov_date+1<sim_days)
						{
						new_event = (struct event_node*)malloc(sizeof( struct event_node ));
						new_event->event_type = 17;
						new_event->animal = fake_animal; //let's see if this works
						new_event->next_node = NULL ;
						add_event_node(event_day,next_non_markov_date+1, new_event) ;	
						} /*Updating waste milk event added*/	
						
				}
			
		/*if the number of new calves reaches more than the number of calves 
		then decalre target_calf_met = 1 and stop keeping the new calves
		essentially, the calf that is just born is the last calf to wean*/
		
				if(List_mng_status[0][0]>herd_size * calf_keep_prop2 && target_calf_met==0)
				{
					target_calf_met = 1;
				//	printf("Target calf secured %d in YEAR %d\n",calendar_day, year);
				//	printf("Calendar is %d", calendar_day);
				//	system("pause")	;
					new_event = (struct event_node*)malloc(sizeof( struct event_node ));
					new_event->event_type = 8; //moving to the next group
					new_event->next_node = NULL ;
					new_event->animal = fake_animal;
				//	printf("Lact number is %d",(int)List_mng_status[id_lact_group][0]);
				//	system("pause");
				//	printf("setting up move date") ;
					if(next_non_markov_date+weaning_wks*7<sim_years*365)
					{
						add_event_node(event_day,next_non_markov_date+weaning_wks*7+1, new_event) ;
						//animals move to the next age group when all the calves get weaned
						
				//	printf("Moving event is %d",next_non_markov_date+weaning_wks*7);
				
					}
			//	printf("last calf to add\n");
				}
				
		/*Now set cull/death/sell rate fro animals just born*/
		
				markov_rate = 
				cull_sell_rate[current_index_cull_sell][1]+
				cull_sell_rate[current_index_cull_sell][2] + 
				mortality[current_index_mortality][1] ;
				animal_node_pointer[current_akey]->sum_markov_rate = markov_rate ;
				//List_mng_status[0][1] = List_mng_status[0][1] + markov_rate ;
			
			/*Add new calve to the pointer*/
				add_animal_group(FarmGroupList,0, animal_node_pointer[current_akey]) ;
				
			/*ADD FIRST HEAT of this calf ASSUMING it's 12MO*/
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
			/*Add an event to change culling/death/sale rate change*/
				new_event = (struct event_node*)malloc(sizeof( struct event_node ));
				new_event->event_type = 3;//change in culling index
				new_event->akey = current_akey;
				new_event->animal = animal_node_pointer[current_akey]; //let's see if this works
				new_event->next_node = NULL ;
					if(next_cull_change_date<sim_years*365)
					{
					add_event_node(event_day,next_cull_change_date, new_event) ;
					}
				
				
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
		 else /*this is male and becomes bobby, but still not enough calves born*/
		 {
			 	(*num_bobby)++;
			 	if(first_bobby_born==0)
			 	{
			 	//arrange calf pick up
				//TASK_31
				new_event = (struct event_node*)malloc(sizeof( struct event_node ));
			    new_event->event_type = 15; //moving to the next group
			    new_event->next_node = NULL ;
			     new_event->animal = fake_animal;
			     last_pick_up = next_non_markov_date + 77 ; //11 weeks from teh first bobby born
				if(next_non_markov_date+bobby_pickup<sim_years*365)
			        {
			        add_event_node(event_day,next_non_markov_date+bobby_pickup, new_event) ;//animals move to the next age group
			
			        }
			        first_bobby_born = 1 ;
				}
		 }
		 /*if the number of lact reaches to the target
		 remove all R2 when the calved cows reached the number of milking animals to keep
		 this happens only once per lactation
		 cows can calve after this but not heifers
		 assume all bobbies from R born on this date*/
		 	if((int)List_mng_status[id_lact_group][0]>=herd_size && target_milk_met==0)
     		{
     		//	printf("Day is %d",calendar_day);
     			//remove all R2 if remaining
     				if(List_mng_status[id_r2_group][0]<0) //this is error
     				{
//     			printf("R2 animals below 0\n");
//     			system("pause");
     			//system("pause");
				 	}
				 	else //we don't need more cows - so heifers expected to calve later than this will be removed
				 	{
			//	 	printf("R2 is %d",(int)List_mng_status[id_r2_group][0]);
				 	current_animal = FarmGroupList[id_r2_group];
				 	while(current_animal!=NULL)
				 	{
				 		(*num_bobby)++;
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
				//error checkinbg: ensure id_dry_group>=0
					if(List_mng_status[id_dry_group][0]<0)
	     			{
	     			printf("Dry animals below 0\n");
	     			//system("pause");
	     			//system("pause");
					 }
			/*Declare the number of animals to milk reached the desired number*/
				 target_milk_met = 1 ;
				// printf("Target milking reached\n");
			//	 system("pause");
			 }
//			 if(farm==39)
//		{
//			printf("Calving done\n");
//		}
     	
} //calving done
		 	
/*=============HEAT==============================================================================================================*/
		 if(current_event->event_type==2 && calendar_day<=PSM+7*mating_period&& current_event->animal->present==1)
		 {

			num_heat++;
			current_animal = current_event->animal;
			if(current_animal->pregnant_status==1)
			{
		//	printf("Heat but pregnant grp is %d day is %d",current_animal->group, next_non_markov_date);
		 //	printf("Calving date is %d",current_animal->calving_date);
		 	if(current_animal->present==0)
		 	{
		 		printf("Ghost heat");
			 }
			// system("pause");
			}
		 		if(((double)rand()/(double)RAND_MAX<=submission_prop ) && 
				 current_animal->pregnant_status==0 &&
				 calendar_day>=PSM)	/*if submission is going to happen*/
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
						
				 } /* Submission did not happen and then next heat event will occur*/
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
//				 		 if(farm==39)
//		{
//			printf("Heat done\n");
//		}	 
			 }
/*===============CHANGE IN CULLING INDEX========================================================================*/
if(current_event->event_type==3&& current_event->animal->present==1)
	{

	current_animal = current_event->animal;
		if(current_animal->present==1)
			{
				  	current_index_cull_sell = current_animal->index_cull_sell;	
		  			pregnant_status = current_animal->pregnant_status;
		// printf("index is %d\n",current_index_cull_sell);
		 /*If it exceeds the limit of row, then come back to row 57*/
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
						/*Set up the next change date*/
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
	
	}
/*===============CHANGE IN MORTALITY INDEX=========================================================*/
if(current_event->event_type==4 && current_event->animal->present==1)
{
	//	 printf("event is %d\n",current_event->event_type);	
	current_animal = current_event->animal;
	if(current_animal->present==1)
		{
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
	    }
//		 		 		 if(farm==39)
//		{
//			printf("culli index done\n");
//		}
	   //  printf("Event 4 done");
}
/*=========Scanning==================================================*/     
if(current_event->event_type==7)
	{

		if(FarmGroupList[id_lact_group]!=NULL)
		{
			current_animal = FarmGroupList[id_lact_group];
			while(current_animal!=NULL)
			{
				next_animal = current_animal->next_node;
				mbovis_status = current_animal->mbovis_status ;
				if(current_animal->pregnant_status==0) /*If animal is empty, they will be culled*/
				{
					current_animal->present = 0;
				
					List_mng_status[id_lact_group][0]--;
					List_mng_status[id_lact_group][mbovis_status+5]--;
//				if(	List_mng_status[id_lact_group][mbovis_status+5]<0)
//				{
//					printf("Warning: status below 0 at scanning adult\n");
//					system("pause");
//				}
				/*Update milk_numbers table*/
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
						if(milk_numbers[bovis_index*3+current_animal->non_colostrum]<0)
						{//error checking
					//	printf("Scanning-2\n");
					//	system("pause");
						}
					}
					else
					{
						milk_numbers[bovis_index*3+2]--; 	
						if(milk_numbers[bovis_index*3+2]<0)
						{//error checking
					//	printf("Scanning-1\n");
					//	system("pause");
						}
					}
			 	remove_animal_group(FarmGroupList,id_lact_group,current_animal);
				(*num_culled)++;
			//	printf("empty culled at scanning\n");	
				}
				current_animal = next_animal ;
			}
		}

		/*Scen R2 too*/
		if(FarmGroupList[id_r2_group]!=NULL)
		{
			current_animal = FarmGroupList[id_r2_group];
			
			while(current_animal!=NULL)
			{
				next_animal = current_animal->next_node;
				if(current_animal->pregnant_status==0)
				{
			
				current_animal->present = 0;
				mbovis_status = current_animal->mbovis_status ;
				List_mng_status[id_r2_group][0]--;
				List_mng_status[id_r2_group][mbovis_status+5]--;

				//TASK_33
				(*num_culled)++;
				remove_animal_group(FarmGroupList,id_r2_group,current_animal);	
				}
				current_animal = next_animal ;
			}
		}

	}

/*===========Move to next subgroup========================================*/
	if(current_event->event_type==8)
	{
	intensity_mc = 0 ; // as all calves get weaned, cow to calf transmission stops happening 	
	
	/*R2 animals should remain here because it means this R2 didn't calve but didn't get removed.
	So just in case clear R2 animal group*/
		if(FarmGroupList[id_r2_group]!=NULL)
		{
			current_animal = FarmGroupList[id_r2_group];
	
			
			while(current_animal!=NULL)
			{
				next_animal = current_animal->next_node;
	
				remove_animal_group(FarmGroupList,id_r2_group,current_animal);
				current_animal->present = 0;
				List_mng_status[id_r2_group][0]--;
				
				(*num_culled)++;
				current_animal = next_animal ;
			}
		//	printf("Now R2 is %d",(int)List_mng_status[id_r2_group][0]);
	//		system("pause");
			FarmGroupList[id_r2_group] = NULL;
			
				for(i=column_Sus; i<= column_Exp; i++)
				{
				List_mng_status[id_r2_group][i] = 0; //all R2 gone
				
				}
			
		}
	/*If R2 is empty then move R1 to R2 group*/
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
			temp_num_animal = (int)List_mng_status[id_r1_group][0]; //R1 heifers
			List_mng_status[id_r2_group][0] = (double)temp_num_animal;
	//		printf("R2 heifer is %d",temp_num_animal);
		/*Now move weaned animals to R1*/	
			if(FarmGroupList[id_calf_weaned]==NULL)
				{
//				printf("There is no animals in calf:YEAR %d\n",year) ;
//				printf("Calf before wean is %d\n",(int)List_mng_status[id_calf_group][0]);
//				printf("Farm is %d Day is %d\n",farm,next_non_markov_date);
//				system("pause");
				}
			else
			{
				if(List_mng_status[id_calf_weaned][0]<=0)
				{
//					printf("FarmGroupList has calves but not in List_mng");
//					system("pause");
				}
				else
				{
					
					temp_num_animal = (int)List_mng_status[id_calf_weaned][0];
					List_mng_status[id_r1_group][0] = (double)temp_num_animal;
				//	printf("R1 heifer is %d",temp_num_animal);
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
			/*transfer List_mng_status information other than the total number of animals*/
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
{

	target_calf_met=0;
	target_milk_met=0;
	first_bobby_born = 0 ;
	dry_season = &var_dry ;
	int num_lact = 0;
	int num_dry = 0;
	int num_culled_nonpregnant = 0;
	int num_pregnant_at_dry = 0;
	for(i=0;i<=column_cli_D; i++)
		{
		milk_numbers[i] = 0 ; //reset to 0, if think about between dried animals transmission need to change this
		}
	
	
		if(FarmGroupList[id_lact_group]==NULL)
		{
			//printf("There are no milking animals A!\n") ;
			//printf("Precisely it's %d\n",(int)List_mng_status[id_dry_group][0]);
			//system("pause");
		}
		else
		{
			current_animal = FarmGroupList[id_lact_group] ;
			if(FarmGroupList[id_dry_group]!=NULL)
			{
			//	printf("There are already dry animals!\n") ;
			//	printf("num is %d",List_mng_status[id_dry_group][0]);
				//next_animal = FarmGroupList[id_dry_group];
				//while(next_animal!=NULL)
				//{
				//	printf("Akey and grp is %lld, %d status %d\n",next_animal->akey,next_animal->group,next_animal->pregnant_status);
				//	printf("Present is %d",next_animal->present);
				//	next_animal = next_animal->next_node;
				//}
			//system("pause");
			}
			else
			{
				while(current_animal!=NULL)
				{
					current_animal->group = id_dry_group;
					//here cull animals if not-pregnant
					next_animal = current_animal->next_node;
					if(current_animal->pregnant_status==0 || current_animal->mastitis_detected==1)//if not pregnant, remove
					{
						//Assumption_4
					//	printf("akey is %lld",current_animal->akey);
						
						current_animal->present = 0 ;
						mbovis_status = current_animal->mbovis_status ;
						remove_animal_group(FarmGroupList,id_lact_group,current_animal);
						(*num_culled)++;
						num_culled_nonpregnant++;
					//	printf("Non-pregnant culled");
					//	system("pause");
					//-1 this group because we need number of remaining lact animals
						List_mng_status[id_lact_group][0]--;
						List_mng_status[id_lact_group][mbovis_status+5]--;
						if(	List_mng_status[id_lact_group][mbovis_status+5]<0)
						{
							printf("Warning: Status below 0 at dryoff\n");
						}
						//also need to change bovis status related number
					}
					else
					{
			//		num_lact++;
					parity = current_animal->parity ;
					//Assumption_5
					current_animal->mastitis_rate = table_CM_rate[parity][0] ;//no mastitis in dryoff
					current_animal->dried= 1 ;	
					}
					current_animal = next_animal ;
				
				}
				current_animal = FarmGroupList[id_lact_group] ;//get the first animal again
				FarmGroupList[id_dry_group] = current_animal; /*All lact animals move to dry animal*/
	
				FarmGroupList[id_lact_group]=NULL;
				temp_num_animal = (int)List_mng_status[id_lact_group][0];
				List_mng_status[id_lact_group][0]=0;
				List_mng_status[id_lact_group][1]=0;
				List_mng_status[id_dry_group][0]=(double)temp_num_animal;
				
				for(i = column_Sus; i<= column_Cli_shed; i++)
					{
					List_mng_status[id_dry_group][i] = List_mng_status[id_lact_group][i];
					List_mng_status[id_lact_group][i] = 0;
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

		//some clinical mstitis will be detected with a certain probability
		current_animal = current_event->animal ; //animal that got treated
		if(current_animal->present==1)
		{
			if(FarmGroupList[id_lact_group]==NULL)
		    {
				printf("There are no milking animals B!\n") ;
				if(current_animal->dried==1)
				{
					printf("This animal was dried\n");
				}
		    }
		    else if(current_animal->dried==1)
		    {
			
				//Assumption_7: CM due to non-bovis cured as soon as dried-off
				current_animal->mastitis_detected = 0;
				current_animal->NB_mastitis_status=0; //CM due to non-bovis cured
				current_animal->mastitis_rate = 0;
			
	
		    }
		    else
		    {
		
				if(current_animal->NB_mastitis_status==0) //if CM was due to bovis
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
		//	    			if(milk_numbers[bovis_index*3+2]<0)
		//					{
		//						printf("T-1\n");
		//						system("pause");
		//					}
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
	
	}
	
	
/*===============EVENT 12: SURVEILLANCE TESTING: USE REAL DATA=============================================================================*/
//TASK_16
//New Task: Add indicator whether mastitis cow was tesed for PCR, then further indicate if any of these mastitis cows tested PCR pos
if(current_event->event_type==12) //testing for MBOVIS 
{

//printf("Test is %d",test_nth);

k = (int)current_event->akey;
//test bulk or dicard
	if(test_schedule[k][btank]==1)//test_schedule is a table with all tests
	{
		shedder = bv_f*milk_numbers[4]+milk_numbers[7]; //effective shedding amount
	//printf("shedder %f\n",shedder);
		total = (milk_numbers[1]+milk_numbers[4]+milk_numbers[7])*av_milk_production;
		if(total>0)
		{
			*density_healthy = shedder/total ;	
		}
		else
		{
			*density_healthy = 0;
		}
		//printf("total %d\n",total);
		test_bulk(0, density_healthy,dilution,PCR_live, test_result_table,farm);
	}
	if(test_schedule[k][bdisc]==1)
	{
		test_bulk(1,density_mbovis_waste_milk,dilution,PCR_live, test_result_table,farm);
	}
/*======Determine the observed sample size for slaughter and live samples.
Get the largest sample size for each sample type (PCR -nomilk or PCR - milk OR ELISA=====*/
//========CALF======================================================
group = 0;	
ELISA_type = test_schedule[k][elisa_type] ;//which ELISA tests?
//Determine the sample size for calf
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
//Determine sample size for adult
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
//Ideally, all required samples come from non-weaned calves and adults
//However, because of stochasticity and not knowing the precise farm manegement
//in some occasions, we may need to sample from weaned calves.
//Below, decide how many samples to take from each management


for(i=0;i<=id_dry_group;i++)
{
//	printf("Determining sample size for group%d\n",i);
	GroupSampleSize[i][slaughter] = 0;//initialise
	GroupSampleSize[i][live] = 0;
	if((target_slaughter_calf +target_live_calf >0) && (i <= id_r2_group)) //if need calf sampling
	{
		if(target_slaughter_calf<=(int)List_mng_status[i][0])
		{
		//	printf("p1 i is %d\n", i);
		//	printf("target_slaughter_calf is %d, size is %f\n",target_slaughter_calf,List_mng_status[i][0]);
			GroupSampleSize[i][slaughter] = target_slaughter_calf; //all calf samples are from id_calf_group
		//	printf("GroupSampleSize is %d\n",GroupSampleSize[i][slaughter]);
			target_slaughter_calf = 0;
			if(target_live_calf <= (int)List_mng_status[i][0] - GroupSampleSize[i][slaughter] )
			{
			GroupSampleSize[i][live] = target_live_calf;
			target_live_calf = 0;
			}
			else
			{
			GroupSampleSize[i][live] =	(int)List_mng_status[i][0]- GroupSampleSize[i][slaughter];
			target_live_calf = target_live_calf - GroupSampleSize[i][live];
			}
		
		}
		else
		{
		//	printf("p2\n");
		//	printf("target_slaughter_calf is %d\n",target_slaughter_calf);
			GroupSampleSize[i][slaughter] = (int)List_mng_status[i][0];
		//		printf("GroupSampleSize is %d\n",GroupSampleSize[i][slaughter]);
			target_slaughter_calf = target_slaughter_calf - (int)List_mng_status[i][0] ;
			GroupSampleSize[i][live] = 0 ;
			//if require_calf_slaughter is <0 then all slaughter sample is collected but not live
			//if >0 then  List_mng_status[i][0] were all used for slaughter and not for live
		}
	}
	else if(i>=id_lact_group && target_adult>0)
	{
		if(FarmGroupList[i]!=NULL)
		{
			if(target_slaughter_adult+target_live_adult<=(int)List_mng_status[i][0])
			{
				GroupSampleSize[i][slaughter] = target_slaughter_adult ;
				GroupSampleSize[i][live] = target_live_adult ;
			}
			else if(target_slaughter_adult <= (int)List_mng_status[i][0])
			{
				GroupSampleSize[i][slaughter] =target_slaughter_adult;
				GroupSampleSize[i][live] = (int)List_mng_status[i][0] - target_slaughter_adult;
			}
			else
			{
				GroupSampleSize[i][slaughter] = (int)List_mng_status[i][0];
				GroupSampleSize[i][live] = 0 ;
			}
			target_adult = 0;
		}
	}

}
/*=========ENDS: Determined how many samples to collect for slaughter and live in each management group======*/

/*=======Do adult and calf sampling as a loop=====================================================*/
for(i=0; i<=id_dry_group; i++)
{

	//TASK_37
N_PPS = 0; //number of remaining previous positive or sick animals: PPS stands for 'previous positive or sick'
N_H = 0; // number of remaining healthy animals
counter_slaughter = 0;//this counts the number of sample collected for slaughter in a given mangement grp
counter_live = 0;//same for live samples
counter_interval_sick = 0;
counter_interval_healthy = 0;
//start sampling
//so this counter needs to be defined only once each in calf and adult sampling
	if(i==id_calf_group)
	{
		
	group = 0 ;
	counter_slaughter_pcr_nomilk = test_schedule[k][slaughter_pcr_nomilk_calf] ;
	counter_slaughter_ELISA = test_schedule[k][slaughter_elisa_calf] ;	
	counter_live_pcr_nomilk = test_schedule[k][live_pcr_nomilk_calf];
	counter_live_ELISA = test_schedule[k][live_elisa_calf] ;
	}
	else if(i==id_lact_group)
	{
	group = 1 ;
	counter_slaughter_pcr_nomilk = test_schedule[k][slaughter_pcr_nomilk_adult] ;
	counter_slaughter_ELISA = test_schedule[k][slaughter_elisa_adult] ;
	counter_live_pcr_nomilk = test_schedule[k][live_pcr_nomilk_adult];
	counter_live_pcr_milk = test_schedule[k][live_pcr_milk];
	counter_live_ELISA = test_schedule[k][live_elisa_adult] ;
	}
	
/////If at least one sample to be taken from this management group
if(GroupSampleSize[i][slaughter]+GroupSampleSize[i][live]>0)//this is not true if all calf before wenaed depleted
{//if sampling from this group
//printf("sampling from %d\n",i);

//Count how many sick animals exist in this group
current_animal = FarmGroupList[i];
	while(current_animal!=NULL)
	{
		if(current_animal->mbovis_elisa_status==1||current_animal->mbovis_pcr_status==1||current_animal->mastitis_detected==1)
		{
		N_PPS++;	
		}
		current_animal=current_animal->next_node;
	}
	N_H = (int)List_mng_status[i][0] - N_PPS ; //Calculated the number of healthy animals

/*=======Scenario 1======================================================================================= =======*/
//Scenario 1: slaughter + live < sick: collecting all samples from sick animals
if(GroupSampleSize[i][slaughter]+ GroupSampleSize[i][live]<= N_PPS)//if sample size from slaughter is less than N_PPS, take all from this group
{
//	printf("Can sample all from sick or previous pos");
	
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
						test(current_animal,test_pcr_slaughter_nomilk,slaughter,group,test_result_table,
						Se_bio, Sp_bio, Se_ID, Sp_ID, PCR_live, PCR_slaughter,farm) ;
					
						current_animal->present = 0;
						}
				
					if(counter_slaughter_ELISA>0 && ELISA_type<=2)
						{
						counter_slaughter_ELISA--;
						test(current_animal,ELISA_type,slaughter,group,test_result_table,
						Se_bio, Sp_bio, Se_ID, Sp_ID, PCR_live, PCR_slaughter,farm) ;
						
						current_animal->present = 0;
						}
					if(current_animal->present == 0) // means this animal was slaughtered
						{
						N_PPS--;//how many PPS remaining
						//Now update List_mng_status and milk_numbers table
						mbovis_status = current_animal->mbovis_status;
						List_mng_status[i][0] --;
			//			if(List_mng_status[i][0]<0)
			//			{
			//				printf("Now group %f goes minus\n",List_mng_status[i][0]);
			//				system("pause");
			//			}
						List_mng_status[i][mbovis_status+5] --;
			//			if(List_mng_status[i][mbovis_status+5]<0)
			//			{
			//				printf("Warning: group %d status negative  at Sc1 A\n",mbovis_status);
			//				system("pause");
			//			}
				
						if(*dry_season==0 && group ==1)
						{
							if(mbovis_status<=1)
							{
							mbovis_status = 0;	
							}
							else
							{
							mbovis_status = mbovis_status - 1 ;
							}
							if(current_animal->mastitis_detected==1)
							{
							milk_numbers[mbovis_status*3+2]--;	
							}
							else
							{
							milk_numbers[mbovis_status*3+current_animal->non_colostrum]--;	
							}
						
					//	if
					//	(milk_numbers[mbovis_status*3+current_animal->non_colostrum]<0)
					//	{
				//			printf("S-4\n");
					//		system("pause");
					//	}
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
						test(current_animal,test_pcr_live_nomilk,live,group,test_result_table,
						Se_bio, Sp_bio, Se_ID, Sp_ID, PCR_live, PCR_slaughter,farm) ;
						}
						if(counter_live_pcr_milk>0)
						{
						counter_live_pcr_milk--;
						test(current_animal,test_pcr_milk,live,group,test_result_table,
						Se_bio, Sp_bio, Se_ID, Sp_ID, PCR_live, PCR_slaughter,farm) ;
						}
						if(counter_live_ELISA>0 && ELISA_type<=2)
						{
						counter_live_ELISA--;
						test(current_animal,ELISA_type,live,group,test_result_table,
						Se_bio, Sp_bio, Se_ID, Sp_ID, PCR_live, PCR_slaughter,farm) ;
						}
						counter_live++;
						current_animal = current_animal->next_node;
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

/*=====================================================================================================================*/

/*=====================Scenario 2: slaughter sample can be sampled from N_PPS but not live samples====================================================================*/
//Need to sample from N_H for slaughter
//this includes slaughter sample is a bit and live sample is a lot, more than N_PPS
else if(GroupSampleSize[i][slaughter] <= N_PPS)
{
		//if sample enough for slaughter while sick still remains, GroupSampleSize[i][live] - N_PPS will be chosen from N_H
		//so if slaughter < sick, we go back to the top then sample all sick and every nth healthy
		
		//	printf("N_PPS is %d N_H is %d\n",N_PPS, N_H);
		//	printf("Slaughter sample is %d\n",GroupSampleSize[i][slaughter]);
		//	printf("live sample is %d\n",GroupSampleSize[i][live]);
		//	printf("slaughter count is %d\n",counter_slaughter);
		//	printf("live count is %d\n",counter_live);
		//take all sick for slaughter
		current_animal = FarmGroupList[i];
		while(counter_slaughter < GroupSampleSize[i][slaughter] && current_animal!=NULL)
		{
			if(current_animal->mbovis_elisa_status==1||current_animal->mbovis_pcr_status==1||current_animal->mastitis_detected==1)
			{
		
						if(counter_slaughter_pcr_nomilk>0)
						{
						counter_slaughter_pcr_nomilk--;
						test(current_animal,test_pcr_slaughter_nomilk,slaughter,group,test_result_table,
						Se_bio, Sp_bio, Se_ID, Sp_ID, PCR_live, PCR_slaughter,farm) ;
						current_animal->present = 0;
						}
						if(counter_slaughter_ELISA>0 && ELISA_type<=2)
						{
						counter_slaughter_ELISA--;
						test(current_animal,ELISA_type,slaughter,group,test_result_table,
						Se_bio, Sp_bio, Se_ID, Sp_ID, PCR_live, PCR_slaughter,farm) ;
						current_animal->present = 0;
						}
						//what do they do when animals found positive?
					if(current_animal->present == 0) // means this animal was slaughtered
					{
					N_PPS--;//how many PPS remaining
				//	printf("N_PPS is %d", N_PPS) ;
					mbovis_status = current_animal->mbovis_status;
					List_mng_status[i][0] --;
		//				if(List_mng_status[i][0]<0)
		//			{
		//				printf("(B) Now group %d goes minus\n",(int)List_mng_status[i][0]);
		//				system("pause");
		//			}
					List_mng_status[i][mbovis_status+5] --;
		//				if(List_mng_status[i][mbovis_status+5]<0)
		//			{
		//				printf("Warning: group %d status negative  at Sc2 B\n",mbovis_status);
		//				system("pause");
		//			}
							if(*dry_season==0 && group ==1)
							{
								if(mbovis_status<=1)
								{
								mbovis_status = 0;	
								}
								else
								{
								mbovis_status = mbovis_status - 1 ;
								}
							if(current_animal->mastitis_detected==1)
							{
							milk_numbers[mbovis_status*3+2]--;	
							}
							else
							{
							milk_numbers[mbovis_status*3+current_animal->non_colostrum]--;	
							}
						
							//	if
						//	(milk_numbers[mbovis_status*3+current_animal->non_colostrum]<0)
						//	{
						//		printf("S-3\n");
						//		system("pause");
						//	}
							}
					animal_to_remove = current_animal;
					current_animal = current_animal->next_node;
					remove_animal_group(FarmGroupList, i, animal_to_remove);	
					counter_slaughter++;		
					}
					else
					{
						break;
					}
			}
			else
			{
			current_animal=current_animal->next_node;
			}	
		}//sampling slaughter done
		//printf("A");
		current_animal = FarmGroupList[i];//reset
		//every sick and 
		counter_interval_healthy = 0;
		sampling_interval_healthy = floor((GroupSampleSize[i][live]-N_PPS)/N_H);
		if(sampling_interval_healthy<1)
		{
		sampling_interval_healthy = 1 ;
		}
		while(counter_live<GroupSampleSize[i][live] && current_animal!=NULL) //round 1 of live: sick animals
		{
			//if sick collect
			if(current_animal->mbovis_elisa_status==1||current_animal->mbovis_pcr_status==1||current_animal->mastitis_detected==1)
			{
						if(counter_live_pcr_nomilk>0)
						{
						counter_live_pcr_nomilk--;
						test(current_animal,test_pcr_live_nomilk,live,group,test_result_table,
						Se_bio, Sp_bio, Se_ID, Sp_ID, PCR_live, PCR_slaughter,farm) ;
						current_animal->tested = 1 ;
						
						}
						if(counter_live_pcr_milk>0)
						{
						counter_live_pcr_milk--;
						test(current_animal,test_pcr_milk,live,group,test_result_table,
						Se_bio, Sp_bio, Se_ID, Sp_ID, PCR_live, PCR_slaughter,farm) ;
						current_animal->tested = 1 ;
						}
						if(counter_live_ELISA>0 && ELISA_type<=2)
						{
						counter_live_ELISA--;
						test(current_animal,ELISA_type,live,group,test_result_table,
						Se_bio, Sp_bio, Se_ID, Sp_ID, PCR_live, PCR_slaughter,farm) ;
						current_animal->tested = 1 ;
						}
						counter_live++;
						current_animal = current_animal->next_node;	
			}
			else
			{
				current_animal = current_animal->next_node ;
			}
		}
		//printf("B");
		//printf("Total is %lf\n",List_mng_status[i][0]);
		//printf("counter_live is %d",counter_live);
		//printf("counter_live_pcr_nomilk is %d\n",counter_live_pcr_nomilk);
		//printf("counter_live_pcr_milk is %d\n",counter_live_pcr_milk);
		//printf("counter_live_ELISA is %d\n",counter_live_ELISA);
		
		current_animal = FarmGroupList[i];//reset
		while(counter_live<GroupSampleSize[i][live] && current_animal!=NULL)//round 2 of live: healthy animals
		{
			if(current_animal->tested==0)
			{
				//printf("not tested\n");
					//if not sick collect every nth
				counter_interval_healthy++;
				if(counter_interval_healthy<sampling_interval_healthy)
				{
				current_animal=current_animal->next_node;
				}	
				else
				{
				//	printf("going to be tested\n");
				counter_interval_healthy = 0;
						if(counter_live_pcr_nomilk>0)
						{
						//	printf("T1\n");
						counter_live_pcr_nomilk--;
						test(current_animal,test_pcr_live_nomilk,live,group,test_result_table,
						Se_bio, Sp_bio, Se_ID, Sp_ID, PCR_live, PCR_slaughter,farm) ;
						current_animal->tested = 1 ;
						
						}
						if(counter_live_pcr_milk>0)
						{
						//	printf("T2\n");
						counter_live_pcr_milk--;
						test(current_animal,test_pcr_milk,live,group,test_result_table,
						Se_bio, Sp_bio, Se_ID, Sp_ID, PCR_live, PCR_slaughter,farm) ;
						current_animal->tested = 1 ;
						}
						if(counter_live_ELISA>0 && ELISA_type<=2)
						{
						//	printf("T3\n");
						counter_live_ELISA--;
						test(current_animal,ELISA_type,live,group,test_result_table,
						Se_bio, Sp_bio, Se_ID, Sp_ID, PCR_live, PCR_slaughter,farm) ;
						current_animal->tested = 1 ;
						}
						//	printf("T4\n");
						counter_live++;
						current_animal->tested = 0;
						current_animal = current_animal->next_node;	//this may become NULL
						if(current_animal==NULL)
						{
						//	printf("End:");
							break;
						}
				//	printf("test done\n");
				}
			}
			else
			{
				current_animal = current_animal->next_node ;
					if(current_animal==NULL)
						{
						//	printf("End:");
							break;
						}
				current_animal->tested = 0;
			}
		} //done of round 2
//printf("C")	;	
}//else if(GroupSampleSize[i][slaughter] <= N_PPS) DONE
/*=================================Scenario 2 done=========================================================*/

/*============Scenario 3: Slaughter can't be sampled from N_PPS, which may be 0============================*/
else if(GroupSampleSize[i][slaughter] > N_PPS) //collect all N_PPS for slaughter then healthy
{
			//ok what if slaughter to ollcet is non-0 but N_PPS is 0?
		//	printf("S2");
		current_animal = FarmGroupList[i];//reset
		//printf("NPS is %d\n",N_PPS);
		//printf("Available size is %d\n",(int)List_mng_status[i][0]);//why no calve though?
		counter_interval_healthy = 0;
		sampling_interval_healthy = floor((GroupSampleSize[i][slaughter]-N_PPS)/N_H);
		if(sampling_interval_healthy<1)
		{
		sampling_interval_healthy = 1 ;
		}
		
		
		//if slaughter > sick, complete all sick remaining slaughter + live will be chosen every nth
		while(current_animal!=NULL && N_PPS>0)//Round 1 of sampling slaughter: sick animals
			{
			
				if(current_animal->mbovis_elisa_status==1||current_animal->mbovis_pcr_status==1||current_animal->mastitis_detected==1)
				{
					
							if(counter_slaughter_pcr_nomilk>0)
							{
							counter_slaughter_pcr_nomilk--;
							test(current_animal,test_pcr_slaughter_nomilk,slaughter,group,test_result_table,
							Se_bio, Sp_bio, Se_ID, Sp_ID, PCR_live, PCR_slaughter,farm) ;
							current_animal->present = 0;
							}
							if(counter_slaughter_ELISA>0 && ELISA_type<=2)
							{
							counter_slaughter_ELISA--;
							test(current_animal,ELISA_type,slaughter,group,test_result_table,
							Se_bio, Sp_bio, Se_ID, Sp_ID, PCR_live, PCR_slaughter,farm) ;
							current_animal->present = 0;
							}
							//what do they do when animals found positive?
						if(current_animal->present == 0) // means this animal was slaughtered
						{
						N_PPS--;//how many PPS remaining
					//	printf("N_PPS is %d", N_PPS) ;
						mbovis_status = current_animal->mbovis_status;
						List_mng_status[i][0] --;
			//				if(List_mng_status[i][0]<0)
			//			{
			//				printf("(C)Now group %f goes minus\n",List_mng_status[i][0]);
			//				system("pause");
			//			}
						List_mng_status[i][mbovis_status+5] --;
			//				if(List_mng_status[i][mbovis_status+5]<0)
			//			{
			//				printf("Warning: group %d status negative at Sc3 A\n",mbovis_status);
			//				system("pause");
			//			}
								if(*dry_season==0 && group ==1)
								{
									if(mbovis_status<=1)
									{
									mbovis_status = 0;	
									}
									else
									{
									mbovis_status = mbovis_status - 1 ;
									}
								if(current_animal->mastitis_detected==1)
								{
								milk_numbers[mbovis_status*3+2]--;	
								}
								else
								{
								milk_numbers[mbovis_status*3+current_animal->non_colostrum]--;	
								}
								
							
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
			
			}//Round 1: slaughter sampling from sick animals done
		//printf("R1 done\n");//collecting sick and previous positive animals done
		
		current_animal = FarmGroupList[i];//reset
			
		while(current_animal!=NULL && counter_slaughter<GroupSampleSize[i][slaughter]) //round 2 of slaughter: healthy animals
		{
		//	printf("a is %d and b is %d\n",a, b) ;
		//printf("today is %d\n",next_non_markov_date);
				counter_interval_healthy++;
				if(counter_interval_healthy<sampling_interval_healthy)
				{
				
				current_animal=current_animal->next_node;
				}	
				else
				{
				counter_interval_healthy = 0;
			//	printf("counter_slaughter is %d\n",counter_slaughter) ;
			//	printf("GS is %d and i is %d\n",GroupSampleSize[i][slaughter],i);
			//	printf("counter_slaughter_pcr_nomilk is %d\n",counter_slaughter_pcr_nomilk);
			//	printf("counter_slaughter_ELISA is %d\n",counter_slaughter_ELISA);
				//both elisa and pcr collected already but still GroupSampleSize says 1 more slaughter needed
		//			if(current_animal==NULL)
		//			{
		//				printf("Animal is null.counter_slaughter is%d",counter_slaughter);
		//				system("pause");
		//				
		//			}	
					
						if(counter_slaughter_pcr_nomilk>0)
						{
						counter_slaughter_pcr_nomilk--;
						test(current_animal,test_pcr_slaughter_nomilk,slaughter,group,test_result_table,
						Se_bio, Sp_bio, Se_ID, Sp_ID, PCR_live, PCR_slaughter,farm) ;
						current_animal->present = 0;
						}
						if(counter_slaughter_ELISA>0 && ELISA_type<=2)
						{
						counter_slaughter_ELISA--;
						test(current_animal,ELISA_type,slaughter,group,test_result_table,
						Se_bio, Sp_bio, Se_ID, Sp_ID, PCR_live, PCR_slaughter,farm) ;
						current_animal->present = 0;
						}
						//what do they do when animals found positive?
					//	printf("Inv1");
					if(current_animal->present == 0) // means this animal was slaughtered
					{
					N_H--;//how many PPS remaining
				//	printf("N_H is %d", N_H) ;
					mbovis_status = current_animal->mbovis_status;
					List_mng_status[i][0] --;
		//				if(List_mng_status[i][0]<0)
		//			{
		//				printf("(D)Now group %f goes minus\n",List_mng_status[i][0]);
		//				system("pause");
		//			}
					List_mng_status[i][mbovis_status+5] --;
		//				if(List_mng_status[i][mbovis_status+5]<0)
		//			{
		//				printf("Warning: group %d status negative at Sc3 B\n",mbovis_status);
		//				system("pause");
		//			}
							if(*dry_season==0 && group ==1)
							{
								if(mbovis_status<=1)
								{
								mbovis_status = 0;	
								}
								else
								{
								mbovis_status = mbovis_status - 1 ;
								}
							if(current_animal->mastitis_detected==1)
							{
							milk_numbers[mbovis_status*3+2]--;	
							}
							else
							{
							milk_numbers[mbovis_status*3+current_animal->non_colostrum]--;	
							}
							//how about mastitis detected ones?
							//well theoretically they should not have CM
						
					
							}
					animal_to_remove = current_animal;
					current_animal = current_animal->next_node;
					remove_animal_group(FarmGroupList, i, animal_to_remove);	
					counter_slaughter++;		
					}
				//	printf("Inv2");
				}
			
		} //Round 2 of slaughter done
		//printf("R2 done\n");
		current_animal = FarmGroupList[i];//reset
		//printf("Target is %d\n",GroupSampleSize[i][live]);
		//printf("Counter is %d\n",counter_live);
		//What if no live sample is colelcted from this group?
		//FInally healthy live sample
		if(GroupSampleSize[i][live]>0)//if live sample is collected
		{
			counter_interval_healthy = 0;
			sampling_interval_healthy = floor((GroupSampleSize[i][live])/N_H);
			if(sampling_interval_healthy<1)
			{
			sampling_interval_healthy = 1 ;
			}
			while(counter_live<GroupSampleSize[i][live] && current_animal!=NULL)
			{
					counter_interval_healthy++;
					if(counter_interval_healthy<sampling_interval_healthy)
					{
					current_animal=current_animal->next_node;
					}	
					else
					{
							if(counter_live_pcr_nomilk>0)
							{
							counter_live_pcr_nomilk--;
							test(current_animal,test_pcr_live_nomilk,live,group,test_result_table,
							Se_bio, Sp_bio, Se_ID, Sp_ID, PCR_live, PCR_slaughter,farm) ;
							}
							if(counter_live_pcr_milk>0)
							{
							counter_live_pcr_milk--;
							test(current_animal,test_pcr_milk,live,group,test_result_table,
							Se_bio, Sp_bio, Se_ID, Sp_ID, PCR_live, PCR_slaughter,farm) ;
							}
							if(counter_live_ELISA>0 && ELISA_type<=2)
							{
							counter_live_ELISA--;
							test(current_animal,ELISA_type,live,group,test_result_table,
							Se_bio, Sp_bio, Se_ID, Sp_ID, PCR_live, PCR_slaughter,farm) ;
						
							}
							counter_live++;
							current_animal = current_animal->next_node;
					}
			}	
		}
		
		//printf("R3 done\n");

}//else if(GroupSampleSize[i][slaughter] > N_PPS)  DONE

//else
//{
//	printf("Missing some scenario!");
//	system("pause");
//}

}//All sampling done


} //for(i=0; i<=id_dry_group; i++) DONE
//printf("Test %d is done\n",test_nth);
	test_nth++ ;
	
} //If test done


/*========EVENT 12: Sampling done=====================================*/

/*===============EVENT 13: Colostrum status change===============*/
	if(current_event->event_type==13) //colostrum animals become non-colostrum (start contributing to milk)
	{
		current_animal = current_event->animal ;
		if(current_animal->present==1)
		{
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
	
				 
	}
/*===============EVENT 13: Colostrum status change ENDS===============*/

/*================EVENT 14: Weaning====================================*/
//TASK_30
	if(current_event->event_type==14)
	{
		current_animal = current_event->animal ;
		if(current_animal->present==1)
		{
				current_grp = current_animal->group ;
		mbovis_status = current_animal->mbovis_status;
		//remove this animal to the next group
		if(current_grp != id_calf_group)
		{
		//	printf("This is not calf!");
		//	system("pause");
		}
		else
		{
			
			List_mng_status[id_calf_group][0] --;
			List_mng_status[id_calf_weaned][0] ++;
//				if(List_mng_status[id_calf_group][0]<0)
//			{
//				printf("Now weaning calf goes minus\n");
//				system("pause");
//			}
			
			List_mng_status[id_calf_group][mbovis_status+column_Sus] --;
			List_mng_status[id_calf_weaned][mbovis_status+column_Sus] ++;
//				if(List_mng_status[id_calf_group][mbovis_status+5]<0)
//			{
//				printf("Warning: group %d status negative at Weaning\n",mbovis_status);
//				printf("Weaned one is %d\n",(int)List_mng_status[id_calf_weaned][mbovis_status+column_Sus]);
//				system("pause");
//			}
			
			//TASK_30
			///change the numebr of mbovis depending on this animal status
			remove_animal_group(FarmGroupList,current_grp,current_animal);
     		current_animal->group = id_calf_weaned ; //change the group
     		add_animal_group(FarmGroupList,id_calf_weaned,current_animal);
		}	
		}
		/*Use for checking
		else
		{
			printf("This animal gone continue\n");
			system("pause");
		}*/

     		
     		
				 
	}
/*===============EVENT 14: Weaning ends=================================*/

/*================EVENT 15: Bobby pick up==================================*/
//TASK_31
if(current_event->event_type==15)
	{
	//	printf("bobby is %d",*num_bobby);
		if((*num_bobby)==0)
		{
		//	printf("No bobby exists");
		}
		else
		{
		//	printf("bobby is %d",*num_bobby);
			var_bobby = 0;
		//	printf("bobby is %d",*num_bobby);
		//	system("pause");
			//(*num_bobby) = 0 ;
		}
		if((next_non_markov_date+bobby_pickup)<sim_days && (next_non_markov_date+bobby_pickup < last_pick_up))
		{
		new_event = (struct event_node*)malloc(sizeof( struct event_node ));
				new_event->event_type = 15;//bobby
			
				new_event->animal = fake_animal; //let's see if this works
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
//	if(farm==39)
//	{
//		printf("Risk event\n");
//		printf("Date is %d\n",next_non_markov_date);
//		printf("Milk is %d\n",(int)List_mng_status[id_lact_group][0]);
//		printf("Dry is %d\n",(int)List_mng_status[id_dry_group][0]);
//		system("pause");
//			}
	//ok it's day 358 and claimed to be milking but no milking cows
	//perhaps all culled?
	//once milk animal hits 0, switch to dry?
	if(*dry_season==0)//milking season
	{
		if((int)List_mng_status[id_lact_group][0]==0)
		{
			group = id_dry_group ; //milking season but no milking animals
		}
		else
		{
				group = id_lact_group ;
		}

    }
    else //dry season
    {
    group = id_dry_group ;	
//    	if((int)List_mng_status[id_dry_group][0]==0)
//		{
//		printf("Dry season but no animals exist for seed\n");
//		system("pause") ;
//		}
	}
	current_animal = FarmGroupList[group];
//	if(farm==39)
//	{
//		printf("Num is %d,dry %d\n",(int)List_mng_status[group][0],*dry_season);
//	}
	
	random_number = rand()%(int)List_mng_status[group][0] + 1;
	while(random_number>0)
	{
	random_number--;
	    if(random_number==0)
	    {
		current_animal->mbovis_status=*init_bovis_status ;
		List_mng_status[group][column_Exp+*init_bovis_status-1]++;
		List_mng_status[group][column_Sus]--;
		break;
		}
	    else
	    {
	    current_animal= current_animal->next_node;	
	    }	
	}
//		if(farm==39)
//	{
//		printf("Risk event done\n");
//	}
   }

/*=================EVENT 16: Risk Event done====================*/


/*=================EVENT 17: UPDATING WASTE MILK VOLUME AND M BOVIS CONATINED*/
/*Updating waste milk volume happens every day during the milking season so put it here rather than
making this as an event*/
if(*dry_season==0 && current_event->event_type==17 && feed_waste_milk == 1)
{
		if(waste_milk_accumulate==0)
		{
		*volume_waste_milk = 0; //no carry over from previous day	
		}
		
		/*As long as pre-weaned calves exist this event happens.
		So every day between the first calving till the last wening
		Purpose here is to calculate intensity_mc: which indicates how large the infection pressure is compared to 
		when a calf receives all the M bovis shed by one clinical M bovis cow
		
		Process:
		1. Calculate total M bovis contained in waste milk tank (effective_mbovis_waste_milk) before adding new waste milk
		2. Update waste milk volume and effective_mbovis_waste_milk, calculate M bovis density per litre 
		3. Calculate the amount of required milk, determine if top up from other source (milk powder or bulk tank) is needed
		4. Compute intensity_mc and pass it
		*/
	
		//Calculate the before-update M bovis amount contained in the waste milk
		effective_mbovis_waste_milk = (*volume_waste_milk) *(*density_mbovis_waste_milk) ; /*Gives the total amount of m bovis in waste milk at this moment*/
		
		/*Update both the total volume of waste milk and density of M bovis in the waste milk*/
		
		/* Who is contributing to waste milk?
		1. Colostrum cow not shedding bovis: milk_numbers[0]
		2. Cows not shedding bovis but detected as CM: milk_numbers[2]
		3. Colostrum cows subclinically shedding bovis: milk_numbers[3]
		4. Cows subclinically shedding bovis and detected as CM (due to non bovis CM): milk_numbers[5]
		5. Cows in colostrum and clinically shedding bovis but not detected: milk_numbers[6]
		6. Cows clinically shedding bovis and detected as CM: milk_numbers[8]
		
		Now among these, who is contributing to M bovis shedding?
		3 (colostrum & subclinical), 4 (subclinical & detected), 5 (colostrum & clinical), 6 (clinical & detected)*/
		*volume_waste_milk = *volume_waste_milk + (milk_numbers[0]+milk_numbers[2]+milk_numbers[3]+milk_numbers[5]+milk_numbers[6]+milk_numbers[8])*av_milk_production;
		effective_mbovis_waste_milk = effective_mbovis_waste_milk + bv_f*(milk_numbers[3] + milk_numbers[5]) + milk_numbers[6] + milk_numbers[8];
		if(*volume_waste_milk>0)
		{
			*density_mbovis_waste_milk = effective_mbovis_waste_milk/(*volume_waste_milk) ;	
		}
		else
		{
			*density_mbovis_waste_milk = 0;
		}
		
		
		
		/*Now if there are any pre-weaned calves, need to calculate intensity_mc.
		Subtract the amount of waste milk consumed and calculate the amount of bulk milk to be used*/
		if(List_mng_status[0][0]>0)
		{
				required_milk = List_mng_status[0][0]*calf_milk_consumption ;
				if(*volume_waste_milk>= required_milk) /*If waste milk is sufficient, the desnity of M bovis in waste milk matters*/
				{
					*volume_waste_milk = *volume_waste_milk - required_milk;
					intensity_mc = (*density_mbovis_waste_milk) * calf_milk_consumption;
				}
				else
				{
					
					if(use_milk_powder==1) //if milk powder is going to be used instead of bulk milk 
					{
						intensity_mc = (*density_mbovis_waste_milk)*(*volume_waste_milk)/List_mng_status[0][0] ; /*desnity times the amount of waste milk distributed to each calf*/
					}
					else /*Insufficient milk comes from bulk milk so calculate how much M bovis conatined in bulk milk*/
					{
						
						healthy_milk = (double)(milk_numbers[1]+milk_numbers[4]+milk_numbers[7])*av_milk_production;	
						
						if(healthy_milk>0)
						{
							shed_healthy = bv_f*milk_numbers[4] + milk_numbers[7];
							topup_milk = required_milk - *volume_waste_milk ; /*insufficient milk amount to top up from bulk milk*/
							intensity_mc = (*density_mbovis_waste_milk)*(*volume_waste_milk)/List_mng_status[0][0] + 
							(topup_milk/List_mng_status[0][0]) * (shed_healthy/healthy_milk) ;	
						}
						else /*If no bulk milk available then milk powder or something is used*/
						{
							intensity_mc = (*density_mbovis_waste_milk)*(*volume_waste_milk)/List_mng_status[0][0] ;
						}
					
						/*The infection pressure wach calf receives is the sum of infection pressure from waste milk and normal milk.
						Density of M bovis in healthy milk is shed_healthy/healthy_milk, and infection pressure of healthy milk can be described as
						(shed_healthy/healthy_milk) multiplied by the amount of healthy milk a calf consumes (=necessary_milk/List_mng_status[0][0]) 
						*/
					}
					*volume_waste_milk = 0; //all waste milk consumed
				}	
		}
		/*Add next waste milk update event*/
		if(next_non_markov_date+1<sim_days)
		{
		new_event = (struct event_node*)malloc(sizeof( struct event_node ));
		new_event->event_type = 17;
		new_event->animal = fake_animal; //let's see if this works
		new_event->next_node = NULL ;
		add_event_node(event_day,next_non_markov_date+1, new_event) ;	
		}

	

}
/*UPDATING WASTE MILK VOLUME AND M BOVIS CONATINED IN IT DONE==========================================================*/
	
/*======Going through all events on this date done=================================================================*/
     	//MOVE to the next event; free the memory of the current event and move to the next one
     	struct event_node *previous_event;
		previous_event = current_event ; //rewire to the next event
	   	current_event = current_event->next_node;
//	   	if (current_event!=NULL)
//		   {
//		 //  	printf("going to free event now") ;
//		   
//		   //printf("next event is %d, %lld", current_event->event_type, current_event->akey);
//	       }
//       	else
//	       {
//	      //	printf("next event is NULL\n") ;
//		   }
	  // printf("now free event");
	   	free(previous_event);///check if this works - previously this was in if(current_event!=NULL) blacket
	  // printf("event was %d",previous_event->event_type);
	  // system("pause");
	   	event_day[next_non_markov_date] = current_event;
	   //printf("event pointing new event\n") ;
     }

	index_column++;	
//	printf("index column is %d\n",index_column);
//	printf("index updated %d",index_column) ;
 	} //LOOP NM1
	today_date = updated_date;
	if((int)today_date==sim_days-1)
	{
		//when it reaches the final day of simulation
	//	printf("Reaches the end of the simulation") ;
		break;
	}
	//printf("Now today date is %lf\n",today_date) ;
	
	
}
//printf("Loop done");
/*===========SIMULATION ENDS=======================================================================================*/
//printf("Total cull is %d\n",*num_culled);
//printf("Total sale is %d\n",*num_sold);
//printf("Total death is %d\n",*num_death);
//system("pause");
//	write_number_animals(NumberAnimalDataFile,NumGrpAnimal,id_bull_group,index_column) ;
//	write_test_result(TestResultFile, test_result_table,nrow_testschedule,ncol_test_table) ;

/*================FREEING MEMORY=============================================================================================*/
for(i = 0; i < id_bull_group+1; i++)
      	  	{
      	  		current_animal = FarmGroupList[i];
      	  		while(current_animal!=NULL)
      	  		{
      	  		///ADD ELISA testing
					if(i==id_dry_group||i==id_lact_group)
					    {
						
							double random = (double)rand()/(double)RAND_MAX;
							if((current_animal->mbovis_sero_status==1 && random<Se_ID) ||
							(current_animal->mbovis_sero_status==0 && random>Sp_ID))
							{
						
						summary_statistics[5] ++;
							}
							
						}	
      	  		previous_animal= current_animal;
      	  		current_animal = current_animal->next_node ;
					free(previous_animal);	
				}
      	  		
     	      FarmGroupList[i] = NULL; // initialise the animal struct
			}
		//	printf("Farm list done\n");
for(i=0;i<sim_days;i++)
            {
            new_event = event_day[i];
            while(new_event!=NULL)
             {
	          previous_event = new_event;
	          new_event = new_event->next_node;
	          free(previous_event);
	         }
			 event_day[i] = NULL;	
            }
         //   printf("event list done\n");

/*===============FREEING DONE=================*/	
//printf("Ite is %d\n",ite);
}//END OF ALL ITERATRIONS - FREE EVENT and animal_node_pointer BEFORE THIS

//printf("up here\n");
//system("pause");

free(animal_node_pointer);
//printf("up here2\n");
//system("pause");
}///simulation for all farm ends
//printf("now calculate\n");
//int sum[10] ;
//int num_tested[10];
//double result[10];
//for(i=0; i<10;i++)
//{
//	sum[i] = 0;
//	num_tested[i] = 0;
//   for(j=0;j< num_total_herd; j++)
//   {
//   //	printf("result is %d\n",test_result_table[j][i]);
//    sum[i]= sum[i]+test_result_table[j][i];	
//    num_tested[i] = num_tested[i] + test_result_table[j][i+10];
//   }
////   if(num_tested[i]==0)
////   {
////   	printf("Warning: Zero tested!\n");
////   	system("pause");
////   }
//  result[i] = sum[i]/(double)num_tested[i] ;
//	 summary_statistics[i] = result[i];
//}
for(i=0;i<4;i++)
{
	summary_statistics[i] = List_mng_status[id_dry_group][column_Sus+i] + List_mng_status[id_lact_group][column_Sus+i] ;
}

//N
summary_statistics[4] = summary_statistics[0]+summary_statistics[1]+summary_statistics[2]+summary_statistics[3];

//PCR+ve
summary_statistics[6] =round((summary_statistics[2]+summary_statistics[3])*PCR_live) ; //crude deterministic - use binomial if need proper estimates


//printf("Calculate Done\n");
for(i=0;i<id_bull_group+1;i++)
{
	free(List_mng_status[i]);
}
free(List_mng_status);
//printf("A\n");
for(i=0;i<num_cull_sell_steps;i++)
{
	free(cull_sell_rate[i]);
}
free(cull_sell_rate);
//printf("B\n");
for(i=0;i<num_mortality_steps;i++)
{
	free(mortality[i]);
}
free(mortality);
//printf("C\n");
//for(i=0;i<id_bull_group+2;i++)
//{
//	free(NumGrpAnimal[i]);
//}
//free(NumGrpAnimal);

for(i=0;i<8;i++)
{
	free(table_CM_rate[i]);
}
free(table_CM_rate);
//printf("D\n");
free(milk_numbers);
//printf("E\n");
for(i=0;i<id_dry_group+1;i++)
{
	free(GroupSampleSize[i]);
}
free(GroupSampleSize);
//printf("F\n");
for(i=0;i<num_total_herd;i++)
{
	free(test_result_table[i]);
}
free(test_result_table);
//printf("G\n");
if(pre_defined_test==1)
{
		for(i=0;i<nrow_testschedule;i++)
		{
		free(test_schedule[i]);
		}
		free(test_schedule);
}

//printf("H\n");

//printf("freeing close\n");
free(animal_to_remove);
free(fake_animal);
//for(i=0;i<10;i++)
//{
//	printf("Summary is %f\n",*summary_statistics[i]);
////	printf("Tested is %d\n",num_tested[i]);
//}

//printf("freeing DONE\n");
} //main ends here

/*=================================================================*/










/*=============================================================================================================================================
**********ALL FUNCTIONS BELOW****************************************************
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
int* num_culled, int* num_sold, int* num_death, double bv_f, double rate_sero_conversion, int id_lact_group2,int id_calf_group2,
int* num_bobby, int* dry_season, int* milk_numbers, struct event_node* event_day[], int* inf_par, double *beta_par, double intensity_mc, int farm)
{
// Calculate the sum of rate of M events
//printf("column_Sus is %d",column_Sus);
//system("pause");
     struct event_node* new_event ;
	double day_to_markov;
	struct animal_node* current_animal;
	struct animal_node* next_animal;
	struct animal_node* previous_animal;
	double current_rate = 0;
	double sum_rate = 0.0;
	double demographic_rate = 0.0 ;
	double total_sero_conversion = 0.0 ;
	double accumu_rate = 0.0 ;
	double r_SE,r_EIs, r_IsIc, r_IcIs, r_IsE;

	double beta ;
	int column_Sus2 = 5;
	double bv_r1 = 1/(double)inf_par[0];
	double bv_r2 = 1/(double)inf_par[1];
	double bv_r3 = 1/(double)inf_par[2];
	double bv_r4 = 1/(double)inf_par[3];
	double beta_mm = (double)beta_par[0];
	double beta_mc = (double)beta_par[1];


	int random_int;
	int i,j;
	int k = 0;
	int pregnant_status ;
	int change ;
	int num_CM_animal_undetected = 0;
	int counter ;
	int mbovis_status;
	int index = 0;
	int current_milk_numbers;


	double disease_transition_sum = 0;

/*================================================================================*/
/*===== Calculate the total rate for Markov events=================================*/
//No need to count the numebr of each status but calculate disease rate
	for(i=0; i<id_lact_group2+1; i++)
	{//for each management group
	//printf("Management grp is %d",i);
	
	current_animal = FarmGroupList[i] ;
	//printf("num in i is %lf",List_mng_status[i][0]);
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
	//	printf("conversion is %lf\n",List_mng_status[i][2]);
	//	system("pause");	
		}
		
		
/*=======Calculating total rate of clinical mastitis=================*/		
	if(i==id_lact_group2 && *dry_season==0)//meaning it's milking season
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
	
		current_animal = current_animal->next_node ;
	}
/*================Calculating Sum of demographic event done==========*/

/*================Calculating sum of Cow to calf transmission=======*/
/*Note that this model does not consider dam to calf transmission in 12 hours after the birth.
If this is going to be added, need to record the number of new calves born on today (which can be done in calving event)
And put some transmission parameters for a direct dam to calf, then compute the total rate of this event*/


if(i==0 && *dry_season==0 && (int)List_mng_status[id_calf_group][0]>0 && feed_waste_milk==1) //wait so in dry period what happens? yes no transmission
		{
	
			if(intensity_mc>0)
			{
				beta = beta_mc*intensity_mc ;/*The infection pressure a calf receives from ingesting contamined milk*/
			
				if(beta<0)
				{
					printf("beta below 0\n");
				}
				if(List_mng_status[i][5]<0)
				{
					printf("Sus below 0 %d\n",(int)List_mng_status[i][5]);
					printf("E is %d\n",(int)List_mng_status[i][6]);
					printf("IS is %d\n",(int)List_mng_status[i][7]);
					printf("Ic is %d\n",(int)List_mng_status[i][8]);
					printf("Whole number %d\n",(int)List_mng_status[i][0]);//so total size is ok but just sus number is wierd
				}
				List_mng_status[i][9] = beta*List_mng_status[i][5];
				//for now calf there is no transition from exposed to infectious
				//Assumption_8
				List_mng_status[i][10] = 0;//exposed to Is
				List_mng_status[i][11] = 0;//Is to Ic
				List_mng_status[i][12] = 0;//Ic to Is
				List_mng_status[i][13] = 0;//Is to E 
				
				List_mng_status[i][14] = List_mng_status[i][9];//for ccalf disease transition is from S to E only
				disease_transition_sum = disease_transition_sum + List_mng_status[i][14] ;
				//List_mng_status[i][10]+List_mng_status[i][11]+List_mng_status[i][12]+
			//	List_mng_status[i][13] ;	
			//TASK_12 DONE		
			}
		
		}
	
/*==========TRANSMISSION FROM COW TO CALVES DONE=======================================================*/

/*==========CALCULATE TRANSMISSION BETWEEN COWS, DISEASE TRANSITION, DETECTION OF NEW CM=====================================*/
		if(i==id_lact_group2 && *dry_season==0)
		{

		
		List_mng_status[i][16] = num_CM_animal_undetected*rate_detect_clinical_mastitis ;
		//Assumption_9
		List_mng_status[i][9] = beta_mm*(milk_numbers[4]*bv_f+milk_numbers[7])*List_mng_status[i][5];
		/*Infection pressure only comes from undetected/untreated cows.
		All susceptible animals whether treated/untreated receive the same infection pressure*/
	
		List_mng_status[i][10] = bv_r1*List_mng_status[i][6];//exposed to Is
		List_mng_status[i][11] = bv_r2*List_mng_status[i][7];//Is to Ic
		List_mng_status[i][12] = bv_r3*List_mng_status[i][8];//Ic to Is
		List_mng_status[i][13] = bv_r4*List_mng_status[i][7];//Is to E 
		
		List_mng_status[i][14] = List_mng_status[i][9]+
		List_mng_status[i][10]+List_mng_status[i][11]+List_mng_status[i][12]+
		List_mng_status[i][13] ;
		disease_transition_sum = disease_transition_sum + List_mng_status[i][14] ;
		} //calculating sum rate for milking group done
/*==========CALCULATE TRANSMISSION BETWEEN COWS, DISEASE TRANSITION, DETECTION OF NEW CM DONE=====================================*/

/*================Calculating sum of m.bovis related Markov sum=======*/
		

		//sum altogether
		//TASK_6 no need to change
		sum_rate = 	sum_rate + List_mng_status[i][1] + List_mng_status[i][2]+ List_mng_status[i][14] +
		List_mng_status[i][15] + List_mng_status[i][16] ;
//		if(List_mng_status[i][1]<0)
//		{
//			printf("Group %d demo is neg",i);
//			system("pause");
//		}
//		if(List_mng_status[i][2]<0)
//		{
//			printf("Group %d sero is neg",i);
//				system("pause");
//		}
//		if(List_mng_status[i][14]<0)
//		{
//			printf("Group %d transition is neg",i);
//				system("pause");
//		}
//		if(List_mng_status[i][15]<0)
//		{
//			printf("Group %d CM is neg",i);
//				system("pause");
//		}
//		if(List_mng_status[i][16]<0)
//		{
//			printf("Group %d Cm detect is neg",i);
//				system("pause");
//		}
		//[1]: demographic sum
		//[2]: seroconversion
		//[14]: disease transition
		//[15]: incidence of CM
		//[16]:detection CM
		
		//printf("sum rate of group %d is %lf\n",i, sum_rate);
	//	printf("demographic is %lf\n", List_mng_status[i][1]);
	//	printf("seroconversion is %lf\n", List_mng_status[i][2]);
	//	printf("transition is %lf\n", List_mng_status[i][14]);
	//	printf("CM incidence is %lf\n", List_mng_status[i][15]);
	//	printf("CM detection is %lf\n", List_mng_status[i][16]);
		
		demographic_rate = demographic_rate + List_mng_status[i][1] ;
		total_sero_conversion = total_sero_conversion + List_mng_status[i][2] ;
//	printf("sum is %lf\n",sum_rate);
	}
	
/*=============================================================================*/
//if(farm==36 && today_date>=745)
//{
//	printf("sum rate is %lf\n",sum_rate);
//}
	


//Now calculate a waiting time
	 double random_value = (double)((double)rand()+1)/((double)RAND_MAX+1);//avoids 0
//	 printf("random value is %f",random_value);
     day_to_markov =(-log(random_value))/sum_rate ; // Waiting time
//     if(day_to_markov<0)
//     {
//     	printf("Error: markov date cannot be negative!") ;
//     	system("pause") ;
//	 }
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
	  
	random_value =  (double)(rand()+1)/((double)RAND_MAX+1)*sum_rate;
	//rand gives a value from 0 to rand_max, so above gives
	//random_value shoud not be 0 because the next if clause becomes false even when demographic_rate = 0
//	printf("random is %lf\n",random_value);
//	printf("demo is %lf\n",demographic_rate);
	//Now determine this value is greater or smaller than sum of demographic rate
	//demographic_rate is always non-zero so I don't really have to worry here
	
/*======================MARKOV is DISEASE EVENT===================*/
	if(random_value-demographic_rate>0)
	{
		//printf("Yes");
	//now do events other than demographic - CM occur + CM detection + MBOVIS transition
//	printf("Non-demographic");
	i = 0;
	accumu_rate = 0.0 ;
	//the problem might is not overwriting the variable
	
	double random_value2 ;

	int index_animal_to_chose = 0;
	struct animal_node* animal_to_chose ;

//	printf("V is %f\n",random_value-List_mng_status[id_lact_group2][15]-demographic_rate) ;
	if(random_value-List_mng_status[id_lact_group2][15]-demographic_rate<=random_value*0.0000000001) //CM incidnece
	/*==============1. CM incidence===================================*/
	{
//if(farm==36 && today_date>=745)
//{
//	printf("CM incidence\n");
//		}		
		//TASK_7 completed
//		if(List_mng_status[id_lact_group2][15]<0.0001)
//		{
//			printf("Zero CM rate\n");
//			
//					printf("Sero is %f\n",total_sero_conversion);
//						printf("Trans Rate is %f\n",List_mng_status[id_lact_group2][14]);
//							printf("CM Rate is %f\n",List_mng_status[id_lact_group2][15]);
//			printf("CM detection Rate is %f\n",List_mng_status[id_lact_group2][16]);
//				printf("Random is %f\n",random_value);
//				printf("Sum is %f\n",sum_rate);
//				printf("Demo is %f\n",demographic_rate);
//				printf("Dry season is %d",*dry_season);
//				printf("sero conversion is %f",rate_sero_conversion);
//			system("pause");
//		}
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
	//	printf("random value A is %f\n",random_value);
		random_value2 =  (double)(rand()+1)/((double)RAND_MAX+1)*List_mng_status[id_lact_group2][15];
		///what happens if random value gets same as the sum of List_mng_status?
	//	printf("random value B is %f\n",random_value);
		while(random_value2>accumu_rate && animal_to_chose!=NULL)
		{
			//TASK_21
		if(animal_to_chose->NB_mastitis_status==0 && animal_to_chose->mastitis_detected==0) //eligible animals is those without CM yet
		{
		 accumu_rate = accumu_rate+ animal_to_chose->mastitis_rate ;
		 if(accumu_rate>=random_value2)	
		 	break ;
		}
		
		animal_to_chose = animal_to_chose->next_node;//yes when it exceeds the rate it's actually next animal
	//	printf("next animal") ;
		}
//		if(animal_to_chose==NULL)
//		{
//			printf("No eligible CM animal\n");
//			system("pause");
//		}
		//above if 
//			if(animal_to_chose->NB_mastitis_status==1)
//			{
//				printf("This animal already has non Bovis CM!");
//				system("pause");
//			}
//			if(animal_to_chose->mastitis_detected==1) //does it matter though? Oh yes because animal won't get CM under treatment
//			{
//				printf("This animal already is detected as CM!");
//				system("pause");
//			}
			
		
			
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
	//	printf("Done\n");
	}
/*====================CM DETECTION========================================*/
	else if(random_value- List_mng_status[id_lact_group2][16]-List_mng_status[id_lact_group2][15]-demographic_rate<=random_value*0.0000000001)
	{
//	if(farm==36 && today_date>=745)
//{
//	printf("CM detection\n");
//		}
//	printf("CM detection");
	animal_to_chose = 	FarmGroupList[id_lact_group2];
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
		//printf("num_CM_animal_undetected is %d",num_CM_animal_undetected);
		index_animal_to_chose = rand()%num_CM_animal_undetected + 1 ;
		//printf("index is %d",index_animal_to_chose);
		counter = 0 ;
				while(counter<index_animal_to_chose)
				{
				if((animal_to_chose->NB_mastitis_status==1||animal_to_chose->mbovis_status==3)
				&& animal_to_chose->mastitis_detected==0) ///CM either non-bovis or bovis and not detected
				{ //if subclinical but not due to M bovis
				counter++;
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
	
	
	else if(random_value- List_mng_status[id_lact_group2][16]-List_mng_status[id_lact_group2][15]-demographic_rate-disease_transition_sum<=random_value*0.00000000001)
/*====================MBOVIS STATUS TRANSITION=============================*/
	{
//				if(farm==36 && today_date>=745)
//{
//	printf("Transition %f\n",disease_transition_sum);
//		}

//	printf("Transition %f\n",disease_transition_sum);

	//i is already 0 in above
	///TASK_15: need to add cow to calf transmission
//	printf("random_value in mbovis is %lf\n",random_value);
//printf("random v2 is %lf",random_value);
//ok problem is likely that when recalculating random value, tiny error occurs and eventhough 
//the subtracted random_value is the same as List_mng_status value, it does not recognise as equal
//alternative is that coming in this section means it's either calf or adult
//get sum and get random number and check if it's calf or adult
i = 0;
	while(i<id_lact_group+1)
	{
		//printf("transition sum is %lf\n",List_mng_status[i][14]);
		accumu_rate = accumu_rate+ List_mng_status[i][14] ;//List_mng_status[i][14] is sum of disease markov in a given mng grp
		i++;
	//	printf("i is %d\n", i) ;
	}
	
	random_value2 = 	((double)(rand()+1)/((double)RAND_MAX+1))*accumu_rate ;
//					if(farm==36 && today_date>=745)
//{
//	printf("random value2 is %f\n",random_value2);
//		}
//	printf("random value2 is %f\n",random_value2);
//	
	if(random_value2>List_mng_status[id_calf_group2][14])//i-1 means it's milkers
	{
	i = id_lact_group + 1;
//					if(farm==36 && today_date>=745)
//{
//printf("Adult transition");	
//		}
//	
	random_value = ((double)(rand()+1)/((double)RAND_MAX+1))* List_mng_status[i-1][14];
//								if(farm==36 && today_date>=745)
//{
//printf("random_value is %f\n",random_value);
//		}
	//********Below obtain how many animals in each category, and get random number below pop size
	//********Then, select this nth animal having this status and change their status
	if(random_value<=List_mng_status[i-1][column_r_SE]+0.000000001*random_value) //S to E
	{
//							if(farm==36 && today_date>=745)
//{
//printf("S to E\n");
//		}
//	
	index_animal_to_chose = rand()%(int)List_mng_status[i-1][column_Sus2]+1 ;
	index = 5;
	change = 1;//change indicates how status changes e.g. from 5 to 6 means 1 plus
	//now jump to the animal linked list and change the status
	//no change in milk_numbers
	
	}
	else if(random_value<List_mng_status[i-1][column_r_SE]+List_mng_status[i-1][column_r_SE+1]+0.000000001*random_value)
	{
//									if(farm==36 && today_date>=745)
//{
//printf("E to Is\n");
//		}
//	
	index_animal_to_chose = rand()%(int)List_mng_status[i-1][column_Sus2+1]+1 ; //E to Is
	index = 6;
	change = 1;
	//now jump to the animal linked list and change the status
	
	}
	else if(random_value<List_mng_status[i-1][column_r_SE]+List_mng_status[i-1][column_r_SE+1]+List_mng_status[i-1][column_r_SE+2]+0.000000001*random_value)
	{
//											if(farm==36 && today_date>=745)
//{
//printf("Is to Ic\n");
//		}
//	
	index_animal_to_chose = rand()%(int)List_mng_status[i-1][column_Sus2+2]+1 ; //Is to Ic
	index = 7;
	change = 1;
	//now jump to the animal linked list and change the status	
	
	}
	else if(random_value<List_mng_status[i-1][column_r_SE]+List_mng_status[i-1][column_r_SE+1]+List_mng_status[i-1][column_r_SE+2]+List_mng_status[i-1][column_r_SE+3]+0.000000001*random_value)
	{
//													if(farm==36 && today_date>=745)
//{
//printf("Ic to IsE\n");
//		}
//	
	index_animal_to_chose = rand()%(int)List_mng_status[i-1][column_Sus2+3]+1 ; //Ic to Is
	index = 8;
	change = -1;
	//now jump to the animal linked list and change the status	
	}
	else if(random_value<=
	List_mng_status[i-1][column_r_SE]+List_mng_status[i-1][column_r_SE+1]+
	List_mng_status[i-1][column_r_SE+2]+List_mng_status[i-1][column_r_SE+3] +List_mng_status[i-1][column_r_SE+4]+0.000000001*random_value )
	{
//															if(farm==36 && today_date>=745)
//{
//printf("Is to E,num is %d\n",(int)List_mng_status[i-1][column_Sus2+2]);
//		}
//	
	index_animal_to_chose = rand()%(int)List_mng_status[i-1][column_Sus2+2]+1 ; //Is to E
	index = 7;
	change = -1;
	//now jump to the animal linked list and change the status	
	}
//	else
//	{
//		printf("Out of range!");
//		system("pause");
//	}

	animal_to_chose = FarmGroupList[i-1];
//	if(animal_to_chose==NULL)
//	{
//		printf("No animal in this group\n");
//		system("pause");
//	}
	counter = 0; //if counter becomes equal to index_animal_to_chose, then that's the animal to choose
	//index_animal_to_chose indicates nth animal to choose
	//index indicates which disease status to choose (e.g. 6 means E staus)
	//index ranges 5 to 8, and status corresponds to 0 to 3
	mbovis_status = index-5 ;//gives current mbovis status
	
	while(counter<index_animal_to_chose)
	{
		if(animal_to_chose->mbovis_status==mbovis_status) 
		{
		counter++;
		}
		if(counter==index_animal_to_chose)
		{
		//TASK_19
	    //Update milk_numbers table
			if(mbovis_status<=1)
			{
				mbovis_status = 0; //E set to 0 for milk_numbers table changing purpose
			}
			else
			{
			mbovis_status = mbovis_status - 1;	
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
//		if(index==0)
//		{
//			printf("Index should not be 0!\n");
//			system("pause");
//		}
		List_mng_status[i-1][index+change]++; // n of (index+change) status animal plus 1
		break;		
		}
		else
		{
		animal_to_chose = animal_to_chose->next_node ; //go to next animal
		}	
	 }//changing disease status done	
    }//if it is milkers done
/*========Calf infectio================================*/
else
    {
//    	if((int)List_mng_status[0][column_Sus2]<=0)
//    	{
//    		printf("No sus, rate is %f\n",List_mng_status[0][14]);
//    		system("pause");
//		}
    	//clearly this is happening when non-sus exists
    	
    	//problem that although calf infection sum being 0 this event keeps happening
//    						if(farm==36 && today_date>=745)
//{
//	printf("Calf infection\n");
//		}
//
 //   printf("CM is %f\n",List_mng_status[id_lact_group2][15]);
//	printf("Detection is %f\n",List_mng_status[id_lact_group2][16]);
//	printf("conversion is %f\n",total_sero_conversion);
//	printf("accumurate is %f\n",accumu_rate);
	//if(total_sero_conversion>0)
	//{
	//	printf("sero over 0\n");
	//	system("pause");
	//}
//	printf("random value2 is %f\n",random_value2);
	i = id_calf_group2 + 1;
//	printf("calf sum is %lf\n",List_mng_status[i-1][14]);
//	printf("column_Sus is %d",column_Sus2);
	//printf("sus num is %lf",List_mng_status[i-1][column_Sus2]);//this number gets minus!
    index_animal_to_chose = rand()%(int)List_mng_status[i-1][column_Sus2]+1 ;//this is getting the number of susceptible
   // printf("index is %d",index_animal_to_chose);
	animal_to_chose = FarmGroupList[i-1];//this is giving null?
//	if(animal_to_chose==NULL)
//	{
//		printf("No animal in calf\n");
//		system("pause");
//	}
	counter = 0;
      while(counter<index_animal_to_chose)
      {
      	 if(animal_to_chose->mbovis_status==0)//if susceptible 
		{
		counter++;
		//printf("counter is %d",counter);
			if(counter==index_animal_to_chose)
		{
		//	printf("Inv3");
		animal_to_chose->mbovis_status=1;
			List_mng_status[id_calf_group2][column_Sus2]--;
//			if((int)List_mng_status[id_calf_group2][column_Sus2]<0)
//			{
//				printf("Warning: Suspectible calf goes below 0!\n");
//				system("pause");
//			}
			List_mng_status[id_calf_group2][column_Sus2+1]++;
			break;
	    }
		}
	
	    else
	    {
	    animal_to_chose = animal_to_chose->next_node ;	
//	    if(animal_to_chose==NULL)
//	    {
//	    	printf("This is end");
//	    	system("pause");
//		}
		}
	  }
	   
	} //if calf done
//printf("transition done\n");
	}
/*===========MBOVIS RELATED DISEASE EVENTS DONE============================*/
		
/*=================SERO-CONVERSION==================================*/


//	else if(random_value <= total_sero_conversion+disease_transition_sum+List_mng_status[id_lact_group2][15]+List_mng_status[id_lact_group2][16]+demographic_rate)
	else
	{
//							if(farm==36 && today_date>=745)
//{
//	printf("Sero-conversion\n");
//		}
	//	printf("Sero-conversion");
	random_value2 =  (double)(rand()+1)/((double)RAND_MAX+1)*total_sero_conversion ;
	if(total_sero_conversion==0)
	{
		printf("No conversion\n");
	//	printf("random b is %f\n",random_value);
      //   printf("Sum is %f\n",total_sero_conversion+disease_transition_sum+List_mng_status[id_lact_group2][15]+List_mng_status[id_lact_group2][16]+demographic_rate);
		system("pause");
	}
	while(random_value2>accumu_rate)
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
//	printf("sero converted\n") ;
	}

/*=================Sero-conversion done===============================*/	
		
	
	
	//event in k-1 happens
	} //end of disease markov
/*=====================MARKOV IS DEMOGRAPHIC EVENT========================================*/
	else
	{
	//	printf("Demographic");
		//now do demographic
		//TASK_19
	i = 0;
	accumu_rate = 0.0 ;
	
	//now choose one of management group for the event
	while(demographic_rate>accumu_rate)
	{
		accumu_rate = accumu_rate+ List_mng_status[i][1] ;
		i++;
		//printf("i is %d\n", i) ;
	} //i is the group at which an event occurs	
//	printf("random and accu are %f, %f\n",random_value,accumu_rate);
//	printf("i is %d",i);
	//printf("group decided");
//	printf("num before is %lf", List_mng_status[i-1][0]);
 //   printf("Demo A\n");
  //  printf("List is %f\n",List_mng_status[i-1][1]);
	random_value = (double)rand()/((double)RAND_MAX+1)* List_mng_status[i-1][1];
	//now decides which animal is going to have an event
	accumu_rate = 0; //reset accumulate, this time use this for sum of rate over animals
	int k = 0;
	//j=0;
	if(FarmGroupList[i-1]==NULL)
	{
	//	printf("There is no animals in this group!") ;
	//	system("pause");
	}
	else                          
	{
	//	printf("Else demo\n");
		current_animal = FarmGroupList[i-1] ;
		accumu_rate = current_animal->sum_markov_rate ;
	//	printf("choose animals");
	while(random_value>accumu_rate)
		{
		//	printf("curernt animal is %d",current_animal->akey) ;
		current_animal = current_animal->next_node ;
		k++;
		if(current_animal==NULL)
		{
			printf("Current animal becomes null");
		//	printf("random V is %f\n",random_value);
		//	printf("accumu_rate is %f\n",accumu_rate);
		//	printf("Total is %f Animal is %d\n",List_mng_status[i-1][1],k);
		//	system("pause");
		}
		accumu_rate = accumu_rate + current_animal->sum_markov_rate ;
		}
	
	//	printf("animal decided");
		//current_animal is the animal that will have an event
	//Choose event
	//printf("Demo B\n");
	random_value = (double)(rand()+1)/((double)RAND_MAX+1)* (current_animal->sum_markov_rate);
	//printf("Demo C\n");
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
	//	printf("this is cull");
		//system("pause");
	}
	else if(random_value<=(cull_sell_rate[current_animal->index_cull_sell][1+2*(current_animal->pregnant_status)]+cull_sell_rate[current_animal->index_cull_sell][1+2*(current_animal->pregnant_status)]))
	{
		//sale
		(*num_sold)++;
	//	printf("this is sold");
		//system("pause");
	}
	else
	{
		//death
		(*num_death)++;
	//	printf("this is death");
		//system("pause");
		
	}
	//printf("i is %d",i);
	List_mng_status[i-1][0]--;
	//printf("num is %lf", List_mng_status[i-1][0]);
	//printf("i is %d",i) ;
//	if(List_mng_status[i-1][0]<0)
//	{
//		printf("Warning! Animals below 0");
//		system("pause");
//	}
	//Need to update List_mng_status-bovis status and milk_numbers table
	mbovis_status = current_animal->mbovis_status ;
	List_mng_status[i-1][mbovis_status+5]-- ;//minus
//		if(List_mng_status[i-1][mbovis_status+5]<0)
//	{
//		printf("Warning! Status went below 0");
//		system("pause");
//	}
	//milk_numbers changes only if it is milking animals
	if(i==id_lact_group2+1)
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

//	printf("event decided");
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
//	printf("Demographic done\n");
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
   // printf("Opening file cull sell");
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
  //  printf("Opening file");
    FILE *temp = fopen(TestScheduleFile,"r"); 
       
    /* DECLARE STORAGE VARIABLES */
    int line_num, id, date;
    int slaughter_PCR_nomilk_calf,slaughter_ELISA_calf;
    int live_PCR_nomilk_calf,live_ELISA_calf;
    int slaughter_PCR_nomilk_adult,slaughter_ELISA_adult;
    int live_PCR_nomilk_adult, live_PCR_milk_adult, live_ELISA_adult ;
    int btank2, bdisc2, ELISA_type ;
   
    
    /* READ LINES OF FARM FILE */
    for(line_num = 0; line_num < nrow_testschedule; line_num++)
      { 
         fscanf(temp, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d, %d,%d,%d,%d"
		 ,&id, &date, 
		 &slaughter_PCR_nomilk_calf,&slaughter_ELISA_calf,
		 &live_PCR_nomilk_calf, &live_ELISA_calf,
		 &slaughter_PCR_nomilk_adult,&slaughter_ELISA_adult,
		 &live_PCR_nomilk_adult,&live_ELISA_adult, &live_PCR_milk_adult,
		 &btank2, &bdisc2, &ELISA_type
		 );
         
             TestSchedule[line_num][0] = id;
             TestSchedule[line_num][1] = date;
             TestSchedule[line_num][2] = slaughter_PCR_nomilk_calf;
             TestSchedule[line_num][3] = slaughter_ELISA_calf;
             
             TestSchedule[line_num][4] = live_PCR_nomilk_calf;
             TestSchedule[line_num][5] = live_ELISA_calf;
             
             TestSchedule[line_num][6] = slaughter_PCR_nomilk_adult;
             TestSchedule[line_num][7] = slaughter_ELISA_adult;
             
             TestSchedule[line_num][8] = live_PCR_nomilk_adult;
             TestSchedule[line_num][9] = live_ELISA_adult;
             TestSchedule[line_num][10] = live_PCR_milk_adult;
             
             
             TestSchedule[line_num][11] = btank2;
             TestSchedule[line_num][12] = bdisc2;
             TestSchedule[line_num][13] = ELISA_type ;
             
      }
      fclose(temp);
}

/*========WRITE DOWN THE NUMBER OF ANIMALS IN EACH MANAGEMENT GROUP IN A GIVEN TIME===================*/
//int write_number_animals(char* NumberAnimalDataFile,double** NumGrpAnimal,int id_bull_group,int n_column_output)
//{
//	FILE *Out = fopen(NumberAnimalDataFile, "w") ;
//	int line_num, col_num;
//	
//	for (line_num = 0 ; line_num < id_bull_group+2; line_num ++)
//	{
//		for (col_num = 0;col_num <n_column_output+1 ; col_num++ )
//		
//		{
//	    fprintf(Out,"%lf,",NumGrpAnimal[line_num][col_num]);
//        }
//        fprintf(Out,"\n");
//    }
//	fclose(Out);
//	return 0;
//}
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
void test(struct animal_node* current_animal,int test_type,int slaughter_live, int group,int** test_result_table,
double Se_bio, double Sp_bio, double Se_ID, double Sp_ID, double PCR_live, double PCR_slaughter, int farm) 
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
	int num_test_vars = 10;
//	int col_mastitis_tested_pcr = num_test_vars*2;
//	int col_mastitis_tested_pos = num_test_vars*2+1;

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
			test_result_table[farm][group*2+ELISA_calf_bio] ++;
			} //detecting true pos
			test_result_table[farm][group*2+ELISA_calf_bio+num_test_vars]++;
		}
		if(test_type!=0)
		{
			Se = Se_ID;
			Sp = Sp_ID;
			if((current_animal->mbovis_sero_status==1 && random<Se) ||
			(current_animal->mbovis_sero_status==0 && random>Sp))
			{
		
			current_animal->mbovis_elisa_status=1 ;
			test_result_table[farm][group*2+ELISA_calf_ID] ++;
			}
			test_result_table[farm][group*2+ELISA_calf_ID+num_test_vars]++;
		}
	
	
	}
	else //ANY PCR
	{
	//First, check if this animal had mastitis
//	if(current_animal->mastitis_detected==1)
//	{
//		if(test_result_table[farm][col_mastitis_tested_pcr]==0)
//		{
//			test_result_table[farm][col_mastitis_tested_pcr] = 1 ;//this round mastitis animal was tested by PCR
//		}
//	}
		if(test_type==3)//live PCR nomilk
		{
			Se = PCR_live;
			if(current_animal->mbovis_status>=1 && random < Se)
			{
			current_animal->mbovis_pcr_status=1 ;
			test_result_table[farm][PCR_calf+group] ++;	
//				if(test_result_table[farm][col_mastitis_tested_pos]==0 && current_animal->mastitis_detected==1)
//					{
//					test_result_table[farm][col_mastitis_tested_pos] = 1 ;//this round mastitis animal was tested Positive by PCR
//					}
			}
			test_result_table[farm][PCR_calf+group+num_test_vars]++;
			
		}
		else if(test_type==4)//slaughter  PCR
		{
			Se = PCR_slaughter;
			if(current_animal->mbovis_status>=1 && random < Se)
			{
			current_animal->mbovis_pcr_status=1 ;
			test_result_table[farm][slaughter_PCR] ++;	
//				if(test_result_table[farm][col_mastitis_tested_pos]==0 && current_animal->mastitis_detected==1)
//					{
//					test_result_table[farm][col_mastitis_tested_pos] = 1 ;//this round mastitis animal was tested Positive by PCR
//					}
			}
			test_result_table[farm][slaughter_PCR+num_test_vars]++;
			
		}
		else //milk PCR
		{
			Se = PCR_live;
			if(current_animal->mbovis_status>=2 && random < Se)
		  	{
		  	current_animal->mbovis_pcr_status=1 ;
			test_result_table[farm][PCR_adult_milk] ++;
//				if(test_result_table[farm][col_mastitis_tested_pos]==0 && current_animal->mastitis_detected==1)
//					{
//					test_result_table[farm][col_mastitis_tested_pos] = 1 ;//this round mastitis animal was tested Positive by PCR
//					}
		  	}
		  	test_result_table[farm][PCR_adult_milk+num_test_vars]++;
		}
	

	
	}
}

/*=========Testing BULK MILK===================================================*/
/*What matters? Desnity of m bovis in milk.
Milk from a clinical cow has density of 1/15.
Dilution factor says how much dilution from 1/15 can be detected by qPCR?*/
void test_bulk(int tank_discard, double* density, double dilution, double PCR_live, int** test_result_table, int farm)
{
	int num_test_vars = 10;
	
	//printf("bulk testing\n");
	double obs_dilution = (*density)*15 ; //Density divided by 1/15
	//printf("Test is %d obs_dilution is %f Dilution is %f\n",tank_discard,obs_dilution,dilution);
	 
	double random = (double)rand()/(double)RAND_MAX;
	if(obs_dilution > dilution) //more than threshold
	{
		if(random < PCR_live)
		{
			test_result_table[farm][8+tank_discard]++;
		//	printf("PCR pos\n");
		}
	}
	test_result_table[farm][8+tank_discard+num_test_vars]++;
}
/*==========Exporting test_result_table==============================================*/
int write_test_result(char* ResultFile,int* test_result_table,int nrow_testschedule,int n_column_output)
{
	FILE *Out = fopen(ResultFile, "w") ;
	//int line_num;
	int col_num;
	
//	for (line_num = 0 ; line_num < nrow_testschedule; line_num ++)
//	{
		for (col_num = 0;col_num <n_column_output ; col_num++ )
		
		{
	    fprintf(Out,"%d,",test_result_table[col_num]);
        }
        fprintf(Out,"\n");
  //  }
	fclose(Out);
	return 0;
}
