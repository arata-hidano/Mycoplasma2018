 /*===========================================================
Mycoplasma project
1. Demographic model - 2018/10
2. Starting date June 1
=============================================================*/

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
#define PSC 62 //Planned start of calving, set as 1 August which is 62 days
#define PSM 144 //Planned start of mating, set as 22 October
int dry_day = 334; //May 1st
int r1_initial_age = 304; //As of June 1st, they are max 304 days old

/*Other parameters*/
#define herd_size 400
#define num_total_herd 1
#define id_calf_group 0
#define id_r1_group 1
#define id_r2_group 2
#define id_lact_group 3
#define id_dry_group 4
#define id_sick_group 5
int const id_bull_group = 6;
#define calf_female_prop 0.5
#define calf_keep_prop 0.3
#define replacement_prop 0.22
#define calv_3weeks 0.6
#define calv_6weeks 0.87
#define calv_9weeks 0.98
#define submission_prop 0.9
#define conception_AI 0.48
#define conception_bull 0.55
#define mating_week_AI 6
#define mating_week_bull 4
#define mating_period 10

#define time_first_heat_min 10  // days until the first oestrus minimum value
#define time_first_heat_max 49  // days until the first Oestrus maximum value
#define interval_heat_min 18  // interval between oestrus events minimum
#define interval_heat_max 24  // interval between oestrus events maximum
#define av_gestation 282
#define error_gestation 10
#define heifer_puberty 365 // reaches puberty and starts heat
#define heifer_puberty_error 30 //margin for reaching puberty
//#define calf_mortality 0.041
//#define R1_mortality 0.017
//#define R2_mortality 0.017
//#define mixed_mortality 0.017
#define weaning_wks 10 
#define sim_years 3
#define num_extra_animal_per_year 150

#define column_cull_empty 1
#define column_sell_empty 2
/*Define tables imported*/
int num_cull_sell_steps = 62 ; //62 steps for culling and death
int num_mortality_steps = 19 ;
int n_column_cull_sell_rate = 5; // days, cull_rate_empty, sell_rate_empty,cull_rate,sell_rate
int n_column_mortality_rate = 2; // day and mortality
int n_column_List_mng_status = 2 ; //expand it later. Now total num of animals in each mng 
int num_column_NumGrpAnimal = 15000;
int temp_num_animal;
// and sum of rates in each mng, but later will include no.sus/exposed/infectious/latent
char CullSellRatesFile[] = "E:/ARATA/Documents/Research/Massey/Postdc_Massey/MBOVIS/Model/CullSellRatesFile2.csv" ;
char MortalityFile[] = "E:/ARATA/Documents/Research/Massey/Postdc_Massey/MBOVIS/Model/mortality.csv" ;
char NumberAnimalDataFile[] = "E:/ARATA/Documents/Research/Massey/Postdc_Massey/MBOVIS/Model/NumberAnimalDataFile.csv" ;
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
int length_animal_pointer = herd_size + num_extra_animal_per_year*sim_years ;

//remaining will give a birth between 63 (9weeks) and 70 (10weeks)

/* Define hospital group related parameters*/



/* STRUCTURE DECLARATIONS */  
struct animal_node {
      long long akey ; // animal id
      int age_day;   /* age in day*/
      int group;         // production type, 0 = calf, 1 = R1, 2 = R2, 3 = Lact, 4 = Dary, 5 = Sick, 6 = bull 
      int sex;          /* sex, female = 1 male = 0*/
      int breed;
      int mbovis_status;    /* disease status. Sus = 0, Exposed = 1, Infectious =2, Latent = 3*/
      int mastitis_status ; //non-mastitis = 0, subclinical mastitis = 1, clinical mastitis = 2
      int pregnant_status; // non-pregnant = 0 preganant to AI = 1 pregnant to bull = 2
      int milking_status; // non-milking age = 0, milking = 1, milking and dried = 2
      int calving_date; //date of calving (Day xx)
      int next_heat_date; //date of the next heat (Day xx)
      int num_births; // record how many births it already gave to
      double bovis_antibody;//MBOVIS antibody
      int present; //whether this animal exists on the farm or not
      long long current_pro_id ;
      int index_cull_sell ;
      int index_mortality ;
      double sum_markov_rate ;
      struct animal_node *previous_node ;/*pointer to the previous node*/
      struct animal_node *next_node; /* pointer to next animal*/
    
   }; 
   
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

double markov_rate ;
int index_column = 0;
int i, j, abc;
int calving_date, next_heat_date,next_non_markov_date ;
double updated_date ;
/*===============================================================================*/
/*Main starts from here*/
int main(void){
srand((unsigned)time(NULL));
int current_akey = 0 ;
double today_date = 0 ;
int year = 0;
int calendar_day = 0;

num_culled = &var_cull;
num_sold = &var_sold;
num_death = &var_death;
struct animal_node* fake_animal;
fake_animal = (struct animal_node*)malloc(sizeof( struct animal_node ));

printf("Starts");	
double** List_mng_status =  (double**)malloc( sizeof(double *) *7); //modify if more than one herd exists	
struct animal_node **animal_node_pointer = malloc( sizeof(struct animal_node*) * length_animal_pointer);
double **cull_sell_rate = (double**)malloc( sizeof(double *) *num_cull_sell_steps );
double **mortality = (double**)malloc( sizeof(double *) *num_mortality_steps );
double **NumGrpAnimal = (double**)malloc( sizeof(double *) *(id_bull_group+2) );
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
	List_mng_status[i] = (double*)malloc( sizeof(double) * n_column_List_mng_status);
	
	for(j = 0; j < n_column_List_mng_status; j++)
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

	      	    
/*Day 0
Add R1/R2/Dry to the linked list*/
//For now, assume no calves and heifers die and farmers only keep 
//as many calves/heifers as needed (i.e. herd size * replacement rate)


/*============ R1 HEIFER========================================================*/
//adding R1 heifer
for(i=0; i< herd_size * calf_keep_prop; i++ )
{
	mng_group = 1 ;
	 //R1 heifers should be between 284 and 304 days old
	current_age =  r1_initial_age - rand()%21 ;
	current_index_cull_sell = 6;
	current_index_mortality = 2 ;
	next_cull_change_date = cull_sell_rate[current_index_cull_sell+1][0] - current_age ;
	next_mortality_change_date = mortality[current_index_mortality+1][0] - current_age ;
		
	animal_node_pointer[current_akey] =malloc(sizeof(struct animal_node)) ;
	animal_node_pointer[current_akey]->akey = current_akey;
	animal_node_pointer[current_akey]->age_day = current_age ;
	animal_node_pointer[current_akey]->index_cull_sell = current_index_cull_sell ;
	animal_node_pointer[current_akey]->index_mortality = current_index_mortality ;
	animal_node_pointer[current_akey]->group = mng_group ; //
	animal_node_pointer[current_akey]->pregnant_status = 0;
	animal_node_pointer[current_akey]->num_births = 0 ;
	animal_node_pointer[current_akey]->present = 1 ;
	animal_node_pointer[current_akey]->previous_node = NULL;
	animal_node_pointer[current_akey]->next_node = NULL;
	
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


}




printf("R1 added");

	
/*==================R2 HEIFER================================================================*/
for(i=0; i< herd_size * calf_keep_prop; i++ )
{
	mng_group = 2 ;
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
	animal_node_pointer[current_akey]->index_cull_sell = current_index_cull_sell ;
	animal_node_pointer[current_akey]->index_mortality = current_index_mortality ;
	animal_node_pointer[current_akey]->group = mng_group ; //
	pregnant_status = 1;
	animal_node_pointer[current_akey]->pregnant_status = pregnant_status;
	animal_node_pointer[current_akey]->num_births = 0 ;
	animal_node_pointer[current_akey]->present = 1 ;
	animal_node_pointer[current_akey]->previous_node = NULL;
	animal_node_pointer[current_akey]->next_node = NULL;
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
}
printf("R2 added");
/*=====================MIXED AGE DRY COWS=============================================*/

int age_indicator ;
double temp_prop ;
int dry_calved = 0;
int dry_culled =0;
for(i=0; i< herd_size * (1-replacement_prop); i++ )
{
	mng_group = 4 ;//dry group when started
	temp_prop = (double)rand()/(double)RAND_MAX ; //get a random value between 0 and 1
	
	animal_node_pointer[current_akey] =malloc(sizeof(struct animal_node)) ;
	animal_node_pointer[current_akey]->group = mng_group ;
	animal_node_pointer[current_akey]->akey = current_akey;
	pregnant_status = 1 ;
	animal_node_pointer[current_akey]->pregnant_status = pregnant_status;
	animal_node_pointer[current_akey]->present = 1 ;
	animal_node_pointer[current_akey]->previous_node = NULL;
	animal_node_pointer[current_akey]->next_node = NULL;
	
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
	//printf("Number is %d",(int)List_mng_status[mng_group][0]);
	/*REPRODUCTION PARAMETERS*/
	animal_node_pointer[current_akey]->num_births = age_indicator - 1 ;// assuming cows have claved every year
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
system("pause");
/*ADD Dry-off event*/

for(i=0;i<sim_years;i++)
{
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
/*============================================================================================
Day procedes
=============================================================================================*/
while(today_date<sim_days)
{
	year = (int)floor(today_date/365);

	
	printf("today is %lf\n",today_date) ;
	
	//visualize_list(event_day,0) ;
/*====START OF DAY========================================================*/
next_non_markov_date = ceil(today_date);
//	if(next_non_markov_date==year*365)
//	{printf("R1 is %d\n",(int)List_mng_status[1][0]);
//	printf("R2 is %d\n",(int)List_mng_status[2][0]);
//	system("pause");
//	}
  while (event_day[next_non_markov_date] == NULL)
            {//loop2
		    next_non_markov_date++; // go to the next day
		    if(next_non_markov_date==sim_days)
		    {
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
  //  printf("NonM is %d",next_non_markov_date);
  if(next_non_markov_date>=sim_days)
  {
  	printf("Break");
  	break;
  }
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
     	if(current_event->event_type==1 && target_milk_met ==0 && current_event->animal->present==1)
     	{
     		
		    current_animal = current_event->animal ;
		    if(current_animal->pregnant_status==0)
		    {
		    	printf("Wait!She is not pregnant!\n");
		    	system("pause");
			}
     		current_grp = current_animal->group ;
     		
     //		printf("event is %d\n",current_event->event_type);
     		if(current_grp==4)
     		{
     		//	printf("Dry cow is %d\n",(int)List_mng_status[id_dry_group][0]);
     		//	dry_calved++;
			 }
			 if(current_animal->present == 0)
			 {
			 	printf("wait! Ghost cow calving!");
			 	printf("akey is %lld and grp is %d",current_animal->akey, current_animal->group);
			 		system("pause");
			 
			 }
			
     		//this animal anyway move into Lactation group[3]
     		
     		//if heifer -1 from heifer, if adult -1 from dry
     
     		List_mng_status[current_grp][0]--;
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
     		
     		List_mng_status[id_lact_group][0]++;//whgether they are heifer or adult, goes to lact
     		//remove R2 if there are 400 calved
		
		 	current_animal->pregnant_status = 0;
     		current_animal->sum_markov_rate = 
					cull_sell_rate[current_animal->index_cull_sell][column_cull_empty+2*pregnant_status]+
					cull_sell_rate[current_animal->index_cull_sell][column_sell_empty+2*pregnant_status] + 
					mortality[current_animal->index_mortality][1] ;
     	
		 	remove_animal_group(FarmGroupList,current_grp,current_animal);
     		current_animal->group = id_lact_group ; //change the group
     		add_animal_group(FarmGroupList,id_lact_group,current_animal);
     		/*ADd First heat event*/
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
			
			
     	if(List_mng_status[0][0]>herd_size * calf_keep_prop)
     	{
     		//do we record the number of bobbied ?
     		//maybe not, but in future calf movement should be considered
     		//printf("Do not keep any more calf!\n");
     		
		 }
		 else if((double)rand()/(double)RAND_MAX<=calf_female_prop) //proportion of female and male 0.5
		 {
		 // Keep at farm new born calves
     	animal_node_pointer[current_akey] =malloc(sizeof(struct animal_node)) ;
		animal_node_pointer[current_akey]->akey = current_akey;	
		animal_node_pointer[current_akey]->age_day = 0; //just born
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
		if(List_mng_status[0][0]==1)
		{
			printf("The first calving date is %d in YEAR %d\n",calendar_day, year);
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
		 //remove all R2 when the calved cows reached 400
		 	if(List_mng_status[id_lact_group][0]==herd_size)
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
				 		current_animal->present = 0;
				 		previous_animal = current_animal;
				 		current_animal = current_animal->next_node;
				 		remove_animal_group(FarmGroupList,id_r2_group,previous_animal);
				 		List_mng_status[id_r2_group][0]--;
				 		(*num_culled)++;
				 		
					 }
				 }
				
				if(List_mng_status[id_dry_group][0]<=0)
     			{
     				printf("No Dry animals remained");
				 }
				 else
				 {
				 	printf("Dry is %d",(int)List_mng_status[id_dry_group][0]);
				 	current_animal = FarmGroupList[id_dry_group];
				 	while(current_animal!=NULL)
				 	{
				 	//	printf("akey is %lld\n", current_animal->akey);
				 	//	if(current_animal->next_node==NULL)
				 	//	{
				 	//		printf("next node NULL");
					//	 }
					//	 else
					//	 {
					//	 	printf("next akey is %lld\n",current_animal->next_node->akey);
					//	 }
				 		current_animal->present = 0;
				 		previous_animal = current_animal;
				 		current_animal = current_animal->next_node;
				 		//printf("now akey is %lld",current_animal->akey);
				 		List_mng_status[id_dry_group][0]--;
				 		(*num_culled)++;
				 	
				 		remove_animal_group(FarmGroupList,id_dry_group,previous_animal);
				 		
				 		
					 }
				//	 if(FarmGroupList[id_dry_group]==NULL)
				//	 {
				//	 	printf("Now all dry animals gone\n");
				//	 }
				 }
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
					 	calving_date = av_gestation+ rand()%error_gestation*2 - error_gestation + next_non_markov_date ;
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
	/*=========Culling==================================================*/     

	/*===========Move to next subgroup========================================*/
	if(current_event->event_type==8)
	{
		target_calf_met=0;
		target_milk_met=0;
	//	printf("Moving Now R2 is %d",(int)List_mng_status[id_r2_group][0]);
	//		printf("Now Lact is %d",(int)List_mng_status[id_lact_group][0]);
	//		system("pause");
	//	printf("event is %d\n",current_event->event_type);
		if(FarmGroupList[2]!=NULL)
		{
			current_animal = FarmGroupList[2];
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
			FarmGroupList[2] = NULL;
		}
	//	if(FarmGroupList[1]==NULL)
	//	{
	//		printf("There is no animals in R1!\n");
	//	}
		if(FarmGroupList[2]==NULL) //now if R2 becomes 0
		{//
		
			current_animal = FarmGroupList[1] ;
			while(current_animal!=NULL)
			{
				current_animal->group = id_r2_group;
				current_animal = current_animal->next_node;
			}
			current_animal = FarmGroupList[1] ;//again get the first animal
			FarmGroupList[2] = current_animal ;
			FarmGroupList[1] = NULL;
			temp_num_animal = List_mng_status[1][0]; //R1 heifers
			List_mng_status[2][0] = temp_num_animal;
	//		printf("R2 heifer is %d",temp_num_animal);
			
			if(FarmGroupList[0]==NULL)
			{
				printf("There is no animals in calf:YEAR %d\n",year) ;
				system("pause");
			}
			else
			{
				if(List_mng_status[0][0]<=0)
				{
					printf("FarmGroupList has calves but not in List_mng");
					system("pause");
				}
				else
				{
					
					temp_num_animal = List_mng_status[0][0];
					List_mng_status[1][0] = temp_num_animal;
					printf("R1 heifer is %d",temp_num_animal);
					List_mng_status[0][0]=0; //all calves gone
					current_animal = FarmGroupList[0] ;
					while(current_animal!=NULL)
					{
					current_animal->group = id_r1_group;
					current_animal = current_animal->next_node;
					}
					current_animal = FarmGroupList[0] ;//get the first animal again
					FarmGroupList[1] = current_animal ;
					FarmGroupList[0] = NULL ;
				}
				
			}
		}
	}
	/*================Drying off=====================================*/
		if(current_event->event_type==9)
	{//cows go from 3 to 4 
	int num_lact = 0;
	int num_dry = 0;
	int num_culled_nonpregnant = 0;
	int num_pregnant_at_dry = 0;
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
				//system("pause");
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
				}
				current_animal = FarmGroupList[id_lact_group] ;//get the first animal again
				FarmGroupList[id_dry_group] = current_animal;
				while(current_animal!=NULL)
				{
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
					if(temp_num_animal<0)
     		{
     			printf("Warning! Number in Dry gets below 0");
     			system("pause");
			 }
			}
		}	
	}
	
	
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
	//printf("Now today date is %lf\n",today_date) ;
	
	
}
printf("Loop done");
/*===========SIMULATION ENDS=============================================*/
printf("Total cull is %d\n",*num_culled);
printf("Total sale is %d\n",*num_sold);
printf("Total death is %d\n",*num_death);
system("pause");
	write_number_animals(NumberAnimalDataFile,NumGrpAnimal,id_bull_group,index_column) ;

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
int* num_culled, int* num_sold, int* num_death)
{
// Calculate the sum of rate of M events
	double day_to_markov;
	struct animal_node* current_animal;
	struct animal_node* next_animal;
	struct animal_node* previous_animal;
	double current_rate;
	double sum_rate = 0.0;
	double accumu_rate = 0.0 ;
	int random_int;
	int i,j;
	int k = 0;
	int pregnant_status ;
	//Calculate MARKOV Date
	for(i=0; i<id_bull_group+1; i++)
	{//for each management group
	//printf("Management grp is %d",i);
	
	current_animal = FarmGroupList[i] ;
	List_mng_status[i][1] = 0; //reset
	while(current_animal!=NULL)
	{
		//printf("summing up rate");
		//visit each linked list and add up rate
		current_rate = current_animal->sum_markov_rate ;
	//	printf("Current rate is %lf",current_rate) ;
		List_mng_status[i][1] = List_mng_status[i][1]+current_rate;
	//	if(current_animal->next_node == NULL)
	//	{
	//		printf("Next node is NULL");
	//		}	
		current_animal = current_animal->next_node ;
	}
//	printf("Group %d is done",i) ;
	sum_rate = 	sum_rate + List_mng_status[i][1] ;
	printf("sum is %lf\n",sum_rate);
	}
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
		i = 0;
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

