 /*===========================================================
Mycoplasma project
1. Demographic model - 2018/10
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
#define submission_prop 0.86
#define conception_AI 0.48
#define conception_bull 0.55
#define mating_week_AI 6
#define mating_week_bull 4
#define PSC 30 //Planned start of calving, set as 1 August which is 30 days
#define PSM 92 //Planned start of mating, set as 1 October
#define time_first_heat_min 10  // days until the first oestrus minimum value
#define time_first_heat_max 49  // days until the first Oestrus maximum value
#define interval_heat_min 18  // interval between oestrus events minimum
#define interval_heat_max 24  // interval between oestrus events maximum
#define av_gestation 282
#define error_gestation 10
//#define calf_mortality 0.041
//#define R1_mortality 0.017
//#define R2_mortality 0.017
//#define mixed_mortality 0.017
#define weaning_wks 13 
#define sim_years 3
#define num_extra_animal_per_year 150

/*Define tables imported*/
int num_cull_sell_steps = 62 ; //62 steps for culling and death
int num_mortality_steps = 19 ;
int n_column_cull_sell_rate = 3; // days, cull_rate, sell_rate
int n_column_mortality_rate = 2; // day and mortality
int n_column_List_mng_status = 2 ; //expand it later. Now total num of animals in each mng 
int num_column_NumGrpAnimal = 3650;
int temp_num_animal;
// and sum of rates in each mng, but later will include no.sus/exposed/infectious/latent
char CullSellRatesFile[] = "E:/ARATA/Documents/Research/Massey/Postdc_Massey/MBOVIS/Model/CullSellRatesFile.csv" ;
char MortalityFile[] = "E:/ARATA/Documents/Research/Massey/Postdc_Massey/MBOVIS/Model/mortality.csv" ;
char NumberAnimalDataFile[] = "E:/ARATA/Documents/Research/Massey/Postdc_Massey/MBOVIS/Model/NumberAnimalDataFile.csv" ;
/*Define general parameters*/



/* Define calf related parameters*/


/* Define heifer related parameters*/ 

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
int dry_day = 304;
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
int current_akey = 0 ;
double today_date = 0 ;
//int fortnight_num = floor(today_date/14) ; //0-364 days, Day 364 becomes week 26, need to make it wek 25
int year = 0;
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
num_culled = &var_cull;
num_sold = &var_sold;
num_death = &var_death;
struct animal_node* fake_animal;
fake_animal = (struct animal_node*)malloc(sizeof( struct animal_node ));

printf("Starts");	
double** List_mng_status =  (double**)malloc( sizeof(double *) *7); //modify if more than one herd exists	
struct animal_node **animal_node_pointer = (struct animal_node**)malloc( sizeof(struct animal_node*) * length_animal_pointer);
double **cull_sell_rate = (double**)malloc( sizeof(double *) *num_cull_sell_steps );
double **mortality = (double**)malloc( sizeof(double *) *num_mortality_steps );
double **NumGrpAnimal = (double**)malloc( sizeof(double *) *(id_bull_group+2) );
/*=======================================================================
Set up linked list for the management group
========================================================================*/
struct animal_node* FarmGroupList[id_bull_group+1];//Just a single pointer, if get confused, check F&Q why
struct animal_node* current_animal;
struct animal_node* next_animal;

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
for(i=0; i< herd_size * replacement_prop; i++ )
{
	mng_group = 1 ;
	 //R1 heifers should be between 314 and 334 days old
	// between 314 and 334
	current_age =  rand()%21+314 ;
	current_index_cull_sell = 6;
	current_index_mortality = 2 ;
	next_cull_change_date = cull_sell_rate[current_index_cull_sell+1][0] - current_age ;
	
	animal_node_pointer[current_akey] =(struct animal_node*)malloc(sizeof(struct animal_node)) ;
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
	cull_sell_rate[current_index_cull_sell][1]+cull_sell_rate[current_index_cull_sell][2] + 
	mortality[current_index_mortality][1] ;
	//then mortality
	List_mng_status[mng_group][1] = List_mng_status[mng_group][1] + markov_rate ;
	
	animal_node_pointer[current_akey]->sum_markov_rate = markov_rate ;
	add_animal_group(FarmGroupList,mng_group, animal_node_pointer[current_akey]) ;
	/*============ADDING EVENTS TO CHANGE CULLING iNDEX==========================*/
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
	current_akey ++;
	List_mng_status[mng_group][0] ++;


}
/*ADD EVENT WHEN R1 MOVES INTO R2 FOR YEAR 1*/
//usually it's when the youngest retained calf becomes 13 weeks old
//but in Year 1 there is no calf, so have to set a date
//that is when youngest R1 reaches 1yr + 13 weeks old
//youngest here is 314 days, means 456 days which is 140 days later
new_event = (struct event_node*)malloc(sizeof( struct event_node ));
new_event->event_type = 8;//moving C to R1 and R1 to R2
new_event->next_node = NULL ;
new_event->animal = fake_animal ;
add_event_node(event_day,140, new_event) ;


printf("R1 added");
/*==================R2 HEIFER================================================================*/
for(i=0; i< herd_size * replacement_prop; i++ )
{
	mng_group = 2 ;
	current_age = rand()%20+322+365;
	//R2 heifers should be between 687 and 706 days old
	if(current_age <= 700)
	{
		current_index_cull_sell = 12;	
	}
	else
	{
		current_index_cull_sell = 13;
	}
	next_cull_change_date = cull_sell_rate[current_index_cull_sell+1][0] - current_age ;
	current_index_mortality = 3 ;
	next_mortality_change_date = mortality[current_index_mortality+1][0] - current_age ;
	
	animal_node_pointer[current_akey] =(struct animal_node*)malloc(sizeof(struct animal_node)) ;
	animal_node_pointer[current_akey]->akey = current_akey;
	animal_node_pointer[current_akey]->age_day = current_age ;
	animal_node_pointer[current_akey]->index_cull_sell = current_index_cull_sell ;
	animal_node_pointer[current_akey]->index_mortality = current_index_mortality ;
	animal_node_pointer[current_akey]->group = mng_group ; //
	animal_node_pointer[current_akey]->pregnant_status = 1;
	animal_node_pointer[current_akey]->num_births = 0 ;
	animal_node_pointer[current_akey]->present = 1 ;
	animal_node_pointer[current_akey]->previous_node = NULL;
	animal_node_pointer[current_akey]->next_node = NULL;
	add_animal_group(FarmGroupList,mng_group, animal_node_pointer[current_akey]) ;
	
	
	markov_rate = 
	cull_sell_rate[current_index_cull_sell][1]+cull_sell_rate[current_index_cull_sell][2] + 
	mortality[current_index_mortality][1] ;
	
	
	List_mng_status[mng_group][1] = List_mng_status[mng_group][1] + markov_rate ;
	
	animal_node_pointer[current_akey]->sum_markov_rate = markov_rate ;
	
	if(i<herd_size * replacement_prop *calv_3weeks)
	{
		calving_date =  today_date + rand()%21 + PSC ; //note today_date = 0
		next_heat_date = calving_date + rand()%interval_first_heat + time_first_heat_min + rand()%interval_heat + interval_heat_min;
	}
	else if(i<herd_size * replacement_prop *calv_6weeks)
	{
		calving_date = today_date + rand()%21 + 21 + PSC ;
		next_heat_date = calving_date + rand()%interval_first_heat + time_first_heat_min + rand()%interval_heat + interval_heat_min;	
	}
	else if(i<herd_size * replacement_prop *calv_9weeks)
	{
		calving_date = today_date + rand()%21 + 42 + PSC ;
		next_heat_date = calving_date + rand()%interval_first_heat + time_first_heat_min + rand()%interval_heat + interval_heat_min;
	}
	else
	{
		calving_date = today_date + rand()%7 + 63 + PSC ;
		next_heat_date = calving_date + rand()%interval_first_heat + time_first_heat_min + rand()%interval_heat + interval_heat_min;
	}
	animal_node_pointer[current_akey]->calving_date = calving_date;
	animal_node_pointer[current_akey]->next_heat_date = next_heat_date ;
	
	/*========== ADDING EVENTS========================================================*/
	//create calving event: calving, heat and mortality and culling index change
	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
	new_event->event_type = 1;//calving
	new_event->akey = current_akey;
	new_event->animal = animal_node_pointer[current_akey]; //let's see if this works
	new_event->next_node = NULL ;
	add_event_node(event_day,calving_date, new_event) ;
	//create heat event
	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
	new_event->event_type = 2;//heat
	new_event->akey = current_akey;
	new_event->animal = animal_node_pointer[current_akey]; //let's see if this works
	new_event->next_node = NULL ;
	add_event_node(event_day,next_heat_date, new_event) ;
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
	
	animal_node_pointer[current_akey] =(struct animal_node*)malloc(sizeof(struct animal_node)) ;
	animal_node_pointer[current_akey]->group = mng_group ;
	animal_node_pointer[current_akey]->akey = current_akey;
	animal_node_pointer[current_akey]->pregnant_status = 1;
	animal_node_pointer[current_akey]->present = 1 ;
	animal_node_pointer[current_akey]->previous_node = NULL;
	animal_node_pointer[current_akey]->next_node = NULL;
	
	add_animal_group(FarmGroupList,mng_group, animal_node_pointer[current_akey]) ;
//printf("Animal added to group") ;
		/*ASSIGN A RANDOM AGE*/
		if(temp_prop<prop_age_2)
		{
		current_age = rand()%20+314+365*2;
		animal_node_pointer[current_akey]->age_day =  current_age ;
		//age between 1044-1063
		current_index_cull_sell = 18 ;
		current_index_mortality = 4 ;
		age_indicator = 2 ;
		age_cat = 3 ;
		//printf("Age 2");
		}
		/*AGE 3/4*/
		else if(temp_prop<sum_age_prop_4)
		{
		age_indicator = rand()%2+3 ;//3 or 4
		current_age = rand()%20+314+365*age_indicator;//1409 and 1428
		animal_node_pointer[current_akey]->age_day = current_age; 
		
		age_cat = 4 ; 
		if(current_age<=1427)
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
		age_indicator = rand()%3+5 ;
		age_cat = 5 ;
		current_age = rand()%20+314+365*age_indicator;
		animal_node_pointer[current_akey]->age_day =  current_age ;
		if(current_age<=2155)
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
		current_age = rand()%20+314+365*age_indicator;
		animal_node_pointer[current_akey]->age_day = current_age ;
		if(current_age<=3247)
		{
		current_index_cull_sell = 55 ;
		current_index_mortality = 17 ;
		}
		else
		{
		current_index_cull_sell = 56 ;
		current_index_mortality = 18 ;
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
	cull_sell_rate[current_index_cull_sell][1]+cull_sell_rate[current_index_cull_sell][2] + 
	mortality[current_index_mortality][1] ;
	animal_node_pointer[current_akey]->sum_markov_rate = markov_rate;
	
	List_mng_status[mng_group][1] = List_mng_status[mng_group][1] + markov_rate ;
	List_mng_status[mng_group][0] ++;
	//printf("Number is %d",(int)List_mng_status[mng_group][0]);
	/*REPRODUCTION PARAMETERS*/
	animal_node_pointer[current_akey]->num_births = age_indicator - 1 ;// assuming cows have claved every year
	if(i<herd_size * replacement_prop *calv_3weeks)
	{
		calving_date =  today_date + rand()%21 + PSC ;
		next_heat_date = calving_date + rand()%interval_first_heat + time_first_heat_min + rand()%interval_heat + interval_heat_min;
	}
	else if(i<herd_size * replacement_prop *calv_6weeks)
	{
		calving_date = today_date + rand()%21 + 21 + PSC ;
		next_heat_date = calving_date + rand()%interval_first_heat + time_first_heat_min + rand()%interval_heat + interval_heat_min;	
	}
	else if(i<herd_size * replacement_prop *calv_9weeks)
	{
		calving_date = today_date + rand()%21 + 42 + PSC ;
		next_heat_date = calving_date + rand()%interval_first_heat + time_first_heat_min + rand()%interval_heat + interval_heat_min;
	}
	else
	{
		calving_date = today_date + rand()%7 + 63 + PSC ;
		next_heat_date = calving_date + rand()%interval_first_heat + time_first_heat_min + rand()%interval_heat + interval_heat_min;
	}
	animal_node_pointer[current_akey]->calving_date = calving_date;
	animal_node_pointer[current_akey]->next_heat_date = next_heat_date ;
	//printf("Adding event");
	/*ADD EVENTS*/
	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
	new_event->event_type = 1;//calving
	new_event->akey = current_akey;
	new_event->animal = animal_node_pointer[current_akey]; //let's see if this works
	new_event->next_node = NULL ;
	add_event_node(event_day,calving_date, new_event) ;
	
	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
	new_event->event_type = 2;//heat
	new_event->akey = current_akey;
	new_event->animal = animal_node_pointer[current_akey]; //let's see if this works
	new_event->next_node = NULL ;
	add_event_node(event_day,next_heat_date, new_event) ;
	
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
printf("Dry is %d",(int)List_mng_status[id_dry_group][0]);
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


/*============================================================================================
Day procedes
=============================================================================================*/
while(today_date<=sim_years*365)
{
//	printf("today is %lf",today_date) ;
	//visualize_list(event_day,0) ;
/*====START OF DAY========================================================*/
next_non_markov_date = ceil(today_date);

  while (event_day[next_non_markov_date] == NULL)
            {//loop2
		    next_non_markov_date++; // go to the next day
	 //   printf("Next non Markov is %d",next_non_markov_date) ;
	 //here I need to cancel events that are associated with animals that are removed
	 if(event_day[next_non_markov_date]!=NULL)
	 {
	 	printf("NonM is %d",next_non_markov_date);
	 	current_event = event_day[next_non_markov_date];
	 	//trial
	 	current_animal = current_event->animal; 	
	 	previous_event = current_event;
	 	while(current_event->animal==NULL)
	 	{
	 		printf("This event does not exist anymore!");
	 		system("pause");
	 		current_event = current_event->next_node;
	 		free(previous_event);
	 		previous_event = current_event;
	 		if(current_event==NULL)
	 		{
	 			printf("current event becomes null");
	 			break;
	 			
			 }
		 }
	 	event_day[next_non_markov_date] = current_event ;
	 	printf("Head updated");
	while(current_event!=NULL)
	{
	//	printf("Type is %d",current_event->event_type) ;
	if(current_event->next_node == NULL)
	{
	//	printf("C");
		break;
	}
	while(current_event->next_node->animal==NULL)
		{
			if(current_event->next_node->next_node!=NULL)
			{
			//	printf("A1\n");
				next_event = current_event->next_node->next_node;
				free(current_event->next_node);
				current_event->next_node = next_event ;
			//	printf("A2\n");
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
    printf("NonM is %d",next_non_markov_date);
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
	printf("Loop started\n");
     current_event = event_day[next_non_markov_date];
     while(current_event!=NULL)
     {
     	//printf("event is not null\n") ;
     	printf("Event is %d", current_event->event_type) ;
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
     	if(current_event->event_type==1)
     	{
		    current_animal = current_event->animal ;
     		current_grp = current_animal->group ;
     		printf("akey is %lld",current_animal->akey);
     //		printf("event is %d\n",current_event->event_type);
     		if(current_grp==4)
     		{
     			printf("Dry cow is %d\n",(int)List_mng_status[id_dry_group][0]);
     			dry_calved++;
			 }
			 if(current_animal->present == 0)
			 {
			 	printf("wait! Ghost cow calving!");
			 	//free(current_event->animal);//this should free but it's not freeing animal
			 	
			//free(animal_node_pointer);
		
			printf("akey is %lld",animal_node_pointer[current_animal->akey]->akey);
		
			 //	current_animal = NULL;
			 	if(current_animal==NULL)//this is not true
			 	{
			 		printf("animal is free");
				 }
			 	if(current_event->animal==NULL) //this is not true because free'd memory is not null
			 	{
			 		printf("Yes nulled");
				 }
			 	printf("akey is %lld",current_event->animal->akey);
			 	printf("akey is %lld",animal_node_pointer[current_animal->akey]->akey);
			 	system("pause");
			 }
			
     		//this animal anyway move into Lactation group[3]
     		
     		//if heifer -1 from heifer, if adult -1 from dry
     
     		List_mng_status[current_grp][0]--;
     			if(current_grp==4)
     		{
     			printf("Dry cow is %d\n",(int)List_mng_status[id_dry_group][0]);
     		}
     		printf("Dry calved cow is %d\n",dry_calved);
     		if(List_mng_status[current_grp][0]<0)
     		{
     			printf("Warning! Number gets below 0 in group %d",current_grp);
     			system("pause");
			 }
     		
     		List_mng_status[id_lact_group][0]++;//whgether they are heifer or adult, goes to lact
     		current_animal->pregnant_status = 0;
     		
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
				add_event_node(event_day,next_heat_date, new_event) ;
			}
			
			
     	if(List_mng_status[0][0]>herd_size * calf_keep_prop)
     	{
     		//do we record the number of bobbied ?
     		//maybe not, but in future calf movement should be considered
     		printf("Do not keep any more calf!\n");
		 }
		 else if((double)rand()/(double)RAND_MAX<=calf_female_prop) //proportion of female and male 0.5
		 {
		 // Keep at farm new born calves
     	animal_node_pointer[current_akey] =(struct animal_node*)malloc(sizeof(struct animal_node)) ;
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
	//	printf("Calf size is %lf\n",List_mng_status[0][0]);
		//if this is the last animal to be kept
		if(List_mng_status[0][0]>=herd_size * calf_keep_prop)
		{
		new_event = (struct event_node*)malloc(sizeof( struct event_node ));
		new_event->event_type = 8; //moving to the next group
		new_event->next_node = NULL ;
		new_event->animal = fake_animal;
	//	printf("setting up move date") ;
		if(next_non_markov_date+91<sim_years*365)
		{
			add_event_node(event_day,(int)next_non_markov_date+91, new_event) ;//animals move to the next age group
		}
	//	printf("last calf to add\n");
		}
		markov_rate = 
		cull_sell_rate[current_index_cull_sell][1]+cull_sell_rate[current_index_cull_sell][2] + 
		mortality[current_index_mortality][1] ;
		animal_node_pointer[current_akey]->sum_markov_rate = markov_rate ;
		//List_mng_status[0][1] = List_mng_status[0][1] + markov_rate ;
	
		add_animal_group(FarmGroupList,0, animal_node_pointer[current_akey]) ;
		
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
     	
		 } //calving done
		 	
	/*=============HEAT============================================*/
		 if(current_event->event_type==2)
		 {
	//	 	printf("event is %d\n",current_event->event_type);
		 	//Heat
		
		 		if((double)rand()/(double)RAND_MAX<=submission_prop )	//submission happen
		 		{
		 			int conception_rate;
		 			if(today_date-365*year-PSM <=7*mating_week_AI)
		 			{
		 				conception_rate = conception_AI ;
					 }
					 else
					 {
					 	conception_rate = conception_bull ;
					 }
		 			if((double)rand()/(double)RAND_MAX <= conception_rate)
					 {
					 	struct animal_node* current_animal;
					 	current_animal = current_event->animal;
					 	current_event->animal->pregnant_status = 1 ;
					 	//add calving event
					 	calving_date = rand()%error_gestation*2 - error_gestation + today_date ;
					 	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
						new_event->event_type = 1;//calving
						new_event->akey = current_event->akey;
						new_event->animal = current_event->animal; //let's see if this works
						new_event->next_node = NULL ;
								if(calving_date<sim_years*365)
			{
			add_event_node(event_day,calving_date, new_event) ;
		}
						
					}
						
				 } // then otherwise set up next heat date
				 else
				 {
				 		next_heat_date = rand()%interval_heat + interval_heat_min ;
				 		new_event = (struct event_node*)malloc(sizeof( struct event_node ));
						new_event->event_type = 2;//calving
						new_event->akey = current_event->akey;
						new_event->animal = current_event->animal; //let's see if this works
						new_event->next_node = NULL ;
						if(calving_date<sim_years*365)
			                  {
			            add_event_node(event_day,next_heat_date, new_event) ;
		                      }
					
				 }
			 }
	/*===============CHANGE IN CULLING INDEX===========================*/
	if(current_event->event_type==3)
		 {
	//	 	printf("event is %d\n",current_event->event_type);
		 	current_animal = current_event->animal;
		  current_index_cull_sell = current_animal->index_cull_sell;	
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
		cull_sell_rate[current_index_cull_sell][1]+cull_sell_rate[current_index_cull_sell][2] + 
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
	if(current_event->event_type==4)
		 {
	//	 printf("event is %d\n",current_event->event_type);	
		 current_animal = current_event->animal;
		 current_index_mortality =current_animal->index_mortality;
		 
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
		cull_sell_rate[current_index_cull_sell][1]+cull_sell_rate[current_index_cull_sell][2] + 
		mortality[current_index_mortality][1] ;	
		//mortality
	new_event = (struct event_node*)malloc(sizeof( struct event_node ));
	new_event->event_type = 4;//change in mortality index
	new_event->akey = current_animal->akey;
	new_event->animal = current_animal; //let's see if this works
	new_event->next_node = NULL ;
		if(current_index_mortality<=17)
		 {
		 	printf("index is %d",current_index_mortality);
		 next_mortality_change_date = next_non_markov_date+
		 mortality[current_index_mortality+1][0] - cull_sell_rate[current_index_mortality][0] ;
		 }
		  else
	     {
	     	printf("index is %d",current_index_mortality);
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
	//	printf("event is %d\n",current_event->event_type);
		if(FarmGroupList[2]!=NULL)
		{
			printf("There are still animals in R2!\n") ;
			//these animals will be culled
			current_animal = FarmGroupList[2];
			while(current_animal!=NULL)
			{
				next_animal = current_animal->next_node;
				remove_animal_group(FarmGroupList,id_r2_group,current_animal);
				free(current_animal);
				List_mng_status[id_r2_group][0]--;
				(*num_culled)++;
				current_animal = next_animal ;
			}
			FarmGroupList[2] = NULL;
		}
		if(FarmGroupList[1]==NULL)
		{
			printf("There is no animals in R1!\n");
		}
		else
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
			printf("R2 heifer is %d",temp_num_animal);
			
			if(FarmGroupList[0]==NULL)
			{
				printf("There is no animals in calf:YEAR %d\n",year) ;
			}
			else
			{
				if(List_mng_status[0][0]<=0)
				{
					printf("FarmGroupList has calves but not in List_mng");
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
						remove_animal_group(FarmGroupList,id_lact_group,current_animal);
						free(current_animal);
						(*num_culled)++;
						List_mng_status[id_lact_group][0]--;
					}
					current_animal = next_animal ;
				}
				current_animal = FarmGroupList[id_lact_group] ;//get the first animal again
				FarmGroupList[id_dry_group] = current_animal;
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
	   printf("now free event");
	   free(previous_event);///check if this works - previously this was in if(current_event!=NULL) blacket
	   printf("event was %d",previous_event->event_type);
	   system("pause");
	   event_day[next_non_markov_date] = current_event;
	   //printf("event pointing new event\n") ;
     }
     printf("day is over\n");
    for(i=0;i<id_bull_group+1;i++)
	{
		//printf("number is %lf",List_mng_status[i][0]);
	//	printf("index column is %d",index_column) ;
	NumGrpAnimal[0][index_column] = today_date ;
	NumGrpAnimal[i+1][index_column] = (int)List_mng_status[i][0] ;
	}
	index_column++;	
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
	//Calculate MARKOV Date
	for(i=0; i<id_bull_group+1; i++)
	{//for each management group
	printf("Management grp is %d",i);
	
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
//	printf("sum is %lf",sum_rate);
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
	 printf("markov date is %lf",day_to_markov);
	 printf("non markov is %d",next_non_markov_date);
/* NOW ASSESS WHETHER MARKOV EVENT OCCURS FIRST OR NOT*/
	if (next_non_markov_date>day_to_markov+today_date)
	{ // if markov comes first, choose markov events 
	   k++;
	   //printf("This is %d th markov",k);
	   today_date = day_to_markov+today_date ;
	   printf("ceil days is %lf",today_date);
	   //system("pause");
	  
	random_value =  (double)(rand()+1)/(double)(RAND_MAX+1)*sum_rate;
		i = 0;
	//now choose one of management group for the event
	while(random_value>accumu_rate)
	{
		accumu_rate = accumu_rate+ List_mng_status[i][1] ;
		i++;
		printf("i is %d\n", i) ;
	} //i is the group at which an event occurs	
	printf("group decided");
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
	//obtain rates for culling, death and sale
	if(random_value<=cull_sell_rate[current_animal->index_cull_sell][1])
	{
		//cull
		(*num_culled)++;
	//	printf("this is cull");
	//	system("pause");
	}
	else if(random_value<=cull_sell_rate[current_animal->index_cull_sell][1]+cull_sell_rate[current_animal->index_cull_sell][2])
	{
		//sale
		(*num_sold)++;
		//printf("this is sold");
	//	system("pause");
	}
	else
	{
		//death
		(*num_death++);
	//	printf("this is death");
	//	system("pause");
		
	}
	printf("i is %d",i);
	List_mng_status[i-1][0]--;
	if(List_mng_status[i-1][0]<0)
	{
		printf("Warning! Animals below 0");
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
		           printf("B ends");
				}	
		else // if next node is null
		        {//Removal Situation B
		   	    prev_animal -> next_node = NULL ;
		   	    printf("A ends");
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
        printf("D ends") ;
		}
		else 
		{//Situation A
	//		printf("E starts") ;
		FarmGroupList[i-1] = NULL ; 
		printf("C ends");
		}
		         
	}
// moving_animal->current_pro_id = -10 ; //means it's sent to slaughterhouse 
 current_animal->previous_node = NULL;
 current_animal->next_node = NULL ;
current_animal->present = 0;
	/*==Removing animal done=========================*/

	printf("free now");
	printf("akey is %lld",current_animal->akey) ;
//	free(current_animal) ; //free memory
	printf("free done");
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
    double cull_rate, sale_rate;
    
    /* READ LINES OF FARM FILE */
    for(line_num = 0; line_num < num_cull_sell_steps; line_num++)
      { 
         fscanf(temp, "%d,%lf,%lf",&day, &cull_rate, &sale_rate);
         
             cull_sell_rate[line_num][0] = day;
             //printf("day is %d",day) ;
             cull_sell_rate[line_num][1] = cull_rate/1000;
             cull_sell_rate[line_num][2] = sale_rate/1000;
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
			 printf("Event is %d", current_node1 -> event_type);
			 // system("pause") ;
		   // }
              
              current_node1 = current_node1 -> next_node; 
              
          }  
       printf("\n");
   }
   
 

}

