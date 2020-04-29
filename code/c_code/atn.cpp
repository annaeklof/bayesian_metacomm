/* Code to generate webs using the allometric trophic network model */ 
 

#include "atn.h"

// *********************** GLOBAL VARIABLES AND INPUT PARAMETERS  *********************************  

//  basic simulation variables:
int webs; 		// total number of webs
int S;          // total number of species
int S_c;        // number consumers (animals, predators)
int S_b;        // number basal species (plants)

//  parameters that determine the network topology:
double zeta_c   = 6;           // log_10 of the mean of the consumer (animal) body masses
double sigma_c	= 3;           // width of the distribution of consumer (animal) body masses
double cutoff_c	= 1e5;         // half relative cutoff for distribution of consumer body masses
double zeta_b	= 5;           // log_10 of the mean of the basal (plant) body masses
double sigma_b	= 3;           // width of the distribution of basal (plant) body masses
double cutoff_b = 1e5;         // half relative cutoff for distribution of basal body masses

double m_p_min = 0;            // minimal and maximal log10 body masses in case of uniform distributions
double m_a_min = 2;
double m_p_max = 6;
double m_a_max = 12;

double cutoff  = 0.01;         // cutoff of the Ricker curve for setting a link between predator and prey
double R_opt   = 100;          // optimal predator-prey body-mass ratio
double g	   = 2;            // width of the Ricker curve (higher g-value -> narrower curve)

double f_herbiv = 0.0;//0.50;  // fraction of species that are strict herbivores
double f_pred 	= 0.0;         // fraction of species that are strict predators

// flags: 
int UNIFORM_MASSES 	= 1;      // 0: body masses drawn from log-normal distribution; 1: log-body masses drawn from uniform distribution (PDEF_STRUCTURE)
int CANNIBALISM    	= 0;      // 0: all cannibalism links are explicitly removed // do we need this here ??? (PDEF_STRUCTURE) 

//~ char* directory 		= getenv("OUTPUTDIR"); 						// diretory for output files
//~ char* inputdir			= getenv("INPUTDIR"); 							// diretory for input files

int seed;

char* C_bash_seed		= getenv("SEED");            					// read in seed from bash
char* C_bash_Sc	    	= getenv("SC");            					    // read in file name from bash
char* C_bash_Sb	    	= getenv("SB");            					    // read in file name from bash
char* name_Sc 		    = getenv("SC");            						// read in web from bash
char* name_Sb 		    = getenv("SB");            						// read in web from bash
char* it 		    	= getenv("ITERATION");            				// read in web from bash

char directory[256]	 = "/homes/jh57masa/github/bayesian_metacomm_round2/data/c_webs/";  			// directory for the output files       							// name completion for output files

//  ********** Main function **********
int main(int argc, char* argv[])
{

	// get parameters from bash script
	seed					= atoi(C_bash_seed); 	
	S_c						= atoi(C_bash_Sc); 							// consumer species 
	S_b						= atoi(C_bash_Sb); 							// plant species
		
	gsl_rng_default_seed = seed;                  						// seed the random number generator 
    gsl_rng *r=gsl_rng_alloc(gsl_rng_default);
	
	S = S_b + S_c; 														// total number of species			
	gsl_matrix *Ap=gsl_matrix_calloc(S,S);                   	  		// adjacency matrix
	gsl_vector *mass=gsl_vector_calloc(S);                     			// mean body masses of the species
	pdef_structure(r,Ap,mass);                                  		// generate (random) food web structure
	output(Ap,mass);													// write food web structure (and bodymasses) to files
	
//  *********** Free memories ***********
    gsl_matrix_free(Ap);
    gsl_vector_free(mass);
    gsl_rng_free(r);

    return 0;
}


//	********** generate the body-mass based network structure **********
static void pdef_structure(gsl_rng *r,gsl_matrix *Ap, gsl_vector *mass)
{
	int i,j,flag=0;
	double Wkeit,sigma_i,zeta_act;
	double temp1,temp2,m_opt,a_max,R;
	double m_crit;

	while(flag==0)
	{
		flag=1;
		gsl_vector_set_all(mass,0);
		gsl_matrix_set_zero(Ap);

//	********* determine body masses *********
		if(UNIFORM_MASSES == 0)
		{
			zeta_act=zeta_c*log(10);			// calculate log(mean) from log10(mean)
			for(i=0;i<S_c;i++)
			{
				temp1=gsl_ran_lognormal(r,zeta_act,sigma_c);	
				if(temp1>exp(zeta_act)/cutoff_c&&temp1<exp(zeta_act)*cutoff_c)	// check whether mass is within certain boundaries
					gsl_vector_set(mass,S_b+i,temp1);		
				else
					i--;
			}
			gsl_sort_vector(mass);

			gsl_vector_view mass_b_vec=gsl_vector_subvector(mass,0,S_b);
			gsl_vector *mass_b=&mass_b_vec.vector;

			zeta_act=zeta_b*log(10);			// calculate log(mean) from log10(mean)
			for(i=0;i<S_b;i++)
			{
				temp1=gsl_ran_lognormal(r,zeta_act,sigma_b);	
				if(temp1>exp(zeta_act)/cutoff_b&&temp1<exp(zeta_act)*cutoff_b)	// check whether mass is within certain boundaries
					gsl_vector_set(mass_b,i,temp1);							
				else
					i--;
			}
			gsl_sort_vector(mass_b);
		}
		else
		{
			for(i = 0; i< S_c; i++)
				gsl_vector_set(mass,S_b+i,pow(10.,gsl_ran_flat(r,m_a_min,m_a_max)));
			gsl_sort_vector(mass);

			gsl_vector_view mass_b_vec=gsl_vector_subvector(mass,0,S_b);
			gsl_vector *mass_b=&mass_b_vec.vector;
			
			for(i = 0; i< S_b; i++)
				gsl_vector_set(mass_b,i,pow(10.,gsl_ran_flat(r,m_p_min,m_p_max)));
			gsl_sort_vector(mass_b);
		}

//	********** fill the adjacency matrix with Ricker attack rates *************
		for(i=0;i<S_c;i++)
		{
			temp1=gsl_vector_get(mass,S_b+i);
			a_max=pow(gsl_vector_get(mass,S_b+i)/R_opt,0.25);

			for(j=0;j<S;j++)
			{
				temp2=gsl_vector_get(mass,j);
				R = temp1/temp2;
				if(pow((R / R_opt) * exp(1 - (R / R_opt)), g) >= cutoff)
					gsl_matrix_set(Ap,S_b+i,j,1);
			}

			if(gsl_rng_uniform(r) < (f_herbiv + f_pred))					// combined probability to mess with this species' links
			{
				if(gsl_rng_uniform(r) < f_herbiv/(f_herbiv+f_pred))		// either make it a strict herbivore...
				{
					if(UNIFORM_MASSES == 0)
						m_crit = pow(10,zeta_b) * cutoff_b * R_opt;
					else
						m_crit = pow(10,m_p_max) * R_opt;
						
					if(gsl_vector_get(mass,S_b+i) < m_crit)
					{				
						for(j = S_b; j < S; j++)
							gsl_matrix_set(Ap,S_b+i,j,0);					// and remove all links from non-plant resources
					}
				}
				else													// ... or make it a strict carnivore 
				{
					for(j = 0; j < S_b; j++)
						gsl_matrix_set(Ap,S_b+i,j,0);					// and remove all links from plant resources
				}
			}

			gsl_vector_view tempp=gsl_matrix_row(Ap,S_b+i);				// reject networks with consumers or predators without prey
			flag=flag*(1-gsl_vector_isnull(&tempp.vector));
		}
		
		
		for(i=0; i<S_b; i++)
		{
			gsl_vector_view tempp = gsl_matrix_column(Ap,i);				// reject networks with uncontrolled basal species
			flag = flag*(1-gsl_vector_isnull(&tempp.vector));
		}
		
	}
	show_matrix(Ap,S,S);

	return;
}


//	********** write adjacency matrix and body masses to file **********
static void output(gsl_matrix *Ap, gsl_vector *mass)
{
    FILE *file1,*file2; 
     
    char web_out[99] = "web";
	strcpy(web_out,directory);
    strcat(web_out,"web");
	strcat(web_out,name_Sb); 
	strcat(web_out,name_Sc);
	strcat(web_out,it);
    strcat(web_out,".out");
     
    file1 = fopen(web_out,"w");
    
// file 1: atn_*.out
    for(int i=0; i<S; i++)
    {
        for(int j=0; j<S; j++)
            fprintf(file1,"%g ",gsl_matrix_get(Ap,i,j));
        fprintf(file1,"\n");
    }
    
    fclose(file1);
 
// file 2: bodymass_*.out
	char mass_out[99] = "bodymass";
    strcpy(mass_out,directory);
    strcat(mass_out,"bodymass"); 
	strcat(mass_out,name_Sb); 
	strcat(mass_out,name_Sc);
	strcat(mass_out,it);
	strcat(mass_out,".out");
    
    file2 = fopen(mass_out,"w");
    
    fprintf(file2, "%s,%s\n", "species", "bodymass","if.basal");
	
	for(int i=0; i<S; i++)
		fprintf(file2,"%d,%.9g,d\n", i, gsl_vector_get(mass,i),1); 
	for(int i=0; i<S; i++)
		fprintf(file2,"%d,%.9g,d\n", i, gsl_vector_get(mass,i),0); 
	
	fclose(file2);
    return;
}



//  ***********************************************************************************************************
//  ********** helper functions that prints a matrix/vector in the standard output (e.g. a terminal) **********
//  ***********************************************************************************************************
/*
1. Function show_matrix prints a matrix with m rows and n cols in the terminal
2. Function show_vector prints a vector with with the length m in the terminal
*/

static void show_matrix(gsl_matrix *M, int rows, int cols)
{

    int i, j;
    for(i=0; i<(rows); i++)
    {
        for(j=0; j<(cols); j++)
        {

            if(j==cols-1)
                printf("%.3g \n",
                       (double)gsl_matrix_get(M, i, j));
            else
                printf("%.3g \t",
                       (double)gsl_matrix_get(M, i, j));
        }
    }
    printf("\n"); 
    return;
}

static void show_vector(gsl_vector *A, int Num)
{
    int i,j;
    for(i = 0; i<Num; i++)
    {
        printf("%g \n",gsl_vector_get(A,i));
    }
    printf("\n"); 
    return;
}

