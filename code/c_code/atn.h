/* Header to "create_adjacency_matrix.h" */ 

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>      // for input/output from or to files
#include <iostream>		// for input/output on terminal
#include <sstream>
#include <vector>
#include <cstdlib>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>				// random number generator
#include <gsl/gsl_randist.h>			// random number distributions
#include <gsl/gsl_blas.h>				// linear algebra routines
#include <gsl/gsl_linalg.h>				// linear algebra
#include <gsl/gsl_sort_vector.h>		// vector operations

static void pdef_structure(gsl_rng *r,gsl_matrix *Ap, gsl_vector *mass);
static void output(gsl_matrix *Ap, gsl_matrix *A, gsl_matrix *H, gsl_vector *mass, gsl_vector* D_Max); 
static void get_dipsersal_dist(gsl_rng *ran, gsl_vector *mass, gsl_vector *D_Max); 

double get_random_parameter(double mean, double sigma, double low_cutoff, double high_cutoff, gsl_rng *r);
static void set_parameters(gsl_rng *r, gsl_matrix *Ap, gsl_matrix *A, gsl_matrix *H, gsl_vector *mass); 

static void show_matrix(gsl_matrix *A, int Num, int N2);
static void show_vector(gsl_vector *A, int Num);
