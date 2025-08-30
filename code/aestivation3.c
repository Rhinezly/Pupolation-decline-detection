#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <time.h>
#include <omp.h>

/////////////////////////////////////////////////
// CALCULATE SUMMARY STATISTICS
/////////////////////////////////////////////////
// ALLELE FREQ
SEXP allele_freq(SEXP input, SEXP ncpu)
{
	int nrow=INTEGER(GET_DIM(input))[0];
	int loci=INTEGER(GET_DIM(input))[1];
	SEXP output=PROTECT(allocVector(REALSXP, loci));
	Rbyte *h_input=RAW(input); double *h_output=REAL(output);
	// OMP STUFF
	omp_set_dynamic(0);
	omp_set_num_threads(asInteger(ncpu));
	#pragma omp parallel for
	for (int i=0; i<loci; i++)
	{
		int temp=0;
		unsigned long int offset=nrow*i;
		for (int j=0; j<nrow; j++)
		{
			temp=temp+(int) h_input[offset+j];
		}
		h_output[i]=(double) temp/nrow;
	}
	UNPROTECT(1);
	return output;
}

// LD r (NOT r2, CALCULATE r^2 IN R SIDE), PHASED
SEXP cal_r(SEXP input, SEXP ncpu)
{
	int nrow=INTEGER(GET_DIM(input))[0];
	int pair=INTEGER(GET_DIM(input))[1]/2;
	SEXP output=PROTECT(allocVector(REALSXP, pair));
	Rbyte *h_input=RAW(input); double *h_output=REAL(output);
	// OMP STUFF
	omp_set_dynamic(0);
	omp_set_num_threads(asInteger(ncpu));
	#pragma omp parallel for
	for (int i=0; i<pair; i++)
	{
		double p[4]={0, 0, 0, 0};
		Rbyte temp=0;
		int offset=2*nrow*i;
		for (int j=0; j<nrow; j++)
		{
			temp=2*h_input[offset+j]+h_input[offset+nrow+j];
			p[temp]=p[temp]+1;
		}
		double D=p[0]*p[3]-p[1]*p[2];
		h_output[i]=(D/sqrt((p[0]+p[1])*(p[2]+p[3])*(p[0]+p[2])*(p[1]+p[3])));
	}
	UNPROTECT(1);
	return output;
}
/////////////////////////////////////////////////
// SIMULATORS
/////////////////////////////////////////////////
// SAMPLE PARENT FROM 0, 1, 2, ..., (N-1). UNIFORMLY. USING R RANDOM NUMBER GENERATOR
int sample_parent(int lower, int N)
{
	double temp=runif(lower, (double) N);
	return floor(temp);
}

// SAMPLE RECOM USING C rand(). RETURN 0 OR 1. NO NEED TO BE TOO PRECISE
// NEED TO SET srand() PER THREAD
int sample_recom(double c)
{
	double temp=(double) rand()/RAND_MAX;
	return temp<c;
}

void start_rand()
{
	rand();
}

// INITIALISATION. A BIT RANDOM START AT AROUND p0, RANDOM LINKAGE (SAMPLING NOISE)
SEXP initialise(SEXP p0, SEXP N, SEXP ncpu)
{
	// OMP STUFF
	omp_set_dynamic(0);
	omp_set_num_threads(asInteger(ncpu));
	// POINTERS
	double *h_p0=REAL(p0); int h_N=asInteger(N); 
	int loci=LENGTH(p0);
	SEXP output=PROTECT(allocMatrix(RAWSXP, 2*h_N, loci));
	Rbyte *h_output=RAW(output);
	// BECAUSE OF THOSE GOMP_barrier THING, IF I SPLIT pragma omp parallel AND pragma omp for INTO TWO LINES
	// I'LL NEED TO REWRITE THE WHOLE THING IN TWO VERSIONS...
	#pragma omp parallel
	{
		// THREAD-PRIVATE RANDOM SEED FOR RECOM
		srand(time(0)+omp_get_thread_num());
		start_rand();
		int offset=0; 
		#pragma omp for
		for (int i=0; i<loci; i++)
		{
			offset=2*h_N*i;
			for (int j=0; j<2*h_N; j++)
			{
				h_output[offset+j]=sample_recom(h_p0[i]);
			}
		}
	}
	UNPROTECT(1);
	return output;
}

// ONE-STEP SEASONAL
SEXP seasonal(SEXP parental, SEXP N, SEXP c, SEXP ncpu)
{
	GetRNGstate();
	// OMP STUFF
	omp_set_dynamic(0);
	omp_set_num_threads(asInteger(ncpu));
	// POPULATION SIZE AT PARENTAL GEN AND AT NEW GEN
	int N_parent=INTEGER(GET_DIM(parental))[0]/2;
	int h_N=asInteger(N); 
	// VECTOR OF RECOM RATES. ALSO DETERMINE HOW MANY PAIRS OF LOCI TO SIM
	double *h_c=REAL(c); int pair=LENGTH(c);
	// MATRIX FOR OFFPSRING AND THE CORRESPONDING POINTERS
	SEXP offspring=PROTECT(allocMatrix(RAWSXP, 2*h_N, 2*pair));
	Rbyte *h_offspring=RAW(offspring); Rbyte *h_parental=RAW(parental);
	// THE TWO PARENTS AS VECTORS, AND THIER POINTERS
	SEXP parent1=PROTECT(allocVector(INTSXP, h_N));
	SEXP parent2=PROTECT(allocVector(INTSXP, h_N));
	int *h_parent1=INTEGER(parent1); int *h_parent2=INTEGER(parent2);
	// SAMPLE PARENTS, SINGLE THREAD. ALL LOCI HAVE THE SAME PARENT
	for (int i=0; i<h_N; i++)
	{
		h_parent1[i]=sample_parent(0, N_parent);
		h_parent2[i]=sample_parent(0, N_parent);
	}
	// REPRODUCTION
	#pragma omp parallel
	{
		// THREAD-PRIVATE RANDOM SEED FOR RECOM
		srand(time(0)+omp_get_thread_num());
		// OTHER THREAD-PRIVATE VARIABLES. MOSTLY ABOUT RECOM
		int p1_l1=0; int p1_l2=0; 
		int p2_l1=0; int p2_l2=0;
		int offset_parental=0; int offset_offspring=0;
		#pragma omp for
		for (int i=0; i<pair; i++)
		{
			// WHICH COLUMNS TO LOOK AT?
			offset_parental=4*N_parent*i;
			offset_offspring=4*h_N*i;
			for (int j=0; j<h_N; j++)
			{
				// CHOOSE ONE CHROMOSOME/GAMETE PER PARENT
				p1_l1=sample_recom(0.5); p1_l2=p1_l1;
				p2_l1=sample_recom(0.5); p2_l2=p2_l1;
				// ANY RECOMBINATION?
				if (sample_recom(h_c[i])==1) {p1_l2=1-p1_l2;}
				if (sample_recom(h_c[i])==1) {p2_l2=1-p2_l2;}
				// PUT INTO CORRECT PLACES. ROW 2*j AND 2*j+1 FORM THE SAME INDIVIDUAL
				// FIRST CHROMOSOME. ROW 2*j
				h_offspring[offset_offspring+2*j]=h_parental[offset_parental+2*h_parent1[j]+p1_l1];
				h_offspring[offset_offspring+2*h_N+2*j]=h_parental[offset_parental+2*N_parent+2*h_parent1[j]+p1_l2];
				// SECOND CHROMOSOME. ROW 2*j+1
				h_offspring[offset_offspring+2*j+1]=h_parental[offset_parental+2*h_parent2[j]+p2_l1];
				h_offspring[offset_offspring+2*h_N+2*j+1]=h_parental[offset_parental+2*N_parent+2*h_parent2[j]+p2_l2];
			}
		}
	}
	PutRNGstate();
	UNPROTECT(3);
	return offspring;
}

// AES AND REP COMPARTMENT SIZES SHIFT PER GENERATION
// WHOLE THING NOT RUN IF AESTIVATING SIZE REMAINS THE SAME
SEXP aestivation_shift(SEXP current_rep, SEXP current_aes, SEXP new_aes_size, SEXP ncpu)
{
	GetRNGstate();
	// OMP STUFF
	omp_set_dynamic(0);
	omp_set_num_threads(asInteger(ncpu));
	// DIM SIZES
	int pair=INTEGER(GET_DIM(current_rep))[1]/2;
	int current_rep_size=INTEGER(GET_DIM(current_rep))[0]/2;
	int current_aes_size=INTEGER(GET_DIM(current_aes))[0]/2;
	int new_rep_size=current_rep_size+current_aes_size-asInteger(new_aes_size);
	// THE TWO MATRICES FOR THE NEW AESTIVATION AND REPRODUCING COMPARTMENTS. AND THE POINTERS. 
	SEXP new_rep=PROTECT(allocMatrix(RAWSXP, 2*new_rep_size, 2*pair));
	SEXP new_aes=PROTECT(allocMatrix(RAWSXP, 2*asInteger(new_aes_size), 2*pair));
	SEXP output=PROTECT(allocVector(VECSXP, 2));
	Rbyte *h_new_rep=RAW(new_rep); Rbyte *h_new_aes=RAW(new_aes);
	Rbyte *h_current_rep=RAW(current_rep); Rbyte *h_current_aes=RAW(current_aes);
	// IF THERE ARE MORE AES INDIVIDUALS
	if (asInteger(new_aes_size)>current_aes_size)
	{
		// FOR EACH LOCUS
		#pragma omp parallel for
		for (int i=0; i<2*pair; i++)
		{
			int offset_current_rep=2*current_rep_size*i;
			int offset_new_rep=2*new_rep_size*i;
			int offset_current_aes=2*current_aes_size*i;
			int offset_new_aes=2*asInteger(new_aes_size)*i;
			// COPY SOME CURRENT REP INDIVIDUALS
			#pragma omp simd
			for (int j=0; j<2*new_rep_size; j++)
			{
				h_new_rep[offset_new_rep+j]=h_current_rep[offset_current_rep+j];
			}
			// COPY ALL CURRENT AES INDIVIDUALS
			#pragma omp simd
			for (int j=0; j<2*current_aes_size; j++)
			{
				h_new_aes[offset_new_aes+j]=h_current_aes[offset_current_aes+j];
			}
			// MOVE REMAINING REP TO THE BOTTOM OF AES
			#pragma omp simd
			for (int j=2*new_rep_size; j<2*current_rep_size; j++)
			{
				h_new_aes[offset_new_aes+2*(current_aes_size-new_rep_size)+j]=h_current_rep[offset_current_rep+j];
			}
		}
	}
	// IF THERE ARE MORE REP INDIVIDUALS
	if (current_aes_size>asInteger(new_aes_size))
	{
		// FOR EACH LOCUS
		#pragma omp parallel for
		for (int i=0; i<2*pair; i++)
		{
			int offset_current_rep=2*current_rep_size*i;
			int offset_new_rep=2*new_rep_size*i;
			int offset_current_aes=2*current_aes_size*i;
			int offset_new_aes=2*asInteger(new_aes_size)*i;
			// COPY ALL CURRENT REP INDIVIDUALS
			#pragma omp simd
			for (int j=0; j<2*current_rep_size; j++)
			{
				h_new_rep[offset_new_rep+j]=h_current_rep[offset_current_rep+j];
			}
			// COPY SOME AES INDIVIDUALS
			#pragma omp simd
			for (int j=0; j<2*asInteger(new_aes_size); j++)
			{
				h_new_aes[offset_new_aes+j]=h_current_aes[offset_current_aes+j];
			}
			// COPY REMAINING AES TO THE BOTTON OF REP
			#pragma omp simd
			for (int j=2*asInteger(new_aes_size); j<2*current_aes_size; j++)
			{
				h_new_rep[offset_new_rep+2*(current_rep_size-asInteger(new_aes_size))+j]=h_current_aes[offset_current_aes+j];
			}
		}
	}
	SET_VECTOR_ELT(output, 0, new_rep);
	SET_VECTOR_ELT(output, 1, new_aes);
	PutRNGstate();
	UNPROTECT(3);
	return output;
}

