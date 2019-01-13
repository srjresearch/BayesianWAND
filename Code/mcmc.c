//header files
#include <stdio.h> 
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


////////////////////////// MCMC scheme for posterior sampling using the Bayesian WAND model \\\\\\\\\\\\\\\\\\\\\\\\\\\\ 


///////////// Global variables /////////////

//data details -- note additional setup is needed within the main function (line 55) to define complete/partial/top rankings
int n = 34; //number of rankings
int K = 30; //maximum number of entities
char datafilename[] = "../Data/nba_data.ssv"; //data file


//mcmc details -- note the total number of iterations performed is burn_in + (nsample*thin)
int burn_in = 5000; //number of burn-in iterations
int nsample = 10000; //number of desired posterior realisations
int thin = 1; //thin factor between iterations




///////////// Functions /////////////

void print_array_int(FILE* file_name, int* array, int size); //print array of ints to file
void print_array_dbl(FILE* file_name, double* array, int size); //print array of doubles to file
int** make_int_mat(int n, int m); //creates a 2D array of integers
void zero_int_mat(int** mat, int dim_1, int dim_2); //zeros an interger matrix of dimension dim_1 by dim_2
double** make_double_mat(int n, int m); //creates a 2D array of doubles
int sample_bern(double U, double prob); //sample from bernoulli distribution
int sample_cum_probs(double U, double* probs, int size); //returns a sample when a cumulative probability array is passed
int vecmax_int(int* b, int size); //returns max value in an array of ints
double vecmax_dbl(double* b, int size); //returns max value in an array of ints
void swap_labels_arr( int* arr, int label_1, int label_2, int* m_vec, int size); //swaps two (cluster) labels in an array of length n, requires the number of times each label is used (contained in m_vec)
void ranker_label_swap(int label_1, int label_2, int* cvec, double** lambda, int* m_vec, int** d,  double* gam); //this function swaps all the required parameters when a ranker label swaps
double loglikeli_oneranking(int** x, double** lambda_dagger, double** z, int* cvec, int** d, int* ni_vec, int ranki, int* K_vec, int* r_vec); //calculates log complete data likelihood of a single ranking
double loglikeli_oneranking_onevec(int** x, double* lambda_dagger, double** z, int* cvec, int** d, int* ni_vec, int ranki, int* K_vec, int* r_vec); //same as above but accepts parameter vector lambda not 2D array
int sample_probs(double U, double* prob, int size); //samples an integer (0,...,size-1) with probaities given by prob



int main()
{
	///////////// Additional data set up to specify complete/partial/top rankings /////////////
	
	int* K_vec = (int*)malloc(n*sizeof(int)); //K_vec[i] holds the number of entities considered by ranker i (K_i)
	for(int i = 0; i < n; i++) 
	{
		K_vec[i] = K; //each ranker considered all entities
	}
	
	int* ni_vec = (int*)malloc(n*sizeof(int)); //ni_vec[i] holds the number of ranks reported by ranker i (n_i)
	for(int i = 0; i < 6; i++) //professional (first 6) rankers report complete rankings
	{
		ni_vec[i] = 30;
	}
	
	for(int i = 6; i < n; i++) //remaining rankers report top-8 rankings
	{
		ni_vec[i] = 8;
	}
	
		
    ///////////// Priors /////////////
    
    double a_lambda = 1; //lambda ~ Ga(a_lambda,b_lambda)
    double b_lambda = 1;
    
	double a_alpha = 3; //alpha ~ Ga(a_alpha,b_alpha)
	double b_alpha = 3;
	
	double a_gamma = 3; //gamma_s ~indep Ga(a_gamma,b_gamma)
	double b_gamma = 3;
    

	double* p_vec = (double*)malloc(n * sizeof(double)); //w_i ~indep Bern(p_i)
	for(int i = 0; i < 6; i++) //professionals
	{
		p_vec[i] = 0.9;
	}
	for(int i = 6; i < 12; i++) //avid fans
	{
		p_vec[i] = 0.7;
	}
	for(int i = 12; i < 18; i++) //fans
	{
		p_vec[i] = 0.5;
	}
	for(int i = 18; i < 25; i++) //infrequent viewers
	{
		p_vec[i] = 0.3;
	}
	for(int i = 25; i < 34; i++) //not interested
	{
		p_vec[i] = 0.1;
	}
	
	
	///////////// Algorithm details /////////////
	
	int m_ranker = 3; //number of auxiliary ranker components in cluster indicator sampling (m^r)
	int m_ent = 3; //number of auxiliary entity components in cluster indicator sampling (m^e)
	
	int weighted_PL = 1; //indicator variable. Set = 1 for Weighted Plackett--Luce model. Set = 0 for standard Plackett--Luce model (w_i = 1).
	
	int fixed_seed = 0; //indicator variable. Set = 1 for a fixed seed. Set = 0 for a random seed.
	
	
	
	///////////// Initalisation /////////////
	
	gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937); //initalise RNG
	if(!fixed_seed)
	{
		gsl_rng_set(r,time(NULL)); //random seed
	}
	double ran_U; //hold random uniform variables
	
	double alpha = gsl_ran_gamma(r,a_alpha,1.0/b_alpha); //alpha - ranker concentration parameter
	
	
	int* cvec = (int*)calloc(n,  sizeof(int)); //holds ranker cluster allocations
	int* m_vec = (int*)calloc(n, sizeof(int)); //keeps track of how many rankers are in each cluster
	double* prior_rank_probs = (double*)malloc(n * sizeof(double)); //holds probailities that rankers are in each cluster
	cvec[0] = 0; //first ranker in cluster 0
	m_vec[0] ++; //1 ranker in cluster 0
	int no_rank_clust = 1; //1 ranker cluster in total
	for(int i = 1; i < n; i++) //sample remaining ranker cluster lables
    {
		for(int q = 0; q < no_rank_clust; q++) //compute (normalised) probabilities
		{
			prior_rank_probs[q] = (double) m_vec[q]/(i + alpha);
		}
		prior_rank_probs[no_rank_clust] = alpha/(i + alpha);
		
		ran_U = gsl_rng_uniform(r); //random U(0,1)
		cvec[i] = sample_probs(ran_U, prior_rank_probs, no_rank_clust + 1); //sample indictor for ranker i
		
		m_vec[cvec[i]] ++; //extra ranker in this cluster
		if(cvec[i] > (no_rank_clust - 1)) //new cluster
					no_rank_clust ++;
	}
	free(prior_rank_probs);
	
	
	double* gam = (double*)malloc(n * sizeof(double)); //gamma[i] holds the entity concentration parameter for ranker cluster i
	for(int i = 0; i < no_rank_clust; i++)
	{
		gam[i] = gsl_ran_gamma(r, a_gamma,1.0/b_gamma);
	}
	
			
	int** d = make_int_mat(n, K); //d[i][j] holds entity cluster indicator for entity j in ranker cluster i
	int* no_ent_clust = (int*)malloc(n * sizeof(int)); // no_ent_clust[i] holds number of entity clusters within ranker cluster i
	double* prior_ent_probs = (double*)malloc(K * sizeof(double)); //holds probailities that entities are in each cluster
	
	for(int i = 0; i < no_rank_clust; i++)
	{
		int* prior_ent_m_vec = (int*)calloc(K, sizeof(int)); //keeps track of how many entities are in each cluster
		d[i][0] = 0; //first entity within cluster 0
		prior_ent_m_vec[0] ++; //1 entity in cluster 0
		no_ent_clust[0] = 1; //1 entity cluster in ranker cluster 0
			
			for(int l = 1; l < K; l++) //sample remaining (entity) cluster labels
			{	
				for(int q = 0; q < no_ent_clust[i]; q++)
				{
					prior_ent_probs[q] = (double) prior_ent_m_vec[q]/(l + gam[i]);
				}
				prior_ent_probs[no_ent_clust[i]] = gam[i]/(l + gam[i]);

				//probabilities don't need normalising by construction
				ran_U = gsl_rng_uniform(r);
				d[i][l] = sample_probs(ran_U, prior_ent_probs, no_ent_clust[i] + 1);
				
				prior_ent_m_vec[d[i][l]] ++; //add 1 to cluster
				
				if(d[i][l] > (no_ent_clust[i] - 1)) //new cluster
					no_ent_clust[i] ++;
			}
		free(prior_ent_m_vec);
	}
	free(prior_ent_probs);
				
	
	double** lambda = make_double_mat(n,K); //holds unique lambda values. lambda[i][j] is the lambda for entity cluster j in ranker cluster i
	for(int i = 0; i < no_rank_clust; i++)
	{
		for(int j = 0; j < no_ent_clust[i]; j++)
		{
			lambda[i][j] = gsl_ran_gamma(r,a_lambda,1.0/b_lambda); //sample lambdas from base
		}
	}
	
	int mix_over_r = weighted_PL; 
	int* r_vec = (int*)malloc(n * sizeof(int)); //r_vec[i] holds binary ranker weight for ranker i (w_i in paper).
	if(mix_over_r) //Weighted Plackett--Luce model
	{
		for(int i = 0; i < n; i++) //initialise r_i from the prior probabilities p_vec
		{
			ran_U = gsl_rng_uniform(r);
			r_vec[i] = sample_bern(ran_U, p_vec[i]);
		}
	}
	else //standard Plackett--Luce model
	{
		for(int i = 0; i < n; i++) //fix equal to one
		{
			r_vec[i] = 1;
		}
	}
	
	
	
	///////////// Read in  data /////////////
	
	int** x = make_int_mat(n,K); //holds data. x[i][j] is ranking i position j.
	FILE *datafile;
	gsl_matrix *data;
	data = gsl_matrix_alloc(n,K);
	datafile=fopen(datafilename,"r");
	gsl_matrix_fscanf(datafile,data);
	fclose(datafile); //close file
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < K; j++)
		{
			x[i][j] = gsl_matrix_get(data,i,j);
		}
	}
	

	///////////// Open output files /////////////
	
	//check if output directories exist, if not create them
	struct stat st = {0};
	if(stat("outputs", &st) == -1)
	{
		mkdir("outputs", 0777);
		printf("Created new output directory as /outputs did not exist\n");
	}
	if(stat("outputs/rankers_lambda", &st) == -1)
	{
		mkdir("outputs/rankers_lambda", 0777);
	}
	
	FILE** lambda_files = malloc(n * sizeof(FILE*)); //output each rankers lambda values
	for(int i = 0; i < n ; i++)
	{
		char myfile[40];
		sprintf(myfile, "outputs/rankers_lambda/r%d_lambda.ssv", i+1);
		lambda_files[i] = fopen(myfile, "w");
	}
	FILE* cvec_out = fopen("outputs/cvec.ssv","w");
	FILE* likeli_out = fopen("outputs/log_comlpete_data_likeli.ssv","w");
	FILE* alpha_out = fopen("outputs/alpha.ssv","w");
	FILE* gamma_out = fopen("outputs/gamma.ssv","w");
	FILE* p_out = fopen("outputs/ranker_weight_probs.ssv","w");
	FILE* rankers_p_out = fopen("outputs/ranker_weights.ssv", "w");
		

	///////////// Variables/arrays used in MCMC /////////////		
			
	int rank_q_minus; //q^- used in sampling ranker allocations
	double** lambda_ranker_aux = make_double_mat(m_ranker, K); //will hold the proposed lambdas for the auxiliary ranker clusters 
	int** d_aux = make_int_mat(m_ranker, K); //holds entity clustering for each auxiliary ranker cluster
	double* gam_aux = (double*)malloc(m_ranker * sizeof(double)); //holds the value of gamma for each auxiliary ranker cluster
	int* aux_ent_m_vec = (int*)malloc(K * sizeof(int)); //keeps track of how many entities are in each cluster within an auxiliary ranker cluster
	double* aux_ent_probs = (double*)malloc(K * sizeof(double)); //holds probabilities for auxiliary clusters
	int ent_q_minus; //q^- used in sampling entity allocations
	int* ent_m_vec = (int*)calloc(K, sizeof(int)); //ent_m_vec[i] holds the number of entities in entity cluster i. Computed  for each ranker cluster 
	int** beta_mat = make_int_mat(n, K); //used in sampling lambda
	double log_likeli_r1, log_likeli_r0; //used in sampling w_i 
	double* ptilde = (double*)malloc(n * sizeof(double)); //ptilde[i] holds probabilty that w_i = 1 given all other quantites, and the data
	
	//introduce and draw latent variables z_ij from  their FCD
	double** z = make_double_mat(n,K); //holds latent variables
	for(int i = 0; i < n; i++)
	{
		if(r_vec[i] == 1)
		{
			double sum = 0;
			for(int q = 0; q < K_vec[i]; q++)
			{
				sum += lambda[cvec[i]][d[cvec[i]][x[i][q]-1]];
			}
			
			for(int j = 0; j < ni_vec[i]; j++)
			{
				z[i][j] = gsl_ran_exponential(r, 1.0/sum);
				sum -= lambda[cvec[i]][d[cvec[i]][x[i][j]-1]];
			}
		}else
		{
			for(int j = 0; j < ni_vec[i]; j++)
			{
				z[i][j] = gsl_ran_exponential(r, 1.0/(K_vec[i]-j));
			}		
		}	
		
	}
   
	//print curret time and total number of iterations to be performed
	time_t rawtime;
	struct tm * timeinfo;
	printf("Starting MCMC for %d iterations\n",  burn_in + (nsample*thin));
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	printf("Current local time and date: %s", asctime(timeinfo));
	
	
///////////// Start MCMC /////////////		
for(int itercount = 0; itercount < ((burn_in + (nsample*thin))+1); itercount++)
{
	
	memset(m_vec, 0, n * sizeof(int)); //zero m_vec
	for(int i = 0; i < n; i++) //compute the number of rankers in each cluster
	{
		m_vec[cvec[i]] ++;
	}
	
	///////////// Sample ranker cluster allocations /////////////
		
	for(int i = 0; i < n; i++)
	{
		int aux_indicator = 0; //we will sample aux_indicator -> m_ranker aux clusters	
		m_vec[cvec[i]] -= 1; //remove 1 from the number in the current cluster. Note: if it was a singleton then m_vec[cvec[i]] == 0. This is n^r_{-i,c}
		
		if(m_vec[cvec[i]] == 0) //singleton
		{
			//relabel this cluster as the fisrt auxillary one.
			for(int j = 0; j < K; j++)
			{
				lambda_ranker_aux[0][j] = lambda[cvec[i]][j];
				d_aux[0][j] = d[cvec[i]][j];
			}		
			gam_aux[0] = gam[cvec[i]];
			
			no_rank_clust = no_rank_clust - 1; //remove one from the total number of ranker clusters
			aux_indicator = 1; //now require m_ranker - 1 aux clusters (first one populated by ranker i's current cluster)
			cvec[i] = -99; //make ranker i's cluster label < 0 so that it doesn't affect the maximum value in cvec
					
		}
		
		rank_q_minus = no_rank_clust; //current number of ranker clusters. This will be the same as the previous iteration unless we are dealing with a singleton
						
		for(int j = aux_indicator; j < m_ranker; j++) //sample auxillary ranker clusters. That is, draw from dirichlet processes via chinese resturant process
		{
			gam_aux[j] = gsl_ran_gamma(r, a_gamma, 1.0/b_gamma); //sample concentration parameter from prior
			
			for(int l = 0; l < K; l++) //sample lambdas from prior
			{
				lambda_ranker_aux[j][l] = gsl_ran_gamma(r, a_lambda, 1.0/b_lambda); 
			}
					
			memset(aux_ent_m_vec, 0, K * sizeof(int));
			
			d_aux[j][0] = 0; //place first entity within cluster 0
			aux_ent_m_vec[0] ++;
			int aux_ent_count = 1; //number of ent clusters
			
			for(int l = 1; l < K; l++) //sample remaining (entity) cluster labels
			{	
				for(int q = 0; q < aux_ent_count; q++)
				{
					aux_ent_probs[q] = (double) aux_ent_m_vec[q]/(l + gam_aux[j]);
				}
				aux_ent_probs[aux_ent_count] = gam_aux[j]/(l + gam_aux[j]);

				//probabilities don't need normalising by construction
				ran_U = gsl_rng_uniform(r);
				d_aux[j][l] = sample_probs(ran_U, aux_ent_probs, aux_ent_count + 1);
				
				aux_ent_m_vec[d_aux[j][l]] ++; //add 1 to cluster
				
				if(d_aux[j][l] > (aux_ent_count - 1)) //new cluster
					aux_ent_count ++;
			}
		}
		
		
		//relabel ranker clusters from 0 to rank_q_minus-1
		int max_ranker = vecmax_int(cvec,n);
		while(max_ranker != (rank_q_minus-1)) //need to relable largest cluster indicator 
		{
			for(int j = 0; j < n; j++) //find empty cluster
			{
				if(m_vec[j] == 0) //cluster j is empty
				{
					ranker_label_swap(j, max_ranker, cvec, lambda,m_vec,d,gam); //swap labels of ranker cluster j and max_ranker
					break;
				}
			}	
			max_ranker = vecmax_int(cvec,n);
		}
		
		int h_ranker = rank_q_minus + m_ranker; //number of possible cluster allocations for ranker i
		
		//compute cluster probabilities
		double* ranker_probs = (double*)malloc(h_ranker * sizeof(double));
		for(int j = 0; j < rank_q_minus; j++) //prob of joining an existing cluster
		{
			cvec[i] = j;
			ranker_probs[j] =  loglikeli_oneranking(x,lambda,z,cvec,d,ni_vec,i,K_vec,r_vec);
			
		}
		for(int j = 0; j < m_ranker; j++) //prob of joining aux cluster
		{
			cvec[i] = j;
			ranker_probs[rank_q_minus + j] = loglikeli_oneranking(x,lambda_ranker_aux,z,cvec,d_aux,ni_vec,i,K_vec,r_vec);
		}
		
	
		//normalise the probailities 		
		double max_ranker_prob = vecmax_dbl(ranker_probs, h_ranker);
		double sum_ranker_probs = 0;
		for(int s = 0; s < rank_q_minus; s++)
		{
			ranker_probs[s] = m_vec[s] * exp(ranker_probs[s] - max_ranker_prob);
			sum_ranker_probs += ranker_probs[s];
		}
		for(int s = rank_q_minus; s < h_ranker; s++)
		{
			ranker_probs[s] = (alpha/m_ranker) * exp(ranker_probs[s] - max_ranker_prob);
			sum_ranker_probs += ranker_probs[s];
		}
		for(int s = 0; s < h_ranker; s++)
		{
			ranker_probs[s] = ranker_probs[s]/sum_ranker_probs;
		}

		ran_U = gsl_rng_uniform(r); //sample U(0,1)
		cvec[i] = sample_probs(ran_U, ranker_probs, h_ranker); //sample a value for cvec[i]
		
		free(ranker_probs);	//free ranker_probs
		
		if(cvec[i] > (rank_q_minus - 1)) //we've sampled an auillary cluster. Need to introduce parameters into lambda/d arrays from the lambda_aux and d_aux arrays
		{
			for(int l = 0; l < K; l++) //move lambda and d
			{
				lambda[rank_q_minus][l] = lambda_ranker_aux[cvec[i] - rank_q_minus][l]; //rank_q_minus is the lowest avalable cluster label as current clusters are labelled 0 to rank_q_minus-1
				d[rank_q_minus][l] = d_aux[cvec[i] - rank_q_minus][l];
			}
			
			gam[rank_q_minus] = gam_aux[cvec[i] - rank_q_minus]; //move gamma
			
			cvec[i] = rank_q_minus; //set cluster label
			no_rank_clust ++; //increase number of ranker clusters by 1
		}
		
		m_vec[cvec[i]] ++; //update m_vec after sampling cluster label
			
	}
	
	
	///////////// Sample alpha /////////////	
	
	double neta = gsl_ran_beta(r, alpha + 1.0, n);
	double frac = ((double) a_alpha + no_rank_clust - 1)/(n * (b_alpha - log(neta)));
	double pi = frac/(1.0 + frac);
	ran_U = gsl_rng_uniform(r);
	if(ran_U < pi) //component 1
	{
		alpha = gsl_ran_gamma(r, a_alpha + no_rank_clust, 1.0/(b_alpha - log(neta)));
	}else  //component 2
	{
		alpha = gsl_ran_gamma(r, a_alpha + no_rank_clust - 1.0, 1.0/(b_alpha - log(neta)));
	}
			

	///////////// Sample entity cluster allocations (within each ranker cluster) /////////////	
	
	//it is usefull to know which rankers are in which cluster
	int* tmp_counts = (int*)calloc(no_rank_clust, sizeof(int));
	int** ranker_allocations = make_int_mat(no_rank_clust, vecmax_int(m_vec, no_rank_clust));
	for(int i = 0; i < n; i++)
	{
		ranker_allocations[cvec[i]][tmp_counts[cvec[i]]] = i;
		tmp_counts[cvec[i]] ++;
	}
	free(tmp_counts);

	for(int s = 0; s < no_rank_clust; s++) //loop over ranker clusters
	{

		//calculate number of entites in each entity cluster for ranker cluster s
		memset(ent_m_vec, 0, K * sizeof(int));
		for(int i = 0; i < K; i++)
		{
			ent_m_vec[d[s][i]] ++;
		}
		
		//calculate number of entity clusters within ranker cluster s
		no_ent_clust[s] = 0; 
		for(int i = 0; i < K; i++) //we can loop over K here as the maximum cluster allocation number is K-1
		{
			if(ent_m_vec[i] > 0)
			{
				no_ent_clust[s] ++;
			}
		}	
		
		//relabel the entity cluster indicators from 0 to no_ent_clust[s]
		int max_ent = vecmax_int(d[s], K); 
		while(max_ent != (no_ent_clust[s]-1))
		{
			for(int l = 0; l < K; l++)
			{
				if(ent_m_vec[l] == 0) //empty cluster
				{
					//swap labels of ent cluster l and max_ent
					swap_labels_arr(d[s] , l , max_ent , ent_m_vec , K);
					
					double tmp_lambda = lambda[s][l]; // swap the lambdas
					lambda[s][l] = lambda[s][max_ent];
					lambda[s][max_ent] = tmp_lambda;
					
					ent_m_vec[l] = ent_m_vec[max_ent]; //swap counts round
					ent_m_vec[max_ent] = 0;
					break;
				}
			}	
			max_ent = vecmax_int(d[s],K);
		}
		
		//sample the cluster labels d[s][i], i = 0,...,K		
		double tmp_single_lambda;
		for(int i = 0; i < K; i++)
		{
			int ent_aux_indicator = 0;
			
			ent_m_vec[d[s][i]] -= 1; //remove 1 from the number in the cluster. Note: if it was a singleton the ent_m_vec[d[s][i]] == 0.
		
			if(ent_m_vec[d[s][i]] == 0) //singleton
			{
				//relabel this cluster as the fisrt auxillary one
				tmp_single_lambda = lambda[s][d[s][i]];
				
				no_ent_clust[s] = no_ent_clust[s] - 1; //remove one from the total number of clusters (its now an auxillary cluster)
				ent_aux_indicator = 1; //now require m_ranker - 1 aux clusters
				d[s][i] = -99; //make this < 0 so that it doesn't affect the maxumum of d[s]
			}
			
			ent_q_minus = no_ent_clust[s];			
			
			//relabel the entity cluster indicators from 0 to no_ent_clust[s]
			int max_ent = vecmax_int(d[s], K); 
			while(max_ent != (ent_q_minus-1))
			{
				for(int l = 0; l < K; l++)
				{
					if(ent_m_vec[l] == 0) //empty cluster
					{
						//swap labels of ent cluster l and max_ent
						swap_labels_arr(d[s] , l , max_ent , ent_m_vec , K);
						
						double tmp_lambda = lambda[s][l]; // swap the lambdas too
						lambda[s][l] = lambda[s][max_ent];
						lambda[s][max_ent] = tmp_lambda;
						
						ent_m_vec[l] = ent_m_vec[max_ent]; //swap counts round
						ent_m_vec[max_ent] = 0; //now an empty cluster
						break;
					}
				}	
				max_ent = vecmax_int(d[s],K);
			}			
			
			
			int h_ent = ent_q_minus + m_ent;
			
			double* lambda_ent_aux = (double*)malloc(h_ent * sizeof(double)); //this will hold the current cluster lambdas, along with the auxillary ones
			for(int j = 0; j < ent_q_minus; j++)
			{
				lambda_ent_aux[j] = lambda[s][j];
			}
			if(ent_aux_indicator == 1) //if we had a singleton we load the value back as the first aux cluster
			{
				lambda_ent_aux[ent_q_minus] = tmp_single_lambda; 
			}
			for(int j = ent_aux_indicator; j < m_ent; j++) //sample aux values from prior
			{
				lambda_ent_aux[ent_q_minus + j] = gsl_ran_gamma(r , a_lambda, b_lambda); 
			}
			
			//compute cluster indicator probabilities
			double* ent_probs = (double*)calloc(h_ent , sizeof(double));
			for(int j = 0; j < h_ent; j++) 
			{
				d[s][i] = j;
				for(int l = 0; l < m_vec[s]; l++) //evalute all rankings in ranker cluster S
				{
					ent_probs[j] += loglikeli_oneranking_onevec(x, lambda_ent_aux, z, cvec, d, ni_vec, ranker_allocations[s][l], K_vec, r_vec);
				}
				
			}

			
			////suitably scale and normalise the probailities 		
			double max_ent_prob = vecmax_dbl(ent_probs, h_ent);
			double sum_ent_probs = 0;
			for(int j = 0; j < ent_q_minus; j++)
			{
				ent_probs[j] = ent_m_vec[j] * exp(ent_probs[j] - max_ent_prob);
				sum_ent_probs += ent_probs[j];
			}
			for(int j = ent_q_minus; j < h_ent; j++)
			{
				ent_probs[j] = (gam[s]/m_ent) * exp(ent_probs[j] - max_ent_prob);
				sum_ent_probs += ent_probs[j];
			}
			for(int j = 0; j < h_ent; j++)
			{
				ent_probs[j] = ent_probs[j]/sum_ent_probs;
			}
			
			ran_U = gsl_rng_uniform(r); 
			d[s][i] = sample_probs(ran_U, ent_probs, h_ent); //sample a value for cvec[s][i]
			
			free(ent_probs); //free ent_probs
			
			if(d[s][i] > (ent_q_minus - 1)) //we've sampled an auillary cluster. Need to introduce parameters into normal arrays from the aux arrays
			{
				lambda[s][ent_q_minus] = lambda_ent_aux[d[s][i]]; //this is the lowest avaiable slot. As current clusters are labbeled 0....rank_q_minus-1	
				d[s][i] = ent_q_minus; 
				no_ent_clust[s] ++;
			}
			
			free(lambda_ent_aux);
			ent_m_vec[d[s][i]] ++; //keep track of number in each entity cluster
		}
	}
	
	free(ranker_allocations); //release ranker_allocations array
	
	///////////// Sample lambda /////////////	
	
	//compute beta_mat. beta[s][t] is the number of times lambda[s][t] corresponds to an entity in an informative ranking 
	zero_int_mat(beta_mat, n, K); //zero beta mat
	for(int i = 0; i < n; i++) 
	{
		if(r_vec[i] == 1)
		{
			for(int j = 0; j < ni_vec[i]; j++)
			{
				beta_mat[cvec[i]][d[cvec[i]][x[i][j]-1]] ++;
			}
		}
	}
	
	for(int s = 0; s < no_rank_clust; s++) //loop over ranker clusters
	{
		for(int t = 0; t < no_ent_clust[s]; t++) //loop over entity clusters within each ranker cluster
		{
			double sum_lambda = 0;
			for(int i = 0; i < n; i++) //compute sum_lambda for FCD
			{
				if(r_vec[i] == 1)
				{
					if(cvec[i] == s)
					{
						int count_lambda = 0;
						for(int m = 0; m < K_vec[i]; m++)
						{
							if(d[cvec[i]][x[i][m]-1] == t)
							{
								count_lambda ++;
							}
						}
						
						for(int j = 0; j < ni_vec[i]; j++)
						{
							sum_lambda += count_lambda*z[i][j];
							if(d[cvec[i]][x[i][j]-1] == t)
							{
								count_lambda -= 1;
							}
							
						}
					}
				}
			}

			lambda[s][t] = gsl_ran_gamma(r, (double) 1.0 + beta_mat[s][t], (double) 1.0/(1.0 + sum_lambda)); //draw from FCD
		}
	}
	
	
	
	///////////// Rescale lambda values /////////////	

	double rho = 0;
	int no_unique = 0;
    for(int s = 0; s < no_rank_clust; s++)
    {
        for(int t = 0; t < no_ent_clust[s]; t++)
        {
			rho += lambda[s][t]; //sum all unique lambda values
			no_unique ++;
        }
    }
 
    double cap_lambda = gsl_ran_gamma(r, (double)no_unique, 1.0); //draw prior realisation of sum of unqiue lambdas
    double ratio = cap_lambda/rho; //scaling ratio
    for(int s = 0; s < no_rank_clust; s++) //loop over all unique lambdas
    {
        for(int t = 0; t < no_ent_clust[s]; t++)
        {
			lambda[s][t] *= ratio; //rescale
        }
    }


	///////////// Sample latent variables Z /////////////	
	
	for (int i = 0; i < n; i++)
	{
		if(r_vec[i] == 1)
		{
			double sum = 0;
			for (int q = 0; q < K_vec[i]; q++) 
			{
				sum += lambda[cvec[i]][d[cvec[i]][x[i][q]-1]];
			}
				
			for(int j = 0; j < ni_vec[i]; j++)
			{
				z[i][j] = gsl_ran_exponential(r, 1.0/sum);
				sum -= lambda[cvec[i]][d[cvec[i]][x[i][j]-1]];
			}
		}else
		{
			for(int j = 0; j < ni_vec[i]; j++)
			{
				z[i][j] = gsl_ran_exponential(r, 1.0/(K_vec[i]-j));
			}
		}		
	}

	
	///////////// Sample ranker weights (w_i in paper) /////////////	
	
	if(mix_over_r) //if using Weighted Plackett--Luce model (if not then no need to sample)
	{
		for(int i = 0; i < n; i++) //loop over rankers
		{
			r_vec[i] = 1;
			log_likeli_r1 = loglikeli_oneranking(x, lambda, z, cvec, d, ni_vec, i, K_vec,r_vec); //log complete data likelihood with w_i = 1
			
			r_vec[i] = 0;
			log_likeli_r0 = loglikeli_oneranking(x, lambda, z, cvec, d, ni_vec, i, K_vec,r_vec); //log complete data likelihood with w_i = 0
			
			ptilde[i] = 1.0/(1.0 + (((1.0-p_vec[i])/p_vec[i])*exp(log_likeli_r0 - log_likeli_r1)));	//probability w_i = 1, given all other quantities and the data
			
			ran_U = gsl_rng_uniform(r); //random U(0,1)
			r_vec[i] = sample_bern(ran_U, ptilde[i]); //sample ranker weight
		}
	}
		
	
	///////////// Sample gamma[s] for each ranker cluster /////////////	
		
	for(int s = 0; s < no_rank_clust; s++) //loop over ranker clusters
	{
		double neta = gsl_ran_beta(r, gam[s] + 1, K);
		double frac = ((double) a_gamma + no_ent_clust[s] - 1)/(K * (b_gamma - log(neta)));
		double pi = frac/(1.0 + frac);
		
		ran_U = gsl_rng_uniform(r); //random U(0,1)
		if(ran_U < pi) //component 1
		{
			gam[s] = gsl_ran_gamma(r, a_gamma + no_ent_clust[s], 1.0/(b_gamma - log(neta)));
		}else //component 2
		{
			gam[s] = gsl_ran_gamma(r, a_gamma + no_ent_clust[s] - 1.0, 1.0/(b_gamma - log(neta)));
		}
		
	}
		
	
	///////////// Print statements /////////////	
	
	if(itercount > burn_in && itercount % thin == 0) //print current realisations of unknown quantities
	{
		//compute log complete data likelihood
		double full_loglikeli_out = 0; 
		for(int i = 0; i < n; i++)
		{	
			full_loglikeli_out += loglikeli_oneranking(x,lambda,z,cvec,d,ni_vec,i,K_vec, r_vec);
		}
		fprintf(likeli_out, "%f \n", full_loglikeli_out); //output
		
		print_array_int(cvec_out, cvec, n); //print ranker cluster indicators
		print_array_dbl(gamma_out, gam, no_rank_clust); //print gamma vec
		print_array_dbl(p_out, ptilde, n); //print probability w_i = 1, given all other quantities and the data
		print_array_int(rankers_p_out, r_vec, n); //print ranker weights w_i
		fprintf(alpha_out, "%f \n", alpha); //print alpha
		
		//print the full parameter vector corresponding to each ranker
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < (K-1); j++)
			{
				fprintf(lambda_files[i], "%f ", log(lambda[cvec[i]][d[cvec[i]][j]]));
			}
			fprintf(lambda_files[i], "%f\n", log(lambda[cvec[i]][d[cvec[i]][K-1]]));
		}	
		
			
	}
	
	if(itercount > 0 && itercount % (thin*1000) == 0) //print iteration and time update to screen
	{
		time(&rawtime);
		timeinfo = localtime(&rawtime);
		printf("Current iteration: %d \t Current local time and date: %s", itercount, asctime(timeinfo));
	}	
	
	
}
	
return 0;
}

void print_array_int(FILE* file_name, int* array, int size)
{	
	for(int i = 0; i < (size-1); i++)
	{
		fprintf(file_name, "%d ", array[i]);
	}
	fprintf(file_name, "%d\n", array[size-1]);
} 

void print_array_dbl(FILE* file_name, double* array, int size)
{	
	for(int i = 0; i < (size-1); i++)
	{
		fprintf(file_name, "%f ", array[i]);
	}
	fprintf(file_name, "%f\n", array[size-1]);
} 



int** make_int_mat(int n, int m)
{
	int** arr = (int**)calloc(n,sizeof(int*));
	for(int i = 0; i < n; i++)
	{
		arr[i] = (int*)calloc(m,sizeof(int));
	}
	return arr;
}

void zero_int_mat(int** mat, int dim_1, int dim_2) 
{
	for(int i = 0; i < dim_1; i++)
	{
		memset(mat[i], 0, dim_2 * sizeof(int));
	}
	
}

double** make_double_mat(int n, int m)
{
	double** arr = (double**)calloc(n,sizeof(double*));
	for(int i = 0; i < n; i++)
	{
		arr[i] = (double*)calloc(m,sizeof(double));
	}
	return arr;
}

int sample_bern(double U, double prob)
{
	if(U < prob)
	{
		return(1);
	}else
	{
		return(0);
	}
}

int sample_cum_probs(double U, double* probs, int size)
{
	int sample;

	for(int i = 0; i < size; i++)
	{
		if(U <= probs[i])
		{
			sample = i;
			break;
		}
	}

	return(sample);
}

int vecmax_int(int* b, int size) 
{
	int max;
	max = b[0];
	for(int i = 1; i < size; i++)
	{
		if(max < b[i])
		{
			max = b[i];
		}
	}

	return(max);
}

double vecmax_dbl(double* b, int size)
{
	double max;
	max = b[0];
	for(int i = 1; i < size; i++)
	{
		if(max < b[i])
		{
			max = b[i];
		}
	}

	return(max);
}


void swap_labels_arr( int* arr, int label_1, int label_2, int* m_vec, int size) //swaps two labels in an array of length n, as long as it knows the number of times each label is used (saved in m_vec)
{	
	if(m_vec[label_1] == 0) 
	{
		if(m_vec[label_2] != 0) //we have m_1 = 0, m_2 > 0 so swap labels in arr from label 2 to label 1
		{
			for(int i = 0; i < size; i++)
			{
				if(arr[i] == label_2)
				{
					arr[i] = label_1;
				}	
			}	
		}
	}
	else //we have m_1 > 0
	{
		if(m_vec[label_2] == 0) //we have m_1 > 0, m_2 = 0 so swap labels in arr from label 1 to label 2
		{
			for(int i = 0; i < size; i++)
			{
				if(arr[i] == label_1)
				{
					arr[i] = label_2;
				}	
			}
		}else //we have m_1 > 0, m_2 > 0, worse case scenario, swap by brute force!
		{
			for(int i = 0; i < size; i++)
			{
				if(arr[i] == label_1)
				{
					arr[i] = -1;
				}
				if(arr[i] == label_2)
				{
					arr[i] = label_1;
				}	
			}
			for(int i = 0; i < size; i++)
			{
				if(arr[i] == -1)
				{
					arr[i] = label_2;
				}
			}
			
		}
		
	}
	
}


void ranker_label_swap(int label_1, int label_2, int* cvec, double** lambda, int* m_vec, int** d, double* gam) //this function swaps all the required parameters when a ranker label swaps
{
	double* tmp_dbl_pt;
	int* tmp_int_pt;
	
	double temp_d;
	int temp_i;
	
	swap_labels_arr(cvec, label_1, label_2, m_vec, n); //swap labels in cvec
		
	tmp_dbl_pt = lambda[label_1]; //swap the rows in lambda_dagger. Pointer
	lambda[label_1] = lambda[label_2];
	lambda[label_2] = tmp_dbl_pt;				
						
	
	temp_i = m_vec[label_1]; //we also swap the numbers in m_vec as they may be neded for the next proposed swap
	m_vec[label_1] = m_vec[label_2];
	m_vec[label_2] = temp_i;
	

	tmp_int_pt = d[label_1]; //swap the rows of d to maintain entity clustering. Pointer
	d[label_1] = d[label_2];
	d[label_2] = tmp_int_pt;

	temp_d = gam[label_1]; //swap the gamma concentration parameter
	gam[label_1] = gam[label_2];
	gam[label_2] = temp_d;
	
} 


double loglikeli_oneranking(int** x, double** lambda_dagger, double** z, int* cvec, int** d, int* ni_vec, int ranki, int* K_vec, int* r_vec)
{
	double likeli = 0;
	
	if(r_vec[ranki] == 1)
	{
		double sum = 0;
		for(int j = 0; j < K_vec[ranki]; j++) //sum all lambdas
		{
			sum += lambda_dagger[cvec[ranki]][d[cvec[ranki]][x[ranki][j]-1]];
		}

		for(int j = 0; j < ni_vec[ranki]; j++)
		{
			likeli += log(lambda_dagger[cvec[ranki]][d[cvec[ranki]][x[ranki][j]-1]]) - z[ranki][j]*sum;
			sum -= lambda_dagger[cvec[ranki]][d[cvec[ranki]][x[ranki][j]-1]]; 
		}	
	}else
	{
		for(int j = 0; j < ni_vec[ranki]; j++) 
		{
			likeli +=  - z[ranki][j]*(K_vec[ranki]-j);
		}
	}
	
	return(likeli);
}


double loglikeli_oneranking_onevec(int** x, double* lambda_dagger, double** z, int* cvec, int** d, int* ni_vec, int ranki, int* K_vec, int* r_vec)
{
	double likeli = 0;
	
	if(r_vec[ranki] == 1)
	{
		double sum = 0;
		for(int j = 0; j < K_vec[ranki]; j++) //sum all lambdas
		{
			sum += lambda_dagger[d[cvec[ranki]][x[ranki][j]-1]];
		}

		for(int j = 0; j < ni_vec[ranki]; j++)
		{
			likeli += log(lambda_dagger[d[cvec[ranki]][x[ranki][j]-1]]) - z[ranki][j]*sum;
			sum -= lambda_dagger[d[cvec[ranki]][x[ranki][j]-1]]; 
		}	
	}else
	{
		for(int j = 0; j < ni_vec[ranki]; j++)
		{
			likeli +=  - z[ranki][j]*(K_vec[ranki]-j);
		}
	}
	
	return(likeli);
}




int sample_probs(double U, double* prob, int size)
{
	int sample;
	long double cum_prob = 0.0;
	
	for(int i = 0; i < size; i++)
	{
		cum_prob += prob[i];
		if(U <= cum_prob)
		{
			sample = i;
			break;
		}
	}

	return(sample);
}

//eof
