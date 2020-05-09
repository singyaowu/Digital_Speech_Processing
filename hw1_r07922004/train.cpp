#include "hmm.h"
#include <math.h>
#include <assert.h>
#include <vector>
#include <string>
#define N 6
#define V 6
#define SEQ 60
using namespace std;

void forward_algoritm(HMM &hmm, double alpha[][N], const char* seq)
{
	int T = strlen(seq);
	for(int i = 0; i < N; ++i)
		alpha[0][i] = hmm.initial[i] * hmm.observation[seq[0] - 'A'][i];
	
	for(int t = 1; t < T; ++t)
	{
		for(int sj = 0; sj < N; ++sj)
		{
			alpha[t][sj] = 0;
			for(int si = 0; si < N; ++si)
				alpha[t][sj] += alpha[t - 1][si] * hmm.transition[si][sj];
			alpha[t][sj] *= hmm.observation[seq[t] - 'A'][sj];
			assert(alpha[t][sj] <=1);
		}
	}
	return;
}
void backward_algoritm(HMM &hmm, double beta[][N], const char* seq)
{
	int T = strlen(seq);
	for(int i = 0; i < N; ++i)
		beta[T - 1][i] = 1;
	
	for(int t = T - 2; t >= 0; --t)
	{
		for(int si = 0; si < N; ++si)
		{
			beta[t][si] = 0;
			for(int sj = 0; sj < N; ++sj)
			{
				beta[t][si] += 
				  hmm.transition[si][sj] * hmm.observation[seq[t+1]- 'A'][sj] * beta[t + 1][sj];
				assert(beta[t][si] <= 1);
			}
		}
	}
	return;
}
void Baum_welch_alogrithm(HMM &hmm, int T, int len_seq, double epsi[][N][N], double gama[][N], double gama_v[][V])
{
	// determine pi
	for(int i = 0; i < N; ++i)
		hmm.initial[i] = gama[0][i] / len_seq;
	
	// determine a
	double gama_sum_n1[N];
	memset(gama_sum_n1, 0, sizeof(gama_sum_n1));

	for(int i = 0; i < N; ++i)
	{
		for(int t = 0; t < T - 1; ++t)
			gama_sum_n1[i] += gama[t][i];
	}
	for(int i = 0; i < N; ++i)
	{
		for(int j = 0; j < N; ++j)
		{
			hmm.transition[i][j] = 0;
			for(int t = 0; t < T - 1; ++t)
				hmm.transition[i][j] += epsi[t][i][j];
			hmm.transition[i][j] /= gama_sum_n1[i];
			assert(hmm.transition[i][j] <= 1 && hmm.transition[i][j] >= 0);
		}
	}
	// determine b
	for(int oi = 0; oi < V; ++oi)
	{
		for(int i = 0; i < N; ++i)
		{
			hmm.observation[oi][i] = gama_v[oi][i] / (gama_sum_n1[i] + gama[T-1][i]);
			//assert(hmm.observation[oi][i] <= 1 && hmm.observation[oi][i] >= 0);
		}
	}

	return;
}
int main(int argc, char* argv[])
{
	//./train iteration model_init.txt seq_model_01.txt model_01.txt
	double alpha[SEQ][N];
	double beta[SEQ][N];
	double epsi[SEQ][N][N];
	double gama[SEQ][N];
	double gama_v[V][N];
	memset(alpha, 0.0, sizeof(alpha));
	memset(beta, 0.0, sizeof(beta));
	memset(epsi, 0.0, sizeof(epsi));
	memset(gama, 0.0, sizeof(gama));
	memset(gama_v, 0.0, sizeof(gama_v));

	assert(argc == 5);
	int iteration = atoi(argv[1]);
	char* models_str = argv[3];
	
	HMM hmm;
	loadHMM( &hmm, argv[2]);
	dumpHMM( stderr, &hmm);

	FILE *fp = fopen( models_str, "r");
	assert(fp != NULL);

	char seq[SEQ];
	vector<string> sequences;
	while(fscanf(fp, "%s ", seq) > 0)
	{
		assert(strlen(seq) <= SEQ);
		sequences.push_back(string(seq));
	}
	fclose(fp);
	
	int len_seq = sequences.size();
	//fprintf(stderr, "len_seq = %d\n", len_seq);

	int T = 0;
	for(int iter = 0; iter < iteration; ++iter)
	{	
		for(int seqi = 0; seqi < len_seq; ++seqi)
		{
			//fprintf(stderr, "iter = %d, seqi = %d\n", iter, seqi);
			strcpy(seq, sequences[seqi].c_str());
			T = strlen(seq);

			forward_algoritm(hmm, alpha, seq);
			backward_algoritm(hmm, beta, seq);
			
			double PO_l = 0;
			for(int si = 0; si < N; ++si)
				PO_l += alpha[0][si] * beta[0][si];
			assert(PO_l <= 1);

			// gama
			for(int t = 0; t < T; ++t)
			{
				for(int i = 0; i < N; ++i)
				{
					double temp = alpha[t][i] * beta[t][i] / PO_l;
					gama[t][i] += temp;
					gama_v[seq[t] - 'A'][i] += temp;
				}
			}
			// epsilon
			for(int t = 0; t < T-1; ++t)
			{
				for(int si = 0; si < N; ++si)
				{
					for(int sj = 0; sj < N; ++sj)
						epsi[t][si][sj] += 
						alpha[t][si] * hmm.transition[si][sj] * hmm.observation[ seq[t+1]-'A'][sj] * beta[t+1][sj] / PO_l;
				}
			}

		}
		// final 
		Baum_welch_alogrithm(hmm, T, len_seq, epsi, gama, gama_v);
		fprintf(stderr, "\nAfter iter %d:\n", iter);
		dumpHMM( stderr, &hmm);
		
		memset(epsi, 0.0, sizeof(epsi));
		memset(gama, 0.0, sizeof(gama));
		memset(gama_v, 0.0, sizeof(gama_v));
		
		// check sum
		double compare1[N];
		double compare2[N];
		memset(compare1, 0.0, sizeof(compare1));
		memset(compare2, 0.0, sizeof(compare2));
		for(int i = 0; i < N; ++i)
		{
			for(int j = 0; j < N ; ++j)
				compare1[i] += hmm.transition[i][j];
		}
		for(int i = 0; i < N; ++i)
		{
			for(int j = 0; j < V ; ++j)
				compare2[i] += hmm.observation[j][i];
		}
		for(int i = 0; i < N; ++i)
			printf("%.10lf ", compare1[i]);
		printf("\n");
		for(int i = 0; i < N; ++i)
			printf("%.10lf ", compare2[i]);
		printf("\n");

	}
	FILE *fp_out = fopen( argv[4], "w");
	assert(fp_out != NULL);
	dumpHMM( fp_out, &hmm);
	fclose(fp_out);

	return 0;
}