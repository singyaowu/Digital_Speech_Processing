#include "hmm.h"
#include <math.h>
#include <assert.h>

#define N 6
#define V 6
#define SEQ 60
using namespace std;

int main(int argc, char* argv[])
{
	// ./test modellist.txt testing_data.txt result.txt
	assert(argc == 4);
	HMM hmms[5];
	int num_model = load_models( argv[1], hmms, 5);
	
	//dump_models( hmms, 5);
	double delta[5][SEQ][N];
	memset(delta, 0, sizeof(delta));

	FILE *fp_data = fopen( argv[2], "r");
	FILE *fp_result = fopen( argv[3], "w");
	assert(fp_data != NULL && fp_result != NULL);

	char seq[SEQ];	
	while(fscanf(fp_data, "%s ", seq) > 0)
	{
		int T = strlen(seq);
		assert(T < SEQ);
		int max_m = 0;
		double final_max = 0;
		for(int mi = 0; mi < 5; ++mi)
		{
			// initailization
			for(int si = 0; si < N; ++si)
			{
				delta[mi][0][si] = 
				hmms[mi].initial[si] * hmms[mi].observation[seq[0] - 'A'][si];
			}
			// recursion
			for(int t = 1; t < T; ++t)
			{
				for(int sj = 0; sj < N; ++sj)
				{
					double max = 0;
					for(int si = 0; si < N; ++si)
					{
						double temp = delta[mi][t-1][si] * hmms[mi].transition[si][sj];
						if(temp > max){
							max = temp;
						}
					}
					delta[mi][t][sj] = max * hmms[mi].observation[seq[t] - 'A'][sj];

					if(t == (T - 1) && delta[mi][t][sj] > final_max){
						final_max = delta[mi][t][sj];
						max_m = mi;
					}
				}
			}
		}// end for each model
		fprintf(fp_result, "%s\t%.10e\n", hmms[max_m].model_name, final_max);
	}
	fclose(fp_data);
	fclose(fp_result);
	return 0;
}


