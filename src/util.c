#include "util.h"
#include "score.h"
#include "BDE.h"

#include <stdlib.h>
#include <stdio.h>
#include <syslog.h>

double get_score(void *buff, int child, int *parents, int size) {
	//return (double)score_calc(buff, child, parents, size);
	return (double)BDE_score(buff, child, parents, size);
}

void util_print_dmatrix(double *c, int sizen, int sizem) {
	FILE *f = fopen("./data/output.txt", "w");
	if (f == NULL) {
		printf("Error opening file\n");
		exit(1);
	}

	for (int i = 0; i < sizen; ++i) {
		for (int j = 0; j < sizem; ++j)
			//printf("%f, ", c[i * sizem + j]);
			fprintf(f, "%f, ", c[i * sizem + j]);
		fprintf(f, "\n");
	}
	fclose(f);
}

void util_print_imatrix(int *c, int size, const char *filename) {
	FILE *f = fopen(filename, "w");
	if (f == NULL) {
		printf("Error opening file\n");
		exit(1);
	}

	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			fprintf(f, "%d", c[i * size + j]);
			if (j < (size-1))
				fprintf(f, ",");
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

void util_errlog(char *s) {
	syslog(LOG_PERROR, "%s\n", s);
	//fprintf(stderr, "%s", s);
}

void util_debuglog(char *s) {
	syslog(LOG_INFO, "%s\n", s);
}

void util_print_score_table(void *buff) {
	int parents[5] = {0,1,2,3,4};
	int four[4] = {0,0,0,0};
	int three[3] = {0,0,0};
	int two[2] = {0,0};
	double score;
	// repeat 5 times;
	for(int i = 0; i < 5; ++i) {
		printf("child: %d =======\n",parents[0]);

		// four set
		for (int j = 0; j < 4; ++j)
			four[j] = parents[j+1];
		printf("parent set %d, {%d, %d, %d, %d} -> ",parents[0], four[0],four[1],four[2],four[3]);
		score = get_score(buff, parents[0], four, 4);
		printf("%f\n",score);

		// three set
		int tmp[4];
		for (int p=0; p<4; ++p) tmp[p] = parents[p+1];
		for (int j = 0; j < 4; ++j) {
			for (int k=0; k<3; ++k) three[k] = tmp[k];
			score = get_score(buff, parents[0], three, 3);
			printf("parent set %d, {%d, %d, %d} -> ",parents[0], three[0],three[1],three[2]);
			printf("%f\n",score);

			int temp = tmp[3];
			for(int m=3; m>0; --m) tmp[m] = tmp[m-1];
			tmp[0] = temp;
		}

		// two set
		printf("parent set %d, {%d, %d} -> ",parents[0],tmp[0],tmp[1]);
		two[0] = tmp[0];
		two[1] = tmp[1];
		score = get_score(buff, parents[0], two, 2);
		printf("%f\n",score);
		//printf("parent set {%d, %d} -> ",tmp[1],tmp[0]);
		//two[0] = tmp[1];
		//two[1] = tmp[0];
		//score = get_score(buff, two, 2);
		//printf("%f\n",score);
		printf("parent set %d, {%d, %d} -> ",parents[0],tmp[1],tmp[2]);
		two[0] = tmp[1];
		two[1] = tmp[2];
		score = get_score(buff, parents[0], two, 2);
		printf("%f\n",score);
		//printf("parent set {%d, %d} -> ",tmp[2],tmp[1]);
		//two[0] = tmp[2];
		//two[1] = tmp[1];
		//score = get_score(buff, two, 2);
		//printf("%f\n",score);
		printf("parent set %d, {%d, %d} -> ",parents[0],tmp[2],tmp[3]);
		two[0] = tmp[2];
		two[1] = tmp[3];
		score = get_score(buff, parents[0], two, 2);
		printf("%f\n",score);
		printf("parent set %d, {%d, %d} -> ",parents[0],tmp[0],tmp[2]);
		two[0] = tmp[0];
		two[1] = tmp[2];
		score = get_score(buff, parents[0], two, 2);
		printf("%f\n",score);
		printf("parent set %d, {%d, %d} -> ",parents[0],tmp[1],tmp[3]);
		two[0] = tmp[1];
		two[1] = tmp[3];
		score = get_score(buff, parents[0], two, 2);
		printf("%f\n",score);
		printf("parent set %d, {%d, %d} -> ",parents[0],tmp[0],tmp[3]);
		two[0] = tmp[0];
		two[1] = tmp[3];
		score = get_score(buff, parents[0], two, 2);
		printf("%f\n",score);

		// one set
		int test[] = { 0 };
		for (int j = 0; j<4; ++j) {
			/*printf("parent set {%d} -> ",tmp[j]); */
			//score = get_score(buff, &tmp[j], 1);
			test[0] = tmp[j];
			printf("parent set %d, {%d} -> ",parents[0], test[0]);
			score = get_score(buff, parents[0], test, 1);
			printf("%f\n",score);
		}

		// null set
		test[0] = i;
		printf("parent set %d, {} -> ", parents[0]);
		score = get_score(buff, parents[0], test, 0);
		printf("%f\n",score);

		//rotate the parent array
		int temp = parents[4];
		for(int j = 4; j > 0; --j) parents[j] = parents[j-1];
		parents[0] = temp;
	}
}
