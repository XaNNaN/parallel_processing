#include <stdio.h>
#include <stdlib.h>
// #include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <omp.h>
#include <time.h>




#define TIME_STEP 0.05
#define FILE_NAME "result_data.txt"
#define IA 10	
#define R 1000
#define C 1e-3
#define U0 100


void data_to_file(double*, int, int);
void config_plot(int, int, int);

// Аргументы ком. строки - исп. ф., кол-во элементов сетки, время моделирования
int main(int argc, char* argv[]){
	int M;
	int N;
	double modeling_time;
	int i, j;
	double *matrix_1, *matrix_2; // две матрицы для старого и нового временных слоёв
	double cur_time = 0; // прошедшее время эксперимента
	FILE *fp;
	double *new_layer, *old, *tmp, *sending_new;
	double my_perfect_const = TIME_STEP / (R * C);
	int myrank, total;
	int rows_per_proc;
	int bot_row, top_row;




	printf("my_perfect_const %f\n", my_perfect_const);
	if (argc < 3){
		printf("Недостаточно аргументов.\n");
		printf("Аргументы ком. строки - исп. ф., кол-во элементов строке, кол-во элементов в столбце, время моделирования.\n");
		return -1;
	}

  	M = atoi(argv[1]);
	N = atoi(argv[2]);
	modeling_time = atof(argv[3]);
	//rows_per_proc = M / total;
	bot_row = 0;
	top_row = M; // индекс верхней строки + 1 
	printf("bot_row - %d, top_row - %d\n", bot_row, top_row);

	if (M % 8 != 0 || N % 8 != 0){
		printf("Количество элементов по каждой из осей сетки должно быть кратно 8.\n");
		return -1;
	}
	 
	omp_set_num_threads(atoi(argv[4]));

	fp = fopen(FILE_NAME, "w+");
	fclose(fp);

	 config_plot(M, N, (double)modeling_time / (double)TIME_STEP + 1 );


	matrix_1 = (double*) calloc(M * N, sizeof(double));
	matrix_2 = (double*) calloc(M * N, sizeof(double));
	sending_new = (double*) calloc(N * rows_per_proc, sizeof(double));
	old = matrix_2;
	new_layer = matrix_1; 


	// инициализация нулевого момента времени, выенесение информации на диск и переход к след слою 
	for (i = 0; i < M; i++){
		for (j = 0; j < N; j++ ){
			old[N*i + j] = 0;
		}
	}

	old[0] = U0;
	old[N-1] = U0;
	old[(M-1) * N] = U0;
	old[M*N -1] = U0; // источники постоянного напряжения по углам

	/*if(!myrank){
		data_to_file(matrix_2, M, N);
	}/*
	printf("my rank before cycle %d\n", myrank);

	MPI_Barrier(MPI_COMM_WORLD); */
	clock_t tStart = clock();
	while (cur_time < modeling_time){
		// расчёт значений для нового временного слоя
		 // нижняя грань

		bot_row = 0;
		top_row = M;

		i = 0;
		for (j = 1; j < N-1; j++){
			new_layer[N * i + j] = (old[N * i + (j-1)] +  // U_left
				old[N * i + (j+1)] +  // U_right
				old[N * (i+1) + j] -  // U_top
				3* old[N * i + j])* my_perfect_const + old[N * i + j];
		}
		new_layer[0] = new_layer[N-1] =  U0;
		bot_row++;
	



		 // верхняя грань
		i = M - 1;
		for (j = 1; j < N-1; j++){
			new_layer[N * i + j] = (old[N * i + (j-1)] +  // U_left
				old[N * (i-1) + j] + 	// U_down
				old[N * i + (j+1)] -  // U_right
				3* old[N * i + j])* my_perfect_const + old[N * i + j];
		}
		new_layer[(M-1) * N] = new_layer[M*N -1] = U0;
		top_row--;
		#pragma omp parallel for private(j)
		for (int i = bot_row; i < top_row; i++ ){
			for (j = 0; j < N; j++){
				if (j == 0) // лево
				{
					new_layer[N * i + j] = (old[N * i + (j-1)] +  // U_down
						old[N * (i+1) + j] +   // U_right
						old[N * i + (j+1)] -   // U_top
						3* old[N * i + j])* my_perfect_const + old[N * i + j];
				}
				if (j == N-1) // право
				{
					new_layer[N * i + j] = (old[N * (i-1) + j] +	// U_left
						old[N * i + (j-1)] + 	// U_down
						old[N * i + (j+1)] -   // U_top
						3* old[N * i + j])* my_perfect_const + old[N * i + j];
				} // centr
				if ( j != 0 && j != N-1 && i != 0 && i != M -1){
					new_layer[N * i + j] = (old[N * i + (j-1)] +  // U_left
						old[N * (i-1) + j] + 	// U_down
						old[N * i + (j+1)] +  // U_right
						old[N * (i+1) + j] -  // U_top
						4* old[N * i + j])* my_perfect_const + old[N * i + j];
				}
			}
		}
		

		new_layer[0] = new_layer[N-1] = new_layer[(M-1) * N] = new_layer[M*N -1] = U0;
		// вынесение значений на жёсткий диск
		// new_to_all();
		for (i = bot_row; i < top_row; i++ ){
			for (j = 0; j < N; j++){
				sending_new[(i - bot_row) * N + j] = new_layer[N * i + j];
			}
		}
		
		
		
		double * ptr = &new_layer[bot_row * N];

		// data_to_file(old, M, N);

 
		tmp = old;
		old = new_layer;
		new_layer = tmp; 
		cur_time = cur_time + TIME_STEP; 
	} 
	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    exit(0);
}



void data_to_file(double* matrix, int M, int N){
	int i = 0;
	int j = 0;
	FILE *fp;

	fp = fopen(FILE_NAME, "a");

	for (i = 0; i < M; i++){
		for (j = 0; j < N; j++){

			fprintf(fp, "%-5.3f  ",  matrix[N * i + j]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp, "\n\n");

	fclose(fp);
}

void config_plot(int M, int N, int time_layers){
	FILE *fp;
	printf("%d\n", time_layers);
	fp = fopen("myplot", "w");
	fprintf(fp, "n = n+1\n");
	fprintf(fp, "splot \"%s\" index n-1 matrix with lines\n", FILE_NAME);
	fprintf(fp, "print n\n");
	fprintf(fp, "pause 0.01\n");
	fprintf(fp, "if (n < %d) reread; else pause -1\n", time_layers);
	fprintf(fp, "\n");
	fclose(fp);

	fp = fopen("cycle", "w");
	fprintf(fp, " n = 0\n");
	fprintf(fp,"set xrange [0:%d]\n", N);
	fprintf(fp,"set yrange [0:%d]\n", M);
	fprintf(fp, "load \"myplot\"\n");

}