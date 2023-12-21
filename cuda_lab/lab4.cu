#include <stdio.h>
#include <stdlib.h>
// #include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>




#define TIME_STEP 0.05
#define FILE_NAME "result_data.txt"
#define IA 10	
#define R 1000
#define C 1e-3
#define U0 100


void data_to_file(double*, int, int);
void config_plot(int, int, int);
__global__ void calculate_new_step(double*, double*, const int, int);

// Аргументы ком. строки - исп. ф., кол-во элементов сетки, время моделирования
int main(int argc, char* argv[]){
	int M, M_gpu;
	int N, N_gpu;
	double modeling_time;
	int i, j;
	double *matrix_1, *matrix_2; // две матрицы для старого и нового временных слоёв
	double cur_time = 0; // прошедшее время эксперимента
	FILE *fp;
	double *new_layer, *old, *tmp;
	double *new_layer_gpu, *old_gpu;
	double my_perfect_const = TIME_STEP / (R * C);
	int bot_row, top_row;
	

	cudaEvent_t start, stop;
	float GPUTime = 0.0f;




	

	cudaEventCreate(&start);
	cudaEventCreate(&stop);


	printf("my_perfect_const %f\n", my_perfect_const);
	if (argc < 3){
		printf("Недостаточно аргументов.\n");
		printf("Аргументы ком. строки - исп. ф., кол-во элементов строке, кол-во элементов в столбце, время моделирования.\n");
		return -1;
	}

  	M = atoi(argv[1]);
	M_gpu = M;
	N = atoi(argv[2]);
	N_gpu = N;
	modeling_time = atof(argv[3]);

	int N_threads = 512;
	int N_blocks;

	if ( (M * N % N_threads) == 0 )
		N_blocks = ( M  / N_threads);
	else
		N_blocks = ( M  / N_threads) + 1;

	dim3 Threads(N_threads);
	dim3 Blocks(N_blocks);

	dim3 block_shape = dim3(32, 32);
	dim3 grid_shape = dim3(max(1.0, ceil((float) M / (float) block_shape.x )),
						   max(1.0, ceil((float) M / (float) block_shape.y )));

	printf("Grid shape: %d %d\n", grid_shape.x, grid_shape.y);
	printf("Grid shape: %d %d\n", block_shape.x, block_shape.y);

	printf("Threads: %d\n", Threads.x);
	printf("Blocks: %d\n", Blocks.x);

	bot_row = 0;
	top_row = M; // индекс верхней строки + 1 
	printf("bot_row - %d, top_row - %d\n", bot_row, top_row);

	if (M % 8 != 0 || N % 8 != 0){
		printf("Количество элементов по каждой из осей сетки должно быть кратно 8.\n");
		return -1;
	}
	 
	fp = fopen(FILE_NAME, "w+");
	fclose(fp);


	config_plot(M, N, (double)modeling_time / (double)TIME_STEP + 1 );


	matrix_1 = (double*) calloc(M * N, sizeof(double));
	matrix_2 = (double*) calloc(M * N, sizeof(double));
	old = matrix_2;
	new_layer = matrix_1; 
	

	// Allocating space for the arrays on the gpu.
	cudaMalloc((void **) &new_layer_gpu, N * M * sizeof(double));
	cudaMalloc((void **) &old_gpu, N * M * sizeof(double));


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
	cudaEventRecord(start, 0);
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

		// bot_row, top_row, new_layer, old, my_perfect_const, N, M
		// int: bot_row, top_row, N, M.
		// double *: new_layer, old. 
		cudaMemcpy(new_layer_gpu, new_layer, N * M *sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(old_gpu, old, N * M * sizeof(double), cudaMemcpyHostToDevice);

		calculate_new_step<<<Blocks, Threads>>>(old_gpu, new_layer_gpu, N_gpu, M_gpu);
		
		cudaMemcpy(new_layer, new_layer_gpu, N * M * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(old, old_gpu, N * M * sizeof(double), cudaMemcpyDeviceToHost);

		new_layer[0] = new_layer[N-1] = new_layer[(M-1) * N] = new_layer[M*N -1] = U0;

		// вынесение значений на жёсткий диск
		data_to_file(old, M, N);

		tmp = old;
		old = new_layer;
		new_layer = tmp; 
		cur_time = cur_time + TIME_STEP; 
	} 
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);

	cudaEventElapsedTime(&GPUTime, start, stop);
	printf("GPU time: %.3f ms or %.0f s\n", GPUTime, GPUTime/1000);

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	free(new_layer);
	free(old);

	cudaFree(new_layer_gpu);
	cudaFree(old_gpu);


    return 0;
}

__global__ void calculate_new_step(double* old_layer, double* new_layer, const int N, int M){
	double my_perfect_const = TIME_STEP / (R * C);

	int i = blockIdx.x * blockDim.x + threadIdx.x;
	// int j = blockIdx.y * blockDim.y + threadIdx.y;
	// i = threadIdx.y;
	// printf("%d\n", i);

	if (i < M) {
		for (int j = 0; j < N; j++){
			if (j == 0) // лево
			{
				new_layer[N * i + j] = (old_layer[N * i + (j-1)] +  // U_down
					old_layer[N * (i+1) + j] +   // U_right
					old_layer[N * i + (j+1)] -   // U_top
					3* old_layer[N * i + j])* my_perfect_const + old_layer[N * i + j];
				
			}
			if (j == N-1) // право
			{
				new_layer[N * i + j] = (old_layer[N * (i-1) + j] +	// U_left
					old_layer[N * i + (j-1)] + 	// U_down
					old_layer[N * i + (j+1)] -   // U_top
					3* old_layer[N * i + j])* my_perfect_const + old_layer[N * i + j];
			} // centr
			if ( j != 0 && j != N-1 && i != 0 && i != M -1){
				new_layer[N * i + j] = (old_layer[N * i + (j-1)] +  // U_left
					old_layer[N * (i-1) + j] + 	// U_down
					old_layer[N * i + (j+1)] +  // U_right
					old_layer[N * (i+1) + j] -  // U_top
					4* old_layer[N * i + j])* my_perfect_const + old_layer[N * i + j];
			}
		}
	}
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