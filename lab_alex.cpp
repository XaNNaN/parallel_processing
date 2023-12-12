#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctime>
// Граничные условия
//Граничные условия первого рода – задается распределение температуры на поверхности (или границе) тела для каждого момента времени
//Граничные условия второго рода – задается значение теплового потока для каждой точки поверхности (или границы) тела в любой момент времени:(производная)
//Граничные условия третьего рода – задается взаимосвязь между потоком тепла за счет теплопроводности от твердой стенки и тепловым потоком из окружающей среды за счет температурного напора
//явная схема по трем узлам находится один узел.
#define L 100 //длина стержня
//#define lambda 45 // коэффициент теплопроводности среды
//#define ct 460 // удельная теплоемкость единицы массы
//#define p 7800 // плотность среды at (lambda/(ct*p))4.2e-6
#define at 4.2 //- коэффициент температуропроводности для стали
#define dx 0.0001 //интервалы между узлами
#define dt 0.1 // интервалы времени
#define time 10.0 //время
const int xx = L / dx; // количество узлов
const int tt = time / dt; // количество промежутков времени
int number_of_threads;
int GU_LEFT;
int GU_RIGHT;
int TEMP1, TEMP2;
//pthread_barrier_t barrier;
struct Data
{
	double* T1; // предыдущее значение температуры на стержне
	double* T2; // нынешнее значение температуры на стержне
	int index; // номер потока
	int num_threads; // количество всех потоков
};
FILE* output; // Файл, в который пишем резльтат работы программы
void thread_work(struct Data* parametr); // Функция потока
void config_file(); // Создаёт файл gnuplot
void value_node(double* T); // запись данных из T[] в файл
void border_conditions(double* T); // Вычисление граничных условий(tt_index)
int main()
{
	printf("Input 3 arguments: number of parallels, BC for left node (1-3), BC for right node (1-3)\n");
	scanf("%d %d %d", &number_of_threads, &GU_LEFT, &GU_RIGHT);
	TEMP1, TEMP2;
	if (GU_LEFT != 3)
	{
		printf("Input temperature for BC of of left node > 0\n");
		while ((!scanf("%d", &TEMP1)) || (TEMP1 < 0))
			printf("Input temperature for BC of of left node > 0\n");
	}
	else TEMP1 = 0;
	if (GU_RIGHT != 3)
	{
		printf("Input temperature for BC of of right node > 0\n");
		while ((!scanf("%d", &TEMP2)) || (TEMP2 < 0))
			printf("Input temperature for BC of of right node > 0\n");
	}
	else TEMP2 = 0;
	config_file();
	output = fopen("plot.dat", "w");
	// Выделяем память под стержень
	double* T1 = (double*)malloc(sizeof(double) * xx);
	double* T2 = (double*)malloc(sizeof(double) * xx);
	// Задаём начальную температуру на стержне
	for (int i = 1; i < xx; i++)
	{
		T1[i] = 0;
	}
	// Задать температуру на границе в начальный момент времени
	border_conditions(T1);
	// Объявляем массив под потоки
	//pthread_t threads[number_of_threads];
	// Выделяем память для данных потоков
	struct Data data;
	// Инициализируем барьер
	//pthread_barrier_init(&barrier, NULL, number_of_threads);
	//инициализация барьера. инициализирует барьер, выделяя необходимую память, устанавливая значения его атрибутов и назначая count "шириной" барьера. В настоящее время атрибуты барьеров не определены поэтому в качестве второго параметра функции pthread_barrier_init следует использовать NULL.

	double begin = clock();
	// Создаём потоки
			data.T1 = T1;
			data.T2 = T2;
			data.index = 0;
			data.num_threads = number_of_threads;
			thread_work(&data);
	// ожидание потоков
	/*for (int i = 0; i < number_of_threads; ++i)
	{
		pthread_join(threads[i], NULL);//закрытие потоков
	}*/
		double end = clock();
	// Вычисляем время работы программы
		double elapsed = (end - begin);
	printf("Threads %3d | Time, s: %.5f\n", number_of_threads, elapsed);
	fclose(output);
	free(T1);
	free(T2);
	return 0;
}
// Основная функция потоков
void thread_work(struct Data* parametr)
{
	struct Data data = *((struct Data*)parametr);
	int numb_nodes = xx;
	int index_left = 0;
	int index_right = xx;
	// вывод значений температуры в начальный момент времени
	if (data.index == 0)
	{
		value_node(data.T1);
	}


	// Решение МКР
	for (int i = 1; i < tt; ++i) // Итерирование по времени
	{
		// Расчёт граничных условий первым потоком

		border_conditions(data.T2);
		//value_node(data.T2);
		// Итерирование по узлам
		omp_set_num_threads(data.num_threads);
		{
			#pragma omp parallel for
			for (int j = 1; j < xx-1; ++j)
			{
				//printf("%d\n", omp_get_num_threads());
				if (j != 0 && j != xx - 1) // Внутренние узлы стержня
				{
					data.T2[j] = at * ((data.T1[j + 1] - 2 * data.T1[j] +
						data.T1[j - 1]) * dt) / pow(dx, 2) + data.T1[j];
				}
			}
		}
		memcpy(data.T1, data.T2, sizeof(double) * xx);
		//pthread_barrier_wait(&barrier); //подключиться к другому потоку и ожидать его завершения; поток, к которому необходимо подключиться, должен быть создан с возможностью подключения. когда все потоки закончат работу выйдет
	}
}
void border_conditions(double* T)
{
	const double C1 = 0.5;
	const double C2 = 5;
	// Расчёт левой границы
	switch (GU_LEFT) {
	case 1:
		T[0] = TEMP1;
		break;
	case 2:
		T[0] = ((double)TEMP1 * dx) / at + T[1];
		break;
	case 3:
		T[0] = (C2 * dx + T[1] * at) / (at - C1 * dx);
		break;
	}
	switch (GU_RIGHT) {
	case 1:
		T[xx - 1] = TEMP2;
		break;
	case 2:
		T[xx - 1] = ((double)TEMP2 * dx) / at + T[xx - 2];
		break;
	case 3:
		T[xx - 1] = (C2 * dx + T[xx - 2] * at) / (at - C1 * dx);
		break;
	}
}
void config_file()
{
	FILE* f = fopen("config", "w");
	fprintf(f, "set terminal gif animate delay 50\nset output 'temp.gif'\nset ticslevel 0\nset key off\n");
	fprintf(f, "set style line 1 \\\nlinecolor rgb '#fc89ac' \\\nlinetype 1 linewidth 2 \\\n\n");
	fprintf(f, "do for [i=0:%d] { \n", tt - 1);
	fprintf(f, "plot 'plot.dat' index i with lines linestyle 1}\n"); //fprintf(f, "pause 0.0005\n}");
	fclose(f);
}
void value_node(double* T)
{
	for (int i = 0; i < xx; ++i)
	{
		fprintf(output, "%d %f\n", i, T[i]);
	}
	fprintf(output, "\n\n");
}