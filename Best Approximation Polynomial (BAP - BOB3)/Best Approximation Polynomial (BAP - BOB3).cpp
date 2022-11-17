#include <iostream>
#include <fstream>
#include <vector>
#include <Dense>


using Eigen::MatrixXd;
using Eigen::VectorXd;

double f(double x) {
	return pow(x - 1, 4) + pow(x - 3, 3) + pow(x - 2, 2);// +100 * sin(100 * x);
	//return sin(x);
}

int main()
{
	double a = -2, b = 2;

	const int K = 4; // количество конечных элементов
	const int N = 4;
	const int L = 160;
	const int M_ = K * N; // количество базисных

	const int M = K * (N - 1) + 1; // количество узлов
	double h = (b - a) / (M - 1);

	//double *mesh = new double [M];
	std::vector<double> mesh(M);
	for (int i = 0; i < M; i++)
	{
		mesh[i] = a + i * h;
	}
	
	const int Mviz = 1920;
	double hviz = (b - a) / (Mviz - 1);
	
	//double *X = new double [Mviz];
	//double *Approx = new double [Mviz];
	std::vector<double> X(Mviz);
	std::vector<double> Approx(Mviz);
	for(int iviz = 0; iviz < Mviz; iviz++)
	{
		X[iviz] = a + iviz * hviz;
	}
	
	///////////////////
	// инициализация //
	///////////////////

	//double *D = new double[M_];
	//double *random_mesh = new double[K * L]; 
	std::vector<double> D(N);
	std::vector<double> random_mesh(K*L);
	for (int ki = 0; ki < K; ki++) 
	{
		if (ki == 0)
		{
			for (int i = 0; i < N; i++)
			{
				D[i] = 1;
				for (int j = 0; j < N; j++)
				{
					if (i != j) {
						D[i] *= (mesh[ki * (N - 1) + i] - mesh[ki * (N - 1) + j]);
					}
				}
			}
		}
		for (int l = 0; l < L; l++)
		{	
			random_mesh[ki * L + l] = mesh[ki * (N - 1)] + (mesh[(ki + 1) * (N - 1)] - mesh[ki * (N - 1)]) * double(rand()) / (RAND_MAX);
		}
	}

	//for (int i = 0; i < M_; i++) {
	//	//std::cout << D[i] << std::endl;
	//}
	//std::array<double, K* L> random_mesh_std;
	//for (int i = 0; i < K * L; i++) {
	//	random_mesh_std[i] = random_mesh[i];
	//}
	//std::sort(random_mesh_std.begin(), random_mesh_std.end());

	//for (auto x : random_mesh_std)
	//{
	//	std::cout << x << std::endl;
	//}

	/////////////////////////
	// Нахождение векторов //
	/////////////////////////

	//double** G = new double* [M];
	//double** F = new double* [M];
	std::vector<std::vector<double>> G(M);
	std::vector<std::vector<double>> F(M);
	for (int i = 0; i < M; i++)
	{
		std::vector<double> tmp1(2 * L);
		std::vector<double> tmp2(2 * L);
		G[i] = tmp1;
		F[i] = tmp2;
		//G[i] = new double[2 * L];
		//F[i] = new double[2 * L];
	}

	for (int ki = 0; ki < K; ki++)
	{
		for (int i = 0; i < N; i++) 
		{
			if (i == N - 1 && ki != K - 1)
			{
				for (int l = 0; l < L; l++)
				{
					G[ki * (N - 1) + i][l] = 1 / D[i];
					G[ki * (N - 1) + i][l + L] = 1 / D[i];
					F[ki * (N - 1) + i][l] = f(random_mesh[ki * L + l]);
					F[ki * (N - 1) + i][l + L] = f(random_mesh[(ki + 1) * L + l]);
					for (int j = 0; j < N; j++)
					{
						if (i != j)
						{
							G[ki * (N - 1) + i][l] *= (random_mesh[ki * L + l] - mesh[ki * (N - 1) + j]);
							G[ki * (N - 1) + i][l + L] *= (random_mesh[(ki + 1) * L + l] - mesh[(ki + 1) * (N - 1) + j]);
						}
					}
				}
			}
			else
			{
				for (int l = 0; l < L; l++)
				{
					G[ki * (N - 1) + i][l] = 1 / D[i];
					F[ki * (N - 1) + i][l] = f(random_mesh[ki * L + l]);
					for (int j = 0; j < N; j++)
					{
						if (i != j)
						{
							G[ki * (N - 1) + i][l] *= (random_mesh[ki * L + l] - mesh[ki * (N - 1) + j]);
						}
					}
				}
			}
		}
	}

	////////////////////////
	// заполнение матрицы //
	////////////////////////

	std::vector<std::vector<double>> LAE(M);
	std::vector<double> right(M);
	for (int i = 0; i < M; i++)
	{
		std::vector<double> tmp1(M);
		LAE[i] = tmp1;
		//G[i] = new double[2 * L];
		//F[i] = new double[2 * L];
	}

	MatrixXd LAE_ = MatrixXd::Zero(M, M);
	VectorXd right_ = VectorXd::Zero(M);


	for (int ki = 0; ki < K; ki++) 
	{
		for (int i = 0; i < N; i++) 
		{
			double scalar_product_ii = 0;
			double scalar_product_fi = 0;
			for (int j = i+1; j < N; j++) 
			{
				double scalar_product_ij = 0;
				for (int l = 0; l < L; l++)
				{
					scalar_product_ij += G[ki * (N - 1) + i][l] * G[ki * (N - 1) + j][l];
				}
				LAE_(ki * (N - 1) + i, ki * (N - 1) + j) = scalar_product_ij;
				LAE_(ki * (N - 1) + j, ki * (N - 1) + i) = scalar_product_ij;
			}
			if (i == N - 1 && ki != K - 1)
			{
				for (int l = 0; l < L; l++)
				{
					scalar_product_ii += G[ki * (N - 1) + i][l] * G[ki * (N - 1) + i][l];
					scalar_product_ii += G[ki * (N - 1) + i][L + l] * G[ki * (N - 1) + i][L + l];
					scalar_product_fi += G[ki * (N - 1) + i][l] * f(random_mesh[ki * L + l]);
					scalar_product_fi += G[ki * (N - 1) + i][L + l] * f(random_mesh[(ki + 1) * L + l]);
				}
			}
			else if (i!=0 || ki==0)
			{
				for (int l = 0; l < L; l++) 
				{
					scalar_product_ii += G[ki * (N - 1) + i][l] * G[ki * (N - 1) + i][l];
					scalar_product_fi += G[ki * (N - 1) + i][l] * f(random_mesh[ki * L + l]);
				}
			}
			else if (i == 0 && ki != 0) 
			{
				continue;
			}
			if (i == N - 1 && ki != K - 1) 
			{
				LAE_(ki * (N - 1) + i, ki * (N - 1) + i) = scalar_product_ii;
				right_(ki * (N - 1) + i) = scalar_product_fi;
			}
			else {

				LAE_(ki* (N - 1) + i, ki* (N - 1) + i) = scalar_product_ii;
				right_(ki* (N - 1) + i) = scalar_product_fi;
			}

			LAE[ki * (N - 1) + i][ki * (N - 1) + i] = scalar_product_ii;
			right[ki * (N - 1) + i] = scalar_product_fi;
		}
	}

	std::cout << LAE_ << std::endl;


	VectorXd c = LAE_.colPivHouseholderQr().solve(right_);

	std::cout << right_ << std::endl;
	std::cout << c << std::endl;

	for (int iviz = 0, ki = 0; iviz < Mviz; iviz++)
	{
		if (iviz == Mviz - 1) {
			ki = K - 1;
		}
		else if (X[iviz] >= mesh[(ki + 1) * (N - 1)])
		{
			ki++;
		}
		Approx[iviz] = 0;
		for (int i = 0; i < N; i++)
		{
			double gviz = 1 / D[i];
			for (int j = 0; j < N; j++) 
			{
				if (i != j)
				{
					gviz *= (X[iviz] - mesh[ki * (N - 1) + j]);
				}
			}
			Approx[iviz] += c(ki * (N - 1) + i) * gviz;
		}

	}

	std::ofstream ParamsFile; // файл параметров - в нем a,b и все нормы
	ParamsFile.open("Params.txt");
	ParamsFile << a << ", " << b /*<< ", "
		<< resREL_1norm << ", " << resREL_2norm << ", " << resREL_infnorm << ", "
		<< resABS_1norm << ", " << resABS_2norm << ", " << resABS_infnorm */ << std::endl;
	ParamsFile.close();

	std::ofstream meshFile, FmeshFile; // файлы точек интерполяции и значения функции в этих точках - для выделения красным
	meshFile.open("mesh.txt");
	FmeshFile.open("Fmesh.txt");
	for (int i = 0; i < M - 1; i++) {
		meshFile << mesh[i] << ", ";
		FmeshFile << f(mesh[i]) << ", ";
	}
	meshFile << mesh[M - 1] << std::endl;
	FmeshFile << f(mesh[M - 1]) << std::endl;
	meshFile.close();
	FmeshFile.close();

	std::ofstream XFile, LFile, FFile; // файлы точек построения графиков
	XFile.open("X.txt"); // Точки X
	LFile.open("L.txt"); // Значения полинома Лагранжа в этих точках
	FFile.open("F.txt"); // Значения функции в этих точках
	for (int i = 0; i < Mviz - 1; i++) {
		XFile << X[i] << ", ";
		FFile << f(X[i]) << ", ";
		LFile << Approx[i] << ", ";
	}
	XFile << X[Mviz - 1] << std::endl;
	FFile << f(X[Mviz - 1]) << std::endl;
	LFile << Approx[Mviz - 1] << std::endl;
	XFile.close();
	LFile.close();

	std::system("python plot.py"); // эта команда вызывает командную строку и включает питоновскую часть задачи
	
	/*delete [] mesh;
	mesh = nullptr;

	delete [] X;
	X = nullptr;

	delete [] Approx;
	Approx = nullptr;

	for (int i = 0; i < M; i++)
	{
		delete[] G[i];
		delete[] F[i];
	}
	delete[] G;
	delete[] F;

	delete [] random_mesh;
	delete [] D;*/
	return 0;
}