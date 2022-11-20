#include <iostream>
#include <fstream>
#include <vector>
#include <Dense>


using Eigen::MatrixXd;
using Eigen::VectorXd;

double f(double x) {
	return pow(x - 1, 4) + pow(x - 3, 3) + pow(x - 2, 2) + 50 + 100 * sin(100 * x);
	//return sin(x);
}

int main()
{
	double a = -2, b = 2; // ������ � ����� �������

	const int K = 100; // ���������� �������� ���������
	const int N = 6; // ���������� ����� �� 1 �������� ��������
	const int L = 40; // ���������� ��������� �����
	const int M_ = K * N; // ���������� �������� (��������) �������

	const int M = K * (N - 1) + 1; // ���������� �������� (��������) �������, ��� �� ���������� �����
	double h = (b - a) / (M - 1); // ��� �� ����������� �����

	double mesh[M]; // ����������� ����� �� �����
	for (int i = 0; i < M; i++)
	{
		mesh[i] = a + i * h; // ��������� ������ �����
	}
	
	const int Mviz = 1920; // ���������� ����� ��� ���������
	double hviz = (b - a) / (Mviz - 1); // ��� ��� ���
	
	double X[Mviz]; // ������ ����� ��� ���������
	double Approx[Mviz]; // �������� ���� � ���� ������
	for(int iviz = 0; iviz < Mviz; iviz++)
	{
		X[iviz] = a + iviz * hviz; // ���������� ������� ����� ��� ���������
	}
	
	///////////////////
	// ������������� //
	///////////////////

	double D[N]; // ������ ������������
				 // ��� ��� ���� �� ����������� �����
				 // �� �� ������ �������� ��������
				 // ����������� ����� �����������
				 // �� ���� �� ���������� ������
				 // N ������������

	double random_mesh[K * L]; // ������ ���� ��������� ����� ��� �������� ���������� ������������
							   // ��� ��� ��� ���������� ���� ����� ������� ������ ����
							   // �� ��� �������� ������� ��� ����� �����
							   // 
	// ��� ����� �������� ������ ������������
	// � ������ ��������� �����
	for (int ki = 0; ki < K; ki++) // ����������� �� �������� ���������
	{
		// ����������� 
		// ki - ����� � 0 �������� ��������� ��������
		if (ki == 0) // ����������� ����� ��������� ������ ���� ���
					 // �������� ����� �������
					 // ����� ��� ����������� 
					 // ������ �� ������ �������� ��������
		{
			for (int i = 0; i < N; i++) // ���������� ������������ �(xi - xj) �� j
			{
				// mesh[ki * (N - 1) + i] - i-�� ����� �������� ��������� ��������

				D[i] = 1;
				for (int j = 0; j < N; j++)
				{
					if (i != j) {
						D[i] *= (mesh[ki * (N - 1) + i] - mesh[ki * (N - 1) + j]);
					}
				}
			}
		}

		for (int l = 0; l < L; l++) // ���������� ����������� �����
		{	
			// ����������� 
			// ki - ����� � 0 �������� ��������� ��������
			// mesh[ki * (N - 1)] - ������ �������� ��������� ��������
			// mesh[(ki + 1) * (N - 1)] - ������ ���������� ��������� ��������
			// random_mesh[ki * L + l] - l-�� ��������� ����� ��������� ��������
			random_mesh[ki * L + l] = mesh[ki * (N - 1)] + (mesh[(ki + 1) * (N - 1)] - mesh[ki * (N - 1)]) * double(rand()) / (RAND_MAX);
		}
	}


	/////////////////////////
	// ���������� �������� //
	/////////////////////////
	
	// G[i] - �������� �� i-�� ������� 
	// �� ���������� ������� �� ��������, 
	// ������� ����� �� M
	// G[i][l] - �������� i-�� ������� �� l-�� ��������� �����
	// �� ���� G[i] - ������ �������� ������� �� L ��������� ������
	// ��� ��� � ��� ���� �������� �� ������� �������� ��������� �������
	// �� ��� G[i], ��� i % N - 1 == 0 (�������� �� ��������)
	// � ��� ���������� 2 �������� ���������
	// ����� ��������� ������� ����� ������� �� 2*L ������
	// � ��������� ��������� ������������
	// �� �������������� �������� ���������

	double G[M][2 * L]; 
	double F[M][2 * L];
		
	for (int ki = 0; ki < K; ki++) // ����������� �� �������� ���������
								   // ����� ����������� ������
								   // N = 4 K = 3
								   // 0 0 0 * 0 0 * 0 0 0
								   // 0 - ������� ������������ �� ������ ����
								   // * - ������� ������������ �� ���� �����
								   // �� ���� * ��������� �������
								   // � ���� ����� ����� ������� �������� ���� �������
								   // � ��������� ������

	{
		for (int i = 0; i < N; i++) 
		{
			if (i == N - 1 && ki != K - 1) // ������� ��������� � * �������, ���� ��� ��������� � �������� ��������
			{
				for (int l = 0; l < L; l++)
				{	
					G[ki * (N - 1) + i][l] = 1 / D[N - 1]; // �������� ������� �� ����� �������� ��������
					G[ki * (N - 1) + i][l + L] = 1 / D[0]; // �������� ������� �� ������ �������� ��������

					// (f,g) ��� ������ ����� � ����
					// ��� ��� �������� f �� 0 �� ���� ������� [a,b]
					// �� ������ g �� 0 ������ �� ����� ��� ���� �������� ���������
					// �� ���������� ������� ������ �� ����������� �������� ���������

					F[ki * (N - 1) + i][l] = f(random_mesh[ki * L + l]); // �������� f �� ����� �������� ��������
					F[ki * (N - 1) + i][l + L] = f(random_mesh[(ki + 1) * L + l]); // �������� f �� ������ �������� ��������

					for (int j = 0; j < N; j++)
					{
						if (j != N - 1) // ��� ������ ��������� �������� ������� ��������� �� ��������� ���� ����� ��������
						{
							G[ki * (N - 1) + i][l] *= (random_mesh[ki * L + l] - mesh[ki * (N - 1) + j]);
							
						}
						if (j != 0) // ��� ������� ��������� �������� ������� ��������� �� ������ ���� ����� ��������
						{
							G[ki * (N - 1) + i][l + L] *= (random_mesh[(ki + 1) * L + l] - mesh[(ki + 1) * (N - 1) + j]);
						}
					}
				}
			}
			else if (i == 0 && ki != 0)
			{
				continue; // ��� ��� ������� ������� ��� ��������� ���������� ��� ��������
			}
			else
			{
				for (int l = 0; l < L; l++)
				{
					// ����������, �� ������ ��� ������, ����������� ��������� ��������
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
	// ���������� ������� //
	////////////////////////

	MatrixXd LAE_ = MatrixXd::Zero(M, M);
	VectorXd right_ = VectorXd::Zero(M);

	// �������� ��� ������ ���� 
	// 0 0 0 * 0 0 * 0 0 0 - �������
	// 0 1 2 3 4 5 6 7 8 9 - �� ������
	// ����� ������� ����� ��������� ���
	// (0,0) (0,1) (0,2) (0,3)   0     0     0     0     0     0 
	// (1,0) (1,1) (1,2) (1,3)   0     0     0     0     0     0
	// (2,0) (2,1) (2,2) (2,3)   0     0     0     0     0     0
	// (3,0) (3,1) (3,2) (3,3) (3,4) (3,5) (3,6)   0     0     0
	//   0     0     0   (4,3) (4,4) (4,5) (4,6)   0     0     0
	//   0     0     0   (5,3) (5,4) (5,5) (5,6)   0     0     0 
	//   0     0     0   (6,3) (6,4) (6,5) (6,6) (6,7) (6,8) (6,9)
	//   0     0     0     0     0     0   (7,6) (7,7) (7,8) (7,9)
	//   0     0     0     0     0     0   (8,6) (8,7) (8,8) (8,9)    
	//   0     0     0     0     0     0   (9,6) (9,7) (9,8) (9,9)
	// 
	//

	for (int ki = 0; ki < K; ki++) // �������� �� �������� ���������
	{
		for (int i = 0; i < N; i++) // ����������� �� �������� ������ ������� ��������� ��������
		{
			double scalar_product_ii = 0; // (gi,gi) i � ��������� ���� �� ������
			double scalar_product_fi = 0; // (f,gi)

			for (int j = i+1; j < N; j++) // ������� ��������� ��������������� ��������
			{
				double scalar_product_ij = 0; // ����������
				
				if (i == 0 && ki != 0)  // ���� ������� �������, �� �� ��������� ������������ ����� ��������� 
										// ������ �� ������ ����� �����
				{
					for (int l = 0; l < L; l++)
					{	// G[ki * (N - 1)][l + L] - ������� �������, ������ ����� �����
						scalar_product_ij += G[ki * (N - 1)][l + L] * G[ki * (N - 1) + j][l];
					}
				}
				else // � ��������� ������ ������� �������� ������� �� ������ L ������
				{
					for (int l = 0; l < L; l++)
					{
						scalar_product_ij += G[ki * (N - 1) + i][l] * G[ki * (N - 1) + j][l];
					}
				}
				LAE_(ki * (N - 1) + i, ki * (N - 1) + j) = scalar_product_ij; // ������� ��� �������� � �������
				LAE_(ki * (N - 1) + j, ki * (N - 1) + i) = scalar_product_ij; // ������� ����� ������������
			}

			// ������ ������ ���������
			// ���� �� ������ ������ ��������� ��������, �� ����� ������ ������� ��� �� ����� ���������
			if ((i > 0 && i < N - 1) || (i == 0 && ki == 0) || (i == N - 1 && ki == K - 1)) 
			{// �� ������� ������� �� ������ L ������
				for (int l = 0; l < L; l++)
				{
					scalar_product_ii += G[ki * (N - 1) + i][l] * G[ki * (N - 1) + i][l];
					scalar_product_fi += G[ki * (N - 1) + i][l] * F[ki * (N - 1) + i][l];
				}
			}
			else 
			// ���� �� � ����� ��������� ��������, �� ���������� ���������
			// �� ������� �� 2 * L ������
			if (i == N - 1 && ki != K - 1)
			{
				for (int l = 0; l < L; l++)
				{
					scalar_product_ii += G[ki * (N - 1) + i][l] * G[ki * (N - 1) + i][l]; 
					scalar_product_ii += G[ki * (N - 1) + i][L + l] * G[ki * (N - 1) + i][L + l];
					scalar_product_fi += G[ki * (N - 1) + i][l] * F[ki * (N - 1) + i][l];
					scalar_product_fi += G[ki * (N - 1) + i][L + l] * F[ki * (N - 1) + i][L + l];
				}
			}
			else if (i == 0 && ki != 0) // �� ��� ��������� ��������� ������� ��������� ������� �� ���������� ��������, �� ����������
			{
				continue;
			}
			else

			LAE_(ki * (N - 1) + i, ki * (N - 1) + i) = scalar_product_ii; // ������� �������� � �������
			right_(ki* (N - 1) + i) = scalar_product_fi; // � � ������ ����� ����
		}
	}

	VectorXd c = LAE_.colPivHouseholderQr().solve(right_); // ������ ����

	for (int iviz = 0, ki = 0; iviz < Mviz; iviz++) // ������� ������������� � ������ ��� �����������
	{
		if (X[iviz] > mesh[(ki + 1) * (N - 1)]) // ������ �������� �� ���� ������, ������� ����� ������� �� ������� �������� ���������
		{
			ki++;
		}
		Approx[iviz] = 0; // �������������� ����� �����
		for (int i = 0; i < N; i++) 
		{
			double gviz = 1 / D[i]; // ������� �������� �������� ������� � �����
			for (int j = 0; j < N; j++) 
			{
				if (i != j)
				{
					gviz *= (X[iviz] - mesh[ki * (N - 1) + j]);
				}
			}
			Approx[iviz] += c(ki * (N - 1) + i) * gviz; // � �������� �� ������ ���������
		}

	}

	std::ofstream ParamsFile; // ���� ���������� - � ��� a,b � ��� �����
	ParamsFile.open("Params.txt");
	ParamsFile << a << ", " << b /*<< ", "
		<< resREL_1norm << ", " << resREL_2norm << ", " << resREL_infnorm << ", "
		<< resABS_1norm << ", " << resABS_2norm << ", " << resABS_infnorm */ << std::endl;
	ParamsFile.close();

	std::ofstream meshFile, FmeshFile; // ����� ����� ������������ � �������� ������� � ���� ������ - ��� ��������� �������
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

	std::ofstream XFile, LFile, FFile; // ����� ����� ���������� ��������
	XFile.open("X.txt"); // ����� X
	LFile.open("L.txt"); // �������� �������� �������� � ���� ������
	FFile.open("F.txt"); // �������� ������� � ���� ������
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

	std::system("python plot.py"); // ��� ������� �������� ��������� ������ � �������� ����������� ����� ������
	
	return 0;
}