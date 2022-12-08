#define _CRT_SECURE_NO_WARNINGS // ��� fopen � visual studio
#include <cmath>
#include <cstdio>
#include <Dense>


using Eigen::MatrixXd;
using Eigen::VectorXd;

double f(double x) 
{
	return sin(x);
}

int main()
{
	double a = -2, b = 2; // ������ � ����� �������

	const int K = 10; // ���������� �������� ���������
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



	// ������ ������������
	// ��� ��� ���� �� ����������� �����, �� �� ������ �������� ��������
	// ����������� ����� �����������, �� ���� �� ���������� ������ N ������������
	double D[N]; 
	for (int i = 0; i < N; i++) // ���������� ������������ �����������
	{
		// mesh[ki * (N - 1) + i] - i-�� ����� �������� ��������� ��������
		D[i] = 1;
		for (int j = 0; j < N; j++)
		{
			if (i != j) {
				D[i] *= (mesh[i] - mesh[j]);
			}
		}
	}


	// ������ ���� ��������� ����� ��� �������� ���������� ������������
	// ��� ��� ��� ���������� ���� ����� ������� ������ ����, �� ��� �������� ������� ��� ����� �����
	double random_mesh[K * L]; 
	srand(time(0)); // �������������� ������������
	for (int ki = 0; ki < K; ki++) // ����������� �� �������� ���������
	{	
		for (int l = 0; l < L; l++) // ��� ������� �������� ������� L �����
		{	
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
	// ������� ����� �� M = K * (N - 1) + 1
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
		
	// ����������� �� �������� ���������
	// ����� ����������� ������
	// N = 4 K = 3
	// 0 0 0 * 0 0 * 0 0 0
	// 0 - ������� ������������ �� ������ ����
	// * - ������� ������������ �� ���� �����
	// �� ���� * ��������� �������
	// � ���� ����� ����� ������� �������� ���� ������� � ��������� ������
	for (int ki = 0; ki < K; ki++) 
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
			else if(i != 0 || ki == 0) // ��� ��������� (i == 0 && ki != 0) - �� �� ������� ��� �������, �� ��������� �� ��� �� ���������� �������
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

	///////////////////////////////////////////////
	// ���������� ������� � ������� ������ ����� //
	///////////////////////////////////////////////

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
			

			LAE_(ki * (N - 1) + i, ki * (N - 1) + i) = scalar_product_ii; // ������� �������� � �������
			right_(ki* (N - 1) + i) = scalar_product_fi; // � � ������ ����� ����
		}
	}

	//////////////////
	// ������� ���� //
	//////////////////

	VectorXd c = LAE_.colPivHouseholderQr().solve(right_); // ������ ����

	/////////////////////////////
	// ���������� ������������ //
	/////////////////////////////
	
	// ���������� ����������� �� ��������� ������
	double er1 = 0, er2 = 0, erinf = 0; 
	// ������������� ���������� = er / fr
	double fr1 = 0, fr2 = 0, frinf = 0;

	for (int il = 0, ki = 0; il < K * L ; il++) // ������� ������������� � ������ ��� �����������
	{
		if (il % L == 0 && il != 0 && il != K * L - 1) // ������ �������� �� ���� ������, ������� ����� ������� �� ������� �������� ���������
		{
			ki++;
		}
		double tmp = 0;
		double ftmp = f(random_mesh[il]);
		for (int i = 0; i < N; i++)
		{
			double g = 1 / D[i]; // ������� �������� �������� ������� � �����
			for (int j = 0; j < N; j++)
			{
				if (i != j)
				{
					g *= (random_mesh[il] - mesh[ki * (N - 1) + j]);
				}
			}
			tmp += c(ki * (N - 1) + i) * g; // � �������� �� ������ ���������
		}
		tmp -= ftmp;
		er1 += abs(tmp);
		er2 += tmp * tmp;
		if (abs(tmp) > erinf) {
			erinf = abs(tmp);
		}

		fr1 += abs(ftmp);
		fr2 += ftmp * ftmp;
		if (abs(ftmp) > frinf) {
			frinf = abs(ftmp);
		}

	}

	er2 = sqrt(er2);
	fr2 = sqrt(fr2);

	// ���������� ����������� �� ����������� ����� � ����� h/100
	double e1 = 0, e2 = 0, einf = 0;
	// ������������� ���������� = e / f
	double f1 = 0, f2 = 0, finf = 0;

	for (int meshi = 0, ki = 0; meshi < M; meshi++) // ����������� �� �������� ���������
	{
		if (meshi > (ki + 1) * (N - 1)) // ������ �������� �� ���� ������, ������� ����� ������� �� ������� �������� ���������
		{
			ki++;
		}
		double h100 = h / 100;
		for (int hi = 0; hi < 100; hi++) 
		{
			double tmp = 0;
			double ftmp = f(mesh[meshi] + hi * h100);
			for (int i = 0; i < N; i++)
			{
				double g = 1 / D[i]; // ������� �������� �������� ������� � �����
				for (int j = 0; j < N; j++)
				{
					if (i != j)
					{
						g *= (mesh[meshi] + hi * h100 - mesh[ki * (N - 1) + j]);
					}
				}
				tmp += c(ki * (N - 1) + i) * g; // � �������� �� ������ ���������
			}

			tmp -= ftmp;
			e1 += abs(tmp);
			e2 += tmp * tmp;
			if (abs(tmp) > einf) {
				einf = abs(tmp);
			}
						
			f1 += abs(ftmp);
			f2 += ftmp * ftmp;
			if (abs(ftmp) > finf) {
				finf = abs(ftmp);
			}
		}
	}




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


	////////////////////
	// ������ � ����� //
	////////////////////


	FILE* ResudialFile; // ���� ������������ - ������� txt ���� � ��������
	ResudialFile = fopen("Resudial.txt", "w");

	fprintf(ResudialFile, "|----------------|----------------|----------------|----------------|\n");
	fprintf(ResudialFile, "|      mesh      |     ||*||1     |     ||*||2     |    ||*||inf    |\n");
	fprintf(ResudialFile, "|----------------|----------------|----------------|----------------|\n");
	fprintf(ResudialFile, "|  h / 100 | abs | %14.8e | %14.8e | %14.8e |\n", e1, e2, einf);
	fprintf(ResudialFile, "|          | rel | %14.8e | %14.8e | %14.8e |\n", e1 / f1, e2 / f2, einf / finf);
	fprintf(ResudialFile, "|----------|-----|----------------|----------------|----------------|\n");
	fprintf(ResudialFile, "|  random  | abs | %14.8e | %14.8e | %14.8e |\n", er1, er2, erinf);
	fprintf(ResudialFile, "|          | rel | %14.8e | %14.8e | %14.8e |\n", er1 / fr1, er2 / fr2, erinf / frinf);
	fprintf(ResudialFile, "|----------|-----|----------------|----------------|----------------|\n");
	fclose(ResudialFile);

	FILE *ParamsFile; // ���� ���������� - � ��� a,b � ��� �����
	ParamsFile = fopen("Params.txt","w");
	fprintf(ParamsFile, "%f, %f", a, b);
	fclose(ParamsFile);

	FILE *meshFile, *FmeshFile; // ����� ����� ������������ � �������� ������� � ���� ������ - ��� ��������� �������
	meshFile = fopen("mesh.txt","w");
	FmeshFile = fopen("Fmesh.txt","w");
	for (int i = 0; i < M - 1; i++) 
	{
		fprintf(meshFile,"%f, ", mesh[i]);
		fprintf(FmeshFile, "%f, ", f(mesh[i]));
	}
	fprintf(meshFile, "%f\n", mesh[M - 1]);
	fprintf(FmeshFile, "%f\n", f(mesh[M - 1]));
	fclose(meshFile);
	fclose(FmeshFile);

	FILE *XFile, *ApproxFile, *FFile; // ����� ����� ���������� ��������
	XFile = fopen("X.txt", "w"); // ����� X
	ApproxFile = fopen("Approx.txt", "w"); // �������� �������� �������� � ���� ������
	FFile = fopen("F.txt","w"); // �������� ������� � ���� ������
	for (int i = 0; i < Mviz - 1; i++) 
	{
		fprintf(XFile, "%f, ", X[i]);
		fprintf(FFile,"%f, ", f(X[i]));
		fprintf(ApproxFile, "%f, ", Approx[i]);
	}
	fprintf(XFile, "%f\n", X[Mviz - 1]);
	fprintf(FFile, "%f\n", f(X[Mviz - 1]));
	fprintf(ApproxFile, "%f\n", Approx[Mviz - 1]);
	fclose(XFile);
	fclose(ApproxFile);
	fclose(FFile);


	FILE *KFile, *KFFile;
	KFile = fopen("K.txt", "w"); // �������� ������� � ���� ������
	KFFile = fopen("KF.txt", "w"); // �������� ������� � ���� ������
	for (int ki = 0; ki < K - 1; ki++) 
	{
		fprintf(KFile, "%f, ", mesh[ki * (N - 1)]);
		fprintf(KFFile, "%f, ", f(mesh[ki * (N - 1)]));
	}
	fprintf(KFile, "%f\n", mesh[M - 1]);
	fprintf(KFFile, "%f\n", f(mesh[M - 1]));
	fclose(KFile);
	fclose(KFFile);


	std::system("python plot.py"); // ��� ������� �������� ��������� ������ � �������� ����������� ����� ������
	std::system("del /s /q Params.txt X.txt F.txt Approx.txt mesh.txt Fmesh.txt K.txt KF.txt"); 
	
	return 0;
}