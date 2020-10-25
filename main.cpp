#include <iostream>
#include <Math.h>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <algorithm>    // std::swap_ranges
#include <vector>       // std::vector
#include <complex>
#include <valarray>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <cstdlib>
#include <iomanip>      // std::setw

using namespace std;



const double PI = 3.141592653589793238460;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;


// Cooley–Tukey FFT (in-place, divide-and-conquer)
// Higher memory requirements and redundancy although more intuitive
void fft(CArray& x)
{
    const size_t N = x.size();
    if (N <= 1) return;

    // divide
    CArray even = x[std::slice(0, N/2, 2)];
    CArray  odd = x[std::slice(1, N/2, 2)];

    // conquer
    fft(even);
    fft(odd);

    // combine
    for (size_t k = 0; k < N/2; ++k)
    {
        Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
        x[k    ] = even[k] + t;
        x[k+N/2] = even[k] - t;
    }
}



int main()
{
	
	ifstream inFile;
	string linha;
	CArray teste;
	int x=0;
	double freq=0.0;

	long int a = 0;
	
	vector<double> x_;
	vector<double> y_;
	vector<double> z_;
	
	
	
	
	inFile.open("1602245833-2715-NAO7856.txt");
	if(inFile.is_open()){
		while(getline(inFile,linha)){
			int temp = linha.size();
			double t;
//			cout << temp << endl;
			string d;
			
			if(linha[0]!='x')
			{
				
				for(int i=0; i<temp; i++)
				{	
					if(linha[i]==',')
					{					
						if(x==0)
						{
							t = stod(d);
							x_.push_back(t);
							cout << a << " " << t << endl;
							x++;
							a = a+1;
							
						}
						else if(x==1)
						{
							t = stod(d);
							y_.push_back(t);
	//						cout << "   " << t;
							x=0;					
						}
						else if(x==2)
						{
							t = stod(d);
							z_.push_back(t);
	//						cout << "   " <<  t << endl;
							x=0;
							break;						
						}
						d = "";
						
					}else if(i==(temp-1)){
						t = stod(d);
						z_.push_back(t);
//						cout << "   " <<  t << endl;
						x=0;
						break;
					}
					//evita valores que nao sejam numéricos
					if(linha[i]!=','){
						
						if(linha[i]==0)
						{
							break;
						}
						d = d + linha[i];
					}
			}	}
		}
	}	
//	cout << a;
////	
//cria as 3 variaveis do tip ocomplexa
	
	ofstream arquivoFFTsaida;
	ofstream arquivoFFTsaida_XYZ;
	//arquivo para debug
//	arquivoFFTsaida.open("FFT OUTPUT");
	arquivoFFTsaida_XYZ.open("FFT X_Y_Z OUTPUT");
	
	Complex eixoX[a];
	Complex eixoY[a];
	Complex eixoZ[a];

	for (int i=0; i<a; i++)
	{
		eixoX[i] = x_[i];
//		cout << i << x_[i] << endl;
	}
	
	CArray dataX(eixoX,a);
	fft(dataX);
	
	for (int i=0; i<a; i++)
	{
		eixoY[i] = y_[i];
//		cout << i << y_[i] << endl;
	}
	
	CArray dataY(eixoY,a);
	fft(dataY);
	
	for (int i=0; i<a; i++)
	{
		eixoZ[i] = z_[i];
//		cout << i << z_[i] << endl;
	}
	
	CArray dataZ(eixoZ,a);
	fft(dataZ);
	
	
//	arquivoFFTsaida << "Complexo FFT eixo X" << setw(30) << "Complexo FFT Eixo Y" << setw(38) << "Complexo FFT Eixo Z" << setw(40)<< "Magnitudes x"<< setw(50)<< "Magnitudes y"<< setw(65)<< "Magnitudes z"<< endl;
	
	std::cout << "fft" << std::endl;
	arquivoFFTsaida_XYZ << "Magnitude (X, Y, Z) e Frequência";
	for (int i = 0; i < a; ++i)
	{
		freq = freq + 1/(2.715);
		
		
		arquivoFFTsaida_XYZ << abs(dataX[i]) <<","<< abs(dataY[i]) << ","<< abs(dataZ[i]) <<","<< abs(dataX[i]) << ","<< freq <<","<< endl;
	    std::cout << dataX[i] << std::endl;
//	    arquivoFFTsaida << dataX[i]<< setw(30) << dataY[i] << setw(38) << dataZ[i] << setw(40)<< abs(dataX[i])<< setw(50)<< abs(dataY[i])<< setw(65)<< abs(dataZ[i])<< endl;
	}
	arquivoFFTsaida.close();

    
    return 0;
}

