#include<iostream>
#include<fstream>
#include<complex>
#include<cmath>
#include<valarray>
#include<cstdlib>

//#include "testfunc.h"
#include <time.h>

using namespace std;

typedef complex<double> dcomp;

//forward declarations
void Gammas(dcomp*, int, bool);
void FFT(dcomp*, int, int, bool);
void DFT(dcomp*, dcomp*, int, bool);
dcomp L2Norm(dcomp*, int);
dcomp ComplexRand();
dcomp Matrix2NormCirculant(dcomp*, dcomp*, int);
void circulantEigenValues(dcomp*, int);
void circulantSystemSolver(dcomp*, dcomp*, int);

dcomp HornerRule(dcomp*,int, dcomp, int);

int main(){
	srand(time(NULL));
	
	//MY_complex::complex ans =testfunc::f1(1.0);
	//cout<<ans.real<<"+ i"<<ans.imag<<endl;
	//define constants.
	dcomp i;
	i=-1;
	i=sqrt(i);
	double pi= 2* asin(1);
	
	clock_t start1, finish1, start2, finish2;
	
	bool IDFT = false; //change for IDFT v DFT
	
	int n =64;
	
	dcomp d[n];
	ifstream infile("in.txt");
	if(infile.is_open()){
		for(int i=0; i<n; i++){
			infile>>d[i];
		}
	}
	
	//int n = 4; //change size for dif d;
	//dcomp d[n];
	/*d[0]=8.0; //deterministic test values
	d[1]=4.0;
	d[2]=7.0;
	d[3]=0.0;
	d[4]=9.0;
	d[5]=3.0;
	d[6]=10.0;
	d[7]=-1.0;
	d[8]=8.0;
	d[9]=4.6732-4.0*i;
	d[10]=8.0;
	d[11]=2.0*i;
	d[12]=8.0;
	d[13]=4.0;
	d[14]=8.0;
	d[15]=3.5;*/
	ofstream myfile2("DemoTest.txt");
	if(myfile2.is_open()){
		//for(int n =4; n<100; n=2*n){
			double max =0.0;
			//for(int k=0; k<n; k++){
				//dcomp originD[n];
				//dcomp d[n];
				//for(int j=0; j<n; j++){
					//d[j]=ComplexRand();
				//}
				
				dcomp a[4]; //circulant matrix test values
				a[0]=1;
				a[1]=4;
				a[2]=3;
				a[3]=2;
				//for(int j=0; j<n; j++)
					//originD[j]=d[j];
				//dcomp OriginalNorm = L2Norm(d,n);
			
				dcomp s[n];
				dcomp l[n];
				//for(int j=0; j<n; j++)
					//s[j]=0.0;
				
				//dcomp max = Matrix2NormCirculant(a, s, n);
				//cout<<"matrix 2norm: "<<max<<endl;
				//circulantEigenValues(a,n);
				//circulantSystemSolver(d, s, n);
				//for(int j=0; j<n; j++)
					//myfile2<<a[j]<<endl;
				
				/*
				dcomp gamma[n]; //for test unitary property.
				dcomp oppGamma[n];
				Gammas(gamma, n, false);
				Gammas(oppGamma, n, true);
				
				for(int j=0; j<n; j++)
					gamma[j]=pow(gamma[j],k)/sqrt(double(n));*/
				
				start1 =clock(); //for timing dft routine.
				//for(int j=0; j<1000; j++){
					DFT(d, s, n, false);
					//DFT(s , l, n, true); 
				//}
				finish1=clock();
			
				
				//FFT
				
				start2=clock(); //for timing fft routine.
				//for(int j=0; j<1000; j++){
					//FFT(d, n, n, false);
				//}
				finish2=clock();
				
				for(int j=0; j<n; j++){
					if(abs(d[j]-l[j])>max)
						max=abs(d[j]-l[j]);
				}
				/*
				for(int j=0; j<n;j++){
					if(j==k && abs(s[j]-1.0)>max)
						max=abs(s[j]-1.0);
					else{
						if(j!=k && abs(s[j])>max)
							max=abs(s[j]);
					}
				}
				
				//FFT(d, s2, n, 1, true);
				/*for(int j=0; j<n; j++)
					myfile2<<oppGamma[j]<<endl;
				myfile2<<endl;
				
				for(int j=0; j<n; j++)
					myfile2<<s[j]<<endl;
				myfile2<<endl;*/
			//}
			
			for(int k=0; k<n;k++)
				s[k]=s[k]/sqrt(double(n));
			
			//for(int k=0; k< n; k++){
				//myfile2<<d[k]<<endl;
				cout<<" "<<HornerRule(s,n, 2.5*M_PI/64.0, 1)<<endl;
				//cout<<HornerRule(d, n, 1, 1)<<endl;
			//}		
				
			//myfile2<<n<<" "<<double(finish1-start1)/CLOCKS_PER_SEC<<endl;
			//myfile2<<"DFT eval time: "<<double(finish1-start1)/CLOCKS_PER_SEC<<endl;
			//myfile2<<"FFT eval time: "<<double(finish2 - start2)/CLOCKS_PER_SEC<<endl;
		//}
		myfile2.close();
	}
	else cout<<"Unable to open file.";
	/*
	//output results to 
	ofstream myfile ("dftout1.txt");
	if (myfile.is_open()){
		for(int j=0; j<n; j++)
			myfile<<originD[j]<<endl;
		myfile.close();
	}
	else cout<<"Unable to open file";
	
	ofstream myfile1 ("fftout.txt");
	if (myfile1.is_open()){
		for(int j=0; j<n; j++)
			myfile1<<d[j]<<endl;
		myfile1.close();
	}
	else cout<<"Unable to open file";*/
	
	return 0;
}


void DFT(dcomp* d, dcomp* s, int n, bool IDFT){
	dcomp i =-1;
	i =sqrt(i);
	double neg=1.0;
	if(IDFT) //swap between IDFT and DFT
		neg =-1.0;
	for(int k =n-1; k>-1; k--){
		dcomp tau = d[n-1];
		for(int j =n-2; j>-1; j--)
			tau = d[j] +exp(M_PI*neg*-2.0* i*double(k)/double(n))*tau;
		s[k]=tau/sqrt(double(n));
	}
}

void FFT(dcomp* x, int n, int originSize, bool IFFT){
	dcomp i =-1;
	i =sqrt(i);
	double neg =1.0;
	if(IFFT)
		neg=-1.0;
	if(n==1){//return once size 1 since then sufficiently broken up
		return;
	}
	else{
		int m = n/2;
		dcomp i =-1;
		i =sqrt(i);
		dcomp even [m];
		dcomp odd[m];
		
		for(int i=0; i<m; i++){ //split based on odd and even indicies
			even[i]=x[i*2];
			odd[i]=x[i*2+1];
		}
		FFT(even, m, n, IFFT);
		FFT(odd, m, n, IFFT);

		for( int j=0; j< m; j++){
			
			x[j] = even[j] +exp(M_PI*neg*-2.0* i*double(j)/double(n))*odd[j];
			x[j+m]=even[j]-exp(M_PI*neg*-2.0* i*double(j)/double(n))*odd[j];
			
			if(originSize==n){ //need to perform the scaling at the end.
				x[j] = x[j]/sqrt(double(n));
				x[j+m] = x[j+m]/sqrt(double(n));
			}
		}
	}
	return;
}

dcomp L2Norm(dcomp* input, int n){
	dcomp sum =0.0;
	
	for(int i=0; i<n; i++){
		sum+=pow(abs(input[i]),2.0);
	}
	
	return sqrt(sum);
	
}

dcomp ComplexRand(){ //produce complex random numbers than are negative and positive randomly
	dcomp i;
	i=-1;
	i=sqrt(i);
	double neg=1.0;
	
	if(rand()%10>5) //hopefully 50% pos, 50% neg on average
		neg=-1.0;
	dcomp raNd = double(rand());
	raNd /= double(rand());
	
	if(rand()%10<5)
		neg=1.0;
	dcomp raNd2 =double(rand())*i;
	raNd2 /= (double(rand()));
	return raNd +raNd2;
}

dcomp Matrix2NormCirculant(dcomp* celements, dcomp* output, int n){ //output should be an array of all zeros, size n.
	dcomp conjCElements [n];
	for(int i=0; i<n; i++)//compute conjugates of circulant matrix elements
		conjCElements[i]=conj(celements[i]);
	
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			if(j-i<0){//if less than lower bound, circle back around to the top. 
				output[j-i+n]+= celements[i]*conjCElements[j];
			}
			else
				output[j-i]+=celements[i]*conjCElements[j];
		}
	}
	
	circulantEigenValues(output, n);
	dcomp max =0.0;
	for(int j=0; j<n; j++){
		if(abs(output[j])>abs(max)){
			max=output[j];
		}
	}
	return sqrt(max);
}

void circulantEigenValues(dcomp* a, int n){ //a holds distinct elements of circulant matrix C_n
	FFT(a, n, n, true);
	for(int i=0; i<n; i++)
		a[i]=sqrt(n)*a[i];
	return;
}

void circulantSystemSolver(dcomp* a, dcomp* b, int n){
	circulantEigenValues(a, n);
	for(int j=0; j<n; j++){//check sufficiently nonsingular
		if(a[j]==0.0){
			cout<<"matrix is singular. Cannot solve system."<<endl;
			return;
		}
	}
	
	FFT(b, n, n, true);
	for(int i=0; i<n;i++){
		b[i]=b[i]/a[i];
	}
	FFT(b,n,n,false);
}

void Gammas(dcomp* gamma, int n, bool IDFT){//for use when checking unitary property. 
	dcomp i;
	i=-1;
	i=sqrt(i);
	
	if(!IDFT){
		for(int j=0; j<n; j++)
			gamma[j]= exp(-2.0 *i * M_PI* double(j)/double(n));
	}
	else{
		for(int j=0; j<n; j++)
			gamma[j]= exp(2.0 *i * M_PI*double(j)/double(n));
	}
	return;
}

dcomp HornerRule(dcomp* d,int n, dcomp x, int k){
	dcomp i;
	i=-1;
	i=sqrt(i);
	
	dcomp tau = d[n-1];
	for(int j =n-2; j>-1; j--)
		tau = d[j] +exp(i*x)*tau;
			
	return tau;
}
