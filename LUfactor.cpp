#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdlib>
#include<vector>
#include<ctime>
using namespace std;

//------------------------begin forward declarations-------------------------
void LUFactor( vector<vector<double> >&, int, int, int*, int*);

void PermuteMatrix(vector<vector<double> >&, int*, int*, int);
void PermuteVector(vector<double>& , int*,int);

vector<double> Solver(vector<vector<double> >& , vector<double>& , int);

vector<vector<double> > LUmult(vector<vector<double> > , int, bool);
//Norms
double Matrix1Norm(vector<vector<double> >, int);
double MatrixMaxNorm( vector<vector<double> > , int);
double MatrixFNorm (vector<vector<double> >, int);
double Vector1Norm(vector<double>, int);
double VectorMaxNorm(vector<double>, int);
double Vector2Norm(vector<double>, int);

//Matrix Generation
void MatrixFromFile(vector<vector<double> >&, int);
//Task1
vector<vector<double> >IntegerDiagMatrix(int, bool);
vector<vector<double> >IntegerAntiDiagMatrix(int, bool);
vector<vector<double> >XMatrixSum(int, bool);
vector<vector<double> >XMatrixOnes(int);
vector<vector<double> >UnitLowerTriangular(int); //for both 4 and 5
vector<vector<double> >UnitLowerTriangularRand(int n);
vector<vector<double> >Tridiagonal(int);
vector<vector<double> >TridiagonalAntiDom(int n);
vector<vector<double> >AlmostULT(int);
vector<vector<double> >SymmetricPositiveDefinite(int);

vector<vector<double> >DiagDom(int);

vector<vector<double> >Transpose(vector<vector<double> >&, int);
vector<double> MatrixVectorMult(vector<vector<double> >&, int, vector<double>&);
vector<double> VectorSubtract(vector<double>, vector<double>, int);
vector<vector<double> >MatrixSubtract(vector<vector<double> >, vector<vector<double> >, int);
double Average(vector<double>, int);
//------------------------------end forward declarations----------------------------

//Main
int main(){
	srand(time(NULL));
	/*//Timing for O(n^3)
	ofstream myfile5("OrderOfComp.txt");
	for(int n=100; n<1000; n+=50){
		cout<<n<<" ";
		
		vector<vector<double> >A = SymmetricPositiveDefinite(n);
		
		//vector<double> b =MatrixVectorMult(A, n, xOriginal);
		
		int P[n-1]; //row permutations indexes
						for(int i=0; i<n-1; i++)
							P[i]=i; //set equal to indexes so that if no swap is performed, applying it will do nothing. 

						int Q[n-1]; //column permutation indexes
						for(int i=0; i<n-1;i++)
							Q[i]=i;
							
		clock_t start1=clock();
		LUFactor(A,n, 0, P, Q); //0= no pivoting, 1 = partial pivoting, 2=complete pivoting
		clock_t finish1=clock();
						
		myfile5<<n<<" "<<double(finish1-start1)/CLOCKS_PER_SEC<<endl;

	}*/

	ofstream myfile("Resid_Size.txt");
	ofstream myfile2("Fact_Err.txt");
	ofstream myfile3("Growth_Factor.txt");
	ofstream myfile4("Error_in_comp.txt");
	ofstream myfile5("Error_in_comp_hist.txt");

	if(myfile2.is_open()&& myfile2.is_open()&&myfile3.is_open()&&myfile4.is_open()){

		for(int n=5; n<50; n=n+5){
			cout<<"N :"<<n<<endl;
			
			vector<double> ResidErr(1000);
			vector<double> GrowthFact(1000);
			vector<double>FactErr(1000);
			vector<double>CompErr(1000);
			
			for(int w =0; w<1000; w++){
				//vector<vector<double> > A(n, vector<double>(n, 0));

				vector<double> xOriginal(n);
				/*for(int i=0; i<n; i++){
					xOriginal[i]=1;
				}*/
				for(int i=0; i<n; i++){
					xOriginal[i]=double(rand()%100)/double(rand()%100+1);
					if(rand()%2==1)
						xOriginal[i]*=-1;
				}
				vector<vector<double> >A = AlmostULT(n);
				
				/*for(int i=0; i<n; i++){
					for(int j=0; j<n; j++)
						cout<<A[i][j]<<" ";
					cout<<endl;
				}*/

				int P[n-1]; //row permutations indexes
				int Q[n-1]; //column permutation indexes
				
				for(int i=0; i<n-1; i++){
					P[i]=rand()%(n-1);
					Q[i]=rand()%(n-1);
				}	
	
				PermuteMatrix(A, P, Q, n);

				/*for(int i=0; i<n; i++){
					for(int j=0; j<n; j++)
						cout<<A[i][j]<<" ";
					cout<<endl;
				}
				*/
				vector<vector<double> > Aoriginal(n, vector<double>(n, 0));
					//2nd A for correctness checking
				for(int i=0; i<n; i++){
					for(int j=0; j<n; j++){
						Aoriginal[i][j]=A[i][j];
					}
				}

				
				vector<double> b =MatrixVectorMult(A, n, xOriginal);
			

					for(int i=0; i<n-1; i++)
						P[i]=i; //set equal to indexes so that if no swap is performed, applying it will do nothing. 

					for(int i=0; i<n-1;i++)
						Q[i]=i;
					/*
					if(k==0){
						myfile<<"Initial A: "<<endl;
						for(int i=0; i<n; i++){
							for(int j=0; j<n; j++){
								myfile<<A[i][j]<<" ";
								if(j!=n-1)
									myfile<<"&";
							}
							myfile<<"\\\\"<<endl;
						}		
					}
					else{
						for(int i=0; i<n; i++){
							for(int j=0; j<n; j++){
								A[i][j]=Aoriginal[i][j];
							}
						}
					}*/

				
				LUFactor(A,n, 0, P, Q); //0= no pivoting, 1 = partial pivoting, 2=complete pivoting
				PermuteMatrix(Aoriginal, P, Q, n);
				PermuteVector(b, P, n);
				vector<double> x = Solver(A, b, n);				
				vector<double> bComp =MatrixVectorMult(Aoriginal, n, x);

				/*if(w==0){
					cout<<"X:"<<endl;
					for(int i=0; i<n; i++)
						cout<<xOriginal[i]<<" ";
					cout<<endl;
					for(int i=0; i<n; i++){
						cout<<x[i]<<" ";
					}
					
					cout<<endl<<"B:"<<endl;
					for(int i=0; i<n; i++){
						cout<<b[i]<<" ";
					}
					cout<<endl;
					for(int i=0; i<n; i++)
						cout<<bComp[i]<<" ";
				}*/
				//residual

				ResidErr[w] = VectorMaxNorm(VectorSubtract(b,bComp, n),n)/VectorMaxNorm(b,n);
				GrowthFact[w]=MatrixMaxNorm(LUmult(A, n, true), n)/MatrixMaxNorm(Aoriginal, n);
				FactErr[w]=MatrixMaxNorm(MatrixSubtract(Aoriginal,LUmult(A,n, false), n),n)/MatrixMaxNorm(Aoriginal,n);
				CompErr[w]=VectorMaxNorm(VectorSubtract(xOriginal, x,n),n)/VectorMaxNorm(xOriginal,n);

					/*
					myfile<<"after lu: "<<(k==0?"no pivoting":k==1?"partial pivoting":"complete pivoting")<<endl;
					for(int i=0; i<n; i++){
						for(int j=0; j<n; j++){
							myfile<<A[i][j]<<" ";
							if(j!=n-1)
								myfile<<"&"; //output directly to latex format
						}
						myfile<<"\\\\"<<endl;
					}*/
					/*
					myfile<<"permutations: "<<endl;
					myfile<<"row-wise (P): ";
					for(int i=0; i<n-1; i++)
						myfile<<P[i]<<" ";
					myfile<<endl<<"column-wise (Q): ";
					for(int i=0; i<n-1; i++)
						myfile<<Q[i]<<" ";*/
					
					//myfile<<endl<<endl<<"growth factors:"<<endl;
					//myfile<<"1-norm: "<<Matrix1Norm(LUmult(A, n, true), n)/Matrix1Norm(Aoriginal, n)<<endl;
					//myfile<<"Max-norm: "<<MatrixMaxNorm(LUmult(A, n, true), n)/MatrixMaxNorm(Aoriginal, n)<<endl;
					//myfile<<"F-norm: "<<MatrixFNorm(LUmult(A, n, true), n)/MatrixFNorm(Aoriginal, n)<<endl;
					//myfile<<endl;
				
				/*vector<vector<double> > M= LUmult(A, n, false);
				cout<<endl<<"final: "<<endl;
				for(int i=0; i<n; i++){
					for(int j=0; j<n; j++)
						cout<<M[i][j]<<" ";
					cout<<endl;
				}
				
				//For outputting the permutations
				/*
				for(int i=0; i<n-1; i++)
					cout<<P[i]<<" ";
				cout<<endl;
				
				for(int i=0; i<n-1; i++)
					cout<<Q[i]<<" ";
				cout<<endl;*/
				/*cout<<"original"<<endl;
				for(int i=0; i<n; i++){
					for(int j=0; j<n; j++)
						cout<<Aoriginal[i][j]<<" ";
					cout<<endl;
				}
				
				//PermuteMatrix(Aoriginal, P, Q, n);
				
				/*PermuteVector(b, P, n);
				for(int i=0; i<n; i++)
					cout<<b[i]<<endl;
				
				vector<double> x = Solver( A, b, n);
				
				cout<<"solution:"<<endl;
				for(int i=0; i<n; i++)
					cout<<x[i]<<endl;
				*/
			}

			
			myfile<<n<<" "<<Average(ResidErr, 1000)<<endl;
			myfile2<<n<<" "<<Average(FactErr, 1000)<<endl;
			myfile3<<n<<" "<<Average(GrowthFact, 1000)<<endl;
			myfile4<<n<<" "<<Average(CompErr, 1000)<<endl;
			
			if(n==45){
				for(int i=0; i<1000; i++){
					myfile5<<CompErr[i]<<endl;
				}
			}
		}
	}
	else cout<<"Failed to open file"<<endl;
	return 0;
}

void LUFactor( vector<vector<double> >&A, int n, int PivotFlag, int* P, int* Q){//P  and Q should be int array, size n-1
	if(PivotFlag!=0 && PivotFlag!=1 && PivotFlag!=2){ //0= no pivoting, 1 = partial pivoting, 2=complete pivoting
		cout<<"Invalid choice of PivotFlag. No pivoting will commence."<<endl;
		PivotFlag=0;
	}
	
	if(PivotFlag==0){
		for(int j=0; j<n-1; j++){
			if(A[j][j]==0){//check for cases when factorization fails. 
				cout<<"pivot element is 0. failed to complete LU factor. "<<endl;
				return;
			}
			
			for(int i=j+1; i<n; i++){//row iterator for active matrix
				A[i][j] /= A[j][j];
				
				for(int k=j+1;k<n;k++){//column iterator for active part of matrix
					A[i][k]= A[j][k]*A[i][j]*(-1.0) +A[i][k];
				}
			}
		}
	}
	
	if(PivotFlag==1){//column-based partial pivoting 
		for(int j=0; j<n-1; j++){
			//find max magnitude in column j
			double max=0.0;
			for(int i=j; i<n; i++){//iterate active matrix only to find swaps
				if(abs(A[i][j])>max){
					P[j]=i;
					max=abs(A[i][j]);
				}
			}
			if(max==0){
				cout<<"All partial pivots 0 at step "<<j<<". failed to compute LU"<<endl;
				return;
			}
			if(max<pow(10, -12))
				cout<<"WARNING: Partial pivots at step "<<j<<" are very small."<<endl;
			
			//Perform row swaps
			for(int i=0; i<n; i++){//column iterator
				double temp = A[j][i];
				A[j][i]= A[P[j]][i];
				A[P[j]][i]=temp;
			}

			//Apply Gauss Transform in-place
			for(int i=j+1; i<n; i++){//row iterator for active matrix
				A[i][j] /= A[j][j];
				
				for(int k=j+1;k<n;k++){//column iterator for active part of matrix
					A[i][k]= A[j][k]*A[i][j]*(-1.0) +A[i][k];
				}
			}
			/*//used with correctness checking
			cout<<"step "<<j<<endl;
			for(int i=0; i<n; i++){
				for(int k=0; k<n; k++)
					cout<<A[i][k]<<" ";
				cout<<endl;
			}*/
		}
	}
	
	if(PivotFlag==2){//complete pivoting
		for(int j=0; j<n-1; j++){
			//find max magnitude in entire active submatrix
			double max=0.0;
			for(int i=j; i<n; i++){//iterate active matrix only to find swaps
				for(int k=j; k<n; k++){
					if(abs(A[i][k])>max){
						P[j]=i;
						Q[j]=k;
						max=abs(A[i][k]);
					}
				}
			}
			if(max==0){
				cout<<"All partial pivots 0 at step "<<j<<". failed to compute LU"<<endl;
				return;
			}
			if(max<pow(10, -12))
				cout<<"WARNING: Partial pivots at step "<<j<<" are very small."<<endl;
			
			//Perform row swaps
			for(int i=0; i<n; i++){//column iterator
				double temp = A[j][i];
				A[j][i]= A[P[j]][i];
				A[P[j]][i]=temp;
			}
			//Perform column swaps
			for(int k=0; k<n;k++){
				double temp =A[k][j];
				A[k][j]=A[k][Q[j]];
				A[k][Q[j]]=temp;
			}

			//Apply Gauss Transform in-place
			for(int i=j+1; i<n; i++){//row iterator for active matrix
				A[i][j] /= A[j][j];
				
				for(int k=j+1;k<n;k++){//column iterator for active part of matrix
					A[i][k]= A[j][k]*A[i][j]*(-1.0) +A[i][k];
				}
			}
		}
	}

}

//Apply Permutations to matrix/vector
void PermuteMatrix(vector<vector<double> >& A, int* P, int* Q, int n){
	for(int i=0; i<n-1; i++){
		//Apply row permutations according to index in P
		for(int j=0; j<n; j++){
			double temp = A[i][j];
			A[i][j]=A[P[i]][j];
			A[P[i]][j]=temp;
		}
		for(int j=0; j<n; j++){
			double temp = A[j][i];
			A[j][i]=A[j][Q[i]];
			A[j][Q[i]]=temp;
		}
	}
}

void PermuteVector(vector<double>& b, int* P, int n){
	for(int i=0; i<n-1; i++){
		double temp= b[i];
			b[i]=b[P[i]];
			b[P[i]]=temp;
	}
}

//Performs forward and backward substitution to solve system
vector<double> Solver(vector<vector<double> >& LU, vector<double>& b, int n){ //assume b is already properly permuted. 
	//forward substitution.
	vector<double> y(n);
	y[0]= b[0];
	for(int i=1; i<n; i++){
		y[i]=b[i];
		for(int j=0; j<i; j++){
			y[i]-=LU[i][j]*y[j];
		}
	}
	/*
	cout<<"y"<<endl;
	for(int i=0; i<n; i++)
		cout<<y[i]<<endl;
	*/
	//backward substitution
	vector<double> x(n);
	x[n-1]=y[n-1]/LU[n-1][n-1];
	for(int i=n-2; i>-1; i--){
		x[i]=y[i];
		for(int j=i+1; j<n;j++){
			x[i]-=LU[i][j]*x[j];
		}
		x[i]/=LU[i][i];
	}
	
	return x;
}


vector<vector<double> > LUmult(vector<vector<double> > LU, int n, bool absVal){
	if(absVal){
		for(int i=0; i<n; i++){
			for(int j=0; j<n; j++)
				LU[i][j]=abs(LU[i][j]);
		}
	}
	
	vector<vector<double> >M(n, vector<double>(n,0));
	
	for(int j=0; j<n; j++){
		for(int i=0; i<n; i++){
			int min=i;
			if(j<i)
				min=j;
			for(int k=0; k<=min; k++){
				if(k==i)
					M[i][j]+=LU[i][j];
				else
					M[i][j]+=LU[k][j]*LU[i][k];
			}
		}
	}
	
	return M;
}

//------------------Norms-----------------
//Matrix Norms
double Matrix1Norm( vector<vector<double> > A , int n){
	double max =0.0;
	for(int i =0; i< n; i++){
		double sum =0.0;
		for(int j=0; j<n; j++){
			sum+=abs(A[j][i]);
		}
		if(sum>max){
			max=sum;
		}
	}
	return max;
}

double MatrixMaxNorm( vector<vector<double> > A , int n){
	double max =0.0;
	for(int i =0; i< n; i++){
		double sum =0.0;
		for(int j=0; j<n; j++){
			sum+=abs(A[i][j]);
		}
		if(sum>max){
			max=sum;
		}
	}
	return max;
}

double MatrixFNorm (vector<vector<double> > A, int n){
	double sum =0.0;
	for(int i=0; i< n; i++){
		for(int j=0; j<n; j++){
			sum+= abs(A[i][j])*abs(A[i][j]);
		}
	}
	return sqrt(sum);
}

//Vector Norms
double Vector1Norm(vector<double> x, int n){
	double sum = 0.0;
	for(int i=0; i<n; i++)
		sum+=abs(x[i]);
	return sum;
}

double VectorMaxNorm(vector<double> x, int n){
	double max =0.0;
	for(int i=0; i<n; i++){
		if(abs(x[i])>max)
			max=abs(x[i]);
	}
	return max;
}

double Vector2Norm(vector<double> x, int n){
	double sum =0.0;
	for(int i=0; i<n; i++)
		sum+=abs(x[i])*abs(x[i]);
	return sqrt(sum);
}


//-------------Matrix Generation-----------------------
void MatrixFromFile(vector<vector<double> >& A, int n){
	ifstream infile("in.txt");
	if(infile.is_open()){
		for(int i=0; i<n; i++){
			for(int j=0; j<n; j++)
				infile>>A[i][j];
		}
	}
	else cout<< "Unable to open file."<<endl;
	infile.close();
}

//Task1
vector<vector<double> >IntegerDiagMatrix(int n, bool reverse){
	vector<vector<double> > A(n, vector<double>(n,0));
	if(!reverse){
		for(int i=0; i<n; i++){
			A[i][i]=double(i)+1.0;
		}
	}
	else{
		for(int i=0; i<n; i++){
			A[i][i]=double(n-i);
		}
	}
	return A;
}

vector<vector<double> >IntegerAntiDiagMatrix(int n, bool reverse){
		vector<vector<double> > A(n, vector<double>(n,0));
	if(!reverse){
		for(int i=0; i<n; i++){
			A[i][n-i-1]=double(i)+1.0;
		}
	}
	else{
		for(int i=0; i<n; i++){
			A[i][n-i-1]=double(n-i);
		}
	}
	return A;
}

vector<vector<double> >XMatrixSum(int n, bool reverse){
	vector<vector<double> > A1 =IntegerDiagMatrix(n, reverse);
	vector<vector<double> >A2 =IntegerAntiDiagMatrix(n, !reverse);
	
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			A1[i][j]+=A2[i][j];
		}
	}
	return A1;
}

vector<vector<double> >XMatrixOnes(int n){
	vector<vector<double> > A (n, vector<double>(n,0));
	
	for(int i=0; i<n; i++){
		A[i][n-i-1]=double(rand()%9+1);
		A[i][i]=double(rand()%9+1);
	}
	
	return A;
}

vector<vector<double> >UnitLowerTriangular(int n){
	vector<vector<double> > A (n, vector<double>(n,0));
	
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			if(i>j){
				A[i][j]=double(i+2-j);
			}
			else if (i==j)
				A[i][j]=2;
			else{}	
		}
	}
	return A;
}

vector<vector<double> >UnitLowerTriangularRand(int n){
	vector<vector<double> > A (n, vector<double>(n,0));
	
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			if(i>=j){
				A[i][j]=double(rand()%10+1)/double(rand()%10+1);
			}
			else{}	
		}
	}
	return A;
}

vector<vector<double> >Tridiagonal(int n){
	vector<vector<double> > A (n, vector<double>(n,0));
	
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			if(j+1==i || i+1==j||(i==j&&i==0)){
				A[i][j]=double(rand()%9+1);
				if(i==j&&i==0)
					A[i][j]*=8.0;
			}
			if(i==j&&i!=0)
				A[i][j]=A[i][j-1]*8.0;
		}
	}
	return A;
}

vector<vector<double> >TridiagonalAntiDom(int n){
	vector<vector<double> > A (n, vector<double>(n,0));
	
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			if(j+1==i||(i==j&&i==0)){
				A[i][j]=double(rand()%100+1);
				//if(i==j&&i==0)
					//A[i][j]/=8.0;
			}
			//if(i==j&&i!=0)
				//A[i][j]=A[i][j-1]/8.0;
			if(i+1==j)
				A[i][j]=double(rand()%100+1)*20;
		}
	}
	return A;
}

vector<vector<double> >AlmostULT(int n){
	vector<vector<double> > A (n, vector<double>(n,0));
	
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			if(i>j){
				A[i][j]=-1.0;
			}
			else if (i==j)
				A[i][j]=1.0;
			else if(j==n-1)
				A[i][j]=1.0;
			else{}	
		}
	}
	
	return A;
}
vector<vector<double> > Transpose(vector<vector<double> >& A, int n){
	vector<vector<double> > L (n, vector<double>(n,0));
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++)
			L[i][j]=A[j][i];
	}
	return L;
}
vector<vector<double> >SymmetricPositiveDefinite(int n){
	vector<vector<double> > L = UnitLowerTriangularRand(n);
	/*
	cout<<"L:"<<endl;
	for(int i=0; i<n; i++){
		for(int k=0; k<n; k++)
			cout<<L[i][k]<<" ";
		cout<<endl;
	}*/
	
	vector<vector<double> >A=Transpose(L, n);

	vector<vector<double> > B(n, vector<double>(n,0));
	for(int k=0;k<n; k++){
		for(int i=0; i<n; i++){
			for(int j=0; j<n; j++){
				B[k][i]+=L[k][j]*A[j][i];
			}
		}
	}

	return B;
}

vector<vector<double> >DiagDom(int n){
	vector<vector<double> > A (n, vector<double>(n,0));
	vector<double>mag(n);
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			A[i][j]=double(rand()%100+1);
			if(rand()%2==1)
				A[i][j]*=-1;
			mag[i]+=abs(A[i][j]);
		}
	}
	for(int i=0; i<n; i++){
		if(A[i][i]<=mag[i])//ensures diagonal dominance by rows
			A[i][i]+=mag[i];
	}
	return A;
}

//utility methods
vector<double> MatrixVectorMult(vector<vector<double> >& A, int n, vector<double>& x){
	vector<double> b(n);
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			b[i]+=A[i][j]*x[j];
		}
	}
	
	return b;
}

vector<double> VectorSubtract(vector<double> a, vector<double>b, int n){
	vector<double> output(n);
	for(int i=0; i<n; i++){
		output[i]=a[i]-b[i];
	}
	return output;
}

vector<vector<double> >MatrixSubtract(vector<vector<double> >A, vector<vector<double> >B, int n){
	vector<vector<double> > output (n, vector<double>(n,0));
	
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			output[i][j]=A[i][j]-B[i][j];
		}
	}
	return output;
}

double Average(vector<double> A, int n){
	double sum=0;
	for(int i=0; i<n; i++){
		sum+=A[i];
	}
	return sum/double(n);
}