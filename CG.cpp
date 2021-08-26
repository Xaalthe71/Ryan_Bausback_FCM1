#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<cstdlib>
#include<ctime>
#include<vector>

using namespace std;
//--------------begin forward declarations-----------------------------------
//Matrix vector product
template <typename type>
vector<type> MVProd(vector<vector<type> >, vector<type>, int);
template <typename type>
vector<type> MVProdBanded(vector<vector<type> >, vector<type>, int, int);

//Iterative methods
template<typename type>
vector<type> SteepestDescent(vector<vector<type> >A, vector<vector<type> >M, vector<type> b, vector<type> x0, vector<double>xStar, vector<double>bStar, int n, bool, ofstream&, int);
template<typename type>
vector<type> ConjugateGradient(vector<vector<type> >A, vector<vector<type> >M, vector<type> b, vector<type> x0, vector<double>xStar, vector<double>bStar, int n, bool banded, ofstream&, int);

//System solvers
template<typename  type>
vector<type> Solver(vector<vector<type> >M, vector<type> r, int n);
template<typename type>
vector<type> LLTsystemSolver(vector<vector<type> >A, vector<type> b, int n, bool Upper);
template<typename type>
vector<type> Residual(vector<type>, vector<vector<type> >, vector<type>, int, bool);
double SolutionError(vector<double>, vector<double>, int);

template<typename type>
vector<type> TridiagonalSPDsolver(vector<vector<type> >A, vector<type> b, int n);
template<typename type>
vector<vector<type> >TridiagCholeskyFactor(vector<vector<type> >A, int n);
template<typename type>
vector<type> DiagonalSolver(vector<vector<type> >A, vector<type> b, int n);

template<typename type>
vector<vector<type> >BlockCholesky(vector<vector<type> >M, int n);
template<typename type>
vector<type> BlockSolver(vector<vector<type> >FlatBlocks, vector<type>b, int n);

//Utility methods
vector<vector<float> > Identity(int);
template<typename type>
type TwoNorm(vector<type>, int);
template<typename type>
type vectorProd(vector<type>, vector<type>, int);
template<typename type>
vector<type> vectorSubtract(vector<type> a, vector<type> b, int n);
vector<vector<float> >TwoDvectorTypecast(vector<vector<double> >, int);
vector<double> MatrixVectorMult(vector<vector<double> > A, vector<double> x, int n);
vector<vector<double> > Transpose(vector<vector<double> >A, int n);
template<typename type>
type ConditionNumber(vector<vector<type> >A, int n);
template<typename type>
type Anorm(vector<vector<type> >A, vector<type>a, int n, bool banded);

//MatrixGenerators
vector<vector<double> >SPDgeneratorRand(int n);
template<typename type>
vector<vector<type> >JacobiPCND(vector<vector<type> >A, int n);
template<typename type>
vector<vector<type> >TriDiagPCND(vector<vector<type> >A, int n);
template<typename type>
vector<vector<type> >SgsPCND(vector<vector<type> >A, int n);
template<typename type>
vector<vector<type> >BlockPCND(vector<vector<type> > A, int n);

template<typename type>
vector<vector<type> >Transformer(vector<vector<type> >Abanded, int n);
template<typename type>
vector<vector<type> >CholeskyFactor(vector<vector<type> >A, int n);
//---------------end forward declarations-----------------------------------

int main(){
	srand(time(NULL));

	ofstream myfile100("withPrecond.txt");
	ofstream myfile101("NoPrecond.txt");
	
	int precond=3; //1= diag, 2=sgs, 3=block, 4=tridiag, 0=no preconditioning
	int n =8;
	
for(int p=0; p<100; p++){
	int k=rand()%6+1;  // k element of [0,6]
	//int k =2;
	cout<<"K: "<<k<<endl;
	
	//vector<vector<double> > A = SPDgeneratorRand(n);
	/*vector<vector<double> >A(n, vector<double>(n,0)); //specific test example for correctness
	A[0][0]=25.0;
	A[1][0]=15.0;
	A[2][0]=-5.0;
	A[0][1]=15.0;
	A[1][1]=18.0;
	A[2][1]=0.0;
	A[0][2]=-5.0;
	A[1][2]=0.0;
	A[2][2]=11.0;*/
	
	//vector<vector<double> >L =CholeskyFactor(A, n); //check SPD with cholesky

	vector<vector<double> >Abanded(n, vector<double>(k+1,0));
	//for(int i=0; i<n; i++)
		//Abanded[i][0]=1; //if identity needed
		
/*	for(int i=0; i<n; i++){ //Toeplitz generator
		Abanded[i][0]=2;
		if(i!=0)
			Abanded[i][1]=-1;
	}*/

	for(int i=0; i<n; i++){
		for(int j=0; j<=((i<k)?i:k); j++){
			Abanded[i][j]=rand()%50+(i+1);
			if(j==0){//assures diagonal dominance
				Abanded[i][j]+=4*n*(rand()%50);
			}
		}
	}

/*	int NondistinctCount=0; //for part1: influence of spectrum thm8.5
	vector<double>disteigen(1,0);
	
	for(int i=0; i<n; i++){
		bool distinct=false;
		for(int j=0; j<disteigen.size(); j++){
			if(Abanded[i][0]==disteigen[j]){
				distinct=true;
			}
		}
		if(!distinct)
			disteigen.push_back(Abanded[i][0]);
	}
	cout<<"distinct eigenvalues: "<<disteigen.size()<<endl;*/

	cout<<"num col:"<<Abanded[0].size()<<endl;
	for(int i=0; i<n; i++){
		for(int j=0; j<=k; j++)
			cout<<Abanded[i][j]<<" ";
		cout<<endl;
	}
	cout<<endl;
	
	vector<vector<double> >L = CholeskyFactor(Transformer(Abanded, n),n); //for checking if actually SPD
	if(L[0][0]==-1)
		continue; // if not SPD, stop computation

	myfile100<<k<<" ";
	myfile101<<k<<" ";
	/*cout<<"L: "<<endl;
	for(int i=0; i<n; i++){
		for(int j=0; j<2; j++)
			cout<<L[i][j]<<" ";
		cout<<endl;
	}
	cout<<endl;*/
	
	vector<double> Xstar(n,0);
	cout<<"X: ";
	for(int i =0; i<n; i++){
		Xstar[i]=double(rand()%10)*(i+1);
		//Xstar[i]=double(i+1);
		cout<<Xstar[i]<<" ";
	}
	cout<<endl;
	
	vector<float> xfloatStar(Xstar.begin(), Xstar.end());
	vector<double> bstar = MVProdBanded(Abanded, Xstar, n, k);
	vector<float> b = MVProdBanded(TwoDvectorTypecast(Abanded,n), xfloatStar, n, k);
	//vector<double> bstar =MVProd(A, Xstar, n);
	
	cout<<"b: ";
	for (int i=0; i<n; i++){
		cout<<bstar[i]<<" ";
	}
	cout<<endl;
	
	
	
	vector<float> x0(n,0);
	cout<<"x0: ";
	for(int i=0; i<n; i++){//generate arbitrary x
		x0[i]=double(rand()%10+1);
		if(rand()%9>5)
			x0[i]*=(-1.0);
		cout<<x0[i]<<" ";
	}
	cout<<endl;
	
	vector<vector<float> >M;
	if(precond==0)
		M=Identity(n);
	else if(precond==1)
		M =JacobiPCND(TwoDvectorTypecast(Abanded,n), n);
	else if(precond==2)
		M=SgsPCND(TwoDvectorTypecast(Abanded,n), n);
	else if(precond==3)
		M=BlockPCND(TwoDvectorTypecast(Abanded,n), n);
	else if(precond==4)
		M=TriDiagPCND(TwoDvectorTypecast(Abanded,n), n);
	else
		cout<<"Invalid preconditioning integer"<<endl;
		
	cout<<endl<<"M: "<<endl;
	for(int i=0; i<n; i++){
		for(int j=0; j<2; j++)
			cout<<M[i][j]<<" ";
		cout<<endl;
	}
		
	vector<float> xfloat =SteepestDescent(TwoDvectorTypecast(Abanded,n), M, b, x0, Xstar, bstar, n, true, myfile100, precond);
	vector<float> xfloat2 =SteepestDescent(TwoDvectorTypecast(Abanded,n), M, b, x0, Xstar, bstar, n, true, myfile101, 0); //0 is without preconditioning
	
/*	
	vector<float>x0Star(Xstar.begin(), Xstar.end());
	
	cout<<endl<<"x0star: ";
	for(int i=0; i<n; i++)
		cout<<x0Star[i]<<" ";
	cout<<endl;
	vector<float> bfloat = MVProdBanded(M, x0Star,n, k);
	cout<<"bfloat: ";
	for(int i=0; i<n; i++)
		cout<<bfloat[i]<<" ";
	cout<<endl;
	vector<float> yfloat = LLTsystemSolver(M, bfloat, n, false);
	vector<float> xfloat = LLTsystemSolver(M, yfloat, n, true);
	//vector<double> xfloat =TridiagonalSPDsolver(L, bstar, n);*/
	cout<<endl<<"Solution:"<<endl;
	for(int i=0; i<n; i++)
		cout<<xfloat[i]<<endl;

}
	myfile100.close();
	myfile101.close();
	return 0;
	
}

template<typename type>
vector<type> MVProd(vector<vector<type> > A, vector<type> x, int n){
	vector<type> b(n);
	for(int i=0; i<n; i++){
		//A only stores the diagonal and subdiagonal elements to not repeat storage of aij=aji
		for(int j=0; j<=i; j++){
			b[i]+=A[i][j]*x[j];
		}
		for(int j=i+1; j<n; j++){
			b[i]+=A[j][i]*x[j];
		}
	}
	
	return b;
}

template<typename type> 
vector<type> MVProdBanded(vector<vector<type>  >A, vector<type> x, int n, int k){ //A is NxK with main and subdiagonals as columns so even though 2d array, it only has O(n) or (k+1)n storage
	vector<type> b(n);
	
	for(int i=0; i<n; i++){
		for(int j=k; j>=0; j--){
			b[i]+=A[i][j]*x[i-j];
		}
		
		if(i+1<n){
			for(int j=i+1,l=1; l<=k && j<n; j++,l++){
				b[i]+=A[j][l]*x[j];
			}
		}
	}
	return b;
}


template<typename type>
vector<type> SteepestDescent(vector<vector<type> >A, vector<vector<type> >M, vector<type> b, vector<type> x0, vector<double>xStar, vector<double>bStar, int n, bool banded, ofstream& myfile100, int precond){
	vector<type> r = Residual(b, A, x0,n, banded);
	for(int i=0; i<n; i++)
		cout<<r[i]<<" ";
	cout<<endl<<endl;
	
	vector<type>z;
	vector<vector<type> >L;
	if(precond==0){
		z = Solver(M, r, n);
	}
	else if(precond==1){
		z=DiagonalSolver(M, r, n);
	}
	else if(precond==3){
		L = BlockCholesky(M, n);
		cout<<"FlatBlocks: "<<endl;
		for(int i=0; i<n/2; i++){
			for(int j=0; j<3; j++)
				cout<<L[i][j]<<" ";
			cout<<endl;
		}
		z=  BlockSolver(L, r, n);
	}
	else if(precond==4){
		L = TridiagCholeskyFactor(M, n); //only factors on the first iteration
		z=TridiagonalSPDsolver(L, r, n);
	}
	else if(precond==2){
		vector<type> y=LLTsystemSolver(M, r, n, false);
		z=LLTsystemSolver(M, y, n, true);
	}
	else{
		cout<<"Invalid preconditioning integer."<<endl;
		return x0;
	}

	cout<<"Z: ";
	for(int i=0; i<n; i++)
		cout<<z[i]<<" ";
	cout<<endl<<endl;
	//for output
	int k=0;
	
	stringstream ss;
	ss<<"SDerror"<<precond<<".txt";
	string filename1 =ss.str();
	ofstream myfile2(filename1.c_str());
	
	stringstream ss2;
	ss2<<"SDresid"<<precond<<".txt";
	string filename2 =ss2.str();
	ofstream myfile(filename2.c_str());
	//ofstream myfile3("SDconv1.txt");
	//ofstream myfile4("SDconv2.txt");
	//bool pos=true;
	
	vector<float>xfloatStar(xStar.begin(), xStar.end());
	double errornew=Anorm(A, vectorSubtract(MVProdBanded(A, b, n, A[0].size()-1), xfloatStar,n),n, banded);
	
	while(TwoNorm(r, n)/TwoNorm(b, n)>0.00000001){
		double errorOld=errornew;
		vector<type> w;
		if(banded)
			w=MVProdBanded(A, z, n, A[0].size()-1);
		else
			w = MVProd(A, z, n);

		type alpha = vectorProd(z, r, n)/vectorProd(z, w, n);
		
		for(int i=0; i<n; i++){
			x0[i]+=alpha*z[i];
			r[i]-=alpha*w[i];
		}
		//cout<<"Hello"<<endl;
		if(precond==0){
			z = Solver(M, r, n);
		}
		else if(precond==1){
			z=DiagonalSolver(M,r,n);
		}
		else if(precond==3){
			z=  BlockSolver(L, r, n);
		}
		else if(precond==4){
			z=TridiagonalSPDsolver(L, r, n);
		}
		else if(precond==2){
			vector<type> y=LLTsystemSolver(M, r, n, false);
			z=LLTsystemSolver(M, y, n, true);
		}
		
		k++;
		vector<double> rStar(r.begin(),r.end());
		vector<double> zStar(z.begin(), z.end());
		vector<double>x0Star(x0.begin(), x0.end());

		myfile<<k<<" "<<TwoNorm(rStar, n)/TwoNorm(bStar, n)<<endl;
		myfile2<<k<<" "<<TwoNorm(vectorSubtract(x0Star, xStar, n), n)/TwoNorm(xStar,n)<<endl;
		
		/*//from Part1:
		errornew =Anorm(A, vectorSubtract(x0, xfloatStar,n), n, banded);
		myfile3<<k<<" "<<errorOld - errornew<<endl;
		//cout<<"condition Number:"<<ConditionNumber(A,n)<<endl;
		if(((ConditionNumber(A,n)-1)/(ConditionNumber(A,n)+1) *errorOld - errornew) <0.0 && pos){*/
		if(TwoNorm(r, n)/TwoNorm(b, n)<0.00000001)
			myfile100<<k<<" "<<TwoNorm(vectorSubtract(x0Star, xStar, n), n)/TwoNorm(xStar,n)<<" "<<TwoNorm(rStar, n)/TwoNorm(bStar, n)<<endl;
			//pos=false;
		
		
	}

	myfile.close();
	myfile2.close();
	//myfile3.close();
	//myfile4.close();
	return x0;
}

template<typename type>
vector<type> ConjugateGradient(vector<vector<type> >A, vector<vector<type> >M, vector<type> b, vector<type> x0, vector<double>xStar, vector<double>bStar, int n, bool banded, ofstream& myfile100, int precond){
	vector<type> r = Residual(b, A, x0,n, banded);
	
	vector<type>z;
	vector<vector<type> >L;
	if(precond==0){
		z = Solver(M, r, n);
	}
	else if(precond==1){
		z=DiagonalSolver(M,r,n);
	}
	else if(precond==3){
		L = BlockCholesky(M,n);
		z=  BlockSolver(L, r, n);
	}
	else if(precond==4){
		L = TridiagCholeskyFactor(M, n); //only factors on the first iteration
		z=TridiagonalSPDsolver(L, r, n);
	}
	else if(precond==2){
		vector<type> y=LLTsystemSolver(M, r, n, false);
		z=LLTsystemSolver(M, y, n, true);
		cout<<"y: ";
		for(int i=0; i<n; i++)
			cout<<y[i]<<" ";
		cout<<endl;
		cout<<"z: ";
		for(int i=0; i<n; i++)
			cout<<z[i]<<" ";
		cout<<endl;
	}
	else{
		cout<<"Invalid preconditioning integer."<<endl;
		return x0;
	}

	vector<type> p =z;
	
	//for output
	int k=0;
	
	
	stringstream ss;
	ss<<"CGerror"<<precond<<".txt";
	string filename1 =ss.str();
	ofstream myfile2(filename1.c_str());
	
	stringstream ss2;
	ss2<<"CGresid"<<precond<<".txt";
	string filename2 =ss2.str();
	ofstream myfile(filename2.c_str());
	
	vector<float>xfloatStar(xStar.begin(), xStar.end());
	//double error0A=Anorm(A, vectorSubtract(MVProdBanded(A, b, n, A[0].size()-1), xfloatStar,n),n, banded);
	//double error0Two=TwoNorm(vectorSubtract(MVProdBanded(A, b, n, A[0].size()-1), xfloatStar,n),n);
	
	double error0Anew;
	double error0Twonew;
	
	while(TwoNorm(r, n)/TwoNorm(b, n)>0.00000001){
		vector<type> w;
		if(banded)
			w=MVProdBanded(A, p, n, A[0].size()-1);
		else
			w = MVProd(A, p, n);
			
		type sigma = vectorProd(z, r, n); //compute seperately cuz need value after z and r updated
		
		type alpha = sigma/vectorProd(p, w, n);
		
		for(int i=0; i<n; i++){
			x0[i]+=alpha*p[i];
			r[i]-=alpha*w[i];
		}
		
		if(precond==0){
			z = Solver(M, r, n);
		}
		else if(precond==1){
			z=DiagonalSolver(M,r,n);
		}
		else if(precond==3){
			z=  BlockSolver(L, r, n);
		}
		else if(precond==4){
			z=TridiagonalSPDsolver(L, r, n);
		}
		else if(precond==2){
			vector<type>y =LLTsystemSolver(M, r, n, false);
			z = LLTsystemSolver(M, y, n, true);
		}
		else{
			cout<<"Invalid preconditioning integer."<<endl;
			return x0;
		}
		
		type Beta = vectorProd(z, r, n)/sigma;
		
		for(int i=0; i<n; i++){
			p[i] = z[i] +Beta*p[i];
		}
		
		
		k++;
		vector<double> rStar(r.begin(),r.end());
		vector<double>x0Star(x0.begin(), x0.end());

		myfile<<k<<" "<<TwoNorm(rStar, n)/TwoNorm(bStar, n)<<endl;
		myfile2<<k<<" "<<TwoNorm(vectorSubtract(x0Star, xStar, n), n)/TwoNorm(xStar,n)<<endl;
		
		//error0Anew=Anorm(A, vectorSubtract(x0, xfloatStar, n),n, banded);
		//error0Twonew=TwoNorm(vectorSubtract(x0, xfloatStar, n),n);
		
		//myfile100<<k<<" "<<error0Twonew<<" "<<error0Two*2*sqrt(ConditionNumber(A,n))*pow((sqrt(ConditionNumber(A,n))-1)/(sqrt(ConditionNumber(A,n))+1),k);
		//myfile100<<" "<<error0Anew<<" "<<error0A*2*pow((sqrt(ConditionNumber(A,n))-1)/(sqrt(ConditionNumber(A,n))+1),k)<<endl;
		
		if(TwoNorm(r, n)/TwoNorm(b, n)<0.00000001)
			myfile100<<k<<" "<<TwoNorm(vectorSubtract(x0Star, xStar, n), n)/TwoNorm(xStar,n)<<" "<<TwoNorm(rStar, n)/TwoNorm(bStar, n)<<endl;
	}
	
	//myfile100<<k-numdist<<" "<<error0Two*2*sqrt(ConditionNumber(A,n))*pow((sqrt(ConditionNumber(A,n))-1)/(sqrt(ConditionNumber(A,n))+1),k)-error0Twonew;
	//myfile100<<" "<<error0A*2*pow((sqrt(ConditionNumber(A,n))-1)/(sqrt(ConditionNumber(A,n))+1),k) - error0Anew<<endl;
	
	return x0;
}

//Method-specific functions

template<typename  type>
vector<type> Solver(vector<vector<type> >M, vector<type> r, int n){
	return r;
}

template<typename type>
vector<type> Residual(vector<type> b, vector<vector<type> > A, vector<type> x, int n, bool banded){
	vector<type> r;

	if(banded)
		r = MVProdBanded(A, x, n, A[0].size()-1);
	else
		r = MVProd(A, x, n);
	cout<<"r: ";

	for(int i=0; i<n; i++)
		r[i]=b[i]-r[i];

	for(int i=0; i<n; i++)
		cout<<r[i]<<" ";
	cout<<endl;
	return r;	
}

double SolutionError(vector<double> x, vector<double> xStar, int n){
	vector<double> subtract;
	
	for(int i =0; i<n; i++){
		subtract[i]= x[i]-xStar[i];
	}
	return TwoNorm(subtract, n)/TwoNorm(xStar, n);
}

template<typename type>
vector<type> TridiagonalSPDsolver(vector<vector<type> >L, vector<type> b, int n){
	//Using the diagonals as columns structure for (k+1)*n storage of Abanded
	//Part1: Compute Factorization

	//vector<vector<type> > U(n, vector<type>(2,0)); //column 1 is main diagonal, column 2 is superdiagonal
/*	//old, cholesky factorization now happens in separate method
	U[0][0]=A[0][0];
	for (int i=0; i<n-1; i++){
		U[i+1][1]=A[i+1][1];
		L[i+1]=A[i+1][1]/U[i][0];
		U[i+1][0]=A[i+1][0]-L[i+1]*A[i+1][1];
	}*/
	
	//Part 2: Solve system using factors L and LT
	vector<type> y(n,0);
	y[0]=b[0]/L[0][0];
	for(int i=1; i<n; i++){
		y[i]=(b[i]-L[i][1]*y[i-1])/L[i][0];
	}

	vector<type> x(n,0);
	x[n-1]=y[n-1]/L[n-1][0];
	for(int i=n-2; i>=0; i--){
		x[i]=(y[i]-L[i+1][1]*x[i+1])/L[i][0];
	}
	return x;
}

template<typename type>
vector<vector<type> >TridiagCholeskyFactor(vector<vector<type> >A, int n){
	vector<vector<type> > L(n, vector<type>(2,0)); //O(n) storage
	if(A[0][0]>0)
		L[0][0] = sqrt(A[0][0]);
	else{
		cout<<"Preconditioner Failed to exist at step 0."<<endl;
		return L;
	}
	
	for(int i=1; i<n; i++){ //O(n) to factor since there are only 2 values to compute per iteration for n iterations
		L[i][1]=A[i][1]/L[i-1][0];
		
		if(A[i][0]-L[i][1] >0){
			L[i][0] = sqrt(A[i][0]-L[i][1]);
		}
		else{
			cout<<"Preconditioner Failed to exist at step: "<<i<<endl;
			return L;
		}
	}
	cout<<"factorization completed."<<endl;
	return L;
}

template<typename type>
vector<type> LLTsystemSolver(vector<vector<type> >A, vector<type> b, int n, bool Upper){
	vector<type> x(n);
	
	if(!Upper){
		x[0]= b[0]/A[0][0];
		for(int i=1; i<n; i++){
			x[i]=b[i];
			for(int j=(i<A[0].size())?i:A[0].size()-1; j>0; j--){
				x[i]-=A[i][j]*x[i-j];
			}
			x[i]/=A[i][0];
		}
	}
	else{ 
		x[n-1]= b[n-1]/A[n-1][0];
		for(int i=n-2; i>=0; i--){
			x[i]=b[i];
			for(int j=n-1; j>i; j--){
				if(j-i <A[0].size())
					x[i]-=A[j][j-i]*x[j];
			}
			x[i]/=A[i][0];
		}
	}
	return x;
}

template<typename type>
vector<type> DiagonalSolver(vector<vector<type> >A, vector<type> b, int n){//assume banded structure of A
	vector<type>x(n);
	
	for(int i=0; i<n; i++){
		x[i]=b[i]/A[i][0];
	}
	return x;
}

template<typename type>
vector<vector<type> >BlockCholesky(vector<vector<type> >M, int n){
	if( n%2!=0){
		cout<<"M is not even. Cannot do Block factor."<<endl;
		return M;
	}
	vector<vector<type> > FlatBlocks(n/2+1, vector<type>(3, 0));
	for(int i=0, j=0; i<n; i=i+2, j++){
		FlatBlocks[j][0]=sqrt(M[i][0]);
		FlatBlocks[j][1]=M[i+1][1]/FlatBlocks[j][0];
		FlatBlocks[j][2]=sqrt(M[i+1][0] - pow(FlatBlocks[j][1], 2));
	}
	return FlatBlocks;
}

template<typename type>
vector<type> BlockSolver(vector<vector<type> >FlatBlocks, vector<type>b, int n){
	vector<type>y(n,0);
	//Forward Solve
	for(int i=0, j=0; i<n; i=i+2, j++){
		y[i] = b[i]/FlatBlocks[j][0];
		y[i+1] =(b[i+1] - y[i]*FlatBlocks[j][1])/FlatBlocks[j][2];
	}
	//Backward Solve
	vector<type>x(n,0);
	for(int i=n-1, j=n/2-1; i>=0; i=i-2, j--){
		x[i] = y[i]/FlatBlocks[j][2];
		x[i-1]=(y[i-1] - x[i]*FlatBlocks[j][1])/FlatBlocks[j][0];
	}
	return x;
}

//Utility Methods
vector<vector<float> > Identity(int n){ //pass this function in terms of n for M for CG/SD w/out preconditioning
	vector<vector<float> >I(n, vector<float>(n,0));
	for(int i=0; i<n; i++)
		I[i][i]=1.0;
	return I;
}

template<typename type>
type TwoNorm(vector<type> x, int n){
	type sum =0.0;
	for(int i=0; i<n; i++)
		sum+=x[i]*x[i];
	
	return sqrt(sum);
}

template<typename type>
type vectorProd(vector<type> a, vector<type> b, int n){
	type prod=0;
	for(int i=0; i<n; i++)
		prod+= a[i]*b[i];
	return prod;
}

template<typename type>
vector<type> vectorSubtract(vector<type> a, vector<type> b, int n){
	vector<type> subtract(n,0);

	for(int i =0; i<n; i++){
		subtract[i]= a[i]-b[i];
	}
	
	return subtract;
}

vector<vector<float> >TwoDvectorTypecast(vector<vector<double> >A, int n){
	vector<vector<float> > Afloat(n, vector<float>(A[0].size(),0));
	
	for (int i=0; i<n; i++){
		for(int j=0; j<A[0].size(); j++)
			Afloat[i][j]=float(A[i][j]);
	}
	return Afloat;
}

//Used for testing correctness only
vector<double> MatrixVectorMult(vector<vector<double> > A, vector<double> x, int n){
	vector<double> b(n);
	for(int i=0; i<n; i++){
		for(int j=0; j<A[0].size(); j++){
			b[i]+=A[i][j]*x[j];
		}
	}
	
	return b;
}

vector<vector<double> > Transpose(vector<vector<double> > A, int n){
	vector<vector<double> > L (n, vector<double>(n,0));
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++)
			L[i][j]=A[j][i];
	}
	return L;
}

template<typename type>
type ConditionNumber(vector<vector<type> >A, int n){//Assume A is the diagonal Lambda matrix of only eigenvalues
	type max =0;
	type min =10000;
	
	for(int i=0; i<n; i++){
		if(A[i][0]>max)
			max=A[i][0];
		if(min>A[i][0])
			min=A[i][0];
	}
	return max/min;
}

template<typename type>
type Anorm(vector<vector<type> >A, vector<type>a, int n, bool banded){
	vector<type> temp;
	if(banded)
		temp = MVProdBanded(A, a, n, A[0].size()-1);
	else
		temp = MVProd(A, a, n);
	
	return vectorProd(temp, a, n);
}

//Matrix Generators

vector<vector<double> >SPDgeneratorRand(int n){
	vector<vector<double> > A (n, vector<double>(n,0));
	
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			if(i>=j){
				A[i][j]=double(rand()%10+1)/double(rand()%10+1);
			}
			if(i==j)
				A[i][j]+=n;
			cout<<A[i][j]<<" ";
		}
		cout<<endl;
	}
	return A;
}

template<typename type>
vector<vector<type> >JacobiPCND(vector<vector<type> >A, int n){
	vector<vector<type> >jacobi(n, vector<type>(1, 0));
	
	for(int i=0; i<n; i++){
		jacobi[i][0]=A[i][0];
	}
	return jacobi;
}

template<typename type>
vector<vector<type> >TriDiagPCND(vector<vector<type> >A, int n){
	vector<vector<type> >TriDiag(n, vector<type>(2,0));
	
	for(int i=0; i<n; i++){
		TriDiag[i][0]=A[i][0];
		TriDiag[i][1]=A[i][1];
	}
	return TriDiag;
}

template<typename type>
vector<vector<type> >BlockPCND(vector<vector<type> > A, int n){
	vector<vector<type> >Block(n, vector<type>(2,0));
	
	for(int i=0; i<n; i++){
		Block[i][0]=A[i][0];
		if(i%2!=0){ //only save odd indexes for the 2x2 blocks
			Block[i][1]=A[i][1];
		}
	}
	
	return Block;
}

template<typename type>
vector<vector<type> >SgsPCND(vector<vector<type> >A, int n){
	vector<vector<type> >SGS(n, vector<type>(A[0].size(), 0));
	
	for(int i=0; i<n; i++){
		for(int j=0, k=i; j<A[0].size() && k<n; j++, k++){
			SGS[k][j]=A[k][j]/sqrt(A[i][0]);
		}
	}
	
	/*cout<<"SGS: "<<endl;
	for(int i=0; i<n;i++){
		for(int j=0; j<A[0].size(); j++)
			cout<<SGS[i][j]<<" ";
		cout<<endl;
	}*/
	return SGS;
}

template<typename type>
vector<vector<type> >Transformer(vector<vector<type> >Abanded, int n){ //expand A for use in dense cholesky 
	vector<vector<type> >A(n, vector<type>(n,0));
	for(int i=0; i<Abanded[0].size(); i++){
		for(int j=i; j<n; j++){
			A[j][j-i]=Abanded[j][i];
		}
	}
	
/*	cout<<"A :"<<endl;
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++)
			cout<<A[i][j]<<" ";
		cout<<endl;
	}*/

	return A;
}

template<typename type>
vector<vector<type> >CholeskyFactor(vector<vector<type> >A, int n){
	vector<vector<type> >L(n, vector<type>(n, 0));
	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < (i + 1); j++) {
			type sum = 0;
			for (int k = 0; k < j; k++)
				sum += L[i][k] * L[j][k];

			if (i == j){
				if(A[i][i] - sum >0){
					L[i][j] = sqrt(A[i][i] - sum);
				}
				else{
					cout<<"Cholesky FAILED. A is NOT SPD."<<endl;
					L[0][0]=-1; //use to check later
					return L;
				}	
			}
			else
				L[i][j] = (1.0 / L[j][j] * (A[i][j] - sum));
		}
	}
	return L; 
}
