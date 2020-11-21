import java.lang.Math;
import java.util.*;
import java.io.*;

public class hw2{
	
	public static void main(String[] args)throws Exception{
	
		int n = 21;
		int meshChoice = 2;
		int ordering = 1;
		
		double[] xi = new double[n+1];
		int rightmost = 22;
		int leftmost = 0;
		double h = ((double)(rightmost-leftmost))/n;
		
		if(meshChoice ==1){
			for(int i =0; i< n+1; i++)
				xi[i]= leftmost + i*h;
		}
		else if(meshChoice==2){	
			for(int i=0; i<=n; i++)
				xi[i] = ((leftmost-rightmost)/2)*(Math.cos(((2*i+1)*Math.PI)/(2*n+2))+1)+leftmost; //change of variables to be altered------------------
		} 
		else{
			for(int i=0; i<=n; i++)
				xi[i]=Math.cos((Math.PI*i)/n);
		}
		
		xi=LejaOrdering(ordering, xi);
		
		double[] xValues = new double[1000];    //generate the z's to evaluate the function at.
		double h1 =((double)(rightmost-leftmost))/1000;

		for(int i=1; i<1000; i++){
			xValues[i]=leftmost+i*h1;
		}

		double[][] DD = DividedDifferences(n, xi, false);
		double[] PnX1= HornerRule(DD, xi, xValues);


		double[][]result = Barycentric2part1(meshChoice, n);
		double[] PnX2 = Barycentric2evaluation(result, xValues, xi);
	
		double[] FX = new double[xValues.length];
		for(int i=0; i<xValues.length; i++)
			FX[i]=function2(xValues[i], n);				//compute F(z)---------------------------------

//Accuracy	
		FileWriter fw = new FileWriter("LagrangeRealvTrueValues.txt");
		for(int i=0;i<xValues.length;i++)
			fw.write(xValues[i]+" "+PnX2[i]+" "+FX[i]+"\n");
		fw.close();
		
		
		FileWriter fw2 = new FileWriter("NewtonRealvTrueValues.txt");
		for(int i=0; i<xValues.length; i++)
			fw2.write(xValues[i]+" "+PnX1[i]+" "+FX[i]+"\n");
		fw2.close();
		
		
		FileWriter fw3 = new FileWriter("PointwiseErrorAccuracyNewton.txt");
		for(int i=0; i<xValues.length; i++)
			fw3.write(xValues[i]+" "+RelErrorAcc(FX[i], PnX1[i])+"\n");
		fw3.close();
		
		FileWriter fw4 = new FileWriter("PointwiseErrorAccuracyLagrange.txt");
		for(int i=0; i<xValues.length; i++)
			fw4.write(xValues[i]+" "+RelErrorAcc(FX[i], PnX2[i])+"\n");
		fw4.close();
		
//Stability
		//repeat everything in single precision
		float[] xiSingle = new float[n+1];
		float hSingle = ((float)(rightmost-leftmost))/n;
		
		if(meshChoice ==1){
			for(int i =0; i< n+1; i++)
				xiSingle[i]= leftmost + i*hSingle;
		}
		else if(meshChoice==2){
			for(int i=0; i<=n; i++)
				xiSingle[i] = (float)Math.cos(((2*i+1)*Math.PI)/(2*n+2)); //change of variables to be altered------------------
		}
		else{
			for(int i=0; i<=n; i++)
				xiSingle[i]=(float)Math.cos((Math.PI*i)/n);
		}
		
		xiSingle = LejaOrderingSingle(ordering, xiSingle);
		
		float[] xValuesSingle = new float[1000];    //generate the z's to evaluate the function at.
		float h1Single =((float)(rightmost-leftmost))/1000;

		for(int i=1; i<1000; i++){
			xValuesSingle[i]=leftmost+i*h1Single;
		}
		
		float[][] DDsingle = DividedDifferencesSingle(n, xiSingle, false);
		float[] PnX1single= HornerRuleSingle(DDsingle, xiSingle, xValuesSingle);
		
		float[][]resultSingle = Barycentric2part1Single(meshChoice, n);
		float[] PnX2single = Barycentric2evaluationSingle(resultSingle, xValuesSingle, xiSingle);
		
		FileWriter fw5 = new FileWriter("PointwiseErrorStabilityNewtonIncreasing.txt");
		for(int i=0; i<xValues.length; i++)
			fw5.write(xValues[i]+" "+RelErrorStab(PnX1[i], PnX1single[i])+"\n");
		fw5.close();
		
		FileWriter fw6 = new FileWriter("PointwiseErrorStabilityLagrange.txt");
		for(int i=0; i<xValues.length; i++)
			fw6.write(xValues[i]+" "+RelErrorStab(PnX2[i], PnX2single[i])+"\n");
		fw6.close();
		
		FileWriter fw7 = new FileWriter("StabilityUpperBound.txt");
		double[] UpperBoundvalues = ConditioningUpperBound(n, result, xValues, xi, meshChoice);
		for(int i=0; i<UpperBoundvalues.length; i++)
			fw7.write(xValues[i]+" "+UpperBoundvalues[i]+"\n");
		fw7.close();
		
		FileWriter fw8 = new FileWriter("Conditioning.txt");
		double[] perturbedPnX2 =Perturber(PnX2, xValues);
		double[] perturbedPnX1 = Perturber(PnX1, xValues);
		double[] perturbedFX = Perturber(FX, xValues);
		fw8.write(MaxNorm(PnX2, perturbedPnX2)+" <? "+LebesgueConstant(n, meshChoice)*MaxNorm(FX, perturbedFX)+"\n");
		fw8.write(MaxNorm(PnX1, perturbedPnX1)+" <? "+LebesgueConstant(n, meshChoice)*MaxNorm(FX, perturbedFX)+"\n");
		fw8.close();
		
		
//Convergence(for function 4 only, comment out otherwise)
		/*FileWriter fw9 = new FileWriter("Barycentric2Convergence.txt");
		for(int N =5; N<501; N++){
			double[][] myResult = Barycentric2part1(meshChoice, N);
			double[] myPnX = Barycentric2evaluation(myResult, xValues, xi);
				
			fw9.write(N+" "+MaxNorm(FX, myPnX)+"\n");
		}
		fw9.close();*/
	}

//Barycentric Form 2	
	public static double[][] Barycentric2part1(int meshChoice, int n)throws Exception{
		double[][] result = new double[n+1][3];  //2d array to be returned, column 0 is Betas, column 1 is function values, column 2 is xi's
		double[] xi = new double[n+1];
		
		//Uniform mesh
		if(meshChoice ==1){
			int rightmost = 1;
			int leftmost = -1;
			double h = ((double)(rightmost-leftmost))/n;
			for(int i =0; i< n+1; i++)
				xi[i]= leftmost + i*h;
			
			double[][] gammaI = new double[n+1][2];
			gammaI =  DividedDifferences(n, xi, true);
			
			for(int i=0; i<n+1; i++)
				result[i][0]=gammaI[i][0];
			
			for(int i=0; i<n+1; i++)
				result[i][2]=xi[i];

		}
		//Chebyshev kind 1
		else if(meshChoice == 2){
			for(int i=0; i<=n; i++){
				xi[i] = 11*Math.cos(((2*i+1)*Math.PI)/(2*n+2))+11; //change of variables to be altered------------------
			}
			for(int i=0; i<=n; i++){
				result[i][0]=11*Math.pow(-1,i)*Math.sin(((2*i+1)*Math.PI)/(2*n+2))+11;  //compute Betas
			}
			
		}
		//Chebyshev kind 2
		else if(meshChoice == 3){
			for(int i=0; i<=n; i++)
				xi[i]=Math.cos((Math.PI*i)/n);
			for(int i=0; i<=n; i++){		//compute Betas
				if(i==0 || i==n)
					result[i][0]=10*Math.pow(-1,i); //change of variables to be altered---------------------------
				else
					result[i][0]=10*Math.pow(-1, i)*0.5;
			}
		}
		else
			System.out.println("Invalid mesh choice");
		
		for(int i=0; i<n+1; i++)
			result[i][2]=xi[i];
		for(int i=0; i<=n; i++)
				result[i][1]=function2(xi[i], n); //compute function values------------------------
		return result;
	}
	
	public static double[] Barycentric2evaluation(double[][] result, double[] xValues , double[] xi){
		double[] interpolation = new double[xValues.length];
		
		for(int k =0; k<xValues.length; k++){
			double sigma = 0;
			double tau =0;
			
			boolean isXValuesanXI =false;   //gives NaN if xi=z so use function value instead. 
			for(int j=0; j<result.length; j++){
				if(xValues[k] == result[j][2])
					isXValuesanXI = true;
			}
			
			if(isXValuesanXI)
				interpolation[k]=function2(xValues[k], xi.length); //compute function values-------------------------------
			else{
				for(int i=0; i<result.length; i++){
					double rho = result[i][0]/(xValues[k]-result[i][2]);
					sigma = sigma + result[i][1]*rho;
					tau=tau+rho;
				}
				interpolation[k]=sigma/tau;
			}
		}

		return interpolation;
	}

//Newton	
	public static double[][] DividedDifferences(int n, double[] xi, boolean flag)throws Exception{
		double[][] result = new double[n+1][2]; //as before column 1 is divided differences, column 2 is function values
		
		double[][]omegas = new double[n+1][n+1];//intialize to 1 so that multiplication can occur
		for(int k=0; k<n+1; k++){
			for(int i=0;i<n+1; i++)
				omegas[k][i]=1;
		}
		
		for(int k=1;k<n+1;k++){
			for(int i=0; i<=k;i++){
				for(int j=0; j<=k; j++){
					if(j!=i)
						omegas[k][i]=omegas[k][i]*(xi[i]-xi[j]);
				}
			}
		}

		double[] gammaI = new double[n+1];
		for(int i=0; i<gammaI.length; i++)  //only need w_{n+1}(xi) to compute gammaI
			gammaI[i]=1/omegas[n][i];

		
		for(int i=0; i<=n; i++){
			result[i][1]=function2(xi[i], n); //compute function values-----------------------------
		}

		//if flag= true, output gammaI's in column 1 of result.
		if(flag == true){
			for(int i=0; i<gammaI.length; i++)
				result[i][0]=gammaI[i];
			
			return result;
		}
		else{
			for(int k=0; k<n+1; k++){
				for(int i=0; i<=k;i++){
					result[k][0]=result[k][0]+result[i][1]/omegas[k][i];
				}
			}
			return result;	
		}
	}
	
	public static double[] HornerRule(double[][]result, double[] xi, double[] xValues){
		double[] interpolation = new double[xValues.length];
		
		for(int k=0; k<xValues.length; k++){
			double s = result[xi.length-1][0];
			for(int i=xi.length-2; i>=0; --i){
				s=s*(xValues[k]-xi[i])+result[i][0];
			}
			interpolation[k]=s;
		}
		
		return interpolation;
	}

//Utility Methods
	public static double[] LejaOrdering(int orderingFlag, double[] xi){
		//increasing
		if(orderingFlag == 1){
			for(int j=0; j<xi.length-1; j++){
				for(int i=0; i<xi.length-j-1; i++){
					if(xi[i+1]<xi[i]){
						double dummy =xi[i];
						xi[i]=xi[i+1];
						xi[i+1]=dummy;
					}
				}
			}
		}
		//decreasing
		else if(orderingFlag == 2){
			for(int j=0; j<xi.length-1; j++){
				for(int i=0; i<xi.length-j-1; i++){
					if(xi[i+1]>xi[i]){
						double dummy =xi[i];
						xi[i]=xi[i+1];
						xi[i+1]=dummy;
					}
				}
			}
		}
		//Leja
		else if(orderingFlag == 3){
			//STEP 1: find max(xi) and set x0=max(xi);
			for(int i=0; i<xi.length; i++){
				if(Math.abs(xi[i])>Math.abs(xi[0])){
					double dummy = xi[0];
					xi[0]=xi[i];
					xi[i]=dummy;
				}
			}
			//STEP 2
			for(int k=1; k<xi.length; k++){
				int nextxposition=k;
				double nextx =0;
				double comparator=0;
				System.out.printf("Iteration: %d %n",k);
				for(int j=k;j<xi.length; j++){
					double product=1;
					for(int i=0; i<k; i++)
						product = product * Math.abs(xi[j]-xi[i]);	
					System.out.println(product);
					System.out.println();
					if(product>comparator){
						nextx=xi[j];
						nextxposition=j;
						comparator=product;
					}
				}
				double dummy=xi[k];
				xi[k]=nextx;
				xi[nextxposition]=dummy;
			}
		}
		else
			System.out.println("Invalid ordering selection.");
		
		return xi;
	}

	public static double MaxNorm(double[]PnX, double[]perturbedPnX){ 
		double max=0;
		for(int i=0; i<PnX.length; i++){
			if(Math.abs(PnX[i]-perturbedPnX[i])>max) //since xValues for each were evaluated at corresponding same indices, comparing against corresponding same x values.
				max=Math.abs(PnX[i]-perturbedPnX[i]);
		}
		return max;
	}
	
	public static double[] Perturber(double[] PnX, double[]xValues){
		double[] perturbedPnX = new double[PnX.length];
		for(int i=0; i<PnX.length; i++)
			perturbedPnX[i]=PnX[i]*0.000000001*xValues[i];
		return perturbedPnX;
	}

//Error and Conditioning	
	public static double RelErrorAcc(double real, double approx){
		return Math.abs(real-approx)/Math.abs(real);
	}
	public static double AbsErrorAcc(double real, double approx){return Math.abs(real-approx);}

	public static double RelErrorStab(double real, float approx){
		return Math.abs(real-approx)/Math.abs(real);
	}
	public static double AbsErrorStab(double real, float approx){return Math.abs(real-approx);}
	
	public static double LebesgueConstant(int n, int meshChoice){
		if(meshChoice == 1)//uniform
			return (Math.pow(2, n+1))/(Math.E*n*Math.log(n));
		else if(meshChoice == 2 || meshChoice == 3) //Chebyshev 1 and 2
			return (2*Math.log(n+1))/(Math.PI)+1;
		else{
			System.out.println("Invalid mesh choice for Lebesgue Constant");
			return -1;
		}
	}
	
	public static double ConditionNumber(int n, double[][]result, double x, double[] xi){
		double numerator=0;
		
		for(int i=0; i<n; i++)
			numerator = numerator+Math.abs(function3(x, i, xi)*result[i][1]);
		
		return numerator/Math.abs(function2(x,n));//change function here------------------------------------------
	}
	
	public static double[] ConditioningUpperBound(int n, double[][]result, double[]xValues, double[]xi, int meshChoice){
		double[] UpperBoundvalues= new double[xValues.length];
		
		for(int i=0; i<xValues.length; i++)
			UpperBoundvalues[i]=(3*n+4)*ConditionNumber(n, result,xValues[i],xi)*0.0000000596+(3*n+2)*LebesgueConstant(n, meshChoice)*0.0000000596;
		
		return UpperBoundvalues;
	}
	
//Double precision functions
	public static double testfunction(double x){return x*x +1;}
	public static double function1(double x){return Math.pow((x-2), 9);}
	public static double function2(double x, int d){
		double result =1;
		for(int i=1; i<=d;i++)
			result= result*(x-i);
		return result;
	}
	public static double function3(double x, int n, double[] xi){
		double numerator=1;
		double denominator=1;
		for(int i=0; i<n; i++)
			numerator = numerator * (x - xi[i]);
		for(int i=0; i<n; i++)
			denominator= denominator * (xi[n]-xi[i]);
		return numerator/denominator;
	}
	public static double function4(double x){return 1/(1+ 25*x*x);}
	
	
	
	
	
//////////////////////////////////////////////////////////////////////////////////////////////////////
	//Single Precision versions of Functions (C++ templates don't exist in java)
/////////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	
	
	

//Barycentric Form 2	
	public static float[][] Barycentric2part1Single(int meshChoice, int n)throws Exception{
		float[][] result = new float[n+1][3];  //2d array to be returned, column 0 is Betas, column 1 is function values, column 2 is xi's
		float[] xi = new float[n+1];
		
		//Uniform mesh
		if(meshChoice ==1){
			int rightmost = 1;
			int leftmost = -1;//change indices here--------------------------------------------
			float h = ((float)(rightmost-leftmost))/n;
			
			for(int i =0; i< n+1; i++)
				xi[i]= leftmost + i*h;
			
			float[][] gammaI = new float[n+1][2];
			gammaI =  DividedDifferencesSingle(n, xi, true);
			
			for(int i=0; i<n+1; i++)
				result[i][0]=gammaI[i][0];
			
			for(int i=0; i<n+1; i++)
				result[i][2]=xi[i];

		}
		//Chebyshev kind 1
		else if(meshChoice == 2){
			for(int i=0; i<=n; i++)
				xi[i] = 11*(float)Math.cos(((2*i+1)*Math.PI)/(2*n+2))+11; //change of variables to be altered---------------------
			for(int i=0; i<=n; i++)
				result[i][0]=11*(float)Math.pow(-1,i)*(float)Math.sin(((2*i+1)*Math.PI)/(2*n+2))+11;  //compute Betas
		}
		//Chebyshev kind 2
		else if(meshChoice == 3){
			for(int i=0; i<=n; i++)
				xi[i]=(float)Math.cos(((float)Math.PI*i)/n);//change of variables to be altered-----------------------------
			for(int i=0; i<=n; i++){		//compute Betas
				if(i==0 || i==n)
					result[i][0]=(float)Math.pow(-1,i);
				else
					result[i][0]=(float)Math.pow(-1, i)* (float)0.5;
			}
		}
		else
			System.out.println("Invalid mesh choice");
		
		for(int i=0; i<n+1; i++)
			result[i][2]=xi[i];
		for(int i=0; i<=n; i++)
				result[i][1]=function2single(xi[i], n); //compute function values-----------------------------------
		return result;
	}
	
	
	public static float[] Barycentric2evaluationSingle(float[][] result, float[] xValues, float[] xi){
		float[] interpolation = new float[xValues.length];
		
		for(int k =0; k<xValues.length; k++){
			float sigma = 0;
			float tau =0;
			
			boolean isXValuesanXI =false;   //gives NaN if xi=z so use function value instead. 
			for(int j=0; j<result.length; j++){
				if(xValues[k] == result[j][2])
					isXValuesanXI = true;
			}
			
			if(isXValuesanXI)
				interpolation[k]=function2single(xValues[k], xi.length); //compute function values-------------------------------
			else{
				for(int i=0; i<result.length; i++){
					float rho = result[i][0]/(xValues[k]-result[i][2]);
					sigma = sigma + result[i][1]*rho;
					tau=tau+rho;
				}
				interpolation[k]=sigma/tau;
			}
		}

		return interpolation;
	}

//Newton	
	public static float[][] DividedDifferencesSingle(int n, float[] xi, boolean flag)throws Exception{
		float[][] result = new float[n+1][2]; //as before column 1 is divided differences, column 2 is function values
		
		float[][]omegas = new float[n+1][n+1];//intialize to 1 so that multiplication can occur
		for(int k=0; k<n+1; k++){
			for(int i=0;i<n+1; i++)
				omegas[k][i]=1;
		}
		
		for(int k=1;k<n+1;k++){
			for(int i=0; i<=k;i++){
				for(int j=0; j<=k; j++){
					if(j!=i)
						omegas[k][i]=omegas[k][i]*(xi[i]-xi[j]);
				}
			}
		}

		float[] gammaI = new float[n+1];
		for(int i=0; i<gammaI.length; i++)  //only need w_{n+1}(xi) to compute gammaI
			gammaI[i]=1/omegas[n][i];

		
		for(int i=0; i<=n; i++){
			result[i][1]=function2single(xi[i], n); //compute function values----------------------------------
		}

		//if flag= true, output gammaI's in column 1 of result.
		if(flag == true){
			for(int i=0; i<gammaI.length; i++)
				result[i][0]=gammaI[i];
			
			return result;
		}
		else{
			for(int k=0; k<n+1; k++){
				for(int i=0; i<=k;i++){
					result[k][0]=result[k][0]+result[i][1]/omegas[k][i];
				}
			}
			return result;	
		}
	}
	
	public static float[] HornerRuleSingle(float[][]result, float[] xi, float[] xValues){
		float[] interpolation = new float[xValues.length];
		
		for(int k=0; k<xValues.length; k++){
			float s = result[xi.length-1][0];
			for(int i=xi.length-2; i>=0; --i){
				s=s*(xValues[k]-xi[i])+result[i][0];
			}
			interpolation[k]=s;
		}
		
		return interpolation;
	}
	
	public static float[] LejaOrderingSingle(int orderingFlag, float[] xi){
		//increasing
		if(orderingFlag == 1){
			for(int j=0; j<xi.length-1; j++){
				for(int i=0; i<xi.length-j-1; i++){
					if(xi[i+1]<xi[i]){
						float dummy =xi[i];
						xi[i]=xi[i+1];
						xi[i+1]=dummy;
					}
				}
			}
		}
		//decreasing
		else if(orderingFlag == 2){
			for(int j=0; j<xi.length-1; j++){
				for(int i=0; i<xi.length-j-1; i++){
					if(xi[i+1]>xi[i]){
						float dummy =xi[i];
						xi[i]=xi[i+1];
						xi[i+1]=dummy;
					}
				}
			}
		}
		//Leja
		else if(orderingFlag == 3){
			//STEP 1: find max(xi) and set x0=max(xi);
			for(int i=0; i<xi.length; i++){
				if(Math.abs(xi[i])>Math.abs(xi[0])){
					float dummy = xi[0];
					xi[0]=xi[i];
					xi[i]=dummy;
				}
			}
			//STEP 2
			for(int k=1; k<xi.length; k++){
				int nextxposition=k;
				float nextx =0;
				float comparator=0;
				System.out.printf("Iteration: %d %n",k);
				for(int j=k;j<xi.length; j++){
					float product=1;
					for(int i=0; i<k; i++)
						product = product * Math.abs(xi[j]-xi[i]);	
					System.out.println(product);
					System.out.println();
					if(product>comparator){
						nextx=xi[j];
						nextxposition=j;
						comparator=product;
					}
				}
				float dummy=xi[k];
				xi[k]=nextx;
				xi[nextxposition]=dummy;
			}
		}
		else
			System.out.println("Invalid ordering selection.");
		
		return xi;
	}
	
//Single Precision functions
	public static float testfunctionsingle(float x){return x*x +1;}
	public static float function1single(float x){return (float)Math.pow((x-2), 9);}
	public static float function2single(float x, int d){
		float result =1;
		for(int i=1; i<=d;i++)
			result= result*(x-i);
		return result;
	}
	public static float function3single(float x, int n, float[] xi){
		float numerator=1;
		float denominator=1;
		for(int i=0; i<n; i++)
			numerator = numerator * (x - xi[i]);
		for(int i=0; i<n; i++)
			denominator= denominator * (xi[n]-xi[i]);
		return numerator/denominator;
	}
	public static float function4single(float x){return 1/(1+ 25*x*x);}
}