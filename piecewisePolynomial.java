import java.lang.Math;
import java.util.*;
import java.io.*;

public class piecewisePolynomial{
	public static void main(String[] args)throws Exception{
		int degree = 1;
		//double a =-1;
		/*//Task 2 file input section
		File inFile = new File("template-program3data.txt");
		Scanner in = new Scanner(inFile);
		
		double[] XI = new double[Integer.parseInt(in.next())+1];

		double[] FI = new double[XI.length];	
		
		int boundaryType = Integer.parseInt(in.next());

		double boundaryLeft = Double.parseDouble(in.next()); //not needed for piecewise

		double boundaryRight = Double.parseDouble(in.next()); //not needed for piecewise

		
		for(int i=0; i<XI.length; i++){
			XI[i]=Double.parseDouble(in.next());
		}
		for(int i=0; i<FI.length; i++){
			FI[i]= Double.parseDouble(in.next());
		}
		
		double a = XI[0];*/
		//for accuracy and convergence evaluation, uncomment to use
		
		double a =-3;
		double b =3;
		FileWriter fw2 = new FileWriter("piecewiseConvergence2Random.txt");
		
		for(int j=1; j<100; j++){
			int numIntervals =j;
			
			double h = (b-a)/(numIntervals*degree);
			double halfH = (b-a)/(2*numIntervals*degree);
			
			double[] XI = new double[degree*numIntervals+1];
			for(int i=0; i<XI.length; i++){
				XI[i]=a+i*h;
			}
			double[] XI2 = new double[2*degree*numIntervals+1];
			for(int i=0; i<XI2.length; i++)
				XI2[i]=a+i*halfH;
			
			double[] FI = new double[degree*numIntervals+1];
			for(int i=0; i<XI.length; i++){
				if(XI[i]>=-2+1/10 && XI[i]<=2-(1/10)){
					FI[i]=1;
				}
				else if(XI[i]>=-2 && XI[i]<-2+(1/10)){
					FI[i]= 10*(XI[i]+2);
				}
				else if(XI[i]>2-(1/10) && XI[i]<=2){
					FI[i]=-10*(XI[i]-2);
				}
				else{
					FI[i]=0;
				}
			}
			double[] FI2 = new double[2*degree*numIntervals+1];
			for(int i=0; i<XI2.length; i++){
				if(XI2[i]>=-2+1/10 && XI2[i]<=2-(1/10)){
					FI2[i]=1;
				}
				else if(XI2[i]>=-2 && XI2[i]<-2+1/10){
					FI2[i]= 10*(XI2[i]+2);
				}
				else if(XI2[i]>2-(1/10) && XI2[i]<=2){
					FI2[i]=-10*(XI2[i]-2);
				}
				else{
					FI2[i]=0;
				}
			}
			
			//evaluation
			//int n = 1000;
			//double newH = (b-a)/n;
			int n =1000;
			double newH = (b-a)/n; //task 2 values
			
			double[] xValues = new double[n];
			for(int i=0;i<n; i++)
				xValues[i]= a+i*newH;
			
			double[] Pi1X = new double[n];
			for(int i=0; i<n; i++)
				Pi1X[i] = piecewisePolynomial(XI, FI, degree, xValues[i]);
			
			double[] Pi1Xderivative = new double[n];
			for(int i=0; i<n; i++)
				Pi1Xderivative[i] = piecewisePolynomialDerivative(XI, FI, degree, xValues[i]);
			
			double[] Pi1X2 = new double[n];
			for(int i=0; i<n; i++)
				Pi1X2[i] = piecewisePolynomial(XI2,FI2, degree, xValues[i]);
			
			FileWriter fw1 = new FileWriter("piecewiseTask2YT1.txt");
			for(int i=0; i<n; i++)
				fw1.write(Pi1X[i]+"\n");
			fw1.close();
			/*
			double[] FT = new double[n];
			for(int i=0; i<n; i++)
				FT[i]=Pi1X[i]+xValues[i]*Pi1Xderivative[i];
			
			FileWriter fw2 = new FileWriter("piecewiseTask2FT1.txt");
			for(int i=0; i<n; i++)
				fw2.write(FT[i]+"\n");
			fw2.close();
			
			double[] GT = new double[n];
			for(int i=0; i<n; i++)
				GT[i]=Math.exp(-xValues[i]*Pi1X[i]);
			
			FileWriter fw3 = new FileWriter("piecewiseTask2GT1.txt");
			for(int i=0; i<n; i++)
				fw3.write(GT[i]+"\n");
			fw3.close();*/

			double[] realF = new double[n];
			for(int i=0; i<n; i++){
				if(xValues[i]>=-2+1/10 && xValues[i]<=2-(1/10)){
					realF[i]=1;
				}
				else if(xValues[i]>=-2 && xValues[i]<-2+1/10){
					realF[i]= 10*(xValues[i]+2);
				}
				else if(xValues[i]>2-(1/10) && xValues[i]<=2){
					realF[i]=-10*(xValues[i]-2);
				}
				else{
					realF[i]=0;
				}
			}
			fw2.write(j+" "+Math.log(MaxNorm(realF, Pi1X)/MaxNorm(realF, Pi1X2))/Math.log(2)+"\n");
		}
		fw2.close();
		
	}

///////////////////////////////Double Precision/////////////////////////////////////

	public static double piecewisePolynomial(double[] XI, double[]FI, int degree, double xValue){//using cardinal basis functions, evaluation combined with creation
		if(degree ==1){
			double Pi1X = 0;
			for(int i=1; i<XI.length;i++){//start at 1 since using i-1 in index
				if(xValue >= XI[i-1] && xValue<=XI[i]){//only consider the case were x is inside interval since 0 otherwise
					Pi1X = (FI[i-1]*(xValue - XI[i])/(XI[i-1]-XI[i]))+(FI[i]*(xValue - XI[i-1])/(XI[i]-XI[i-1]));
				}
			}
			return Pi1X;
		}
		else if(degree ==2){
			double Pi2X =0;
			for(int i=2; i<XI.length;i++){//start at 2 since using i-2 in index
				if(xValue >= XI[i-2] && xValue <=XI[i]){
					Pi2X = FI[i-2]*((xValue-XI[i-1])*(xValue-XI[i]))/((XI[i-2]-XI[i-1])*(XI[i-2]-XI[i])) + 
								FI[i-1]*((xValue-XI[i-2])*(xValue-XI[i]))/((XI[i-1]-XI[i-2])*(XI[i-1]-XI[i])) +
									FI[i]*((xValue-XI[i-2])*(xValue-XI[i-1]))/((XI[i]-XI[i-2])*(XI[i]-XI[i-1]));
				}
			}
			return Pi2X;
		}
		else{
			System.out.println("Invalid degree input.");
			return 0;
		}
	}
	
	//for use with derivative calculation in Task2 
	public static double piecewisePolynomialDerivative(double[] XI, double[]FI, int degree, double xValue){
		if(degree ==1){
			double Pi1X = 0;
			for(int i=1; i<XI.length;i++){//start at 1 since using i-1 in index
				if(xValue >= XI[i-1] && xValue<=XI[i]){//only consider the case were x is inside interval since 0 otherwise
					Pi1X = (FI[i-1]/(XI[i-1]-XI[i]))+(FI[i]/(XI[i]-XI[i-1]));
				}
			}
			return Pi1X;
		}
		else if(degree ==2){
			double Pi2X =0;
			for(int i=2; i<XI.length;i++){//start at 2 since using i-2 in index
				if(xValue >= XI[i-2] && xValue <=XI[i]){
					Pi2X = FI[i-2]*(xValue-XI[i-1])/((XI[i-2]-XI[i-1])*(XI[i-2]-XI[i])) + 
							FI[i-2]*(xValue-XI[i])/((XI[i-2]-XI[i-1])*(XI[i-2]-XI[i]))+
								FI[i-1]*(xValue-XI[i])/((XI[i-1]-XI[i-2])*(XI[i-1]-XI[i])) +
									FI[i-1]*(xValue-XI[i-2])/((XI[i-1]-XI[i-2])*(XI[i-1]-XI[i]))+
										FI[i]*(xValue-XI[i-2])/((XI[i]-XI[i-2])*(XI[i]-XI[i-1]))+
											FI[i]*(xValue-XI[i-1])/((XI[i]-XI[i-2])*(XI[i]-XI[i-1]));
				}
			}
			return Pi2X;
		}
		else{
			System.out.println("Invalid degree input.");
			return 0;
		}
	}
	
	public static double MaxNorm(double[]PnX, double[]perturbedPnX){ 
		double max=0;
		for(int i=0; i<PnX.length; i++){
			if(Math.abs(PnX[i]-perturbedPnX[i])>max) //since xValues for each were evaluated at corresponding same indices, comparing against corresponding same x values.
				max=Math.abs(PnX[i]-perturbedPnX[i]);
		}
		return max;
	}
///////////////////////////Single Precision/////////////////////////////////////


	public static float piecewisePolynomialSingle(float[] XI, float[]FI, int degree, float xValue){//using cardinal basis functions, evaluation combined with creation
		if(degree ==1){
			float Pi1X = 0;
			for(int i=1; i<XI.length;i++){//start at 1 since using i-1 in index
				if(xValue >= XI[i-1] && xValue<=XI[i]){//only consider the case were x is inside interval since 0 otherwise
					Pi1X = (FI[i-1]*(xValue - XI[i])/(XI[i-1]-XI[i]))+(FI[i]*(xValue - XI[i-1])/(XI[i]-XI[i-1]));
				}
			}
			return Pi1X;
		}
		else if(degree ==2){
			float Pi2X =0;
			for(int i=2; i<XI.length;i++){//start at 2 since using i-2 in index
				if(xValue >= XI[i-2] && xValue <=XI[i]){
					Pi2X = FI[i-2]*((xValue-XI[i-1])*(xValue-XI[i]))/((XI[i-2]-XI[i-1])*(XI[i-2]-XI[i])) + 
								FI[i-1]*((xValue-XI[i-2])*(xValue-XI[i]))/((XI[i-1]-XI[i-2])*(XI[i-1]-XI[i])) +
									FI[i]*((xValue-XI[i-2])*(xValue-XI[i-1]))/((XI[i]-XI[i-2])*(XI[i]-XI[i-1]));
				}
			}
			return Pi2X;
		}
		else{
			System.out.println("Invalid degree input.");
			return 0;
		}
	}
	
	//for use with derivative calculation in Task2 
	public static float piecewisePolynomialDerivativeSingle(float[] XI, float[]FI, int degree, float xValue){
		if(degree ==1){
			float Pi1X = 0;
			for(int i=1; i<XI.length;i++){//start at 1 since using i-1 in index
				if(xValue >= XI[i-1] && xValue<=XI[i]){//only consider the case were x is inside interval since 0 otherwise
					Pi1X = (FI[i-1]/(XI[i-1]-XI[i]))+(FI[i]/(XI[i]-XI[i-1]));
				}
			}
			return Pi1X;
		}
		else if(degree ==2){
			float Pi2X =0;
			for(int i=2; i<XI.length;i++){//start at 2 since using i-2 in index
				if(xValue >= XI[i-2] && xValue <=XI[i]){
					Pi2X = FI[i-2]*(xValue-XI[i-1])/((XI[i-2]-XI[i-1])*(XI[i-2]-XI[i])) + 
							FI[i-2]*(xValue-XI[i])/((XI[i-2]-XI[i-1])*(XI[i-2]-XI[i]))+
								FI[i-1]*(xValue-XI[i])/((XI[i-1]-XI[i-2])*(XI[i-1]-XI[i])) +
									FI[i-1]*(xValue-XI[i-2])/((XI[i-1]-XI[i-2])*(XI[i-1]-XI[i]))+
										FI[i]*(xValue-XI[i-2])/((XI[i]-XI[i-2])*(XI[i]-XI[i-1]))+
											FI[i]*(xValue-XI[i-1])/((XI[i]-XI[i-2])*(XI[i]-XI[i-1]));
				}
			}
			return Pi2X;
		}
		else{
			System.out.println("Invalid degree input.");
			return 0;
		}
	}
	
	public static float MaxNormSingle(float[]PnX, float[]perturbedPnX){ 
		float max=0;
		for(int i=0; i<PnX.length; i++){
			if(Math.abs(PnX[i]-perturbedPnX[i])>max) //since xValues for each were evaluated at corresponding same indices, comparing against corresponding same x values.
				max=Math.abs(PnX[i]-perturbedPnX[i]);
		}
		return max;
	}

}
