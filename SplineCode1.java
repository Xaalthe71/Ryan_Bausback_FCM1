import java.lang.Math;
import java.util.*;
import java.io.*;

public class SplineCode1{
	public static void main(String[] args)throws Exception{
		//To use hard-coded values instead of accepting from file, uncomment below and re-comment file input.
		/*File inFile = new File("template-program3data.txt");
		Scanner in = new Scanner(inFile);
		
		double[] XI = new double[Integer.parseInt(in.next())+1];

		double[] FI = new double[XI.length];	
		
		int boundaryType = Integer.parseInt(in.next());

		double boundaryLeft = Double.parseDouble(in.next());

		double boundaryRight = Double.parseDouble(in.next());

		
		for(int i=0; i<XI.length; i++){
			XI[i]=Double.parseDouble(in.next());
		}
		for(int i=0; i<FI.length; i++){
			FI[i]= Double.parseDouble(in.next());
		}*/
		
		int boundaryType =2;  //1 if s'0 and s'n given, 2 if s"0 and s"n given
		double boundaryLeft=0;
		double boundaryRight=0; //if not using input from file
		
		double a =-2;
		double b =2;
		
		//For use with accuracy assessment (uncomment to use)
		FileWriter fw2 = new FileWriter("SplineRatioConvergence2.txt"); 
		
		for(int n=5; n<=1000; n=n*2){			
			//double[] XI = {0.5, 1, 2, 4, 5, 10, 15, 20};
			//double[] FI = {0.04, 0.05, 0.0682, 0.0801, 0.0940, 0.0981, 0.0912, 0.0857};
			/*
			double[] xValues = new double[10000];
			for(int i=0; i<xValues.length; i++){
				xValues[i]= a + i*h;
			}*/
			double h = (b-a)/(n-1);

			double[] XI = new double[n];
			for(int i=0; i<XI.length; i++){
				XI[i]=a+i*h;
			}
			
			double[] FI = new double[n];
			for(int i=0; i<XI.length; i++){
				//FI[i]=XI[i]*XI[i]*XI[i];
				//FI[i]=90/(95+76*XI[i]*XI[i]);//test function 
				if(XI[i]>-1 && XI[i]<1)
					FI[i]= Math.pow(Math.E, -1/(1-XI[i]*XI[i]));
			}
			//evaluation
			int N = 100000;
			double newH = (b-a)/(N-1);
			double[] xValues = new double[N];
			for(int i=0;i<N; i++)
				xValues[i]= a+i*newH;
			
			double[] RealF = new double[N];
			for(int i=0; i<N; i++){
				//RealF[i]=xValues[i]*xValues[i]*xValues[i];
				//RealF[i]= 90/(95+76*xValues[i]*xValues[i]);
				if(xValues[i]>-1 && xValues[i]<1)
					RealF[i]= Math.pow(Math.E, -1/(1-xValues[i]*xValues[i]));
			}
			double[] PX = splineCode1(XI, FI, boundaryLeft, boundaryRight,xValues, boundaryType, false);
			//double[] YprimeT = splineCode1(XI, FI, boundaryLeft, boundaryRight, xValues, boundaryType, true);
			
			//for estimating convergence
			double halfH =h/2;
			
			double[] XI2 = new double[n*2-1];
			for(int i=0; i<XI2.length;i++){
				XI2[i]=a+i*halfH;
			}
			
			
			double[] FI2= new double[n*2-1];
			for(int i=0; i<FI2.length; i++){
				//FI2[i]=XI2[i]*XI2[i]*XI2[i];
				//FI2[i]=90/(95+76*XI2[i]*XI2[i]);
				if(XI2[i]>-1 && XI2[i]<1)
					FI2[i]= Math.pow(Math.E, -1/(1-XI2[i]*XI2[i]));
			}
			double[] PX2 = splineCode1(XI2, FI2, boundaryLeft, boundaryRight, xValues, boundaryType,false);
			//used for checking correctness of output
			System.out.println(MaxNorm(RealF, PX));
			System.out.println(MaxNorm(RealF,PX2));
			System.out.println(MaxNorm(RealF, PX)/MaxNorm(RealF,PX2));
			System.out.println(Math.log(MaxNorm(RealF, PX)/MaxNorm(RealF,PX2)));
			System.out.println();
			
			fw2.write(n+" "+(Math.log(MaxNorm(RealF, PX)/MaxNorm(RealF,PX2))/Math.log(2))+"\n");
			/*
			double[] FT = new double[xValues.length];
			double[] GT = new double[xValues.length];
			
			for(int i=0; i< FT.length; i++)
				FT[i]= YT[i] + xValues[i]*YprimeT[i];
			for(int i=0; i<GT.length; i++)
				GT[i]= Math.exp(-xValues[i]*YT[i]);
			
			
			FileWriter fw1 = new FileWriter("SplineTask2YT2.txt");
			for(int i=0; i<xValues.length; i++)
				fw1.write(YT[i]+"\n");
			fw1.close();
			
			FileWriter fw3 = new FileWriter("SplineTask2FT2.txt");
			for(int i =0; i<xValues.length; i++)
				fw3.write(FT[i]+"\n");
			fw3.close();
			
			FileWriter fw4 = new FileWriter("SplineTask2GT2.txt");
			for(int i=0; i<xValues.length; i++)
				fw4.write(GT[i]+"\n");
			fw4.close();*/
		}
		fw2.close();
	}
	
	
//////////////////////////////DOUBLE PRECISION/////////////////////////////////
	public static double[] splineCode1(double[] XI, double[] FI, double boundaryLeft, double boundaryRight, double[] xValues, int boundaryType, boolean derivative)throws Exception{

		double[] HI = new double[XI.length]; //compute all the interval lengths for nonuniform
		for(int i=1; i<XI.length; i++){
			HI[i] = XI[i]-XI[i-1];
		}
		
		double[] MuI = new double[XI.length-1]; //compute parameters for the matrix
		for(int i=1; i<XI.length-1; i++){
			MuI[i]= HI[i]/(HI[i]+HI[i+1]);
		}
		double[] LambdaI = new double[XI.length-1];
		for(int i=1; i<XI.length-1; i++)
			LambdaI[i]=HI[i+1]/(HI[i]+HI[i+1]);
		
		
		double[] DI = new double[XI.length-1]; //compute vector for RHS
		for(int i=1; i<XI.length-1; i++)
			DI[i]= (6/(HI[i]+HI[i+1]))*((FI[i+1]-FI[i])/HI[i+1] - (FI[i]-FI[i-1])/HI[i]);
		
		double[] mainDiag= new double[DI.length];
		for(int i=1; i<mainDiag.length; i++){
			mainDiag[i]=2;
		}
		if(boundaryType==1){//adjust to get a tridiagonal matrix in the same form as for s"i
			DI[1]=DI[1]- 3*MuI[1]/HI[1]*((FI[1]-FI[0])/HI[1] - boundaryLeft);
			DI[XI.length-2]=DI[XI.length-2] -3*LambdaI[XI.length-2]/HI[XI.length-1]*(boundaryRight-(FI[FI.length-1]-FI[FI.length-2])/HI[HI.length-1]);
			
			mainDiag[1]=mainDiag[1]-MuI[1]/2;
			mainDiag[mainDiag.length-1]=mainDiag[mainDiag.length-1] - LambdaI[LambdaI.length-1]/2;
		}
		
		if(boundaryType==2){
			DI[1]=DI[1]-MuI[1]*boundaryLeft;//adjust the RHS vector to get tridiagonal matrix
			DI[XI.length-2]=DI[XI.length-2]-LambdaI[XI.length-2]*boundaryRight;
		}
		
		double[] sValues = TriDiagMatrixSolver(mainDiag, MuI, LambdaI, DI);//solve tridiagonal matri
		
		if(boundaryType ==1){
			sValues[0]=3/HI[1] * ((FI[1]-FI[0])/HI[1] - boundaryLeft) - 0.5*sValues[1];
			sValues[sValues.length-1] = 3/HI[HI.length-1] * (boundaryRight-(FI[FI.length-1]-FI[FI.length-2])/HI[HI.length-1]) -0.5*sValues[sValues.length-2];
		}
		
		if(boundaryType==2){
			sValues[0]=boundaryLeft;
			sValues[sValues.length-1]=boundaryRight;
		}
		
		//for(int i=0; i<sValues.length; i++)
			//System.out.println(sValues[i]);
		
		double[] gamma = new double[XI.length];
		double[] gammaSquiggle = new double[XI.length];
		
		for(int i=1; i<XI.length;i++)
			gammaSquiggle[i-1]=FI[i-1]-(sValues[i-1]*HI[i]*HI[i])/6;
		
		for(int i=1; i< XI.length; i++)
			gamma[i-1]=(FI[i]-FI[i-1])/HI[i] -(HI[i]*(sValues[i]-sValues[i-1]))/6;
		
		
		double[] result = new double[xValues.length];
		if(!derivative){
			//evaluate the xValues for splines
			for(int i=0; i<xValues.length; i++){
				for(int j=1; j<XI.length; j++){
					if(xValues[i]>=XI[j-1] && xValues[i]<=XI[j]){
						result[i]=sValues[j-1]*Math.pow((XI[j]-xValues[i]), 3)/(6*HI[j]) +sValues[j]*Math.pow((xValues[i]-XI[j-1]), 3)/(6*HI[j]) +gamma[j-1]*(xValues[i]-XI[j-1]) + gammaSquiggle[j-1];
					}
				}
			}
		}
		else{
			for(int i=0; i<xValues.length; i++){
				for(int j=1; j<XI.length; j++){
					if(xValues[i]>=XI[j-1] && xValues[i]<=XI[j]){
						result[i]=-sValues[j-1]*Math.pow((XI[j]-xValues[i]), 2)/(2*HI[j]) +sValues[j]*Math.pow((xValues[i]-XI[j-1]), 2)/(2*HI[j]) +gamma[j-1];
					}
				}
			}
		}
		return result;
	}
	
	public static double[] TriDiagMatrixSolver(double[] a, double[] b, double[] c, double[] f){//uses thomas's algorithm form pg 96 of textbook
		int n = a.length;

		double[] X = new double[n+1];
		
		double[] gamma = new double[n];
		gamma[1]=1/a[1];
		
		for(int i=2; i<n; i++)
			gamma[i]= 1/(a[i]-b[i]*gamma[i-1]*c[i-1]);
		
		double[] Y = new double[n];
		Y[1]=gamma[1]*f[1];
		for(int i=2; i<n; i++)
			Y[i]=gamma[i]*(f[i]-b[i]*Y[i-1]);
		
		X[n-1]=Y[n-1];
		
		for(int i=n-2;i>0; i=i-1)
			X[i]=Y[i]-gamma[i]*c[i]*X[i+1];
		
		return X;
	}
	
	public static double MaxNorm(double[]PnX, double[]perturbedPnX){ 
		double max=0;
		for(int i=0; i<PnX.length; i++){
			if(Math.abs(PnX[i]-perturbedPnX[i])>max) //since xValues for each were evaluated at corresponding same indices, comparing against corresponding same x values.
				max=Math.abs(PnX[i]-perturbedPnX[i]);
		}
		return max;
	}
	
	
	/*
/////////////////////////////////////////Single PRECISION///////////////////////////////////

	public static float[] splineCode1Single(float[] XI, float[] FI, float boundaryLeft, float boundaryRight, float[] xValues, int boundaryType, boolean derivative)throws Exception{

		float[] HI = new float[XI.length]; //compute all the interval lengths for nonuniform
		for(int i=1; i<XI.length; i++){
			HI[i] = XI[i]-XI[i-1];
		}
		
		float[] MuI = new float[XI.length-1]; //compute parameters for the matrix
		for(int i=1; i<XI.length-1; i++){
			MuI[i]= HI[i]/(HI[i]+HI[i+1]);
		}
		float[] LambdaI = new float[XI.length-1];
		for(int i=1; i<XI.length-1; i++)
			LambdaI[i]=HI[i+1]/(HI[i]+HI[i+1]);
		
		
		float[] DI = new float[XI.length-1]; //compute vector for RHS
		for(int i=1; i<XI.length-1; i++)
			DI[i]= (6/(HI[i]+HI[i+1]))*((FI[i+1]-FI[i])/HI[i+1] - (FI[i]-FI[i-1])/HI[i]);
		
		float[] mainDiag= new float[DI.length];
		for(int i=1; i<mainDiag.length; i++){
			mainDiag[i]=2;
		}
		if(boundaryType==1){//adjust to get a tridiagonal matrix in the same form as for s"i
			DI[1]=DI[1]- 3*MuI[1]/HI[1]*((FI[1]-FI[0])/HI[1] - boundaryLeft);
			DI[XI.length-2]=DI[XI.length-2] -3*LambdaI[XI.length-2]/HI[XI.length-1]*(boundaryRight-(FI[FI.length-1]-FI[FI.length-2])/HI[HI.length-1]);
			
			mainDiag[1]=mainDiag[1]-MuI[1]/2;
			mainDiag[mainDiag.length-1]=mainDiag[mainDiag.length-1] - LambdaI[LambdaI.length-1]/2;
		}
		
		if(boundaryType==2){
			DI[1]=DI[1]-MuI[1]*boundaryLeft;//adjust the RHS vector to get tridiagonal matrix
			DI[XI.length-2]=DI[XI.length-2]-LambdaI[XI.length-2]*boundaryRight;
		}
		
		float[] sValues = TriDiagMatrixSolverSingle(mainDiag, MuI, LambdaI, DI);//solve tridiagonal matri
		
		if(boundaryType ==1){
			sValues[0]=3/HI[1] * ((FI[1]-FI[0])/HI[1] - boundaryLeft) - 0.5*sValues[1];
			sValues[sValues.length-1] = 3/HI[HI.length-1] * (boundaryRight-(FI[FI.length-1]-FI[FI.length-2])/HI[HI.length-1]) -0.5*sValues[sValues.length-2];
		}
		
		if(boundaryType==2){
			sValues[0]=boundaryLeft;
			sValues[sValues.length-1]=boundaryRight;
		}
		
		//for(int i=0; i<sValues.length; i++)
			//System.out.println(sValues[i]);
		
		float[] gamma = new float[XI.length];
		float[] gammaSquiggle = new float[XI.length];
		
		for(int i=1; i<XI.length;i++)
			gammaSquiggle[i-1]=FI[i-1]-(sValues[i-1]*HI[i]*HI[i])/6;
		
		for(int i=1; i< XI.length; i++)
			gamma[i-1]=(FI[i]-FI[i-1])/HI[i] -(HI[i]*(sValues[i]-sValues[i-1]))/6;
		
		
		float[] result = new float[xValues.length];
		if(!derivative){
			//evaluate the xValues for splines
			for(int i=0; i<xValues.length; i++){
				for(int j=1; j<XI.length; j++){
					if(xValues[i]>=XI[j-1] && xValues[i]<=XI[j]){
						result[i]=sValues[j-1]*Math.pow((XI[j]-xValues[i]), 3)/(6*HI[j]) +sValues[j]*Math.pow((xValues[i]-XI[j-1]), 3)/(6*HI[j]) +gamma[j-1]*(xValues[i]-XI[j-1]) + gammaSquiggle[j-1];
					}
				}
			}
		}
		else{
			for(int i=0; i<xValues.length; i++){
				for(int j=1; j<XI.length; j++){
					if(xValues[i]>=XI[j-1] && xValues[i]<=XI[j]){
						result[i]=-sValues[j-1]*Math.pow((XI[j]-xValues[i]), 2)/(2*HI[j]) +sValues[j]*Math.pow((xValues[i]-XI[j-1]), 2)/(2*HI[j]) +gamma[j-1];
					}
				}
			}
		}
		return result;
	}
	
	public static float[] TriDiagMatrixSolverSingle(float[] a, float[] b, float[] c, float[] f){//uses thomas's algorithm form pg 96 of textbook
		int n = a.length;

		float[] X = new float[n+1];
		
		float[] gamma = new float[n];
		gamma[1]=1/a[1];
		
		for(int i=2; i<n; i++)
			gamma[i]= 1/(a[i]-b[i]*gamma[i-1]*c[i-1]);
		
		float[] Y = new float[n];
		Y[1]=gamma[1]*f[1];
		for(int i=2; i<n; i++)
			Y[i]=gamma[i]*(f[i]-b[i]*Y[i-1]);
		
		X[n-1]=Y[n-1];
		
		for(int i=n-2;i>0; i=i-1)
			X[i]=Y[i]-gamma[i]*c[i]*X[i+1];
		
		return X;
	}
	
	public static float MaxNormSingle(float[]PnX, float[]perturbedPnX){ 
		float max=0;
		for(int i=0; i<PnX.length; i++){
			if(Math.abs(PnX[i]-perturbedPnX[i])>max) //since xValues for each were evaluated at corresponding same indices, comparing against corresponding same x values.
				max=Math.abs(PnX[i]-perturbedPnX[i]);
		}
		return max;
	}*/
}