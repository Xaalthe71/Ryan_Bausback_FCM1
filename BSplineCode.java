import java.lang.Math;
import java.util.*;
import java.io.*;

public class BSplineCode{
	public static void main(String[] args)throws Exception{
		/* //Uncomment if accepting input from file instead of using hard-coded values (then recomment the hard coded values in main)
		File inFile = new File("template-program3data.txt");
		Scanner in = new Scanner(inFile);
		
		double[] XI = new double[Integer.parseInt(in.next())+1];

		double[] FI = new double[XI.length];	
		
		int boundaryType = Integer.parseInt(in.next());

		double boundaryLeft = Double.parseDouble(in.next());

		double boundaryRight = Double.parseDouble(in.next());
		*/
		/*
		for(int i=0; i<XI.length; i++){
			XI[i]=Double.parseDouble(in.next());
		}
		for(int i=0; i<FI.length; i++){
			FI[i]= Double.parseDouble(in.next());
		}*/
		
		
		
		
		int boundaryType =2; //1 if s'0 and s'n given, 2 if s"0 and s"n given
		double boundaryLeft = 0;
		double boundaryRight= 0;
		
		double a =-2;
		double b =2;

		FileWriter fw2 = new FileWriter("BSplineRatioConvergence.txt");
		
		for(int n=5; n<1000; n=n*2){
			double h = (b-a)/(n-1);

			double[] XI = new double[n]; //uniform mesh ONLY
			for(int i=0; i<XI.length; i++){
				XI[i]=a+i*h;
				//System.out.println(XI[i]);
			}
			System.out.println();
			double[] FI = new double[n];
			for(int i=0; i<XI.length; i++){
				//FI[i]= XI[i]*XI[i]*XI[i];//test function 
				//FI[i]=9/(3+71*XI[i]*XI[i]);
				if(XI[i]>-1 && XI[i]<1)
					FI[i]= Math.pow(Math.E, -1/(1-XI[i]*XI[i]));
			}
			double halfH=h/2;
			
			double[] XI2 = new double[n*2-1];
			for(int i=0; i<XI2.length;i++){
				XI2[i]=a+i*halfH;
				//System.out.println(XI2[i]);
			}
			double[] FI2= new double[n*2-1];
			for(int i=0; i<FI2.length; i++){
				//FI2[i]=XI2[i]*XI2[i]*XI2[i];
				//FI2[i]=9/(3+71*XI2[i]*XI2[i]);
				if(XI2[i]>-1 && XI2[i]<1)
					FI2[i]= Math.pow(Math.E, -1/(1-XI2[i]*XI2[i]));
			}
			//evaluation
			int N = 100000;
			double newH = (b-a)/(N-1);
			double[] xValues = new double[N];
			for(int i=0;i<N; i++)
				xValues[i]= a+i*newH;
			
			double[] RealF = new double[N];
			for(int i=0; i<N; i++){
				//RealF[i]= xValues[i]*xValues[i]*xValues[i];
				//RealF[i]=9/(3+71*xValues[i]*xValues[i]);
				if(xValues[i]>-1 && xValues[i]<1)
					RealF[i]= Math.pow(Math.E, -1/(1-xValues[i]*xValues[i]));
			}
			double[] PX = Bspline(XI, FI, boundaryLeft, boundaryRight,xValues, boundaryType);
			System.out.println(MaxNorm(RealF, PX));
			double[] PX2 = Bspline(XI2, FI2, boundaryLeft, boundaryRight, xValues, boundaryType);
			
			FileWriter fw1 = new FileWriter("BSplineTest1.txt");
			for(int i=0; i<N; i++)
				fw1.write(xValues[i]+" "+PX[i]+"\n");
			fw1.close();
			
			fw2.write(n+" "+(Math.log(MaxNorm(RealF, PX)/MaxNorm(RealF,PX2))/Math.log(2))+"\n");
		}
		fw2.close();
	}

////////////////////////////////DOUBLE PRECISION////////////////////////////////////
	public static double[] Bspline(double[] XI, double[] FI, double boundaryLeft, double boundaryRight, double[] xValues, int boundaryType){
		double h = XI[1]-XI[0]; //assuming uniform spacing so all intervals should be the same length
		
		//set up the tridiagonal matrix
		double[] mainDiag = new double[XI.length+1];
		for(int i=1; i<mainDiag.length; i++){
			mainDiag[i]=4;
		}
		
		double[] upperDiag = new double[XI.length];
		for(int i=1; i<upperDiag.length; i++){
			upperDiag[i]=1;
		}

		double[] lowerDiag = new double[XI.length+1];
		for(int i=2; i<lowerDiag.length; i++)
			lowerDiag[i]=1;
		
		double[] RHS = new double[FI.length+2];
		for(int i=1; i<RHS.length-1; i++)
			RHS[i]=FI[i-1];
		
		if(boundaryType==1){//adjust matrix based on boundary condition
			RHS[1]= (3.0/h)* RHS[1]+boundaryLeft;
			RHS[RHS.length-2]=(-3.0/h) * RHS[RHS.length-2]+boundaryRight;
			
			mainDiag[1]=mainDiag[1]*3/h;
			mainDiag[mainDiag.length-1]=mainDiag[mainDiag.length-1]*-3/h;
			
			upperDiag[1]=6/h;
			lowerDiag[lowerDiag.length-1]=-6/h;
		}
		if(boundaryType==2){
			RHS[1]=(-6.0/(h*h))*RHS[1]+boundaryLeft;
			RHS[RHS.length-2]=(-6.0/(h*h))*RHS[RHS.length-2]+boundaryRight;
			
			mainDiag[1]=(-36.0/(h*h));
			mainDiag[mainDiag.length-1]=(-36.0/(h*h));
			
			upperDiag[1]=0;
			lowerDiag[lowerDiag.length-1]=0;
		}
	
		double[] aVector=TriDiagMatrixSolver(mainDiag, lowerDiag, upperDiag, RHS);
		
		//aVector adjustment is indepedent of boundary condition as a_{-1} and a_{n+1} are present in normal function.
		aVector[0]=FI[0]-4*aVector[1]-aVector[2];
		aVector[aVector.length-1]=FI[FI.length-1]-4*aVector[aVector.length-2]-aVector[aVector.length-3];
		
		//evaluate the Bspline for each z 
		double[] result=new double[xValues.length];
		for(int i=0; i<xValues.length; i++){
			for(int j=0; j<XI.length-1; j++){
				if(xValues[i]>=XI[j] && xValues[i]<=XI[j+1]){
					result[i]=aVector[j]*(1/(Math.pow(h, 3))*Math.pow((XI[j+1]-xValues[i]),3))+ //a[i-1]*B[i-1]
									aVector[j+1]*1/(Math.pow(h,3))*(Math.pow(h,3)+3*h*h*(XI[j+1]-xValues[i])+3*h*Math.pow((XI[j+1]-xValues[i]),2)-3*Math.pow((XI[j+1]-xValues[i]),3))+
										aVector[j+2]*1/(Math.pow(h,3))*(Math.pow(h,3)+3*h*h*(xValues[i]-XI[j])+3*h*Math.pow((xValues[i]-XI[j]),2)-3*Math.pow((xValues[i]-XI[j]),3))+
											aVector[j+3]*1/(Math.pow(h, 3))*Math.pow((xValues[i]-XI[j]),3);
											
					
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
/////////////////////////////SIngle Precision////////////////////////////////////////////
	public static float[] Bspline(float[] XI, float[] FI, float boundaryLeft, float boundaryRight, float[] xValues, int boundaryType){
		float h = XI[1]-XI[0]; //assuming uniform spacing so all intervals should be the same length
		
		//set up the tridiagonal matrix
		float[] mainDiag = new float[XI.length+1];
		for(int i=1; i<mainDiag.length; i++){
			mainDiag[i]=4;
		}
		
		float[] upperDiag = new float[XI.length];
		for(int i=1; i<upperDiag.length; i++){
			upperDiag[i]=1;
		}

		float[] lowerDiag = new float[XI.length+1];
		for(int i=2; i<lowerDiag.length; i++)
			lowerDiag[i]=1;
		
		float[] RHS = new float[FI.length+2];
		for(int i=1; i<RHS.length-1; i++)
			RHS[i]=FI[i-1];
		
		if(boundaryType==1){//adjust matrix based on boundary condition
			RHS[1]= (3.0/h)* RHS[1]+boundaryLeft;
			RHS[RHS.length-2]=(-3.0/h) * RHS[RHS.length-2]+boundaryRight;
			
			mainDiag[1]=mainDiag[1]*3/h;
			mainDiag[mainDiag.length-1]=mainDiag[mainDiag.length-1]*-3/h;
			
			upperDiag[1]=6/h;
			lowerDiag[lowerDiag.length-1]=-6/h;
		}
		if(boundaryType==2){
			RHS[1]=(-6.0/(h*h))*RHS[1]+boundaryLeft;
			RHS[RHS.length-2]=(-6.0/(h*h))*RHS[RHS.length-2]+boundaryRight;
			
			mainDiag[1]=(-36.0/(h*h));
			mainDiag[mainDiag.length-1]=(-36.0/(h*h));
			
			upperDiag[1]=0;
			lowerDiag[lowerDiag.length-1]=0;
		}
	
		float[] aVector=TriDiagMatrixSolverSingle(mainDiag, lowerDiag, upperDiag, RHS);
		
		//aVector adjustment is indepedent of boundary condition as a_{-1} and a_{n+1} are present in normal function.
		aVector[0]=FI[0]-4*aVector[1]-aVector[2];
		aVector[aVector.length-1]=FI[FI.length-1]-4*aVector[aVector.length-2]-aVector[aVector.length-3];
		
		//evaluate the Bspline for each z 
		float[] result=new float[xValues.length];
		for(int i=0; i<xValues.length; i++){
			for(int j=0; j<XI.length-1; j++){
				if(xValues[i]>=XI[j] && xValues[i]<=XI[j+1]){
					result[i]=aVector[j]*(1/(Math.pow(h, 3))*Math.pow((XI[j+1]-xValues[i]),3))+ //a[i-1]*B[i-1]
									aVector[j+1]*1/(Math.pow(h,3))*(Math.pow(h,3)+3*h*h*(XI[j+1]-xValues[i])+3*h*Math.pow((XI[j+1]-xValues[i]),2)-3*Math.pow((XI[j+1]-xValues[i]),3))+
										aVector[j+2]*1/(Math.pow(h,3))*(Math.pow(h,3)+3*h*h*(xValues[i]-XI[j])+3*h*Math.pow((xValues[i]-XI[j]),2)-3*Math.pow((xValues[i]-XI[j]),3))+
											aVector[j+3]*1/(Math.pow(h, 3))*Math.pow((xValues[i]-XI[j]),3);
											
					
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