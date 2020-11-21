import java.lang.Math;
import java.util.*;
import java.io.*;

public class CompTrapRule{
	public static void main(String[] args)throws Exception{
		
		//FileWriter fw1 = new FileWriter("ErrorEstimatesTrapex.txt");
		FileWriter fw2 = new FileWriter("ExactErrorsTrap4.txt");
		//FileWriter fw3 = new FileWriter("IntegralValueTrapex.txt");
		
		double error=1;
		double Im =1;
		double threshold= Math.pow(10, -4);
	
		int m=2;
		double a=0.1;
		double b=2.5;
		
		//double exactIntegralValue=1.71828182;
		//double exactIntegralValue = 2.0/3.0; //need .0 to get floats, otherwise does integer division
		//double exactIntegralValue = Math.pow(Math.E, 3)-1;//19.085553692;
		//double exactIntegralValue = 0.5*(-1+Math.pow(Math.E, Math.sqrt(3.0)/2.0));//0.688721337618;
		//double exactIntegralValue = Math.log(Math.cosh(1)/Math.cosh(2));
		//double exactIntegralValue = -1/(2*Math.PI*Math.PI);
		double exactIntegralValue = 3.12 +Math.log(2.5/0.1);//6.338875824868;
		
		while(Math.abs((exactIntegralValue-Im))>threshold){
			double Hm = (b-a)/m;
			double H2m = (b-a)/(2*m);
			double[] XI = new double[m+1];
			for(int i=0; i<=m; i++){
				XI[i]=a+i*Hm;
			}
			
			double[] XI2 = new double[m]; //the finer mesh points, since halving mesh, take average of the 2 points on either side.
			for(int i=0; i<m; i++){
				XI2[i]=(XI[i]+XI[i+1])/2;
			}
			double[] FI = new double[m+1];
			for(int i=0; i<=m;i++){
				//FI[i]=XI[i]*XI[i]; //test functions(uncomment respectively to use)
				//FI[i]=Math.pow(Math.E, XI[i]);
				//FI[i]=Math.pow(Math.E, Math.sin(2*XI[i]))*Math.cos(2*XI[i]);
				//FI[i]=Math.tanh(XI[i]);
				//FI[i]=XI[i]*Math.cos(2*Math.PI*XI[i]);
				FI[i]=XI[i] + (1/XI[i]);
			}
			Im = CompTrapRule(Hm, FI);
			
			System.out.println(Im);
			
			
			double[] FInew = new double[m];
			for(int i=0; i<m; i++){
				//FInew[i]=XI2[i]*XI2[i];//test function
				//FInew[i]=Math.pow(Math.E, XI2[i]);
				//FInew[i]=Math.pow(Math.E, Math.sin(2*XI2[i]))*Math.cos(2*XI2[i]);
				//FInew[i]=Math.tanh(XI2[i]);
				//FInew[i]=XI2[i]*Math.cos(2*Math.PI*XI2[i]);
				FInew[i]=XI2[i] + (1/XI2[i]);
			}	
			
			double I2m = TrapRuleFineMesh(Hm, Im, FInew);
			System.out.println(I2m);
			System.out.println();
			
			
			error = FineError(Im, I2m);
			
			//fw1.write(H2m+" "+error+"\n");
			fw2.write(H2m+" "+(exactIntegralValue-Im)+"\n");
			//fw3.write(I2m+"\n");
			
			m = 2*m;
		}
		
		///fw1.close();
		fw2.close();
		//fw3.close();
	}
	
	public static double CompTrapRule(double Hm, double[] FI){
		double sum=0;
		for(int i=1; i<=FI.length-2;i++){
			sum = sum + FI[i];
		}
		return Hm/2*(FI[0]+FI[FI.length-1]+2*sum);
	}
	
	public static double TrapRuleFineMesh(double Hm, double Im, double[]FInew){
		double sum=0;
		for(int i=0; i<FInew.length; i++)
			sum= sum + FInew[i];
			
		return 0.5*(Im + Hm*sum);
	}
	
	public static double CoarseError(double Im, double I2m){return 4/3 * (I2m-Im);}
	public static double FineError(double Im, double I2m){return (I2m - Im)/3;}

}