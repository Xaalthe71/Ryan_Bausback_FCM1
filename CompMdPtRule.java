import java.lang.Math;
import java.util.*;
import java.io.*;

public class CompMdPtRule{
	public static void main(String[] args)throws Exception{
		//FileWriter fw1 = new FileWriter("ErrorEstimatesMdptex.txt");
		FileWriter fw2 = new FileWriter("ExactErrorsMdpt8.txt");
		//FileWriter fw3 = new FileWriter("IntegralValueMdptex.txt");
		
		double error =1;
		double exacterror=1;
		double threshold = Math.pow(10, -8);
		
		int m =2;
		double a=0;
		double b=3.0;
		
		//double exactIntegralValue=1.71828182;
		//double exactIntegralValue = 2.0/3.0;
		double exactIntegralValue =Math.pow(Math.E, 3)-1;//19.085553692;
		//double exactIntegralValue = 0.5*(-1+Math.pow(Math.E, Math.sqrt(3.0)/2.0));//0.688721337618;
		//double exactIntegralValue = Math.log(Math.cosh(1)/Math.cosh(2));
		//double exactIntegralValue = -1/(2*Math.PI*Math.PI);
		//double exactIntegralValue =3.12 +Math.log(2.5/0.1); //6.338875824868;
		
		while(Math.abs(exacterror)>threshold){
			
			double Hm = (b-a)/m;
			double thirdHm = (b-a)/(3*m);
			
			double[] Ai = new double[m+1];
			for(int i=0; i<=m; i++){
				Ai[i]=a+i*Hm;
			}
			System.out.println();
			double[] XI = new double[m];
			for(int i=0; i<m; i++){
				XI[i]=(Ai[i]+Ai[i+1])/2;
			}
			
			double[] FI = new double[m];
			for(int i=0; i<m;i++){
				FI[i]=2*XI[i];
				//FI[i]=XI[i]*XI[i]; //test functions (uncomment to use)
				FI[i]=Math.pow(Math.E, XI[i]);
				//FI[i]=Math.pow(Math.E, Math.sin(2*XI[i]))*Math.cos(2*XI[i]);
				//FI[i]=Math.tanh(XI[i]);
				//FI[i]=XI[i]*Math.cos(2*Math.PI*XI[i]);
				//FI[i]=XI[i] + (1/XI[i]);
			}
			
			double[] XI3 = new double[2*m];
			for(int i=0;i<m; i++){
				XI3[2*i]= XI[i] - thirdHm;
				XI3[2*i+1]= XI[i] + thirdHm;
			}
			
			double[] FInew = new double[2*m];
			for(int i=0; i<2*m; i++){
				FInew[i]=2*XI3[i];
				//FInew[i]=XI3[i]*XI3[i];
				FInew[i]=Math.pow(Math.E, XI3[i]);
				//FInew[i]=Math.pow(Math.E, Math.sin(2*XI3[i]))*Math.cos(2*XI3[i]);
				//FInew[i]=Math.tanh(XI3[i]);
				//FInew[i]=XI3[i]*Math.cos(2*Math.PI*XI3[i]);
				//FInew[i]=XI3[i] + (1/XI3[i]);
			}
			double Im = CompMidpointRule(Hm, FI);
			double I3m = MdptRuleFineMesh(thirdHm, Im, FInew);
			System.out.println(Im);
			System.out.println(I3m);
			System.out.println();
			
			error = FineError(Im, I3m);
			exacterror =(exactIntegralValue-Im);
			
			//fw1.write(thirdHm+" "+error+"\n");
			fw2.write(thirdHm+" "+exacterror+"\n");
			//fw3.write(I3m+"\n");
			
			m=3*m;
		}
		
		//fw1.close();
		fw2.close();
		//fw3.close();
	}
	
	public static double CompMidpointRule(double Hm, double[] FI){
		double sum=0;
		for(int i=0; i<FI.length;i++){
			sum = sum + FI[i];
		}
		return Hm*sum;
	}
	
	public static double MdptRuleFineMesh(double thirdHm, double Im, double[]FInew){
		double sum =0; 
		for(int i=0; i<FInew.length; i++)
			sum = sum + FInew[i];
			
		return Im/3 + thirdHm*sum;
	}
	
	public static double CoarseError(double Im, double I3m){return 9/8 * (I3m-Im);}
	public static double FineError(double Im, double I3m){return (I3m - Im)/8;}
}