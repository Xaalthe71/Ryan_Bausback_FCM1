import java.lang.Math;
import java.security.SecureRandom;
import java.io.*;
import java.util.Scanner;

public class StabilityandConditioning1{
	public static void main(String[] args)throws Exception{
		double floatmachineEpsilon = 0.0000000596; //stored as double since java automatically assumes decimals are double.

		File inFile = new File("input3.txt"); //change number in input_.txt and output_.txt to get different sum computation
		Scanner in = new Scanner(inFile);
		
		double[] XI_DP = new double[Integer.parseInt(in.nextLine())];
		for(int i=0; i< XI_DP.length;i++)
			XI_DP[i]=Double.parseDouble(in.nextLine());
		
	//used to compute the Leibniz Sum to test correctness.
		/*double[] XI_DP = new double[(int)Math.pow(2,20)];
		for(int i=0; i<XI_DP.length; i++)
			XI_DP[i]=Math.pow(-1,i)/(2*i+1);*/
		
		float[] XI_SP = new float[XI_DP.length];
		for(int i=0;i<XI_DP.length;i++)
			XI_SP[i] = (float) XI_DP[i];
	
	//summation algoritm computation
		float S_SPfan = SinglePrecisionFanIn(XI_SP);
		float S_SPaccumS = SinglePrecisionAccumulation(XI_SP);
		float S_SPaccumD = DoublePrecisionAccumulation(XI_SP);
		double exactSum = ExactSum(XI_DP);
		
		double conditionNumber = ConditionNumber(XI_DP);
		
	//stability bound computation
		double accumstabilitybound = floatmachineEpsilon*XI_DP.length*conditionNumber;
		
		boolean sizechecker=false;  
		int k=0;									//for use in determining p(n) for the fan-in stability bound
		while(!sizechecker){
			if(XI_DP.length <= Math.pow(2, k))
				sizechecker=true;
			else
				k++;
		}
		double faninstabilitybound = floatmachineEpsilon*k*conditionNumber;
		
	//conditioning bound computation
		double[] XI_DPperturbed = new double[XI_DP.length];
		double perturbation = 0.000001; 
		for(int i=0; i<XI_DP.length; i++)
			XI_DPperturbed[i]=XI_DP[i]*(1+perturbation); 
		
		double normD =0;
		for(int i=0; i<XI_DP.length; i++)
			normD=normD + Math.abs(XI_DP[i]);
		
		double normE =0;
		for(int i=0; i<XI_DP.length; i++)
			normE = normE + Math.abs(XI_DP[i]*perturbation); //since perturbation is relative, normE must be relative.
		
		double conditioningBound = (conditionNumber * normE)/normD;
		
	//data output
		FileWriter fw = new FileWriter("output3.txt");
		fw.write("Results: \n\nNumber of Elements added: "+XI_DP.length+"\n\nExact Sum: \n"+exactSum+"\n\nSingle Precision Fan In:\n"+S_SPfan+"\n\nError: \n"+AnalyzeError(exactSum,S_SPfan)+"\n\n\n");
		fw.write("Single Precision Accumulation:\n"+S_SPaccumS+"\n\nError:\n"+AnalyzeError(exactSum, S_SPaccumS)+"\n\n");
		fw.write("Double Precision Accumulation:\n"+S_SPaccumD+"\n\nError:\n"+AnalyzeError(exactSum, S_SPaccumD)+"\n\n");
		fw.write("Condition Number:\n"+conditionNumber+"\n\nSingle Accumulation Stability Bound: \n"+accumstabilitybound+"\n\n");
		
		if(AnalyzeError(exactSum, S_SPaccumS)<accumstabilitybound)
			fw.write("Simple accumulation is weakly stable\n\n");
		
		fw.write("Binary Fan-In Tree Stability Bound: \n"+faninstabilitybound+"\n\n");
		
		if(AnalyzeError(exactSum, S_SPfan)<faninstabilitybound)
			fw.write("Binary Fan In is weakly stable.\n\n");
		
		fw.write("Conditioning Error: \n"+AnalyzeError(exactSum, ExactSum(XI_DPperturbed))+"\n\nConditioning Bound: \n"+conditioningBound+"\n\n");
		if(AnalyzeError(exactSum, ExactSum(XI_DPperturbed))<conditioningBound)
			fw.write("Less than conditioning bound\n\n");
		
		//fw.write("True sum: "+(Math.PI/4+"\n\n");
		//fw.write("Relative Error from True Sum: "+AnalyzeError((Math.PI/4),exactSum));
		fw.close();
	}
	
//Summation Algorithms
	public static float SinglePrecisionAccumulation(float[] XI_SP)throws Exception{
		float S_SP=0;
		FileWriter fw = new FileWriter("SinglePrecisionAccum.txt"); //for outputting raw data to create graphs
		for(int i=0;i<XI_SP.length;i++){
			S_SP=S_SP+XI_SP[i];
			fw.write(S_SP+"\n");
		}
		fw.close();
		return S_SP;
	}

	public static float DoublePrecisionAccumulation(float[] XI_SP){
		double S_DP=0;
		for(int i=0; i<XI_SP.length; i++)
			S_DP=S_DP+XI_SP[i];
		float S_SP = (float) S_DP;
		return S_SP;
	}

	public static float SinglePrecisionFanIn(float[] XI_SP){
		boolean sizechecker =false;
		int i=0;
		while(!sizechecker){
			if(XI_SP.length <= Math.pow(2, i))
				sizechecker=true;
			else
				i++;
		}
		float[]S_SPprelim = new float[(int)Math.pow(2,i)]; //since i is int, no loss of data from typecast
		
		for(int j=0;j<S_SPprelim.length;j++) //compute only on 2^i indices with 0s in empty slots for simplicity
			S_SPprelim[j]=0;
		
		for(int j=0;j<XI_SP.length;j++)
			S_SPprelim[j]=XI_SP[j];
		
		int indicator = S_SPprelim.length;
		for(int k=0; k<i; k++){
			for(int j=0;indicator/2-j>0;j++){ //add first to last and save in first position, progressing inwards toward center
				S_SPprelim[j]=S_SPprelim[j]+S_SPprelim[indicator-1-j];
			}
			indicator=indicator/2;
		}
		return S_SPprelim[0];
	}
	
	public static double ExactSum(double[] XI_DP)throws Exception{
		double S_DP=0;
		FileWriter fw = new FileWriter("ExactSum.txt"); //output to file to graph sum's progression to justify choosing it for experiments
		for(int i=0; i<XI_DP.length;i++){
			S_DP=S_DP+XI_DP[i];
			fw.write(S_DP+"\n");
		}
		fw.close();
		return S_DP;
	}
	
//Analysis Tools
	public static double AnalyzeError(double S_DP, double S_SP){
		return (Math.abs(S_SP - S_DP) / Math.abs(S_DP));
	}
	
	public static double ConditionNumber(double[] XI_DP)throws Exception{
		double numerator =0;
		for(int i=0; i<XI_DP.length; i++)
			numerator = numerator + Math.abs(XI_DP[i]);
		
		double denominator = Math.abs(ExactSum(XI_DP));
		
		return numerator/denominator;
	}
}
