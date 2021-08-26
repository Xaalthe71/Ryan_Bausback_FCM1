import java.lang.Math;
import java.util.*;
import java.io.*;

public class explicitRK{
	public static void main(String[] args)throws Exception{
		//convergence
		/*FileWriter fwconv1 = new FileWriter("FEconvIVP3.txt");
		FileWriter fwconv2 = new FileWriter("RK2convIVP3MDPT.txt");
		FileWriter fwconv2_1 = new FileWriter("RK2convIVP3TRAP.txt");
		FileWriter fwconv3 = new FileWriter("RK3convIVP3.txt");
		FileWriter fwconv3_1 = new FileWriter("RK3convIVP3m1.txt");
		FileWriter fwconv4 = new FileWriter("RK4convIVP3.txt");
		
		/*FileWriter fwconv1_1 = new FileWriter("FEconvIVP2_y2.txt");
		FileWriter fwconv2_2 = new FileWriter("RK2convIVP2_y2MDPT.txt");
		FileWriter fwconv2_3 = new FileWriter("RK2convIVP2_y2TRAP.txt");
		FileWriter fwconv3_3 = new FileWriter("RK3convIVP2_y2.txt");
		FileWriter fwconv4_4 = new FileWriter("RK4convIVP2_y2.txt");*/
		
		double tn_1=0;
		
		//changeable hardcoded values
		double tn = 5;
		double N = 50; //stages in between each time step.
		double y0 = 1; //for 2nd IVP, y0 will be the initial condition for y2, since y1(0)=0 always. 
		double lambda =-3; //for the 2nd IVP, omega will be named lambda
		int IVP =1;		// 1 is growth, 3 is driven system, 2 is top half of oscillator, 4 is bottom half of oscillator
		
		//for(double y0=-10; y0<=10; y0=y0+0.1){ //change indicator initial based on variable iterating
			double h =(tn-tn_1)/N; //stepsize
			
			
			double[][] ynEuler = ForwardEuler(N, h, y0, tn_1, lambda, IVP);
			double[][] ynRK2MDPT = TwoStageRK(N, h, y0, 0.5, tn_1, lambda, IVP);
			double[][] ynRK2TRAP = TwoStageRK(N, h, y0, 1, tn_1, lambda, IVP);
			double[][] ynRK3 = ThreeStageRK(N, h, y0, 0.5, tn_1, lambda, IVP);
			double[][] ynRK3m1 = ThreeStageRK(N, h, y0, 1, tn_1, lambda, IVP);
			double[][] ynRK4 = FourStageRK(N, h, y0, tn_1, lambda, IVP);
			
			double exactIVP1 =y0*Math.pow(Math.E, lambda*tn);
			double exactIVP2_y1 = y0*Math.sin(lambda*tn);
			double exactIVP2_y2 = y0*Math.cos(lambda*tn);
			double exactIVP3 = (y0 - 0)*Math.pow(Math.E, lambda*tn) + tn*tn;
			//double exactIVP3 = (y0 - 0)*Math.pow(Math.E, lambda*tn) + Math.sin(tn);
			//double exactIVP3 = Math.cos(tn);
			
			if(N==50){ //only need for 1 value
				//graphs against the t-values
				FileWriter fw1 = new FileWriter("FEsolutionIVP2.txt");
				for(int i=0; i<=N; i++)
					fw1.write(tn_1+i*h+" "+ynEuler[i][0]+"\n");
				fw1.close();
				
				FileWriter fw2 = new FileWriter("RK2solutionIVP3MDPT.txt");
				for(int i=0; i<=N; i++)
					fw2.write(tn_1+i*h+" "+ynRK2MDPT[i][0]+"\n");
				fw2.close();
				
				FileWriter fw2_1 = new FileWriter("RK2solutionIVP3TRAP.txt");
				for(int i=0; i<=N; i++)
					fw2_1.write(tn_1+i*h+" "+ynRK2TRAP[i][0]+"\n");
				fw2_1.close();
				
				FileWriter fw3 = new FileWriter("RK3solutionIVP2.txt");
				for(int i=0; i<=N; i++)
					fw3.write(tn_1+i*h+" "+ynRK3[i][0]+"\n");
				fw3.close();
				
				FileWriter fw4 = new FileWriter("RK4solutionIVP2.txt");
				for(int i=0; i<=N; i++)
					fw4.write(tn_1+i*h+" "+ynRK4[i][0]+"\n");
				fw4.close();
				
				if(IVP==2){//graphs against the t-values
					FileWriter fw1_2 = new FileWriter("FEsolutionIVP2y2.txt");
					for(int i=0; i<=N; i++)
						fw1_2.write(tn_1+i*h+" "+ynEuler[i][1]+"\n");
					fw1_2.close();
					
					FileWriter fw2_2 = new FileWriter("RK2solutionIVP2MDPTy2.txt");
					for(int i=0; i<=N; i++)
						fw2_2.write(tn_1+i*h+" "+ynRK2MDPT[i][1]+"\n");
					fw2_2.close();
					
					FileWriter fw2_3 = new FileWriter("RK2solutionIVP2TRAPy2.txt");
					for(int i=0; i<=N; i++)
						fw2_3.write(tn_1+i*h+" "+ynRK2TRAP[i][1]+"\n");
					fw2_3.close();
					
					FileWriter fw3_2 = new FileWriter("RK3solutionIVP2y2.txt");
					for(int i=0; i<=N; i++)
						fw3_2.write(tn_1+i*h+" "+ynRK3[i][1]+"\n");
					fw3_2.close();
					
					FileWriter fw4_2 = new FileWriter("RK4solutionIVP2y2.txt");
					for(int i=0; i<=N; i++)
						fw4_2.write(tn_1+i*h+" "+ynRK4[i][1]+"\n");
					fw4_2.close();
				}
			}
			
			
		//Error computation
			/*FileWriter fwErr = new FileWriter("GlobalErrroIVP1.txt");
			//IVP1
			double exactIVP1 =y0*Math.pow(Math.E, lambda*tn);
			fwErr.write(Math.abs(ynEuler[(int)N]-exactIVP1)+"\n"+Math.abs(ynRK2MDPT[(int)N]-exactIVP1)+"\n"+Math.abs(ynRK2TRAP[(int)N]-exactIVP1)+"\n"+Math.abs(ynRK3[(int)N]-exactIVP1)+"\n"+Math.abs(ynRK4[(int)N]-exactIVP1));
			fwErr.close();*/
			
		
		//convegence
			/*fwconv1.write(lambda+" "+absStability(h, lambda, 1)+"\n");
			fwconv2.write(lambda+" "+absStability(h, lambda, 2)+"\n");
			fwconv2_1.write(lambda+" "+absStability(h, lambda, 2)+"\n");
			fwconv3.write(lambda+" "+absStability(h, lambda, 3)+"\n");
			fwconv4.write(lambda+" "+absStability(h, lambda, 4)+"\n");
			
			fwconv1_1.write(lambda+" "+absStability(h, lambda, 1)+"\n");
			fwconv2_2.write(lambda+" "+absStability(h, lambda, 2)+"\n");
			fwconv2_3.write(lambda+" "+absStability(h, lambda, 2)+"\n");
			fwconv3_3.write(lambda+" "+absStability(h, lambda, 3)+"\n");
			fwconv4_4.write(lambda+" "+absStability(h, lambda, 4)+"\n");*/
		
			/*fwconv1.write(y0+" "+Math.abs(ynEuler[(int)N][0]-exactIVP3)+"\n");
			fwconv2.write(y0+" "+Math.abs(ynRK2MDPT[(int)N][0]-exactIVP3)+"\n");
			fwconv2_1.write(y0+" "+Math.abs(ynRK2TRAP[(int)N][0]-exactIVP3)+"\n");
			fwconv3.write(y0+" "+Math.abs(ynRK3[(int)N][0]-exactIVP3)+"\n");
			fwconv3_1.write(y0+" "+Math.abs(ynRK3m1[(int)N][0]-exactIVP3)+"\n");
			fwconv4.write(y0+" "+Math.abs(ynRK4[(int)N][0]-exactIVP3)+"\n");
			
			/*fwconv1_1.write(y0+" "+Math.abs(ynEuler[(int)N][1]-exactIVP2_y2)+"\n");
			fwconv2_2.write(y0+" "+Math.abs(ynRK2MDPT[(int)N][1]-exactIVP2_y2)+"\n");
			fwconv2_3.write(y0+" "+Math.abs(ynRK2TRAP[(int)N][1]-exactIVP2_y2)+"\n");
			fwconv3_3.write(y0+" "+Math.abs(ynRK3[(int)N][1]-exactIVP2_y2)+"\n");
			fwconv4_4.write(y0+" "+Math.abs(ynRK4[(int)N][1]-exactIVP2_y2)+"\n");*/
		//}
		
		/*fwconv1.close();
		fwconv2.close();
		fwconv2_1.close();
		fwconv3.close();
		fwconv3_1.close();
		fwconv4.close();*/
		
		/*fwconv1_1.close();
		fwconv2_2.close();
		fwconv2_3.close();
		fwconv3_3.close();
		fwconv4_4.close();*/
	}
	
//Explicit RK methods
	public static double[][] ForwardEuler(double N, double h, double y0, double tn_1, double lambda, int IVP ){
		double[][] yn = new double[(int) N+1][2]; //IVP 1 & 3 will only use column 1.
		
		if(IVP ==1){
			yn[0][0]=y0;
			
			for(int i=1; i<=N; i++)
				yn[i][0] = yn[i-1][0] + h*f1(yn[i-1][0], tn_1+(i-1)*h, lambda); 
		}
		else if(IVP ==2){
			yn[0][0]=0;
			yn[0][1]=y0;
			
			for(int i=1; i<=N; i++){
				yn[i][0] = yn[i-1][0] + h*f2_1(yn[i-1][1], tn_1+(i-1)*h, lambda);  //y1
				yn[i][1] = yn[i-1][1] + h*f2_2(yn[i-1][0], tn_1+(i-1)*h, lambda);  //y2
			}
		}
		else if(IVP ==3){
			yn[0][0]=y0;
			
			for(int i=1; i<=N; i++)
				yn[i][0] = yn[i-1][0] + h*f3(yn[i-1][0], tn_1+(i-1)*h, lambda); 
		}
		return yn;
	}
	
	public static double[][] TwoStageRK(double N, double h, double y0, double mu, double tn_1, double lambda, int IVP){
		double[][] yn = new double[(int)N+1][2];

		if(IVP ==1){
			yn[0][0]=y0;
			for(int i=1; i<=N; i++){
				yn[i][0]= yn[i-1][0];
				double functionVal = f1(yn[i-1][0], tn_1, lambda); 
				yn[i][0]= yn[i][0]+ h*(1.0-(1.0/(2.0*mu)))*functionVal;
				functionVal = f1(yn[i-1][0]+h*mu*functionVal, tn_1+(i-1)*h+mu*h, lambda); 
				yn[i][0]= yn[i][0]+h*(1.0/(2.0*mu))*functionVal;
			}
		}
		else if(IVP ==2){
			yn[0][0]=0;
			yn[0][1]=y0;
			
			for(int i=1; i<=N; i++){ //direct application of the output of computing the matrix-vector products by hand
				yn[i][0] =yn[i-1][0]+ h*(1.0-(1.0/(2.0*mu)))*(lambda*yn[i-1][1])+ h*(1.0/(2.0*mu))*(lambda*(yn[i-1][1]-h*mu*lambda*yn[i-1][0]));
				yn[i][1] =yn[i-1][1]+ h*(1.0-(1.0/(2.0*mu)))*(-1*lambda*yn[i-1][0])+ h*(1.0/(2.0*mu))*(-1*lambda*(yn[i-1][0]+h*mu*lambda*yn[i-1][1]));
			}
		}
		else if(IVP==3){
			yn[0][0]=y0;
			for(int i=1; i<=N; i++){
				yn[i][0]= yn[i-1][0];
				double functionVal = f3(yn[i-1][0], tn_1+(i-1)*h, lambda); 
				yn[i][0]= yn[i][0]+ h*(1.0-(1.0/(2.0*mu)))*functionVal;
				functionVal = f3(yn[i-1][0]+h*mu*functionVal, tn_1+(i-1)*h+mu*h, lambda);
				yn[i][0]= yn[i][0]+h*(1.0/(2.0*mu))*functionVal;
			}
		}
		return yn;
	}
	
	public static double[][] ThreeStageRK(double N, double h, double y0, double mu, double tn_1, double lambda, int IVP){
		double[][] yn = new double[(int)N+1][2];

		if(IVP==1){
			yn[0][0]=y0;
			
			for(int i=1; i<=N; i++){
				yn[i][0]=yn[i-1][0];
				//f1 hat
				double fVal1 = f1(yn[i-1][0], tn_1+(i-1)*h, lambda); 
				yn[i][0]= yn[i][0]+ h*0.25*fVal1;
				//f2 hat
				double fVal2 = f1(yn[i-1][0]+h*(2.0/3.0)*fVal1, tn_1+(i-1)*h+(2.0/3.0)*h, lambda); 
				yn[i][0]= yn[i][0]+h*(0.75-mu)*fVal2;
				//f3 hat
				double fVal3 = f1(yn[i-1][0]+h*((2.0/3.0)-(1.0/(4.0*mu)))*fVal1 +h*(1.0/(4.0*mu))*fVal2, tn_1+(i-1)*h+(2.0/3.0)*h, lambda);
				yn[i][0] = yn[i][0] + h*mu*fVal3;
			}
		}
		else if(IVP==2){
			yn[0][0]=0;
			yn[0][1]=y0;
			
			for(int i=1; i<=N; i++){
				yn[i][0]=yn[i-1][0];
				yn[i][1]=yn[i-1][1];
				
				//f1 hat
				double fVal1 = f2_1(yn[i-1][1], tn_1+(i-1)*h, lambda); 
				double fVal1_2 = f2_2(yn[i-1][0], tn_1+(i-1)*h, lambda); 
				
				yn[i][0]= yn[i][0]+ h*0.25*fVal1;
				yn[i][1]= yn[i][1]+ h*0.25*fVal1_2;
				
				//f2 hat
				double fVal2 = f2_1(yn[i-1][1]+h*(2.0/3.0)*fVal1_2, tn_1+(i-1)*h+(2.0/3.0)*h, lambda); 
				double fVal2_2 = f2_2(yn[i-1][0]+h*(2.0/3.0)*fVal1, tn_1+(i-1)*h+(2.0/3.0)*h, lambda); 
				
				yn[i][0]= yn[i][0]+h*(0.75-mu)*fVal2;
				yn[i][1]= yn[i][1]+h*(0.75-mu)*fVal2_2;
				
				//f3 hat
				double fVal3 = f2_1(yn[i-1][1]+h*((2.0/3.0)-(1.0/(4.0*mu)))*fVal1_2 +h*(1.0/(4.0*mu))*fVal2_2, tn_1+(i-1)*h+(2.0/3.0)*h, lambda);
				double fVal3_2 = f2_2(yn[i-1][0]+h*((2.0/3.0)-(1.0/(4.0*mu)))*fVal1 +h*(1.0/(4.0*mu))*fVal2, tn_1+(i-1)*h+(2.0/3.0)*h, lambda);
				
				yn[i][0] = yn[i][0] + h*mu*fVal3;
				yn[i][1] = yn[i][1] + h*mu*fVal3_2;
			}
		}
		else if(IVP==3){
			yn[0][0]=y0;
			
			for(int i=1; i<=N; i++){
				yn[i][0]=yn[i-1][0];
				//f1 hat
				double fVal1 = f3(yn[i-1][0], tn_1+(i-1)*h, lambda); 
				yn[i][0]= yn[i][0]+ h*0.25*fVal1;
				//f2 hat
				double fVal2 = f3(yn[i-1][0]+h*(2.0/3.0)*fVal1, tn_1+(i-1)*h+(2.0/3.0)*h, lambda); 
				yn[i][0]= yn[i][0]+h*(0.75-mu)*fVal2;
				//f3 hat
				double fVal3 = f3(yn[i-1][0]+h*((2.0/3.0)-(1.0/(4.0*mu)))*fVal1 +h*(1.0/(4.0*mu))*fVal2, tn_1+(i-1)*h+(2.0/3.0)*h, lambda);
				yn[i][0] = yn[i][0] + h*mu*fVal3;
			}
		}
		return yn;
	}
	
	public static double[][] FourStageRK(double N, double h, double y0, double tn_1, double lambda, int IVP){
		double[][] yn = new double[(int)N+1][2];

		if(IVP == 1){
			yn[0][0]=y0;
			
			for(int i=1; i<=N; i++){
				yn[i][0]=yn[i-1][0];
				//f1 hat
				double fVal = f1(yn[i-1][0], tn_1+(i-1)*h, lambda); 
				yn[i][0]= yn[i][0]+ (1.0/6.0)*h*fVal;
				//f2 hat
				fVal = f1(yn[i-1][0]+h*0.5*fVal, ((tn_1+i*h)+tn_1+(i-1)*h)/2, lambda); 
				yn[i][0]= yn[i][0]+ h*(1.0/3.0)*fVal; 
				//f3 hat
				fVal = f1(yn[i-1][0]+0.5*h*fVal, ((tn_1+i*h)+tn_1+(i-1)*h)/2, lambda); 
				yn[i][0] = yn[i][0] +h*(1.0/3.0)*fVal;
				//f4 hat
				fVal= f1(yn[i-1][0]+h*fVal, (tn_1+i*h), lambda); 
				yn[i][0] = yn[i][0] + h*(1.0/6.0)*fVal;
			}
		}
		else if(IVP==2){
			yn[0][0]=0;
			yn[0][1]=y0;
			
			for(int i=1; i<=N; i++){
			//y1
				yn[i][0]=yn[i-1][0];
				yn[i][1]=yn[i-1][1];
				
				//f1 hat
				double fVal = f2_1(yn[i-1][1], tn_1+(i-1)*h, lambda); 
				double fVal2 = f2_2(yn[i-1][0], tn_1+(i-1)*h, lambda);
				
				yn[i][0]= yn[i][0]+ (1.0/6.0)*h*fVal;
				yn[i][1]= yn[i][1]+ (1.0/6.0)*h*fVal2;
				
				//f2 hat
				fVal = f2_1(yn[i-1][1]+h*0.5*fVal2, ((tn_1+N*h)+tn_1+(i-1)*h)/2, lambda); 
				fVal2 = f2_2(yn[i-1][0]+h*0.5*fVal, ((tn_1+N*h)+tn_1+(i-1)*h)/2, lambda);
				
				yn[i][0]= yn[i][0]+ h*(1.0/3.0)*fVal;
				yn[i][1]= yn[i][1]+ h*(1.0/3.0)*fVal2; 
				
				//f3 hat
				fVal = f2_1(yn[i-1][1]+0.5*h*fVal2, ((tn_1+i*h)+tn_1+(i-1)*h)/2, lambda); 
				fVal2 = f2_2(yn[i-1][0]+0.5*h*fVal, ((tn_1+i*h)+tn_1+(i-1)*h)/2, lambda); 
				
				yn[i][0] = yn[i][0] +h*(1.0/3.0)*fVal;
				yn[i][1] = yn[i][1] +h*(1.0/3.0)*fVal2;
				
				//f4 hat
				fVal= f2_1(yn[i-1][1]+h*fVal2, (tn_1+i*h), lambda);
				fVal2= f2_2(yn[i-1][0]+h*fVal, (tn_1+i*h), lambda); 
				
				yn[i][0] = yn[i][0] + h*(1.0/6.0)*fVal; 
				yn[i][1] = yn[i][1] + h*(1.0/6.0)*fVal2;
			}
		}
		else if(IVP ==3){
			yn[0][0]=y0;
			
			for(int i=1; i<=N; i++){
				yn[i][0]=yn[i-1][0];
				//f1 hat
				double fVal = f3(yn[i-1][0], tn_1+(i-1)*h, lambda); 
				yn[i][0]= yn[i][0]+ (1.0/6.0)*h*fVal;
				//f2 hat
				fVal = f3(yn[i-1][0]+h*0.5*fVal, ((tn_1+i*h)+(tn_1+(i-1)*h))/2, lambda); 
				yn[i][0]= yn[i][0]+ h*(1.0/3.0)*fVal; 
				//f3 hat
				fVal = f3(yn[i-1][0]+0.5*h*fVal, ((tn_1+i*h)+(tn_1+(i-1)*h))/2, lambda); 
				yn[i][0] = yn[i][0] +h*(1.0/3.0)*fVal;
				//f4 hat
				fVal= f3(yn[i-1][0]+h*fVal, (tn_1+i*h), lambda); 
				yn[i][0] = yn[i][0] + h*(1.0/6.0)*fVal;
			}
		}
		return yn;
	}

//IVP functions
	public static double f1(double y, double tn_1, double lambda){return lambda*tn_1*y;}
	public static double f2_1(double y, double tn_1, double omega){return omega*y;}
	public static double f2_2(double y, double tn_1, double omega){return -1*omega*y;}
	public static double f3(double y, double tn_1,double lambda){return lambda*y + Math.pow(Math.E, tn_1);}				//change F(t) and F'(t) here
	
	
//Absolute Stability 
	public static double absStability(double h, double lambda, int order){
		double Rz =1 + h*lambda;
		if(order>=2)
			Rz = Rz + Math.pow(h*lambda, 2)/2.0;
		if(order >=3)
			Rz=  Rz + Math.pow(h*lambda, 3)/6.0;
		if(order == 4)
			Rz = Rz + Math.pow(h*lambda, 4)/24.0;
		return Math.abs(Rz);
	}
	
}