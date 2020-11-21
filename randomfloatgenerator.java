import java.io.*;
import java.lang.Math;
import java.security.SecureRandom;

public class randomfloatgenerator{
	public static void main(String[] args)throws Exception{
		FileWriter fw = new FileWriter("input3.txt");
		double value =0;
		SecureRandom rn = new SecureRandom();
		
		int numElements = (int) Math.pow(2, 20); //change this value to get different lengths of sums
		fw.write(numElements + "\n");
		for(int i=1; i<=numElements; i++){
			value = Math.pow(-1, rn.nextInt(2)) * (rn.nextDouble()+rn.nextInt(3));
			fw.write(value + "\n");
		}
		fw.close();
	}
}