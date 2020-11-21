import java.io.*;
import java.lang.Math;

public class rationalsumgenerator{
	public static void main(String[] args)throws Exception{
		FileWriter fw = new FileWriter("input2.txt");
		double value =0;
		
		int numElements = 10000; //change this value to get different lengths of sums
		fw.write(numElements + "\n");
		for(int i=1; i<=numElements; i++){
			value = Math.pow(-1, i) * (1/Math.pow(i, 2));
			fw.write(value + "\n");
		}
		fw.close();
	}
}