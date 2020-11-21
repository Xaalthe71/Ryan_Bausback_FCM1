import java.io.*;
import java.lang.Math;

public class SinElementsGenerator1{
	public static void main(String[] args)throws Exception{
		FileWriter fw = new FileWriter("input1.txt");
		double value =0;
		
		int numElements = (int)Math.pow(2, 20); //change this value to get different length sums
		fw.write(numElements + "\n");
		for(int i=1; i<=numElements; i++){
			value = i * Math.sin(1/i + i/2);
			fw.write(value + "\n");
		}
		fw.close();
	}
}