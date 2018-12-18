package signal;

import java.io.File;
import java.util.Arrays;

import signal.FourierTransform.VectorCondition;

public class testing {

	public static void main(String[] args) {
		double[] data = new double[16];
		data[0] = 26160.0;
		data[1] = 19011.0;
		data[2] = 18757.0;
		data[3] = 18405.0;
		data[4] = 17888.0;
		data[5] = 14720.0;
		data[6] = 14285.0;
		data[7] = 17018.0;
		data[8] = 18014.0;
		data[9] = 17119.0;
		data[10] = 16400.0;
		data[11] = 17497.0;
		data[12] = 17846.0;
		data[13] = 15700.0;
		data[14] = 17636.0;
		data[15] = 17181.0;

		ComplexNumber[] test = FourierTransform.makeComplex(data);
		for (ComplexNumber c : test) {
			System.out.println(c);
		}
		System.out.println();

		System.out.println("Test run ");
		ComplexNumber[] ans = FourierTransform.FFT(test, 1);
		for (ComplexNumber cn : ans) {
			System.out.println(cn);
		}

		System.out.println("\nSunspot trial");
		File f = new File("SUNSPOTdataAnnually.txt");
		ComplexNumber[] sunTest = FourierTransform.readInData(f);

		for (ComplexNumber c : sunTest) {
			System.out.println(c);
		}
		System.out.println();
		
		//sunTest = Arrays.copyOf(sunTest, 256);
		System.out.println(sunTest.length + " " + Math.sqrt(sunTest.length));
		
		System.out.println("\nFFT trail run on sunspot data");
		ComplexNumber[] ans2 = FourierTransform.FFT(sunTest, 1, VectorCondition.CUT_OFF_DATA);
		for (ComplexNumber c : ans2) {
			System.out.println(c);
		}
		System.out.println();
		
		
		
	}

}
