package signal;

public class sinuTesting {

	public static void main(String[] args) {
				
		Sinusoid h = new Sinusoid(1, 2, 0, 0, 256);
		Sinusoid g2 = new Sinusoid(1, 5, 0, 0, 256);
		Sinusoid g3 = new Sinusoid(1, 8, 0, 0, 256);
		h.computeWave();
		g2.computeWave();
		g3.computeWave();
		h.addData(g2);
		h.addData(g3);
		ComplexNumber[] t1 = h.FFT();
		
		double[] p = FourierTransform.computePSD(t1);
		FourierTransform.graphPSD(p, 4);
		h.showImg();
		ComplexNumber[] H = h.FFT();
		System.out.println(g3);
		System.out.println(h);
		ComplexNumber.printArray(H);
		double[] ampH = new double[H.length];
		for (int i = 0; i < ampH.length; i++) {
			ampH[i] = H[i].mod();
		}
		System.out.println("Amplitude");
		for (double amp : ampH)
			System.out.println(amp);
		
		System.out.println("\n\nImaginary Compliment then multiply H and H* (PSD?)");
		ComplexNumber[] conjH = new ComplexNumber[H.length];
		for (int k = 0; k < conjH.length; k++) {
			conjH[k] = H[k].conjugate();
		}
		ComplexNumber[] HHstar = new ComplexNumber[H.length];
		for (int hh = 0; hh < HHstar.length; hh++) {
			HHstar[hh] = H[hh].multiply(conjH[hh]);
		}
		ComplexNumber.printArray(HHstar);
		
		System.out.println("\nFiltration test");
		
		ComplexNumber.printArray(conjH);
		ComplexNumber[] HI = FourierTransform.highPassFilter(conjH, 8);
		ComplexNumber.printArray(HI);
		
		ComplexNumber[] shouldBeFeqEight = FourierTransform.FFT(HI, -1);
		double[] realsbfe = ComplexNumber.getRealVals(shouldBeFeqEight);
		ComplexNumber.printArray(shouldBeFeqEight);
		for (int ore = 0; ore < realsbfe.length; ore++) {
			System.out.println(ore + " " + realsbfe[ore]);
		}
	}

}
