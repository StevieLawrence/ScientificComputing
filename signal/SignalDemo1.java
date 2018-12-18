package signal;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Scanner;

import signal.FourierTransform.FILTER;

public class SignalDemo1 {

	public static void main(String[] args) {
		double phase = 0, bias = 0, amp = 0, f = 0;
		int N = 512;

		Sinusoid f3 = new Sinusoid(amp, f, phase, bias, N);
		f3.computeWave();
		// make f3 with 3 terms
		int numberOfTerms = 3;
		for (int k = 1; k <= numberOfTerms; k++) {
			f = 2 * k - 1;
			amp = 1 / f;
			Sinusoid term = new Sinusoid(amp, f, phase, bias, N);
			term.computeWave();
			f3.addData(term);
		}
		// f3.showImg("f3");

		// make g3 with 3 terms
		phase = 0;
		bias = 0;
		amp = 0;
		f = 0;
		Sinusoid g3 = new Sinusoid(amp, f, phase, bias, N);
		g3.computeWave();
		for (int k = 1; k <= numberOfTerms; k++) {
			f = 2 * k;
			amp = 1 / f;
			Sinusoid term = new Sinusoid(amp, f, phase, bias, N);
			term.computeWave();
			g3.addData(term);
		}
		// g3.showImg("g3");

		// make f10 with 10 terms
		numberOfTerms = 10;
		phase = 0;
		bias = 0;
		amp = 0;
		f = 0;
		Sinusoid f10 = new Sinusoid(amp, f, phase, bias, N);
		f10.computeWave();
		for (int k = 1; k <= numberOfTerms; k++) {
			f = 2 * k - 1;
			amp = 1 / f;
			Sinusoid term = new Sinusoid(amp, f, phase, bias, N);
			term.computeWave();
			f10.addData(term);
		}
		// f10.showImg("f10");

		// make g10 with 10 terms
		phase = 0;
		bias = 0;
		amp = 0;
		f = 0;
		Sinusoid g10 = new Sinusoid(amp, f, phase, bias, N);
		g10.computeWave();
		for (int k = 1; k <= numberOfTerms; k++) {
			f = 2 * k;
			amp = 1 / f;
			Sinusoid term = new Sinusoid(amp, f, phase, bias, N);
			term.computeWave();
			g10.addData(term);
		}
		// g10.showImg("g10");

		// make f50 with 50 terms
		numberOfTerms = 50;
		phase = 0;
		bias = 0;
		amp = 0;
		f = 0;
		Sinusoid f50 = new Sinusoid(amp, f, phase, bias, N);
		f50.computeWave();
		for (int k = 1; k <= numberOfTerms; k++) {
			f = 2 * k - 1;
			amp = 1 / f;
			Sinusoid term = new Sinusoid(amp, f, phase, bias, N);
			term.computeWave();
			f50.addData(term);
		}
		// f50.showImg("f50");

		// make g50 with 50 terms
		phase = 0;
		bias = 0;
		amp = 0;
		f = 0;
		Sinusoid g50 = new Sinusoid(amp, f, phase, bias, N);
		g50.computeWave();
		for (int k = 1; k <= numberOfTerms; k++) {
			f = 2 * k;
			amp = 1 / f;
			Sinusoid term = new Sinusoid(amp, f, phase, bias, N);
			term.computeWave();
			g50.addData(term);
		}
		// g50.showImg("g50");

		// really big f
		numberOfTerms = 10000;
		phase = 0;
		bias = 0;
		amp = 0;
		f = 0;
		Sinusoid fbig = new Sinusoid(amp, f, phase, bias, N);
		fbig.computeWave();
		for (int k = 1; k <= numberOfTerms; k++) {
			f = 2 * k - 1;
			amp = 1 / f;
			Sinusoid term = new Sinusoid(amp, f, phase, bias, N);
			term.computeWave();
			fbig.addData(term);
		}
		// fbig.showImg("f10000");

		// really big g
		phase = 0;
		bias = 0;
		amp = 0;
		f = 0;
		Sinusoid gbig = new Sinusoid(amp, f, phase, bias, N);
		gbig.computeWave();
		for (int k = 1; k <= numberOfTerms; k++) {
			f = 2 * k;
			amp = 1 / f;
			Sinusoid term = new Sinusoid(amp, f, phase, bias, N);
			term.computeWave();
			gbig.addData(term);
		}
		// gbig.showImg("g10000");

		/*
		 * Sinusoid.makeGraph(f3.getDataValues(), f10.getDataValues() ,
		 * f50.getDataValues()); Sinusoid.makeGraph(g3.getDataValues(),
		 * g10.getDataValues(), g50.getDataValues());
		 * Sinusoid.makeGraph(f50.getDataValues(), fbig.getDataValues());
		 * Sinusoid.makeGraph(g50.getDataValues(), gbig.getDataValues());
		 * FourierTransform.graphPSD(f50.FFT(), 20);
		 * FourierTransform.graphPSD(g50.FFT(), 15);
		 * FourierTransform.graphPSD(fbig.FFT(), 20);
		 * FourierTransform.graphPSD(gbig.FFT(), 15);
		 */

		// ************** Problem 2 *************
		Sinusoid v1 = new Sinusoid(1, 13, 0, 0, 512);
		Sinusoid v2 = new Sinusoid(1, 31, 0, 0, 512);
		Sinusoid xt = new Sinusoid(0, 0, 0, 0, 512); // add
		Sinusoid yt = new Sinusoid(0, 0, 0, 0, 512); // multiply

		v1.computeWave();
		v2.computeWave();
		xt.computeWave();
		yt.computeWave();
		xt.addData(v1);
		xt.addData(v2);
		yt.addData(v2);
		yt.multiplyData(v1);
		
		/*
		 * xt.showImg(); yt.showImg(); FourierTransform.graphPSD(xt.FFT(), 40);
		 * FourierTransform.graphPSD(yt.FFT(), 60);
		 */

		// ************** Problem 3 ****************

		Sinusoid whiteLight = new Sinusoid(0, 0, 0, 0, 256);
		whiteLight.getDataValues()[1] = 1;
		// FourierTransform.graphPSD(whiteLight.FFT(), 256);
		Sinusoid whiteLight2 = new Sinusoid(0, 0, 0, 0, 256);
		whiteLight2.getDataValues()[92] = 1;
		// FourierTransform.graphPSD(whiteLight2.FFT(), 256);
		Sinusoid whiteLight3 = new Sinusoid(0, 0, 0, 0, 256);
		whiteLight3.getDataValues()[226] = 1;
		// FourierTransform.graphPSD(whiteLight3.FFT(), 256);

		Sinusoid sample = new Sinusoid(1, 10, 0, 0, 256);
		sample.computeWave();
		// sample.showImg();
		ComplexNumber[] sampleFFT = sample.FFT();
		// ComplexNumber.printArray(sampleFFT);
		// FourierTransform.graphPSD(sampleFFT, 100);

		Sinusoid sample2 = new Sinusoid(1, 10, 0.33, 0, 256);
		sample2.computeWave();
		// sample2.showImg();
		ComplexNumber[] sample2FFT = sample2.FFT();
		// System.out.println("Sample2");
		// ComplexNumber.printArray(sample2FFT);
		// FourierTransform.graphPSD(sample2FFT, 100);

		// ********* Problem 4 ************

		// lowpass

		Sinusoid f7 = new Sinusoid(0, 0, 0, 0, N);
		f7.computeWave();
		amp = 0;
		f = 0;
		phase = 0;
		bias = 0;
		for (int k = 1; k <= 7; k++) {
			f = 2 * k - 1;
			amp = 1 / f;
			Sinusoid term = new Sinusoid(amp, f, phase, bias, N);
			term.computeWave();
			f7.addData(term);
		}
		// f7.showImg();
		// FourierTransform.graphPSD(f7.FFT(), 20);
		double[] filter = Arrays.copyOf(f50.getDataValues(), f50.N);
		// odds 1, 3, 5, 7, 9, 11, 13 = 7 lowest freqs
		double[] lowFiltrate = FourierTransform.filter(filter, FILTER.LOW, 13);
		// Sinusoid.makeGraph(lowFiltrate);

		// high pass

		Sinusoid fHigh43 = new Sinusoid(0, 0, 0, 0, 512);
		fHigh43.computeWave();
		amp = 0;
		f = 0;
		phase = 0;
		bias = 0;
		for (int k = 8; k <= 50; k++) { // freq 15 to 99
			f = 2 * k - 1;
			amp = 1 / f;
			Sinusoid term = new Sinusoid(amp, f, phase, bias, N);
			term.computeWave();
			fHigh43.addData(term);
		}
		// fHigh43.showImg();
		// FourierTransform.graphPSD(fHigh43.FFT(), 102);

		filter = Arrays.copyOf(f50.getDataValues(), f50.N);
		double[] highFiltrate = FourierTransform.filter(filter, FILTER.HIGH, 14);
		// Sinusoid.makeGraph(highFiltrate);
		// ComplexNumber[] highFFT =
		// FourierTransform.FFT(FourierTransform.makeComplex(highFiltrate), 1);
		// FourierTransform.graphPSD(highFFT, 106);

		// band

		Sinusoid f5to8 = new Sinusoid(0, 0, 0, 0, 512);
		f5to8.computeWave();
		amp = 0;
		f = 0;
		phase = 0;
		bias = 0;
		for (int k = 5; k <= 8; k++) { // freq 9 to 15
			f = 2 * k - 1;
			amp = 1 / f;
			Sinusoid term = new Sinusoid(amp, f, phase, bias, N);
			term.computeWave();
			f5to8.addData(term);
		}
		// f5to8.showImg();
		// FourierTransform.graphPSD(f5to8.FFT(), 30);

		filter = Arrays.copyOf(f50.getDataValues(), f50.N);
		double[] bandFiltrate = FourierTransform.filter(filter, FILTER.BAND, 9, 15);
		//Sinusoid.makeGraph(bandFiltrate);

		// notch
		Sinusoid fno5to8 = new Sinusoid(0, 0, 0, 0, 512);
		fno5to8.computeWave();
		amp = 0;
		f = 0;
		phase = 0;
		bias = 0;
		for (int k = 1; k <= 4; k++) { // freq 1 to 7
			f = 2 * k - 1;
			amp = 1 / f;
			Sinusoid term = new Sinusoid(amp, f, phase, bias, N);
			term.computeWave();
			fno5to8.addData(term);
		}

		for (int k = 9; k <= 50; k++) { // freq 17 to 99
			f = 2 * k - 1;
			amp = 1 / f;
			Sinusoid term = new Sinusoid(amp, f, phase, bias, N);
			term.computeWave();
			fno5to8.addData(term);
		}
		// fno5to8.showImg();
		// FourierTransform.graphPSD(fno5to8.FFT(), 25);

		filter = Arrays.copyOf(f50.getDataValues(), f50.N);
		double[] notchFiltrate = FourierTransform.filter(filter, FILTER.NOTCH, 9, 15);
		// Sinusoid.makeGraph(notchFiltrate);

		// ********* Problem 5 ************
		File toneAFile = new File("tonedataA1.txt");
		ComplexNumber[] toneA = FourierTransform.readInData(toneAFile);
		File toneBFile = new File("tonedataB1.txt");
		ComplexNumber[] toneB = FourierTransform.readInData(toneBFile);
		ComplexNumber[] toneAFFT = FourierTransform.FFT(toneA, 1);
		ComplexNumber[] toneBFFT = FourierTransform.FFT(toneB, 1);
		// FourierTransform.graphPSD(toneAFFT, 180);
		// FourierTransform.graphPSD(toneBFFT, 130);

		System.out.println("Tone A k1 = 65, k2 = 152, fs = 44.1kHz");
		double fk1 = ((65 * 44100) / toneA.length);
		double fk2 = ((152 * 44100 / toneA.length));
		System.out.println("f65 = " + fk1 + " f152 = " + fk2);
		System.out.println("Keypad key = A");
		System.out.println();

		System.out.println("Tone B k1 = 72, k2 = 112, fs = 44.1Hz");
		double fbk1 = ((72 * 44100) / toneB.length);
		double fbk2 = ((112 * 44100 / toneB.length));
		System.out.println("f72 = " + fbk1 + " f112 = " + fbk2);
		System.out.println("Keypad key = 4");
		System.out.println();

		// ************ Problem 6 ********************
		File rangeTest = new File("rangeTestDataSpring2018.txt");
		Scanner sc = null;
		try {
			sc = new Scanner(rangeTest);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		sc.nextLine();
		int i = 0;
		int k = 0;
		double[] rawPulse = new double[1024];
		double[] echoSamp = new double[1024];
		while (sc.hasNext()) {
			String line = sc.nextLine();
			if (line.equals("") || line.equals("1024 signal samples"))
				continue;
			if (i < 52)
				rawPulse[i] = Double.parseDouble(line);
			else {
				echoSamp[k] = Double.parseDouble(line);
				k++;
			}
			i++;
		}
		
		// Sinusoid.makeGraph(echoSamp);
		// Sinusoid.makeGraph(rawPulse);
		ComplexNumber[] echoSample = FourierTransform.makeComplex(echoSamp);
		ComplexNumber[] sampPulse = FourierTransform.makeComplex(rawPulse);
		// perform FFT on sample and pulse
		echoSample = FourierTransform.FFT(echoSample, 1);
		sampPulse = FourierTransform.FFT(sampPulse, 1);
		ComplexNumber[] sampPulseConj = Arrays.copyOf(sampPulse, sampPulse.length);
		// get the conjugate of P
		for (int p = 0; p < sampPulse.length; p++) {
			sampPulseConj[p] = sampPulseConj[p].conjugate();
		}

		ComplexNumber[] sampPulseProduct = new ComplexNumber[sampPulse.length];
		for (int x = 0; x < sampPulseProduct.length; x++) {
			sampPulseProduct[x] = echoSample[x].multiply(sampPulseConj[x]);
		}
		sampPulseProduct = FourierTransform.FFT(sampPulseProduct, -1);
		// ComplexNumber.printArray(sampPulseProduct);
		double[] rCorr = ComplexNumber.getRealVals(sampPulseProduct);
		// Sinusoid.makeGraph(rCorr);

		// find peak
		double secPeak = 0;
		int secTime = 0;
		double peak = rCorr[0];
		int timeSteps = 0;
		for (int index = 0; index < rCorr.length; index++) {
			if (peak < rCorr[index]) {
				secPeak = peak;
				secTime = timeSteps;
				peak = rCorr[index];
				timeSteps = index;
			}
		}
		System.out.println("k = " + timeSteps + " peak = " + peak);
		System.out.println("k2 = " + secTime + " peak2 = " + secPeak);
		double timeInterval = 1 / 50000.0;
		System.out.println("fs = 50khz so deltaT = " + timeInterval);
		double timeDelay = timeInterval * timeSteps;
		System.out.println("time delay = time interval * timesteps to peak = " + timeDelay);
		timeDelay = timeDelay + 0.002;
		System.out.println("Dead time = 0.002s so add it -> " + timeDelay);
		double dist = 1500 * timeDelay;
		System.out.println("Distance = r * t = 1500m/s * " + timeDelay + "s = " + dist + "m");
		System.out.println("Want distance to reflector not round trip distances so devide by 2 = " + dist / 2 + "m");

		// convolution
		ComplexNumber[] sixPointP = new ComplexNumber[256];
		for (int p = 0; p < sixPointP.length; p++) {
			if (p < 6)
				sixPointP[p] = new ComplexNumber(1.0 / 6.0, 0);
			else
				sixPointP[p] = new ComplexNumber(0, 0);
		}
		
		sixPointP = FourierTransform.FFT(sixPointP, 1);
		ComplexNumber[] origSamp = Arrays.copyOf(echoSample, 256); // already FFT it
		ComplexNumber[] ux = new ComplexNumber[256];
		for (int ii = 0; ii < origSamp.length; ii++) {
			ux[ii] = origSamp[ii].multiply(sixPointP[ii]);
		}

		ux = FourierTransform.FFT(ux, -1);
		double[] smoothedSamp = ComplexNumber.getRealVals(ux);
		double[] original = Arrays.copyOf(echoSamp, 256);
		
		// Sinusoid.makeGraph(original);
		// Sinusoid.makeGraph(smoothedSamp);
		// Sinusoid.makeGraph(original, smoothedSamp);
	}

}
