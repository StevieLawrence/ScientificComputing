package signal;

import java.awt.Color;

public class SignalImageDemo {

	public static void main(String[] args) {

		// make the image black
		Picture testSamp = new Picture(512, 512);
		for (int i = 0; i < testSamp.height(); i++) {
			for (int j = 0; j < testSamp.width(); j++) {
				testSamp.set(i, j, new Color(0, 0, 0));
			}
		}

		// create the white rectangular region R
		for (int j = 220; j < 220 + 110; j++) {
			for (int i = 180; i < 180 + 140; i++) {
				testSamp.set(j, i, new Color(255, 255, 255));
			}
		}

		// erase part of the side to make the c-shape
		// 180 + 25 = 205 205 + 90 = 295
		for (int row = 205; row < 295; row++) {
			// 330 - 30 = 300
			for (int col = 300; col < 330; col++) {
				testSamp.set(col, row, new Color(0, 0, 0));
			}
		}
		testSamp.show();

		// make the test pulse image black
		Picture testPulse = new Picture(512, 512);
		for (int i = 0; i < testPulse.height(); i++) {
			for (int j = 0; j < testPulse.width(); j++) {
				testPulse.set(i, j, new Color(0, 0, 0));
			}
		}

		// create top horizontal rectangle
		for (int row = 0; row < 15; row++) {
			for (int col = 0; col < 30; col++) {
				testPulse.set(col, row, new Color(255, 255, 255));
			}
		}

		// create tall mid-section rectangle
		for (int row = 15; row < 105; row++) {
			for (int col = 0; col < 15; col++) {
				testPulse.set(col, row, new Color(255, 255, 255));
			}
		}

		// create bottom rectangle
		for (int row = 105; row < 120; row++) {
			for (int col = 0; col < 30; col++) {
				testPulse.set(col, row, new Color(255, 255, 255));
				;
			}
		}

		testPulse.show();

		Color[][] sample = testSamp.getColorArray();
		Color[][] pulse = testPulse.getColorArray();

		// convert to complex numbers
		ComplexNumber[][] sampleData = new ComplexNumber[sample.length][sample[0].length];
		ComplexNumber[][] pulseData = new ComplexNumber[pulse.length][pulse[0].length];
		Color white = new Color(255, 255, 255);
		for (int row = 0; row < sampleData.length; row++) {
			for (int col = 0; col < sampleData[0].length; col++) {
				if (sample[row][col].equals(white))
					sampleData[row][col] = new ComplexNumber(255, 0);
				else
					sampleData[row][col] = new ComplexNumber(0, 0);
			}
		}

		for (int row = 0; row < sampleData.length; row++) {
			for (int col = 0; col < sampleData[0].length; col++) {
				if (pulse[row][col].equals(white))
					pulseData[row][col] = new ComplexNumber(255, 0);
				else
					pulseData[row][col] = new ComplexNumber(0, 0);
			}
		}

		// perfrom 2dFFT
		ComplexNumber[][] sampleDataF = FourierTransform.twoDFFT(sampleData, 1);
		ComplexNumber[][] pulseDataF = FourierTransform.twoDFFT(pulseData, 1);

		// get the conjugate of the pulse
		for (int row = 0; row < pulseDataF.length; row++) {
			for (int col = 0; col < pulseDataF[0].length; col++) {
				pulseDataF[row][col] = pulseDataF[row][col].conjugate();
			}
		}

		// multiply them together
		ComplexNumber[][] sampTimesPConj = new ComplexNumber[sampleDataF.length][sampleDataF[0].length];
		for (int row = 0; row < sampTimesPConj.length; row++) {
			for (int col = 0; col < sampTimesPConj[0].length; col++) {
				sampTimesPConj[row][col] = sampleDataF[row][col].multiply(pulseDataF[row][col]);
			}
		}

		// perform inverse FFT to get correlation
		sampTimesPConj = FourierTransform.twoDFFT(sampTimesPConj, -1);
		double[][] correlationRaw = new double[sampTimesPConj.length][sampTimesPConj[0].length];
		double minimum = 0;
		double value;
		for (int row = 0; row < correlationRaw.length; row++) {
			for (int col = 0; col < correlationRaw[0].length; col++) {
				value = sampTimesPConj[row][col].getReal();
				correlationRaw[row][col] = value;
				if (value < minimum)
					minimum = value;
			}
		}

		// scale by logarithmic
		double[][] correlationLog = new double[correlationRaw.length][correlationRaw[0].length];
		double val;
		double logMax = 0;
		double logMinNotZero = 200;
		for (int row = 0; row < correlationLog.length; row++) {
			for (int col = 0; col < correlationLog[0].length; col++) {
				// minimum translates data to appropriate positive domain for log
				val = Math.log(correlationRaw[row][col] + minimum);
				if (val > logMax) {
					logMax = val;
				}
				if (val < 0)
					val = 0;
				if (val == Double.NaN)
					val = 0;
				if (val > 0 && logMinNotZero > val )
					logMinNotZero = val;
				correlationLog[row][col] = val;
			}
		}

		// scale to linear
		double[][] corrLinear = new double[correlationLog.length][correlationLog[0].length];
		
		double valueRange = logMax - logMinNotZero;
		double slope = 255.0 / valueRange;
		for (int row = 0; row < corrLinear.length; row++) {
			for (int col = 0; col < corrLinear[0].length; col++) {
				value = correlationLog[row][col] - logMinNotZero;
				if (value < 0 )
					corrLinear[row][col] = 0;
				else 
					corrLinear[row][col] = (slope * value);
			}
		}

		Picture twoDCorrelation = new Picture(512, 512);
		int corr;
		
		for (int row = 0; row < twoDCorrelation.height(); row++) {
			for (int col = 0; col < twoDCorrelation.width(); col++) {
				corr = (int) corrLinear[row][col];
				twoDCorrelation.set(col, row, new Color(corr, corr, corr));
			}

		}

		twoDCorrelation.show();

		Picture heatMap = new Picture(512, 512);
		int correlation;
		for (int row = 0; row < twoDCorrelation.height(); row++) {
			for (int col = 0; col < twoDCorrelation.width(); col++) {
				correlation = (int) corrLinear[row][col];
				if (correlation >= 0.95 * 255)
					heatMap.set(col, row, new Color(204, 0, 0));
				
				else if (correlation >= (0.9 * 255))
					heatMap.set(col, row, Color.RED);
				
				else if (correlation > 0.8 * 255)
					heatMap.set(col, row, new Color(255, 102, 102));

				else if (correlation >= (0.6 * 255))
					heatMap.set(col, row, Color.PINK);

				else if (correlation >= (0.4 * 255))
					heatMap.set(col, row, Color.orange);
				
				else if (correlation > 0)
					heatMap.set(col, row, Color.YELLOW);

				else
					heatMap.set(col, row, Color.BLUE);
			}
		}

		heatMap.show();
	}

}
