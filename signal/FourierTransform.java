package signal;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

import javax.swing.JFrame;
import javax.swing.JPanel;

public class FourierTransform {

	public enum VectorCondition {
		LEAVE_ALONE, PAD_DATA, CUT_OFF_DATA
	}

	/**
	 * The Fast Fourier Transform (FFT), computes the Discrete Fourier Transform
	 * (DFT) much faster. The DFT is on order O(N^2), while FFT is O(N*log(N)). The
	 * DFT can be implemented much faster if N (the number of data samples) is a
	 * power of two. Consider an NxN Matrix F (that is the Vandermonde matrix) which
	 * transforms a vector of our samples to the corresponding Fourier coefficients.
	 * It can be decomposed into smaller DFT problems based on the basic principle
	 * that from the original DFT the even sums and the odd sums can be separated
	 * out saving a factor of 2. This process can be repeated recursively which is
	 * where the logN term comes from in the big O complexity. This particular
	 * method also allows for the inverse Fourier transform to be computed based on
	 * direction. If direction is a positive one it transforms the function from a
	 * time domain to a frequency domain. if the direction is negative then it
	 * transforms the Fourier frequency coefficients to time data.
	 * 
	 * @param z
	 *            the data vector of real numbers.
	 * @param dir
	 *            direction of the transform
	 * @return the complex vector transformed to its transformation pair
	 */
	public static ComplexNumber[] FFT(ComplexNumber[] data, int direction, VectorCondition vCond) {
		data = conditionVector(data, vCond);
		int N = data.length; // number of data samples
		// direction
		int dir;
		if (direction == 1 || direction == -1)
			dir = direction;
		else
			dir = 1;
		double theta = (-2 * Math.PI * dir) / N; // how much around the circle sectioned off into N slices

		int r = N / 2; // reverse

		ComplexNumber fourierCoeff;
		for (int i = 1; i < N;) { // N - 1 instead of N?
			fourierCoeff = new ComplexNumber(Math.cos(i * theta), Math.sin(i * theta)); // initialize coeff
			// for each compute the geometric progression of theta
			for (int k = 0; k < N;) {
				ComplexNumber u = new ComplexNumber(1, 0); // scaler
				for (int m = 0; m < r; m++) { // r - 1 instead of r?
					ComplexNumber t = data[k + m].sub(data[k + m + r]);
					data[k + m] = data[k + m].add(data[k + m + r]);
					data[k + m + r] = t.multiply(u);
					u = fourierCoeff.multiply(u);
				}
				k = k + (2 * r);
			}
			// for FFT want power of 2
			i = 2 * i;
			r = r / 2;
		}
		// re-arange the results
		for (int i = 0; i < N; i++) { //
			r = i;
			int k = 0;
			for (int m = 1; m < N;) { //
				k = (2 * k) + (r % 2);
				r = r / 2;
				m = 2 * m;
			}
			if (k > i) {
				ComplexNumber t = data[i];
				data[i] = data[k];
				data[k] = t;
			}
		}
		// then it is the inverse FFT
		if (dir < 0) {
			ComplexNumber n = new ComplexNumber(N, 0);
			for (int i = 0; i < N; i++)
				data[i] = data[i].divide(n);
		}
		return data;

	}

	public static ComplexNumber[] FFT(ComplexNumber[] data, int direction) {
		return FFT(data, direction, VectorCondition.LEAVE_ALONE);
	}

	public static ComplexNumber[] conditionVector(ComplexNumber[] data, VectorCondition v) {
		switch (v) {
		case LEAVE_ALONE:
			return data;
		case PAD_DATA: {
			int N = data.length;
			while ((N & (N - 1)) != 0) { // bitwise test
				N++;
			}
			ComplexNumber[] newData = new ComplexNumber[N];
			for (int b = 0; b < N; b++) {
				if (b < data.length)
					newData[b] = data[b];
				else
					newData[b] = new ComplexNumber(0, 0);
			}
			return newData;
		}
		case CUT_OFF_DATA: {
			int N = data.length;
			while ((N & (N - 1)) != 0) {
				N--;
			}
			ComplexNumber[] d = Arrays.copyOf(data, N);
			return d;
		}
		default:
			throw new IllegalArgumentException("Can't yet handle " + v);
		}
	}

	public static ComplexNumber[] makeComplex(double[] data) {
		ComplexNumber[] c = new ComplexNumber[data.length];
		for (int i = 0; i < c.length; i++) {
			c[i] = new ComplexNumber(data[i], 0);
		}
		return c;
	}

	public static ComplexNumber[] readInData(File file) {
		Scanner sc = null;
		ArrayList<Double> readin = new ArrayList<>();
		try {
			sc = new Scanner(file);
			while (sc.hasNext()) {
				String line = sc.nextLine();
				String[] form = line.split(" ");
				String l = form[form.length - 1];
				double d = Double.parseDouble(l);
				readin.add(d);
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} finally {
			sc.close();
		}
		double[] data = new double[readin.size()];
		for (int k = 0; k < readin.size(); k++) {
			data[k] = readin.get(k);
		}

		ComplexNumber[] compData = makeComplex(data);
		return compData;
	}

	public static ComplexNumber[] highPassFilter(ComplexNumber[] fData, double cutOffFreq) {
		ComplexNumber[] filtrate = new ComplexNumber[fData.length];
		filtrate[0] = new ComplexNumber(0, 0);
		for (int i = 1; i <= fData.length / 2; i++) {
			if (i < cutOffFreq) {
				filtrate[i] = new ComplexNumber(0, 0);
				filtrate[fData.length - i] = new ComplexNumber(0, 0);
			} else {
				filtrate[i] = fData[i];
				filtrate[fData.length - i] = fData[fData.length - i];
			}
		}
		return filtrate;
	}

	public static ComplexNumber[] lowPassFilter(ComplexNumber[] fData, double cutOffFreq) {
		ComplexNumber[] filtrate = new ComplexNumber[fData.length];
		filtrate[0] = fData[0];
		for (int i = 1; i <= fData.length / 2; i++) {
			if (i > cutOffFreq) {
				filtrate[i] = new ComplexNumber(0, 0);
				filtrate[fData.length - i] = new ComplexNumber(0, 0);
			} else {
				filtrate[i] = fData[i];
				filtrate[fData.length - i] = fData[fData.length - i];
			}
		}
		return filtrate;
	}

	public static ComplexNumber[] bandPassFilter(ComplexNumber[] fData, double f1, double f2) {
		ComplexNumber[] filtrate = new ComplexNumber[fData.length];
		filtrate[0] = new ComplexNumber(0, 0);
		for (int i = 1; i <= fData.length / 2; i++) {
			if (i < f1 || i > f2) {
				filtrate[i] = new ComplexNumber(0, 0);
				filtrate[fData.length - i] = new ComplexNumber(0, 0);
			} else {
				filtrate[i] = fData[i];
				filtrate[fData.length - i] = fData[fData.length - i];
			}
		}
		return filtrate;
	}

	public static ComplexNumber[] notchPassFilter(ComplexNumber[] fData, double f1, double f2) {
		ComplexNumber[] filtrate = new ComplexNumber[fData.length];
		filtrate[0] = fData[0];
		for (int i = 1; i <= fData.length / 2; i++) {
			if (i >= f1 && i <= f2) {
				filtrate[i] = new ComplexNumber(0, 0);
				filtrate[fData.length - i] = new ComplexNumber(0, 0);
			} else {
				filtrate[i] = fData[i];
				filtrate[fData.length - i] = fData[fData.length - i];
			}
		}
		return filtrate;
	}

	public enum FILTER {
		LOW, HIGH, BAND, NOTCH
	};

	public static double[] filter(double[] filter, FILTER filterType, double cutoff) {
		ComplexNumber[] fData = FourierTransform.makeComplex(filter);
		fData = FourierTransform.FFT(fData, 1);
		if (filterType == FILTER.HIGH) {
			fData = FourierTransform.highPassFilter(fData, cutoff);
			fData = FourierTransform.FFT(fData, -1);
			return ComplexNumber.getRealVals(fData);
		} else if (filterType == FILTER.LOW) {
			fData = FourierTransform.lowPassFilter(fData, cutoff);
			fData = FourierTransform.FFT(fData, -1);
			return ComplexNumber.getRealVals(fData);
		} else {
			throw new IllegalArgumentException("Missing second frequency");
		}
	}

	public static double[] filter(double[] filter, FILTER filterType, double f1, double f2) {
		ComplexNumber[] fData = FourierTransform.makeComplex(filter);
		fData = FourierTransform.FFT(fData, 1);
		if (filterType == FILTER.BAND) {
			fData = FourierTransform.bandPassFilter(fData, f1, f2);
			fData = FourierTransform.FFT(fData, -1);
			return ComplexNumber.getRealVals(fData);

		} else if (filterType == FILTER.NOTCH) {
			fData = FourierTransform.notchPassFilter(fData, f1, f2);
			fData = FourierTransform.FFT(fData, -1);
			return ComplexNumber.getRealVals(fData);

		} else {
			throw new IllegalArgumentException("To many parameters given for filter type");
		}
	}

	public static double[] computePSD(ComplexNumber[] data) {
		double[] psd = new double[data.length / 2];
		for (int i = 0; i < (data.length / 2); i++) {
			psd[i] = Math.pow(data[i].mod(), 2);
		}
		return psd;
	}

	public static void graphPSD(double[] psd, int span) {
		JFrame f = new JFrame();
		int s;
		if (span < 0 || span > psd.length) {
			s = psd.length;
		} else
			s = span;
		JPanel jp = new JPanel() {
			@Override
			public void paintComponent(Graphics g) {
				super.paintComponent(g);
				drawPSD(psd, g, s);
			}
		};
		f.setSize(new Dimension(1500, 1500));
		jp.setSize(new Dimension(1500, 1500));
		f.add(jp);
		f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		f.setVisible(true);
	}

	private static void drawPSD(double[] psd, Graphics g, int span) {
		int HEIGHT, WIDTH;
		HEIGHT = WIDTH = 1500;
		BufferedImage waveImg = new BufferedImage(HEIGHT, WIDTH, BufferedImage.TYPE_INT_RGB);
		int xAxis, yAxis;
		yAxis = (int) (0.9 * HEIGHT);
		xAxis = (int) (0.9 * WIDTH);
		Graphics2D g2 = (Graphics2D) waveImg.getGraphics();
		g2.setColor(Color.white);
		g2.fillRect(0, 0, WIDTH, HEIGHT);

		// normalize
		double maxHeight = psd[0];
		for (int i = 0; i < psd.length; i++) {
			if (psd[i] > maxHeight)
				maxHeight = psd[i];
		}

		double[] nrmlized = new double[psd.length];
		for (int j = 0; j < psd.length; j++) {
			nrmlized[j] = psd[j] / maxHeight;
		}

		// draw the axis
		int xCord = WIDTH - xAxis;
		g2.setColor(Color.black);
		g2.setStroke(new BasicStroke(5)); // set line thickness
		g2.drawLine(xCord, yAxis, WIDTH, yAxis); // x-axis
		g2.drawLine(xCord, yAxis, xCord, 0); // y-axis

		// draw markers
		g2.setColor(Color.red);
		g2.setFont(new Font("TimesRoman", Font.PLAIN, 48));
		g2.drawString("Y = " + maxHeight, 20, 40);

		// draw the wave
		g2.setColor(Color.blue);
		double t1, t2;
		int y1, y2;
		double space;
		space = (1 / (double) span) * (double) xAxis;
		// System.out.println("PSD");
		for (int k = 1; k < span; k++) {
			t1 = nrmlized[k - 1];
			t2 = nrmlized[k];
			y2 = (yAxis - (int) (t2 * yAxis));
			g2.drawLine(xCord, yAxis, (int) (xCord + space), yAxis);
			g2.drawLine((int) (xCord + space), yAxis, (int) (xCord + space), y2);
			if (y2 < yAxis) {
				g2.drawString("" + k, (int) (xCord + space) - 15, yAxis + 50);
				// System.out.println("k = " + k + " y2 = " + y2 + " t2 = " + t2);
			}
			xCord = (int) (xCord + space);
		}
		g.drawImage(waveImg, 0, 0, null);
	}

	public static void graphPSD(ComplexNumber[] data, int span) {
		double[] psd = computePSD(data);
		/*
		 * for (int i = 0; i < psd.length; i++) { System.out.println(i + " " + psd[i]);
		 * }
		 */
		graphPSD(psd, span);
	}

	public static ComplexNumber[][] twoDFFT(ComplexNumber[][] data, int direction) {
		int n = data.length; // rows
		int m = data[0].length; // columns

		// perform column 1dFFT
		for (int k = 0; k < m; k++) {
			// get the column
			ComplexNumber[] col = new ComplexNumber[n];
			for (int r = 0; r < n; r++) {
				col[r] = data[r][k];
			}
			col = FourierTransform.FFT(col, direction);

			// put the info back into data
			for (int j = 0; j < n; j++) {
				data[j][k] = col[j];
			}
		}

		// perform row 1dFFT on the intermediate
		for (int i = 0; i < m; i++) {
			data[i] = FourierTransform.FFT(data[i], direction);
		}
		return data;
	}

}
