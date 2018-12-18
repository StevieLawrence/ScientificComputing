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
import java.util.Scanner;

import javax.swing.JFrame;
import javax.swing.JPanel;

public class Sinusoid {

	double amplitude, bias, frequency, phase, startTime, endTime;
	private double dt;
	int N;
	private int dynamicVal;
	private double[] dataValues, contingentVals;
	boolean isDynamic;

	public Sinusoid(double amp, double freq, double phase, double bias, int N) {
		amplitude = amp;
		frequency = freq;
		this.phase = phase;
		this.bias = bias;
		this.N = N;
		startTime = 0;
		endTime = 1;
		dt = (endTime - startTime) / N;
		dataValues = new double[N];
		dynamicVal = -1;
		isDynamic = false;

	}

	public void computeWave() {
		int sample = 0;
		for (double i = startTime; i < endTime; i += dt) {
			double x = amplitude * Math.sin(2 * Math.PI * frequency * (i - phase)) + bias;
			dataValues[sample] = x;
			sample++;
		}
		isDynamic = false;
		contingentVals = null;
	}

	public ComplexNumber[] FFT() {
		ComplexNumber[] d = FourierTransform.makeComplex(dataValues);
		return FourierTransform.FFT(d, 1);
	}

	public ComplexNumber[] invFFT() {
		ComplexNumber[] d = FourierTransform.makeComplex(dataValues);
		return FourierTransform.FFT(d, -1);
	}

	public void addData(Sinusoid s) {
		if (N != s.getN()) {
			System.out.println("Sinusoids have different sample amounts");
		} else {
			double[] sd = s.getDataValues();
			for (int i = 0; i < N; i++) {
				dataValues[i] = dataValues[i] + sd[i];
			}

		}
	}

	public void multiplyData(Sinusoid s) {
		if (N != s.getN()) {
			System.out.println("Sinusoids have different sample amounts");
		} else {
			double[] sd = s.getDataValues();
			for (int i = 0; i < N; i++) {
				dataValues[i] = dataValues[i] * sd[i];
			}
		}
	}

	public double[] getDataValues() {
		return dataValues;
	}

	private void setDataValues(double[] d) {
		dataValues = d;
	}

	public void beginDynamicWave() {
		contingentVals = new double[N];
		dynamicVal = 0;
	}

	public boolean isDynamComputing() {
		if (dynamicVal < 0)
			return true;
		else
			return false;
	}

	public boolean isDyanmic() {
		return isDynamic;
	}

	public void calcDynamicVal(double amp, double freq, double phase, double bias) {
		if (dynamicVal < 0) {
			throw new IllegalArgumentException("Error: Need to begin new dynamic wave computation");
		} else {
			double t = dynamicVal * dt;
			double x = amp * Math.sin(2 * Math.PI * (t - phase)) - bias;
			contingentVals[dynamicVal] = x;
			dynamicVal++;
			if (dynamicVal == N) {
				dynamicVal = -1;
			}
		}
	}

	public void setDynamic() {
		if (contingentVals != null) {
			dataValues = contingentVals;
			isDynamic = true;
		} else
			System.out.println("Need to make dynamic session");
	}

	private void draw(Graphics g) {
		int range = 5;
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
		double maxHeight = dataValues[0], smallestHeight = dataValues[0];
		for (int i = 0; i < N; i++) {
			if (dataValues[i] > maxHeight)
				maxHeight = dataValues[i];
			if (dataValues[i] < smallestHeight)
				smallestHeight = dataValues[i];
		}
		while (Math.abs(maxHeight) > range || Math.abs(smallestHeight) > range)
			range += 5;
		double[] nrmlized = new double[N];
		for (int j = 0; j < N; j++) {
			nrmlized[j] = dataValues[j] / range;
		}

		// draw the axis
		int xCord = WIDTH - xAxis;
		g2.setColor(Color.black);
		g2.setStroke(new BasicStroke(5)); // set line thickness
		g2.drawLine(xCord, yAxis / 2, WIDTH, yAxis / 2); // x-axis
		g2.drawLine(xCord, yAxis, xCord, 0); // y-axis

		// draw markers
		g2.setColor(Color.red);
		g2.setFont(new Font("TimesRoman", Font.PLAIN, 48));
		g2.drawString("Y = " + range, 20, 40);

		// draw the wave
		g2.setColor(Color.blue);
		double t1, t2;
		int y1, y2;
		double space = dt * xAxis;
		for (int k = 1; k < N; k++) {
			t1 = nrmlized[k - 1];
			t2 = nrmlized[k];
			y1 = (yAxis / 2) - (int) (t1 * (yAxis / 2));
			y2 = (yAxis / 2) - (int) (t2 * (yAxis / 2));
			g2.drawLine(xCord, y1, (int) (xCord + space), y2);
			// System.out.println(xCord + " " + y1 + " " + (int) (xCord + space) + " " +
			// y2);
			xCord = (int) (xCord + space);
		}
		g2.setColor(Color.red);
		g2.drawString("X = 1", xCord, yAxis / 2);

		g.drawImage(waveImg, 0, 0, null);
	}

	public void showImg(String title) {
		JFrame f = new JFrame();
		f.setTitle(title);
		JPanel jp = new JPanel() {
			@Override
			public void paintComponent(Graphics g) {
				super.paintComponent(g);
				draw(g);
			}
		};
		f.setSize(new Dimension(1500, 1500));
		jp.setSize(new Dimension(1500, 1500));
		f.add(jp);
		f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		f.setVisible(true);
	}

	public void showImg() {
		showImg("");
	}

	public void setStartTime(double t) {
		startTime = t;
		dt = (endTime - startTime) / N;
	}

	public void setEndTime(double t) {
		endTime = t;
		dt = (endTime - startTime) / N;
	}

	public double getAmplitude() {
		return amplitude;
	}

	public void setAmplitude(double amplitude) {
		this.amplitude = amplitude;
	}

	public double getBias() {
		return bias;
	}

	public void setBias(double bias) {
		this.bias = bias;
	}

	public double getFrequency() {
		return frequency;
	}

	public void setFrequency(double frequency) {
		this.frequency = frequency;
	}

	public double getPhase() {
		return phase;
	}

	public void setPhase(double phase) {
		this.phase = phase;
	}

	public int getN() {
		return N;
	}

	public void setN(int n) {
		N = n;
		dt = (endTime - startTime) / N;
		dataValues = new double[N];
	}

	public double getStartTime() {
		return startTime;
	}

	public double getEndTime() {
		return endTime;
	}

	public static void makeGraph(double[] dataVals) {
		Sinusoid s = new Sinusoid(0, 0, 0, 0, dataVals.length);
		s.dataValues = dataVals;
		s.showImg();

	}

	public static void makeGraph(double[] dataVals1, double[] dataVals2) {
		JFrame f = new JFrame();
		JPanel jp = new JPanel() {
			@Override
			public void paintComponent(Graphics g) {
				super.paintComponent(g);
				int range = 5;
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
				double maxHeight = dataVals1[0], smallestHeight = dataVals1[0];

				for (int i = 0; i < dataVals1.length; i++) {
					if (dataVals1[i] > maxHeight)
						maxHeight = dataVals1[i];
					if (dataVals2[i] > maxHeight)
						maxHeight = dataVals2[i];
					if (dataVals1[i] < smallestHeight)
						smallestHeight = dataVals1[i];
					if (dataVals2[i] < smallestHeight)
						smallestHeight = dataVals2[i];
				}
				while (Math.abs(maxHeight) > range || Math.abs(smallestHeight) > range)
					range += 5;
				double[] nrmlized1 = new double[dataVals1.length];
				double[] nrmlized2 = new double[dataVals2.length];
				for (int j = 0; j < dataVals1.length; j++) {
					nrmlized1[j] = dataVals1[j] / range;
					nrmlized2[j] = dataVals2[j] / range;
				}

				// draw the axis
				int xCord = WIDTH - xAxis;
				g2.setColor(Color.black);
				g2.setStroke(new BasicStroke(5)); // set line thickness
				g2.drawLine(xCord, yAxis / 2, WIDTH, yAxis / 2); // x-axis
				g2.drawLine(xCord, yAxis, xCord, 0); // y-axis

				// draw markers
				g2.setColor(Color.red);
				g2.setFont(new Font("TimesRoman", Font.PLAIN, 48));
				g2.drawString("Y = " + range, 20, 40);

				// draw the waves
				double t1, t2;
				int y1, y2;
				double space = (1 / (double) dataVals1.length) * (double) xAxis;
				for (int k = 1; k < dataVals1.length; k++) {
					g2.setColor(Color.blue);
					t1 = nrmlized1[k - 1];
					t2 = nrmlized1[k];
					y1 = (yAxis / 2) - (int) (t1 * (yAxis / 2));
					y2 = (yAxis / 2) - (int) (t2 * (yAxis / 2));
					g2.drawLine(xCord, y1, (int) (xCord + space), y2);

					g2.setColor(Color.green);
					t1 = nrmlized2[k - 1];
					t2 = nrmlized2[k];
					y1 = (yAxis / 2) - (int) (t1 * (yAxis / 2));
					y2 = (yAxis / 2) - (int) (t2 * (yAxis / 2));
					g2.drawLine(xCord, y1, (int) (xCord + space), y2);

					xCord = (int) (xCord + space);
				}
				g2.setColor(Color.red);
				g2.drawString("X = 1", xCord, yAxis / 2);

				g.drawImage(waveImg, 0, 0, null);
			}
		};
		f.setSize(new Dimension(1500, 1500));
		jp.setSize(new Dimension(1500, 1500));
		f.add(jp);
		f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		f.setVisible(true);
	}

	public static void makeGraph(double[] dataVals1, double[] dataVals2, double[] dataVals3) {
		JFrame f = new JFrame();
		JPanel jp = new JPanel() {
			@Override
			public void paintComponent(Graphics g) {
				super.paintComponent(g);
				int range = 5;
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
				double maxHeight = dataVals1[0], smallestHeight = dataVals1[0];

				for (int i = 0; i < dataVals1.length; i++) {
					if (dataVals1[i] > maxHeight)
						maxHeight = dataVals1[i];
					if (dataVals2[i] > maxHeight)
						maxHeight = dataVals2[i];
					if (dataVals3[i] > maxHeight)
						maxHeight = dataVals3[i];
					if (dataVals1[i] < smallestHeight)
						smallestHeight = dataVals1[i];
					if (dataVals2[i] < smallestHeight)
						smallestHeight = dataVals2[i];
					if (dataVals3[i] < smallestHeight)
						smallestHeight = dataVals3[i];
				}
				while (Math.abs(maxHeight) > range || Math.abs(smallestHeight) > range)
					range += 5;
				double[] nrmlized1 = new double[dataVals1.length];
				double[] nrmlized2 = new double[dataVals2.length];
				double[] nrmlized3 = new double[dataVals3.length];
				for (int j = 0; j < dataVals1.length; j++) {
					nrmlized1[j] = dataVals1[j] / range;
					nrmlized2[j] = dataVals2[j] / range;
					nrmlized3[j] = dataVals3[j] / range;
				}

				// draw the axis
				int xCord = WIDTH - xAxis;
				g2.setColor(Color.black);
				g2.setStroke(new BasicStroke(5)); // set line thickness
				g2.drawLine(xCord, yAxis / 2, WIDTH, yAxis / 2); // x-axis
				g2.drawLine(xCord, yAxis, xCord, 0); // y-axis

				// draw markers
				g2.setColor(Color.red);
				g2.setFont(new Font("TimesRoman", Font.PLAIN, 48));
				g2.drawString("Y = " + range, 20, 40);

				// draw the waves
				double t1, t2;
				int y1, y2;

				double space = (1 / (double) dataVals1.length) * (double) xAxis;
				for (int k = 1; k < dataVals1.length; k++) {
					g2.setColor(Color.blue);
					t1 = nrmlized1[k - 1];
					t2 = nrmlized1[k];
					y1 = (yAxis / 2) - (int) (t1 * (yAxis / 2));
					y2 = (yAxis / 2) - (int) (t2 * (yAxis / 2));
					g2.drawLine(xCord, y1, (int) (xCord + space), y2);

					g2.setColor(Color.green);
					t1 = nrmlized2[k - 1];
					t2 = nrmlized2[k];
					y1 = (yAxis / 2) - (int) (t1 * (yAxis / 2));
					y2 = (yAxis / 2) - (int) (t2 * (yAxis / 2));
					g2.drawLine(xCord, y1, (int) (xCord + space), y2);

					g2.setColor(Color.MAGENTA);
					t1 = nrmlized3[k - 1];
					t2 = nrmlized3[k];
					y1 = (yAxis / 2) - (int) (t1 * (yAxis / 2));
					y2 = (yAxis / 2) - (int) (t2 * (yAxis / 2));
					g2.drawLine(xCord, y1, (int) (xCord + space), y2);

					xCord = (int) (xCord + space);
				}
				g2.setColor(Color.red);
				g2.drawString("X = 1", xCord, yAxis / 2);

				g.drawImage(waveImg, 0, 0, null);
			}
		};
		f.setSize(new Dimension(1500, 1500));
		jp.setSize(new Dimension(1500, 1500));
		f.add(jp);
		f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		f.setVisible(true);
	}

	public String toString() {
		String s = "Time | x\n";
		for (int i = 0; i < N; i++) {
			s += (startTime + (dt * i)) + " | " + dataValues[i] + "\n";
		}

		return s;
	}

}
