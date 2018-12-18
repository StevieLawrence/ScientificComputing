package TSP;

import java.awt.BasicStroke;
import java.awt.Canvas;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.FileWriter;
import java.io.IOException;

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JPanel;


public class Histogram {

	double[] histogram, nrmHist;
	double peak, maxBound, minBound, binWidth;
	BufferedImage histImage;
	int numOfElements;

	/**
	 * A histogram is a sequence of bins and frequencies for values that fall within those bins 
	 * @param binAmount
	 * @param minBound
	 * @param maxBound
	 */
	public Histogram(int binAmount, double minBound, double maxBound) {
		histogram = new double[binAmount];
		if (maxBound == 0)
			this.maxBound = 1;
		else
			this.maxBound = maxBound;
		this.minBound = minBound;
		numOfElements = 0;
		binWidth = (maxBound - minBound) / binAmount;
		nrmHist = new double[histogram.length];
	}
	
	/**
	 * increment the count of the bin that the value falls in
	 * @param value
	 */
	public void add(double value) {
		double dx = value - minBound;
		histogram[(int) (dx / binWidth)]++;
		numOfElements++;
	}
	
	/**
	 * bin with the largest count
	 * @return
	 */
	public double getPeak() {
		peak = 1; // for normalization method to avoid division by zero
		for (double r : histogram) {
			if (r > peak) {
				peak = r;
			}
		}
		return peak;
	}
	
	/**
	 * divide all the bins by the peak so that the largest bin coresponds to 1
	 */
	public void normalize() {
		int index = 0;
		peak = getPeak();
		for (double freq : histogram) {
			/*
			 * normalize the image basedon the peak and then multiply it by the y-axis range
			 * to re-scale the histogram to fit on the screen
			 */
			double newFreq = (double) freq / peak;
			nrmHist[index] = newFreq;
			index++;
		}
	}
	
	/**
	 * Draw a picture of the histogram
	 * @param g
	 */
	private void draw(Graphics g) {
		int HEIGHT, WIDTH;
		HEIGHT = WIDTH = 600;
		int yRange, xRange;
		histImage = new BufferedImage(HEIGHT, WIDTH, BufferedImage.TYPE_INT_RGB);
		yRange = (int) (0.9 * HEIGHT);
		xRange = (int) (0.9 * WIDTH);
		Graphics2D g2 = (Graphics2D) histImage.getGraphics();
		
		g2.setColor(Color.white);
		g2.fillRect(0, 0, WIDTH, HEIGHT);

		// draw the axis
		int xCord = WIDTH - xRange;
		g2.setColor(Color.black);
		g2.setStroke(new BasicStroke(5)); // set line thickness
		g2.drawLine(xCord, yRange, WIDTH, yRange); // x-axis
		g2.drawLine(xCord, yRange, xCord, 0); // y-axis
		
		// normalize
		normalize();
		
		/*
		 * draw the scaled histogram. Re-scale the histogram in relation to the x-axis.
		 * Evenly spaced objects of a known total width can go into a larger known width
		 * by using the formula (outerWidth - widthOfObjects) / (widthOfObjects + 1).
		 */
		// get the space length
		double spaceLength = (double) (xRange - nrmHist.length) / (nrmHist.length + 1);

		int k = 0, y1Cord, y2Cord, x2Cord;
		double x = xCord; // since spacing is a double, x needs to be a double
		for (int i = 0; i < nrmHist.length - 1; i++) {
			k = i + 1;
			// top left is (0, 0) so need to flip the y values
			y1Cord = yRange - (int) (nrmHist[i] * yRange);
			y2Cord = yRange - (int) (nrmHist[k] * yRange);
			xCord = (int) x;
			x2Cord = (int) (x + spaceLength + 1);
			g2.drawLine(xCord, y1Cord, x2Cord, y2Cord);
			x += spaceLength + 1;
		}
		
		g.drawImage(histImage, 0, 0, null);
	}

	public BufferedImage getHistImage() {
		return histImage;
	}

	public int numberOfElements() {
		return numOfElements;
	}
	
	/**
	 * Write the counts to a file
	 * @param filename
	 */
	public void logInfo(String filename) {
		normalize();
		try {
			FileWriter writer = new FileWriter(filename);
			for (int index = 0; index < nrmHist.length; index++) {
				String data = ""+ index + "," + nrmHist[index] + "\n";
				writer.write(data);
			}
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/**
	 * Display the histogram image
	 */
	public void showHistImage() {
		JFrame f = new JFrame();
		JPanel jp = new JPanel(){
            @Override
            public void paintComponent(Graphics g) {
                super.paintComponent(g);
                draw(g);
            }
        };
        f.setSize(new Dimension(700,700));
        jp.setSize(new Dimension(605, 605));
		f.add(jp);
		f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		f.setVisible(true);
	}
	
}
