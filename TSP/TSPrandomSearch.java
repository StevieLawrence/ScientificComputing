package TSP;

import java.awt.Graphics;
import java.util.Arrays;
import java.util.Random;

import javax.swing.JFrame;

public class TSPrandomSearch {
	
	public static void randomSearch(DistanceMatrix dm, int sampleSize) {
		int BINS = 100;
		DistanceMatrix dMatrix = dm;
		double summation, sumSqrd, minTripLength, maxTripLength;
		int[] minTrip, maxTrip;
		maxTrip = new int[1];
		minTrip = new int[1];
		minTripLength = dMatrix.getUpperBound();
		maxTripLength = dMatrix.getLowerBound();
		summation = 0;
		sumSqrd = 0;
		double minRange = dm.getLowerBound();
		double maxRange = dm.getUpperBound();
				
		// initialize the city positions	
		int[] cities = new int[dMatrix.getRowSize()];
		for (int k = 0; k < cities.length; k++) {
			cities[k] = k;
		}
		
		// initialize the histogram
		Histogram hist = new Histogram(BINS, minRange, maxRange);
		
		for (int i = 0; i < sampleSize; i++) {
			double tLength = dMatrix.calcTripLength(cities);
			
			// greater than current max?
			if (tLength > maxTripLength) {
				maxTripLength = tLength;
				maxTrip = Arrays.copyOf(cities, cities.length);
			}
			
			// smaller than current min?
			if (tLength < minTripLength) {
				minTripLength = tLength;
				minTrip = Arrays.copyOf(cities, cities.length);
			}
			
			// accumulate
			summation += tLength;
			sumSqrd += tLength * tLength;
			
			// increment Histogram
			hist.add(tLength);
			
			// get new trip
			cities = randomShuffle(cities);
		}
		
		// calc the mean
		double mean = summation / sampleSize;
		
		// calc standard deviation    sqrt((sumOfSquares - (sumsquared / size)) - size)
		double stdDev = Math.sqrt((sumSqrd - ((summation * summation) / sampleSize)) / sampleSize);
		
		// print the stats
		System.out.println("The max trip and trip length are:");
		System.out.print("Max Trip: ");
		dMatrix.encodeLabels(maxTrip);
		System.out.println("Trip length: " + maxTripLength);
		
		System.out.println("\nThe min trip and trip length are:");
		System.out.print("Min Trip: ");
		dMatrix.encodeLabels(minTrip);
		System.out.println("Trip length: " + minTripLength);
		
		System.out.println("\nThe mean is: " + mean);
		
		System.out.println("\nThe standard deviation is: " + stdDev);
		
		hist.logInfo("randomSearchData.csv");
		hist.showHistImage();
		
	}
	
	private static int[] randomShuffle(int[] positions) {
		Random rand = new Random();
		int temp;
		
		// random shuffle
		for (int i = 0; i < positions.length; i++) {
			int r = rand.nextInt(positions.length);
			temp = positions[i];
			positions[i] = positions[r];
			positions[r] = temp;
			
		}
		return positions;
	}
	
}
