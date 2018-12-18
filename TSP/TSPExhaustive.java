package TSP;

import java.util.Arrays;

public class TSPExhaustive {

	public static void exhaustive(DistanceMatrix dm) {
		// initialize
		int BINS = 100;
		double minRange = dm.getLowerBound();
		double maxRange = dm.getUpperBound();
		double minTripLength = dm.getUpperBound();
		double maxTripLength = dm.getLowerBound();
		int[] minTrip, maxTrip, currentTrip;
		Histogram hist = new Histogram(BINS, minRange, maxRange);
		double summation = 0;
		double sumSqrd = 0;
		int size = dm.getColSize();
		minTrip = new int[1];
		maxTrip = new int[1];
		
		PermutationGenerator pg = new PermutationGenerator(size, 1); // nail down one, decrease by factor size

		while (pg.hasNext()) {
			// get the next permutation
			currentTrip = pg.next();
			double tLength = dm.calcTripLength(currentTrip);

			if (tLength < minTripLength) {
				minTripLength = tLength;
				minTrip = Arrays.copyOf(currentTrip, currentTrip.length);
			} else if (tLength > maxTripLength) {
				maxTripLength = tLength;
				maxTrip = Arrays.copyOf(currentTrip, currentTrip.length);
			}
			summation += tLength;
			sumSqrd += tLength * tLength;
			hist.add(tLength);
		}
		
		// since decrease by a factor due to nail down the population size is (size - 1)!
		double N = 1;
		for (double i = 1; i <= size - 1; i++) {
			N *= i;
		}
		
		System.out.println("N is " + N);
		
		// calc the mean
		double mean = summation / N;
		
		// calc standard deviation sqrt((sumOfSquares - (sumsquared / size)) - size)
		double stdDev = Math.sqrt((sumSqrd - ((summation * summation) / N)) / N);

		// print the stats
		System.out.println("\nThe max trip and trip length are:");
		System.out.print("Max Trip: ");
		dm.encodeLabels(maxTrip);
		System.out.println("Trip length: " + maxTripLength);

		System.out.println("\nThe min trip and trip length are:");
		System.out.print("Min Trip: ");
		dm.encodeLabels(minTrip);
		System.out.println("Trip length: " + minTripLength);

		System.out.println("\nThe mean is: " + mean);

		System.out.println("\nThe standard deviation is: " + stdDev);

		hist.logInfo("bruteForce.csv");
		hist.showHistImage();

	}

}
