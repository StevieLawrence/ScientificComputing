package TSP;

import java.util.Arrays;
import java.util.Random;

public class TSPsimulatedAnnealing {
	static double temperature, initTemp;
	static double time = 1;

	public static void simmulatedAnnealing(DistanceMatrix dm, double temp, int numOfTries) {
		int BINS = 100;
		DistanceMatrix dMatrix = dm;
		double summation, sumSqrd, minTripLength, maxTripLength;
		int[] minTrip, maxTrip;
		double minRange = dm.getLowerBound();
		double maxRange = dm.getUpperBound();
		Random rndm = new Random();
		int[] nextCandidate;

		// initialize the city positions
		minTrip = new int[dMatrix.getRowSize()];
		for (int k = 0; k < minTrip.length; k++) {
			minTrip[k] = k;
		}
		minTrip = randomShuffle(minTrip);

		// initialize the histogram
		Histogram hist = new Histogram(BINS, minRange, maxRange);

		// handle x^0
		double tLength = dMatrix.calcTripLength(minTrip);
		maxTrip = Arrays.copyOf(minTrip, minTrip.length);
		maxTripLength = tLength;
		minTripLength = tLength;
		summation = tLength;
		sumSqrd = tLength * tLength;

		// initialize the temperature and initial temperature
		if (temp < 1 && temp > 0) {
			temperature = temp;
			initTemp = temp;
		} else {
			temperature = 0.5;
			initTemp = 0.5;
		}

		for (int i = 0; i < numOfTries - 1;) {

			// generate next candidate by changing only a little
			// initialze next candidate solution
			nextCandidate = Arrays.copyOf(minTrip, minTrip.length);
			nextCandidate = randomSwap(nextCandidate);
			double newTLength = dMatrix.calcTripLength(nextCandidate);

			// If the candidate tour is better than the existing tour, accept it as the new
			// tour.
			double dx = newTLength - minTripLength;
			double deltaF = (minTripLength + dx) - minTripLength;

			if (deltaF < 0) {
				if (newTLength < minTripLength) {
					minTrip = Arrays.copyOf(nextCandidate, nextCandidate.length);
					minTripLength = newTLength;
				}
				summation += newTLength;
				sumSqrd += newTLength * newTLength;
				hist.add(newTLength);
				i++;
			} else {
				// if it is worse still maybe accept it
				double pDeltaF = prob(deltaF);

				// select a value z from uniform random distribution from 0 to 1/kt (higher the
				// temperature the smaller the interval)
				double z = rndm.nextDouble() * (1 / (Math.sqrt(2 * Math.PI) * temperature));

				if (pDeltaF > z) {
					if (newTLength < minTripLength) {
						minTrip = Arrays.copyOf(nextCandidate, nextCandidate.length);
						minTripLength = newTLength;
					}
					if (newTLength > maxTripLength) {
						maxTripLength = newTLength;
						maxTrip = Arrays.copyOf(nextCandidate, nextCandidate.length);
					}
					summation += newTLength;
					sumSqrd += newTLength * newTLength;
					hist.add(newTLength);
					i++;
				} else {
					if (newTLength > maxTripLength) {
						maxTripLength = newTLength;
						maxTrip = Arrays.copyOf(nextCandidate, nextCandidate.length);
					}
					// don't increment i, keep generating new candidates until one is accepted this
					// way
				}
			}
			updateTemp();
		}

		// calc the mean
		double mean = summation / numOfTries;

		// calc standard deviation sqrt((sumOfSquares - (sumsquared / size)) - size)
		double stdDev = Math.sqrt((sumSqrd - ((summation * summation) / numOfTries)) / numOfTries);

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

		hist.logInfo("simAnnealing.csv");
		hist.showHistImage();

	}
 
	// probability function
	private static double prob(double df) {
		double oneOverKt = 1 / (temperature * Math.sqrt(2.0 * Math.PI));
		double exponent = -df / (temperature * Math.sqrt(2.0 * Math.PI));
		double px = oneOverKt * Math.exp(exponent);
		return px;
	}

	private static void updateTemp() {
		// reciprocal of logtime
		double coolingInit = initTemp / (1 + Math.log(1 + time));
		time++;
		double coolingFin = initTemp / (1 + Math.log(1 + time));

		// make proportional to temperature, coolingInit/coolingFin = temperature / x,
		// set the new temperature to x
		temperature = temperature / (coolingInit / coolingFin);

	}
	
	// change a little bit
	private static int[] randomSwap(int[] positions) {
		Random rand = new Random();

		int city1 = rand.nextInt(positions.length);
		int city2 = rand.nextInt(positions.length);
		int temp;
		temp = positions[city1];
		positions[city1] = positions[city2];
		positions[city2] = temp;

		return positions;
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
