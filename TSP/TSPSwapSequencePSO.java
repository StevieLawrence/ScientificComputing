package TSP;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class TSPSwapSequencePSO {
	static double minTripLength, maxTripLength, summation, sumSqrd;
	static double alpha, beta;
	static int[] minTrip, maxTrip;
	static Particle[] swarm;
	static int numOfCities;
	static Histogram hist;

	/**
	 * Population based optimization algorithm modeled after swarms of
	 * animal/insects in nature
	 * 
	 * @param dm
	 *            - the distance matrix having all the distances from city to city
	 * @param iterations
	 *            - number of cycles the algorithm is run
	 * @param size
	 *            - The number of particles in the swarm
	 * @param a
	 *            - alpha is the probability that ALL swap operators from (pBest(i)
	 *            - x(i)) will be added to velocity
	 * @param b
	 *            - beta is the probability that ALL swap operators from (gBest(i) -
	 *            x(i)) will be add to velocity
	 */
	public static void feed(DistanceMatrix dm, int iterations, int size, double a, double b) {
		numOfCities = dm.getColSize();
		swarm = new Particle[size];
		summation = 0;
		sumSqrd = 0;
		alpha = a;
		beta = b;
		hist = new Histogram(100, dm.getLowerBound(), dm.getUpperBound());
		initSwarm(dm);

		for (int i = 0; i < iterations; i++) {
			for (int p = 0; p < swarm.length; p++) {
				swarm[p].updatePosition();
			}
			updateVariables(dm);
		}

		double sampSize = size * iterations;
		// calc the mean
		double mean = summation / sampSize;

		// calc standard deviation sqrt((sumOfSquares - (sumsquared / size)) - size)
		double stdDev = Math.sqrt((sumSqrd - ((summation * summation) / sampSize)) / sampSize);

		// print the stats
		System.out.println("The max trip and trip length are:");
		System.out.print("Max Trip: ");
		dm.encodeLabels(maxTrip);
		System.out.println("Trip length: " + maxTripLength);

		System.out.println("\nThe min trip and trip length are:");
		System.out.print("Min Trip: ");
		dm.encodeLabels(minTrip);
		System.out.println("Trip length: " + minTripLength);

		System.out.println("\nThe mean is: " + mean);

		System.out.println("\nThe standard deviation is: " + stdDev);

		hist.logInfo("particleSwarmTSP.csv");
		hist.showHistImage();

	}

	// Generate a random permutation of cities
	public static int[] genRandomTour(int numOfCities) {
		int[] tour = new int[numOfCities];
		for (int i = 0; i < numOfCities; i++) {
			tour[i] = i;
		}
		return randomShuffle(tour);
	}

	// shuffle the cities
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

	// Create the swarm intializing the tours randomly
	private static void initSwarm(DistanceMatrix dm) {
		minTrip = new int[numOfCities];
		maxTrip = new int[numOfCities];
		maxTripLength = dm.getLowerBound();
		minTripLength = dm.getUpperBound();
		for (int i = 0; i < swarm.length; i++) {
			swarm[i] = new Particle();
			swarm[i].setFitness(dm.calcTripLength(swarm[i].getTrip()));
			swarm[i].setBestFitness(swarm[i].getFitness());

			if (swarm[i].getFitness() < minTripLength) {
				minTrip = Arrays.copyOf(swarm[i].getTrip(), swarm[i].getTrip().length);
				minTripLength = swarm[i].getFitness();
			}

			if (swarm[i].getFitness() > maxTripLength) {
				maxTrip = Arrays.copyOf(swarm[i].getTrip(), swarm[i].getTrip().length);
				maxTripLength = swarm[i].getFitness();
			}

			summation += swarm[i].getFitness();
			sumSqrd += swarm[i].getFitness() * swarm[i].getFitness();
			hist.add(swarm[i].getFitness());
		}
	}

	// Tracking variables including the best and worst tours seen so far
	private static void updateVariables(DistanceMatrix dm) {
		for (int i = 0; i < swarm.length; i++) {
			swarm[i].setFitness(dm.calcTripLength(swarm[i].getTrip()));

			if (swarm[i].getFitness() < swarm[i].getBestFitness()) { // update personal best
				swarm[i].setBestTrip(swarm[i].getTrip());
				swarm[i].setBestFitness(swarm[i].getFitness());
			}

			if (swarm[i].getFitness() < minTripLength) { // update global best
				minTrip = Arrays.copyOf(swarm[i].getTrip(), swarm[i].getTrip().length);
				minTripLength = swarm[i].getFitness();
			}

			if (swarm[i].getFitness() > maxTripLength) { // update worst solution
				maxTrip = Arrays.copyOf(swarm[i].getTrip(), swarm[i].getTrip().length);
				maxTripLength = swarm[i].getFitness();
			}

			summation += swarm[i].getFitness();
			sumSqrd += swarm[i].getFitness() * swarm[i].getFitness();
			hist.add(swarm[i].getFitness());
		}
	}

	/**
	 * A single permutation solution representing a tour of cities
	 */
	private static class Particle {
		int[] trip;
		int[] bestTrip;
		double fitness, bestFitness;
		ArrayList<int[]> velocity;

		Particle() {
			trip = genRandomTour(numOfCities);
			bestTrip = Arrays.copyOf(trip, trip.length);
			velocity = new ArrayList<int[]>();
		}

		/**
		 * Swap Sequence PSO uses the concept of swap operators and swap sequences to
		 * represent velocity, while position is the permutation of cities. A swap
		 * operator represents a single swap between two cities. A swap sequence is a
		 * collection of swap operators. A basic swap sequence is the minimal amount of
		 * swap operators needed to get generate a particular permutation.
		 */
		void updatePosition() {
			/*
			 * First calculate the new velocity with Vt = Vt-1 + alpha*r1*(Pid-Xt) +
			 * beta*r2*(Pgd - Xt) where Vt-1 is the previous iterations velocity alpha * r1
			 * * (Pid-Xt) is the cognitive component beta * r2 * (Pgd - Xt) is the social
			 * component
			 */
			Random rand = new Random();

			// find all swap operators to calculate (pbest - x(t-1))
			ArrayList<int[]> pBestMinX = new ArrayList<int[]>();
			for (int i = 0; i < numOfCities; i++) {
				if (trip[i] != bestTrip[i]) {
					int[] swapOperator = new int[2];
					swapOperator[0] = i;
					swapOperator[1] = index(bestTrip, trip[i], i);
					pBestMinX.add(swapOperator);
				}
			}

			// find all swap operators to calculate (gBest - x(t-1))
			ArrayList<int[]> gBestMinX = new ArrayList<int[]>();
			for (int i = 0; i < numOfCities; i++) {
				if (trip[i] != minTrip[i]) {
					int[] swapOperator = new int[2];
					swapOperator[0] = i;
					swapOperator[1] = index(minTrip, trip[i], i);
					gBestMinX.add(swapOperator);
				}
			}

			// input in Vt-1
			for (int[] swapOp : velocity) {
				swap(swapOp);
			}
			velocity.clear();

			// alpha is the probability that all pBests swap operators will occur and
			// likewise for gbest and beta
			for (int[] swapOp : pBestMinX) {
				if (rand.nextDouble() <= alpha) {
					velocity.add(swapOp);
				}
			}

			for (int[] swapOp : gBestMinX) {
				if (rand.nextDouble() <= beta) {
					velocity.add(swapOp);
				}
			}

			for (int[] swapOp : velocity) {
				swap(swapOp);
			}
		}

		void swap(int[] swapOp) {
			int temp = trip[swapOp[0]];
			trip[swapOp[0]] = trip[swapOp[1]];
			trip[swapOp[1]] = temp;
		}

		int index(int[] tour, int value, int defaultVal) {
			for (int i = 0; i < tour.length; i++) {
				if (tour[i] == value) {
					return i;
				}
			}
			return defaultVal;
		}

		void setFitness(double fitness) {
			this.fitness = fitness;
		}

		double getFitness() {
			return fitness;
		}

		void setBestFitness(double bestFitness) {
			this.bestFitness = bestFitness;
		}

		double getBestFitness() {
			return bestFitness;
		}

		int[] getTrip() {
			return trip;
		}

		void setBestTrip(int[] bestTrip) {
			this.bestTrip = bestTrip;
		}
	}
}
