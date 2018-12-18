package TSP;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.ArrayList;

public class TSPgeneticAlg {

	static int[][] population;
	static int mutationIntensity;
	final static int POP_SIZE = 50;
	static DistanceMatrix dMatrix;
	static Random rand;
	static double minTripLength, maxTripLength, summation, sumSqrd, avgFitness;
	static int[] minTrip, maxTrip;
	static Histogram hist;
	

	public static void evovleTSP(DistanceMatrix dm, int generationCycles, int muIntensity) {
		int BINS = 100;
		double minRange = dm.getLowerBound();
		double maxRange = dm.getUpperBound();

		if (muIntensity < 13 && muIntensity > 0)
			mutationIntensity = 12 - muIntensity;
		else
			mutationIntensity = 5;
		dMatrix = dm;
		rand = new Random();

		// randomly create a population of solutions and initialize variables
		createPop();
		minTrip = Arrays.copyOf(population[0], population[0].length);
		maxTrip = Arrays.copyOf(population[0], population[0].length);
		minTripLength = dMatrix.getUpperBound();
		maxTripLength = dMatrix.getLowerBound();
		
		/*avgFitness = 0;
		for (int i = 0; i < POP_SIZE; i++) {
			avgFitness += dMatrix.calcTripLength(population[i]);
		}
		avgFitness /= POP_SIZE;*/

		// initialize histogram
		hist = new Histogram(BINS, minRange, maxRange);

		// loop until "end time"
		for (int i = 0; i < generationCycles; i++) {

			// measure the pop
		   calcVar();

			// select two Parents to generate offspring
			int[][] newPop = new int[POP_SIZE][dMatrix.getColSize()];
	
			for (int p = 0; p < population.length; p++) {
				// choose the parents
				int[] parent1 = tournamentSelect(population);
				int[] parent2 = tournamentSelect(population);

				// perform crossover
				 //int[] child = recombination(parent1, parent2);
				int[] child = crossover(parent1, parent2);

				// mutate with some probability
				if (rand.nextInt(mutationIntensity) == 1) {
					child = mutate(child);
				}

				// add the child to the new population
				newPop[p] = child;
			}

			// update population
			population = newPop;

		}

		double sampSize = POP_SIZE * generationCycles;
		// calc the mean
		double mean = summation / sampSize;

		// calc standard deviation sqrt((sumOfSquares - (sumsquared / size)) - size)
		double stdDev = Math.sqrt((sumSqrd - ((summation * summation) / sampSize)) / sampSize);

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

		hist.logInfo("geneticAlg.csv");
		hist.showHistImage();

	}
	// make the initial pop
	private static void createPop() {
		population = new int[POP_SIZE][dMatrix.getRowSize()];

		for (int i = 0; i < POP_SIZE; i++) {
			population[i] = immigrate();
		}

	}
	
	private static void calcVar() {
		for (int i = 0; i < POP_SIZE; i++) {
			double tLength = dMatrix.calcTripLength(population[i]);
			if (tLength < minTripLength) {
				minTripLength = tLength;
				minTrip = Arrays.copyOf(population[i], population[i].length);
			}
			if (tLength > maxTripLength) {
				maxTripLength = tLength;
				maxTrip = Arrays.copyOf(population[i], population[i].length);
			}
			// update variables
			summation += tLength;
			sumSqrd += tLength * tLength;
			hist.add(tLength);
			
		}
		
	}
	
	/**
	 * Randomly select from to population and then have them duke it by choosing the winner 
	 * based on its fitness measure
	 * @param pop
	 * @return
	 */
	private static int[] tournamentSelect(int[][] pop) {
		int lim = rand.nextInt(POP_SIZE);
		double winner = dMatrix.getUpperBound();
		int index = 1;
		for (int i = 1; i <= lim; i ++) {
			double tLength = dMatrix.calcTripLength(pop[i]);
			if (tLength < winner) {
				winner = tLength;
				index = i;

			}
		}
		return pop[index];
	}
	
	

	// two point cross over
	public static int[] recombination(int[] parent1, int[] parent2) {
		Random rand = new Random();
		int[] child = new int[parent1.length];
		HashSet<Integer> cities = new HashSet<Integer>();
		HashSet<Integer> pos = new HashSet<Integer>();

		int start = rand.nextInt(parent1.length); // inclusive
		int end = rand.nextInt(parent1.length); // exclusive (no particular reason just easier)
		while (end == start) {
			end = rand.nextInt(parent1.length);
		}

		// splice out the choosen section
		int current = start;
		String s = "";
		while (current != end) {
			pos.add(current);
			cities.add(parent1[current]);
			s += parent1[current] + " ";
			if (current + 1 == parent1.length)
				current = 0;
			else
				current++;
		}
		s = s.substring(0, s.length() - 1); // remove end space

		/*
		 * // reverse the order with some probability if (rand.nextInt(2) == 0) { String
		 * s2 = ""; for (String c : s.split(" ")) { s2 = " " + c + s2; } s =
		 * s2.substring(1); }
		 */

		// add the splice to the child array
		String[] positions = s.split(" ");
		int size = positions.length;
		int i = 0;
		current = start;
		while (i < size) {

			String position = positions[i];
			child[current] = Integer.parseInt(position);
			if (current + 1 == parent1.length)
				current = 0;
			else
				current++;
			i++;
		}

		// populate the rest of the positions with cities from parent 2 that haven't
		// been used yet

		for (int j = 0; j < parent1.length; j++) {
			if (pos.contains(j)) {
				continue;
			}

			int k = 0;
			boolean flag = true;
			while (flag) {
				if (!cities.contains(parent2[k])) {
					break;
				}
				k++;
			}
			// add the chosen city to the city set and the child
			child[j] = parent2[k];
			cities.add(parent2[k]);
		}
		return child;
	}
	
	// simple single point cross
	public static int[] crossover(int[] parent1, int[] parent2) {
		HashSet<Integer> cities = new HashSet<Integer>();
		int[] child = new int[parent1.length];

		int mid = rand.nextInt(parent1.length);

		for (int i = 0; i < mid; i++) {
			child[i] = parent1[i];
			cities.add(parent1[i]);
		}

		for (int j = 0; j < parent1.length; j++) {
			if (j < mid)
				continue;

			int k = 0;
			while (true) {
				if (!cities.contains(parent2[k])) {
					break;
				}
				k++;
			}
			// add the chosen city to the city set and the child
			child[j] = parent2[k];
			cities.add(parent2[k]);
		}

		return child;
	}

	public static int[] mutate(int[] trip) {
		int city1 = rand.nextInt(trip.length);
		int city2 = rand.nextInt(trip.length);
		int temp;
		temp = trip[city1];
		trip[city1] = trip[city2];
		trip[city2] = temp;

		return trip;
	}

	public static int[] immigrate() {

		int[] trip = new int[dMatrix.getRowSize()];
		for (int k = 0; k < trip.length; k++) {
			trip[k] = k;
		}

		int temp;
		for (int i = 0; i < trip.length; i++) {
			int r = rand.nextInt(trip.length);
			temp = trip[i];
			trip[i] = trip[r];
			trip[r] = temp;

		}
		return trip;

	}

	private static void printArray(int[] cities, String mes) {
		System.out.print(mes);
		String s = "";
		for (int c : cities) {
			s += c + " ";
		}
		System.out.println(s);
	}

}
