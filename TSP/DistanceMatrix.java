package TSP;

import java.util.ArrayList;

import linAlg.Matrix;

public class DistanceMatrix extends Matrix {
	double maxDistance;
	double minDistance;

	final static String ALPHA = "ABCDEFGHIJKLMNOPQRSTUVWYXZ";

	/**
	 * Matrix that contains the distances between points, it is a symmetrical matrix
	 * 
	 * @param xValues
	 * @param yValues
	 * @param n
	 */
	public DistanceMatrix(ArrayList<Double> xValues, ArrayList<Double> yValues, int n) {
		super(n, n);
		maxDistance = 0;
		minDistance = 9999999999999.0;

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				double x1 = xValues.get(i);
				double y1 = yValues.get(i);
				double x2 = xValues.get(j);
				double y2 = yValues.get(j);
				double distance = Math.hypot(x2 - x1, y2 - y1); // calculate the distance
				A[i][j] = distance;

				// keep track of the large and smallest distances
				if (distance > maxDistance)
					maxDistance = distance;

				if (distance != 0 && distance < minDistance)
					minDistance = distance;
			}
		}

	}

	// number of positions times the max distance
	public double getUpperBound() {
		return A.length * maxDistance;
	}

	// number of positions times the min distance
	public double getLowerBound() {
		return A.length * minDistance;
	}

	// instead of integers use the alphabet to label the positions
	public void encodeLabels(int[] positions) {
		String cities = ALPHA.substring(0, A.length);
		char[] cityLabels = cities.toCharArray();
		char[] trip = new char[cityLabels.length];

		for (int i = 0; i < positions.length; i++) {
			int pos = positions[i];
			trip[i] = cityLabels[pos];
		}

		printArray(trip);

	}

	// for testing
	private void printArray(char[] cities, String mes) {
		System.out.print(mes);
		String s = "";
		for (char c : cities) {
			s += c + " ";
		}
		System.out.println(s);
	}

	private void printArray(char[] cities) {
		printArray(cities, "");
	}

	/**
	 * The value inside the distance matrix
	 * 
	 * @param city1
	 * @param city2
	 * @return
	 */
	public double calc_Distance(int city1, int city2) {
		return A[city1][city2];
	}

	/**
	 * given an array of positions, calculates the total distances to visit all of
	 * those positions and back
	 * 
	 * @param trip
	 * @return
	 */
	public double calcTripLength(int[] trip) {

		double totalDistance = 0;
		int second;
		for (int first = 0; first < trip.length - 1; first++) {
			second = first + 1;
			int city1 = trip[first];
			int city2 = trip[second];
			totalDistance += calc_Distance(city1, city2);
		}
		totalDistance += calc_Distance(trip[trip.length - 1], trip[0]); // rounds back to the start location
		return totalDistance;
	}

	// for testing
	public static void main(String[] args) {
		ArrayList<Double> xv = new ArrayList<Double>();
		ArrayList<Double> yv = new ArrayList<Double>();
		xv.add(0.5);
		yv.add(1.0);
		xv.add(2.0);
		yv.add(0.5);
		xv.add(2.0);
		yv.add(4.0);
		xv.add(10.0);
		yv.add(10.0);
		DistanceMatrix dm = new DistanceMatrix(xv, yv, 4);
		System.out.println(dm);
		System.out.println(dm.maxDistance);
		System.out.println(dm.minDistance);
		System.out.println(dm.getUpperBound());
		System.out.println(dm.getLowerBound());

		int[] test = { 0, 3, 1, 2 };
		dm.encodeLabels(test);
		int[] test2 = { 3, 1, 2, 0 };
		dm.encodeLabels(test2);
		System.out.println("Trip lenght of test 1: " + dm.calcTripLength(test));
		System.out.println("Trip length of test 2: " + dm.calcTripLength(test2));

	}

}
