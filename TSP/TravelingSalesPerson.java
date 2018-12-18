package TSP;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

public class TravelingSalesPerson {
	
	static DistanceMatrix distMatrix;
	
	public static void main(String[] args) {
		Scanner sc = new Scanner(System.in);
		System.out.println("Please input filename for TSP data: ");
		String filename = sc.nextLine();
		File tspData = new File(filename);
		setUpTSP(tspData);
		System.out.println("Which TSP whould you like to run? Enter corresponding number:\n"
				+ "1 bruteforce\n"
				+ "2 random search\n"
				+ "3 genetic algorithms\n"
				+ "4 simulated annealing\n"
				+ "5 particle swarm optimization");
		int choice = sc.nextInt();
		
		System.out.println("The distance matrix is:");
		System.out.println(distMatrix);
		System.out.println("The upper bound is: " + distMatrix.getUpperBound());
		System.out.println("The lower bound is: " + distMatrix.getLowerBound() + "\n");
		runTSP(choice, distMatrix);
		sc.close();
		

	}
	/**
	 * Get the tsp data and make the distance matrix
	 * @param positions
	 */
	public static void setUpTSP(File positions) {
		try {
			Scanner sc2 = new Scanner(positions);
			int n = 0;
			ArrayList<Double> xVals = new ArrayList<Double>();
			ArrayList<Double> yVals = new ArrayList<Double>();
			while (sc2.hasNextLine()) {
				n++;
				String line = sc2.nextLine();
				String[] data = line.split("\t");
				double x = Double.parseDouble(data[0]);
				xVals.add(x);
				double y = Double.parseDouble(data[1]);
				yVals.add(y);
			}
			
			distMatrix = new DistanceMatrix(xVals, yVals, n);
			sc2.close();
			
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
	}
	
	
	/**
	 * Run the search the user selected
	 * @param choice
	 * @param dMatrix
	 */
	public static void runTSP(int choice, DistanceMatrix dMatrix) {
		double start, stop;
		switch(choice) {
		case 1:
			start = System.currentTimeMillis();
			TSPExhaustive.exhaustive(dMatrix);
			stop = System.currentTimeMillis();
			break;
		case 2: 
			start = System.currentTimeMillis();
			TSPrandomSearch.randomSearch(dMatrix, 1000000);
			stop = System.currentTimeMillis();
			break;
		case 3:
			start = System.currentTimeMillis();
			TSPgeneticAlg.evovleTSP(dMatrix, 50000, 5);
			stop = System.currentTimeMillis();
			break;
		case 4:
			start = System.currentTimeMillis();
			TSPsimulatedAnnealing.simmulatedAnnealing(dMatrix, 0.6, 50000); 
			stop = System.currentTimeMillis();
			break;
		case 5:
			start = System.currentTimeMillis();
			TSPSwapSequencePSO.feed(dMatrix, 20000, 100, 0.4, 0.4);
			stop = System.currentTimeMillis();
			break;
		default:
			System.out.println("invalid choice decided to run random search");
			start = System.currentTimeMillis();
			TSPrandomSearch.randomSearch(dMatrix, 1000000);
			stop = System.currentTimeMillis();
		}
		double timeMilli = stop - start;
		double timeSec = timeMilli / 1000;
		int minutes = (int) timeSec / 60;
		int seconds = (int) (timeSec % 60);
		System.out.println("\nThe process took " + minutes + " minutes and " + seconds + " seconds.");
		System.out.println("\n Time in millseconds: " + timeMilli);
				
		
		
	}

}
