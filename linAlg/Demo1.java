package linAlg;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

public class Demo1 {

	public static void main(String[] args) {
		File f = new File("2018 Spring Project 1 data.txt");
		Matrix[] class1 = new Matrix[110];
		Matrix[] class2 = new Matrix[110];
		double x1, y1, x2, y2;
		int u = 0;
		try {
			Scanner sc = new Scanner(f);
			sc.nextLine();
			sc.nextLine();
			while (sc.hasNext()) {
				String line = sc.nextLine();
				String[] data = line.split("\t");
				x1 = Double.parseDouble(data[0]);
				y1 = Double.parseDouble(data[1]);
				x2 = Double.parseDouble(data[2]);
				y2 = Double.parseDouble(data[3]);
				Matrix p1 = new Matrix(2, 1);
				p1.setEntry(1, 1, x1);
				p1.setEntry(2, 1, y1);
				class1[u] = p1;
				Matrix p2 = new Matrix(2, 1);
				p2.setEntry(1, 1, x2);
				p2.setEntry(2, 1, y2);
				class2[u] = p2;
				u++;
			}
			
			Matrix mean1 = new Matrix(2,1);
			Matrix mean2 = new Matrix(2,1);
			int k;
			for (k = 0; k < class1.length; k++) {
				mean1.add(class1[k]);
				mean2.add(class2[k]);
			}
			double scaler = 1 / (double) k;
			
			mean1.scale(scaler);
			mean2.scale(scaler);
			System.out.println("1.) ----------------------");
			System.out.println("The mean vectors are:\n");
			System.out.println(mean1);
			System.out.println(mean2);
			
			// number 2
			
			// accumalators
			Matrix coVarC1 = new Matrix(2,2);
			Matrix coVarC2 = new Matrix(2,2);
			
			for (int i = 0; i < class1.length; i++){
				//subtract away the mean vectors
				Matrix devC1 = Matrix.getDiff(class1[i], mean1);
				Matrix devC2 = Matrix.getDiff(class2[i], mean2);
				
				// get the transpose matrix
				Matrix c1DevTran = Matrix.getTranspose(devC1);
				Matrix c2DevTran = Matrix.getTranspose(devC2);
				
				// square them by multiplying the deviation matrix by its transpose matrix
				devC1.times(c1DevTran);
				devC2.times(c2DevTran);
				
				coVarC1.add(devC1);  // sum them
				coVarC2.add(devC2);
			}
				
			// multiply by scaler
			coVarC1.scale(scaler);
			coVarC2.scale(scaler);
			
			System.out.println("2.) -----------------------------");
			System.out.println("The Covariant matrices are:\n");
			System.out.println(coVarC1);
			System.out.println(coVarC2);
	
			// number 3
			double class1Det = coVarC1.determinant();
			double class2Det = coVarC2.determinant();
			System.out.println("3.) ------------------------------");
			System.out.println("The determinants are:\n");
			System.out.println(class1Det);
			System.out.println(class2Det);
			System.out.println();
			
			// number 4
			Matrix invertC1 = Matrix.getInverse(coVarC1);
			Matrix invertC2 = Matrix.getInverse(coVarC2);
			System.out.println("4.) ------------------------------");
			System.out.println("The inverse Matrices are:\n");
			System.out.println(invertC1);
			System.out.println(invertC2);
			
			// number 5 
			// g(x) = (-1/2) (x-u)^T (covariance^-1) (x-u) (-1/2)(ln(Det)) + ln(probability)
			// First Calculate (x - u) matrix and its transpose
			// then calculate (x-u)^t * invertC * (x - u)
			// retrieve the value from the resulting 1x1 matrix (labeled D below)
			// Plug it into g(x) = (-1/2) * D - (1/2) * (ln(CoDet)) + ln(0.5)
			
			// number 6
			System.out.println("6.) ------------------------------");
			Matrix m1 = new Matrix(mean1);
			Matrix m2 = new Matrix(mean2);
			// g1 function (mean1)
			Matrix g1DevM1 = Matrix.getDiff(m1, mean1);
			Matrix g2DevM1 = Matrix.getDiff(m1, mean2);
			
			// x - m1 trans
			Matrix g1TranM1 = Matrix.getTranspose(g1DevM1);
			Matrix g2TranM1 = Matrix.getTranspose(g2DevM1);
			
			// m1 D
			Matrix g1DM1 = Matrix.getMult(g1TranM1, invertC1);
			g1DM1.times(g1DevM1);
			Matrix g2DM1 = Matrix.getMult(g2TranM1, invertC2);
			g2DM1.times(g2DevM1);
			
			// get D value
			double g1D1 = g1DM1.getEntry(1,1);
			double g2D1 = g2DM1.getEntry(1, 1);
			
			// input into function
			double g1M1 = (-0.5 * g1D1) - (0.5 * Math.log(class1Det)) + Math.log(0.5);
			double g2M1 = (-0.5 * g2D1) - (0.5 * Math.log(class2Det)) + Math.log(0.5);
			
			if (g1M1 > g2M1) {
				System.out.println("Point m1 goes into class 1 (Expected for mean 1 to go into class 1)");
				System.out.println("g1(m1) = " + g1M1);
				System.out.println("g2(m1) = " + g2M1);
			} else {
				System.out.println("Weird result m1 was put into class 2");
			}
			
			// now for m2 
			Matrix g1DevM2 = Matrix.getDiff(m2, mean1);
			Matrix g2DevM2 = Matrix.getDiff(m2, mean2);
			
			// x - m2 trans
			Matrix g1TranM2 = Matrix.getTranspose(g1DevM2);
			Matrix g2TranM2 = Matrix.getTranspose(g2DevM2);
			
			// m1 D
			Matrix g1DM2 = Matrix.getMult(g1TranM2, invertC1);
			g1DM2.times(g1DevM2);
			Matrix g2DM2 = Matrix.getMult(g2TranM2, invertC2);
			g2DM2.times(g2DevM2);
			
			// get D value
			double g1D2 = g1DM2.getEntry(1,1);
			double g2D2 = g2DM2.getEntry(1, 1);
			
			// input into function
			
			double g1M2 =  (-0.5 * g1D2) - (0.5 * Math.log(class1Det)) + Math.log(0.5);
			double g2M2 =  (-0.5 * g2D2) - (0.5 * Math.log(class2Det)) + Math.log(0.5);
			
			if (g2M2 > g1M2) {
				System.out.println("Point m2 goes into class 2 (Expected for mean 2 to go into class 2)");
				System.out.println("g2(m2) = " + g2M2);
				System.out.println("g1(m2) = " + g1M2);
			} else {
				System.out.println("Weird result m2 was put into class 1");
			}
			
			// number 7
			
			// class 1
			ArrayList<double[]> missClassed1 = new ArrayList<double[]>();
			for (int i = 0; i < class1.length; i ++) {
				// x - u
								
				Matrix g1Dev = Matrix.getDiff(class1[i], mean1);
				Matrix g2Dev = Matrix.getDiff(class1[i], mean2);				
				
				// transpose
				Matrix g1Tran = Matrix.getTranspose(g1Dev);
				Matrix g2Tran = Matrix.getTranspose(g2Dev);
				
				// D matrix
				Matrix matD1 = Matrix.getMult(g1Tran, invertC1);
				matD1.times(g1Dev);
				Matrix matD2 = Matrix.getMult(g2Tran, invertC2);
				matD2.times(g2Dev);
				
				// D value
				double D1 = matD1.getEntry(1, 1);
				double D2 = matD2.getEntry(1, 1);
				
				// distinction function
				double g1Value = (-0.5 * D1) - (0.5 * Math.log(class1Det)) + Math.log(0.5);
				double g2Value = (-0.5 * D2) - (0.5 * Math.log(class2Det)) + Math.log(0.5);
				
				// classify
				if (g2Value > g1Value) {
					double[] miss = new double[4];
					miss[0] = class1[i].getEntry(1,1);
					miss[1] = class1[i].getEntry(2, 1);
					miss[2] = g1Value;
					miss[3] = g2Value;
					missClassed1.add(miss);
				}
			}
			System.out.println();
			System.out.println("7a.) -----------------------------");
			System.out.println(missClassed1.size() + " points misclassed from class 1.\n");
			for (double[] msPoint : missClassed1) {
				System.out.println("(x1 = " + msPoint[0] + " , y1 = " + msPoint[1] + ") ");
				System.out.println("g1(x) = " + msPoint[2]);
				System.out.println("g2(x) = " + msPoint[3]);
				System.out.println();
			}
			
			// class 2
			ArrayList<double[]> missClassed2 = new ArrayList<double[]>();
			for (int i = 0; i < class2.length; i ++) {
				// x - u				
				Matrix g1Deviation = Matrix.getDiff(class2[i], mean1);
				Matrix g2Deviation = Matrix.getDiff(class2[i], mean2);				
				
				// transpose
				Matrix g1Transpose = Matrix.getTranspose(g1Deviation);
				Matrix g2Transpose = Matrix.getTranspose(g2Deviation);
				
				// D matrix
				Matrix matD3 = Matrix.getMult(g1Transpose, invertC1);
				matD3.times(g1Deviation);
				Matrix matD4 = Matrix.getMult(g2Transpose, invertC2);
				matD4.times(g2Deviation);
				
				// D value
				double D3 = matD3.getEntry(1, 1);
				double D4 = matD4.getEntry(1, 1);
				
				// distinction function
				double g1Val = (-0.5 * D3) - (0.5 * Math.log(class1Det)) + Math.log(0.5);
				double g2Val = (-0.5 * D4) - (0.5 * Math.log(class2Det)) + Math.log(0.5);
				
				// classify
				if (g2Val < g1Val) {
					double[] miss = new double[4];
					miss[0] = class2[i].getEntry(1,1);
					miss[1] = class2[i].getEntry(2, 1);
					miss[2] = g1Val;
					miss[3] = g2Val;
					missClassed2.add(miss);
				}
			}
			System.out.println();
			System.out.println("7b.)--------------------------------");
			System.out.println(missClassed2.size() + " points misclassed from class 2.\n");
			for (double[] msPoint : missClassed2) {
				System.out.println("(x1 = " + msPoint[0] + " , y1 = " + msPoint[1] + ") ");
				System.out.println("g1(x) = " + msPoint[2]);
				System.out.println("g2(x) = " + msPoint[3]);
				System.out.println();
			}
			
			// number 8
			System.out.println("8.) --------------------------------");
			double epsi = 0.01;
			ArrayList<Matrix> boundary = new ArrayList<Matrix>();
			FileWriter writer = new FileWriter("src//linAlg//data.txt");
			for (double y = -6; y < 9; y += 0.01) {
				for (double x = -3.5; x < 4; x += 0.1) {
					Matrix xVector = new Matrix(2, 1);
					xVector.setEntry(1, 1, x);
					xVector.setEntry(2, 1, y);
					
					Matrix g1Diff = Matrix.getDiff(xVector, mean1);
					Matrix g2Diff = Matrix.getDiff(xVector, mean2);
					
					Matrix g1T = Matrix.getTranspose(g1Diff);
					Matrix g2T = Matrix.getTranspose(g2Diff);
					
					Matrix g1MatD = Matrix.getMult(g1T, invertC1);
					g1MatD.times(g1Diff);
					Matrix g2MatD = Matrix.getMult(g2T, invertC2);
					g2MatD.times(g2Diff);
					
					double g1D = g1MatD.getEntry(1, 1);
					double g2D = g2MatD.getEntry(1, 1);
					
					double g1 = (-0.5 * g1D) - (0.5 * Math.log(class1Det));
					double g2 = (-0.5 * g2D) - (0.5 * Math.log(class2Det));
					
					double value = Math.abs(g1 - g2);
					if (value < epsi) {
						//System.out.println(xVector);
						boundary.add(xVector);
						String data = "" + xVector.getEntry(1, 1) + "\t" + xVector.getEntry(2, 1) + "\n";
						try{
							writer.write(data);
							writer.flush();
						}
						catch(Exception e)
						{
							e.printStackTrace();
						}
					}
				}
			}
			writer.close();
			System.out.println("total Boundary points = " + boundary.size());
			
			// number 9
			System.out.println("9a.) ---------------------------------");
			/*double[][] linSysData = {{2, 1, -1, -1, 1, 0, -1, -1},
					                 {1, 0, 2, 0, -1, -2, 2, 2},
					                 {0, -2, 5, 4, -1, 0, 3, 1},
					                 {1, 1, -7, 3, 2, 1, -1, 0},
					                 {1, 1, 2, 3, -2, 2, 2, 9},
					                 {0, -3, -2, 2, 0, 2, 4, -5},
					                 {-2, 5, -1, 1, 1, 3, 0, -2},
					                 {1, 0, 1, 1, 0, 2, 1, 1}};
					                 */
			double[][] linSysData = {{0, 1, 3, -1, 1, 0, -1, -1},
			                         {5, 0, 2, 0, -1, 3, 1, 1},
			                         {2, -2, 2, -1, -1, 2, 3, 1},
			                         {1, 1, 0, 3, 2, 1, -1, 0},
			                         {4, 1, 2, 3, -2, 2, 2, 1},
			                         {-1, -3, -2, 2, 0, 2, 4, 1},
			                         {3, 5, -1, 1, 1, 3, 0, -2},
			                         {1, 0, 1, 1, 0, 2, 2, 1}};

			Matrix linSys = new Matrix(linSysData);
			System.out.println("For the linear system:");
			System.out.println(linSys);
			
			//double[][]linSysAData = {{1}, {-1}, {2}, {-2}, {3}, {-3}, {4}, {-4}};
			double[][]linSysAData = {{1}, {2}, {2}, {-2}, {1}, {7}, {14}, {6}};
			Matrix linb = new Matrix(linSysAData);
			System.out.println("and b of:");
			System.out.println(linb);
			
			Matrix ans = linSys.solve(linb, true);
			System.out.println("The solution matrix is:");
			System.out.println(ans);
			
			double linSysDet = linSys.determinant();
			System.out.println("The Determinant is: " + linSysDet);
			System.out.println();
			
			// 9b
			System.out.println("9b.) --------------------------");
			Matrix linInv = Matrix.getInverse(linSys);
			System.out.println("The inverse matrix is: ");
			System.out.println(linInv);
			
			//Matrix ttt = Matrix.getMult(linSys, linInv);
			//System.out.println(ttt);
			
			System.out.println("9c.) --------------------------");
			double linInvDet = linInv.determinant();
			System.out.println("The determinant of the inverse matrix is: " + linInvDet);
			System.out.println("1 / det(A) = " + 1/linSysDet);
			System.out.println("det(A) * det(A^-1) = " + linSysDet * linInvDet);
			System.out.println();
			
			// 9d
			System.out.println("9d.) --------------------------");  //A-1 * A * x = A-1 * b
			Matrix checkMat = Matrix.getMult(linInv, linb);
			System.out.println("Should be the solution matrix:");
			System.out.println(checkMat);
			
			// number 10
			System.out.println("10.) ---------------------------");
			double conditioning = linSys.infCond();
			System.out.println("The Condition number is: " + conditioning);
			
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
