/**
 * 
 */
package linAlg;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Scanner;

/**
 * @author dmx
 *
 */
public class Demo2 {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		File file = new File("eigendataS2018.txt");
		try {
			Matrix[] vectors = new Matrix[202];
			double x,y;
			int u = 0;
			Scanner sc = new Scanner(file);
			sc.nextLine();
			sc.nextLine();
			while (sc.hasNext()) {
				String line = sc.nextLine();
				String[] data = line.split("\t");
				x = Double.parseDouble(data[0]);
				y = Double.parseDouble(data[1]);
				Matrix point = new Matrix(2,1);
				point.setEntry(1, 1, x);
				point.setEntry(2, 1, y);
				vectors[u] = point;
				u++;
			}
			// find the mean
			Matrix mean = new Matrix(2,1);
			for (int k = 0; k < vectors.length; k++) {
				mean.add(vectors[k]);
			}
			double scaler = 1.0 / (double) vectors.length;
			mean.scale(scaler);
			
			System.out.println("ai.) ----------------------");
			System.out.println("The mean vector is: ");
			System.out.println(mean);
			
			// find the covariance matrix
			Matrix covariance = new Matrix(2,2);
			for (int v = 0; v < vectors.length; v++) {
				//subtract away the mean vectors
				Matrix deviation = Matrix.getDiff(vectors[v], mean);
				
				// get the transpose matrix
				Matrix tran = Matrix.getTranspose(deviation);
				
				// square them by multiplying the deviation matrix by its transpose matrix (2x1 and 1x2 gives 2x2)
				deviation.times(tran);
				
				// accumalate 
				covariance.add(deviation);
			}
			// scale it
			covariance.scale(scaler);
			System.out.println("\nThe Covariant matrix is: ");
			System.out.println(covariance);
			
			System.out.println("aii.) ----------------------");
			double trace = covariance.trace();
			System.out.println("The trace of the covariance matrix = " + trace);
			
			System.out.println("\naiii.) ----------------------");
			double determinant = covariance.determinant();
			System.out.println("The determinant of the covariance matrix = " + determinant);

			System.out.println("\naiv.) ----------------------");
			System.out.println("LeVier's method give characteristic polynomial:");
			Polynomial levierCoeff = covariance.LeVerrier();
			System.out.println(levierCoeff);
			
			System.out.println("\nThe companion Matrix is:");
			Matrix covarCompanion = new Matrix(levierCoeff);
			System.out.println(covarCompanion);
			
			System.out.println("Sovling the characteristic equation with quadratic equation:");
			double[] coQuadEigVals = levierCoeff.quadratic();
			System.out.println("eigen values are: ");
			for (int i = 0; i < coQuadEigVals.length; i++) {
				System.out.println("EigenValue" + (i+1) + " = " + coQuadEigVals[i]);
			}
			
			System.out.println("\nUsing rational root theorem (non-integer so it should fail):");
			System.out.println(levierCoeff.rationalRootThm());
			
			System.out.println("\nPower Method");
			Matrix covarInit = new Matrix(2,1);
			covarInit.setEntry(1,1, 1);
			covarInit.setEntry(2, 1, 1);
			double coPowerEig1 = covariance.powerMethod(covarInit, 0.0000001, 100);
			System.out.println("Largest Eigenvalue is: " + coPowerEig1);
			
			System.out.println("\n inverse Powermethod");
			double coInvPow = covariance.inversePowerMethod(covarInit, 0.0000001, 100);
			System.out.println("The smallest Eigenvalue is: " + coInvPow);
			
			System.out.println("\nPolynomial deflation");
			Polynomial coPolyDeflate = levierCoeff.syntheticRootDiv(coPowerEig1);
			System.out.println(coPolyDeflate);
			
			System.out.println("\nQR Method Check");
			covariance.qr(0.0000001, 100);	
			
			System.out.println("\nav.) ----------------------");
			
			// find the charecteristic equation
			Polynomial compCoeff = covarCompanion.LeVerrier();
			System.out.println(compCoeff);
			
			// find its eigen vectors using the power methods, Don't have to do deflation cause there is only 2 eigenvalue/vectors
			System.out.println(covarCompanion);
			System.out.println("Eigenvector 1");
			double compEig1 = covarCompanion.powerMethod(covarInit, 0.0000001, 100); // larger eigenvalue
			System.out.println("Check (should be 0.741556): " + compEig1); 
			
			System.out.println("\n Eigenvector 2");
			covarInit.setEntry(1, 1, -1);
			double compEig2 = covarCompanion.inversePowerMethod(covarInit, 0.0000001, 100); // smaller eigenvalue
			System.out.println("Check (should be 0.0401256): " + compEig2);
			
			System.out.println("Eigenvector 1 is (from solving by hand):");
			Matrix coEigVect1 = new Matrix(2,1);
			coEigVect1.setEntry(1, 1, 0.886854);
			coEigVect1.setEntry(2, 1, 0.462049);
			System.out.println(coEigVect1);
			
			System.out.println("Eigenvector 2 is (from solving by hand):");
			Matrix coEigVect2 = new Matrix(2,1);
			coEigVect2.setEntry(1, 1, -0.462049);
			coEigVect2.setEntry(2, 1, 0.886854);
			System.out.println(coEigVect2);
			
			System.out.println("Translating the eigenvectors by the mean:");
			System.out.println(mean);
			coEigVect1.add(mean);
			coEigVect2.add(mean);
			
			System.out.println("Eigenvector 1 translated is:");
			System.out.println(coEigVect1);
			System.out.println("Eigenvector 2 translated is:");
			System.out.println(coEigVect2);
			
			System.out.println("2.) ********************************");
			System.out.println("Problem 2\n");
			double[] pXcoeff = {30, -139, -1689, 4903, -2733, -756};
			Polynomial px = new Polynomial(pXcoeff);
			System.out.println("For polynomial:");
			System.out.println(px);
			
			System.out.println("\nThe monic poynomial is:");
			px.makeMonic();
			System.out.println(px);
			Matrix pxMat2 = new Matrix(px);
			
			int continueDeflate = px.deg - 1;
			double[] eigenVals = new double[5];
			for (int k = 0; k < continueDeflate; k++) {
				System.out.println(String.format("\nfor iteration %d the polynomial is:", k+1));
				System.out.println(px);

				System.out.println(String.format("\nThe Companion matrix for iteration %d is:", k+1));
				Matrix pxMat = new Matrix(px);
				System.out.println(pxMat);

				System.out.println(String.format("The characteristic equajtion is (should spit back the same polynomial) for iteration %d is:",k+1));
				Polynomial charEq = pxMat.LeVerrier();
				System.out.println(charEq);

				Matrix initVect = new Matrix(pxMat.getRowCount(), 1);
				for (int i = 1; i <= initVect.getRowCount(); i++) {
					initVect.setEntry(i, 1, 1);
				}

				double pxEigen = pxMat.powerMethod(initVect, 0.0000000000000000000000001, 1000);
				System.out.println("The Eigenvalue for iteration " + (k+1) + " is: " + pxEigen);
				eigenVals[k] = pxEigen;
				
				px = px.syntheticRootDiv(pxEigen);
			}
			System.out.println("\nThe final root is: ");
			System.out.println(px);
			eigenVals[4] = -1 * px.coefficients[px.coefficients.length - 1];
			System.out.println("\n All of the eigenvalues are:");
			double summ = 0;
			for (double ei : eigenVals) {
				summ += ei;
				System.out.print(ei + " ");
			}System.out.println();
			
			System.out.println("\nTrace check");
			System.out.println(String.format("The added sum of the eigen values is %f and trace of original companion matrix is %f", summ, pxMat2.trace()));
			
			System.out.println("\n QR method check");
			pxMat2.qr(0.0000000000000000000000001, 1000);
			
		}catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
