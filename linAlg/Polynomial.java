package linAlg;

import java.util.HashSet;

public class Polynomial {
	double[] coefficients;
	int deg;

	public Polynomial(double[] coeff, int deg) {
		coefficients = coeff;
		this.deg = deg;
	}

	public Polynomial(double[] coeff) {
		this(coeff, coeff.length - 1);
	}
	
	public Polynomial(Polynomial poly) {
		deg = poly.deg;
		coefficients = poly.coefficients.clone();
	}

	/**
	 * The quadratic formula
	 * 
	 * @param a
	 * @param b
	 * @param c
	 * @return
	 */
	public static double[] quadratic(double a, double b, double c) {
		double r1, r2, discriminant;

		discriminant = (b * b) - (4 * a * c);

		if (discriminant > 0) {
			r1 = (-b + Math.sqrt(discriminant)) / (2 * a);
			r2 = (-b - Math.sqrt(discriminant)) / (2 * a);
			System.out.println(String.format("Roots are real and distinct: root1 = %f, root2 = %f", r1, r2));

		} else if (discriminant == 0) {
			r1 = r2 = (-b + Math.sqrt(discriminant)) / (2 * a);
			System.out.println(String.format("Roots are real and repeated: root1 = root2 = %f", r1));

		} else {
			r1 = r2 = Double.NaN;
			System.out.println("Roots are imaginary ie NaN");
		}
		double[] roots = { r1, r2 };
		return roots;
	}

	public double[] quadratic() {
		if (deg == 2) {
			double r1, r2, discriminant;
			double a = coefficients[0];
			double b = coefficients[1];
			double c = coefficients[2];

			discriminant = (b * b) - (4 * a * c);

			if (discriminant > 0) {
				r1 = (-b + Math.sqrt(discriminant)) / (2 * a);
				r2 = (-b - Math.sqrt(discriminant)) / (2 * a);
				System.out.println(String.format("Roots are real and distinct: root1 = %f, root2 = %f", r1, r2));

			} else if (discriminant == 0) {
				r1 = r2 = (-b + Math.sqrt(discriminant)) / (2 * a);
				System.out.println(String.format("Roots are real and repeated: root1 = root2 = %f", r1));

			} else {
				r1 = r2 = Double.NaN;
				System.out.println("Roots are imaginary ie NaN");
			}
			double[] roots = { r1, r2 };
			return roots;
		} else {
			System.out.println("Polynomial not degree of 2");
			return null;
		}
	}

	/**
	 * Evaluates a given input to the polynomial
	 * 
	 * @param input
	 * @return
	 */
	public double solve(double input) {
		double summation = 0;
		int power = deg;
		for (double co : coefficients) {
			double val = Math.pow(input, power);
			val *= co;
			summation += val;
			power--;
		}
		return summation;
	}

	public static double solve(Polynomial poly, double input) {

		double summation = 0;
		int power = poly.deg;
		for (double co : poly.coefficients) {
			double val = Math.pow(input, power);
			val *= co;
			summation += val;
			power--;
		}
		return summation;
	}

	/**
	 * The remainder of the division of a polynomial by a linear polynomial such as
	 * (x - a) is equal to plugging a into the equation. (roots should have
	 * remainders of zero)
	 * 
	 * @param poly
	 * @param root
	 * @param remainder
	 * @return
	 */
	public static boolean remainderThm(Polynomial poly, double root, double remainder) {
		double ans = poly.solve(root);
		return ans == remainder;
	}

	private static boolean nearZero(double value, double epsilon) {
		if (value > Math.abs(epsilon) || value < -1 * Math.abs(epsilon)) {
			return false;
		}else
			return true;
	}

	/**
	 * Finds all of the possible rational roots of the number given polynomial
	 * coefficients. Finds all of the factors of the constant and the leading
	 * coefficient. The quotients if dividing the constant's factors with the
	 * leading coefficient's factors are the possible rational roots of the
	 * polynomial. (a polynomial can have irrational roots as well) Then remainder
	 * Theorem can be applied to determine if the canidate is a rational root of the
	 * polynomial.
	 * 
	 * @param poly
	 * @return
	 */
	public HashSet<Double> rationalRootThm() {
		HashSet<Double> rationalRoots = new HashSet<Double>();
		HashSet<Double> constFactors = new HashSet<Double>();
		HashSet<Double> coeffFactors = new HashSet<Double>();

		double constant = coefficients[coefficients.length - 1];
		double leadingCoeff = coefficients[0];

		// get the factors of the constant
		for (int i = 1; i <= Math.abs(constant); i++) {
			if (constant % i == 0) {
				constFactors.add((double) i);
			}
		}

		// get the factors of the leading coefficient
		for (int i = 1; i <= Math.abs(leadingCoeff); i++) {
			if (leadingCoeff % i == 0) {
				coeffFactors.add((double) i);
			}
		}

		// make sure neither is empty
		if (constFactors.isEmpty() || coeffFactors.isEmpty()) {
			System.out.println("Numbers dont have factors");
		} else {
			// add all of the potential canidates
			for (double con : constFactors) {
				for (double coef : coeffFactors) {
					double possibleRoot = con / coef;
					rationalRoots.add(possibleRoot);
					rationalRoots.add(-1.0 * possibleRoot);
				}
			}
			System.out.println("Potential candidates: " + rationalRoots);
			HashSet<Double> removeSet = new HashSet<Double>();
			// find the roots from the candidates
			for (double candidate : rationalRoots) {
				boolean zeroTest = Polynomial.remainderThm(this, candidate, 0);
				if (!zeroTest) {
					removeSet.add(candidate);
				}
			}
			// remove all of the invalid candidates
			rationalRoots.removeAll(removeSet);
		}

		return rationalRoots;
	}

	public static HashSet<Double> rationalRootThm(Polynomial poly) {
		HashSet<Double> rationalRoots = poly.rationalRootThm();
		return rationalRoots;
	}

	public Polynomial syntheticRootDiv(double root) {
		double EPSILON = 0.0001; // error allowance for comparing floating point numbers

		// Use remainder theorem to check if it is a root
		double rootCheck = solve(root);
		if (!Polynomial.nearZero(rootCheck, EPSILON)) {
			System.out.println("Given a non-root for the polynomial");
			return null;
		} else {
			Polynomial deflatedPoly;
			
			double [] divisidand = coefficients;
			double[] coeff = new double[coefficients.length - 1];
			// loop through the divisidand
			double newCoeff = 0;
			// ignore the last one
			for (int i = 0; i < divisidand.length - 1; i++) {
				newCoeff += divisidand[i];
				coeff[i] = newCoeff;
				newCoeff *= root;
			}
			
			deflatedPoly = new Polynomial(coeff);
			return deflatedPoly;
		}
	}
	
	public boolean isMonic() {
		if (coefficients[0] == 1)
			return true;
		else
			return false;
	}
	
	public void makeMonic() {
		if (isMonic()) {
			return;
		} else {
			double leadingCo = coefficients[0];
			for (int i = 0; i < coefficients.length; i++) {
				coefficients[i] = coefficients[i] / leadingCo;
			}
		}
	}
	
	public String toString() {
		String result = "";
		for (int i = 0; i < coefficients.length; i++) {
			result += coefficients[i] + "x^" + (deg - i) + "  ";
		}
		return result;
	}

}
