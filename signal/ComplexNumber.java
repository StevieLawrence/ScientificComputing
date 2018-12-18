package signal;

public class ComplexNumber {

	double real, imaginary;

	public ComplexNumber(double real, double imaginary) {
		this.real = real;
		this.imaginary = imaginary;
	}

	public ComplexNumber() {
		real = 0;
		imaginary = 0;
	}

	public static ComplexNumber add(ComplexNumber c1, ComplexNumber c2) {
		return new ComplexNumber(c1.real + c2.real, c1.imaginary + c2.imaginary);
	}

	public ComplexNumber add(ComplexNumber c) {
		return new ComplexNumber(real + c.real, imaginary + c.imaginary);
	}

	public void addTo(ComplexNumber c) {
		real = real + c.real;
		imaginary = imaginary + c.imaginary;
	}

	public static ComplexNumber sub(ComplexNumber c1, ComplexNumber c2) {
		return new ComplexNumber(c1.real - c2.real, c1.imaginary - c2.imaginary);
	}

	public ComplexNumber sub(ComplexNumber c) {
		return new ComplexNumber(real - c.real, imaginary - c.imaginary);
	}

	public void minus(ComplexNumber c) {
		real = real - c.real;
		imaginary = imaginary - c.imaginary;
	}

	public static ComplexNumber multiply(ComplexNumber c1, ComplexNumber c2) {
		return new ComplexNumber(c1.real * c2.real - c1.imaginary * c2.imaginary,
				c1.real * c2.imaginary + c1.imaginary * c2.real);
	}

	public ComplexNumber multiply(ComplexNumber c) {
		return new ComplexNumber(real * c.real - imaginary * c.imaginary, real * c.imaginary + imaginary * c.real);
	}

	public void scale(ComplexNumber c) {
		real = real * c.real - imaginary * c.imaginary;
		imaginary = real * c.imaginary + imaginary * c.real;
	}
	
	public void scale(double x) {
		real = real * x;
		imaginary = imaginary * x;
	}

	/**
	 * Returns the reciprocal of this complex number.
	 *
	 * @return the complex number whose value is {@code (1 / this)}
	 */
	public ComplexNumber reciprocal() {
		double denom = real * real + imaginary * imaginary;
		if (denom != 0)
			return new ComplexNumber(real / denom, -imaginary / denom);
		else 
			throw new java.lang.ArithmeticException("/ by zero");
	}
	
	public static ComplexNumber divison(ComplexNumber c1, ComplexNumber c2) {
		return c1.multiply(c2.reciprocal());
	}
	
	public ComplexNumber divide(ComplexNumber c) {
		return this.multiply(c.reciprocal());
	}
	
	public void divideBy(ComplexNumber c) {
		this.scale(c.reciprocal());
	}

	/**
	 * Modulus of this Complex number (the distance from the origin in polar
	 * coordinates).
	 * 
	 * @return |z| where z is this Complex number.
	 */
	public double mod() {
		return Math.sqrt(real * real + imaginary * imaginary);
		
	}

	/**
	 * Argument of this Complex number (the angle in radians with the x-axis in
	 * polar coordinates).
	 * 
	 * @return arg(z) where z is this Complex number.
	 */
	public double argument() {
		return Math.atan2(imaginary, real);
	}
	
	public double getReal() {
		return real;
	}
	
	public double getImaginary() {
		return imaginary;
	}
	
	public String toString() {
		return "" + real + " + " + imaginary + "i"; 
	}
	
	 public ComplexNumber exp() {
	        return new ComplexNumber(Math.exp(real) * Math.cos(imaginary), Math.exp(real) * Math.sin(imaginary));
	    }
	 
	
	    public ComplexNumber sin() {
	        return new ComplexNumber(Math.sin(real) * Math.cosh(imaginary), Math.cos(real) * Math.sinh(imaginary));
	    }

	    /**
	     * Returns the complex cosine of this complex number.
	     *
	     * @return the complex cosine of this complex number
	     */
	    public ComplexNumber cos() {
	        return new ComplexNumber(Math.cos(real) * Math.cosh(imaginary), -Math.sin(real) * Math.sinh(imaginary));
	    }

	    /**
	     * Returns the complex tangent of this complex number.
	     *
	     * @return the complex tangent of this complex number
	     */
	    public ComplexNumber tan() {
	        return sin().divide(cos());
	    }
	    
	    public double expForm() {
	    	return mod() * Math.exp(argument());
	    }
	    
	    public static void printArray(ComplexNumber[] cn) {
	    	for (int i = 0; i < cn.length; i++) {
	    		System.out.println("k = " + i + " | " + cn[i]);
	    	}System.out.println();
	    }
	    
	    public ComplexNumber conjugate() {
	    	return new ComplexNumber(real, -1 * imaginary);
	    }
	    
	    public static double[] getRealVals(ComplexNumber[] complexInfo) {
			double[] d = new double[complexInfo.length];
			for (int i = 0; i < complexInfo.length; i++) {
				d[i] = complexInfo[i].getReal();
			}
			return d;
		}
	    
	    public static double[] getImgVals(ComplexNumber[] complexInfo) {
			double[] d = new double[complexInfo.length];
			for (int i = 0; i < complexInfo.length; i++) {
				d[i] = complexInfo[i].getImaginary();
			}
			return d;
		}
	 

}
