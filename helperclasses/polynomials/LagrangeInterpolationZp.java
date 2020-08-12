package helperclasses.polynomials;
import java.io.File;
import com.herumi.mcl.Fr;
import com.herumi.mcl.Mcl;

public class LagrangeInterpolationZp {
	
	//multiplies two polynomials
	//input: coefficients of 1st polynomial, coefficients of second polynomial
	//output: coefficients of resultant polynomial
	//ex. (x - 1) *  (x - 3) = x^2 - 4x + 3  --> multiply([1, -1], [1, -3]) = [1, -4, 3]
	public static Fr[] multiply(Fr[] c1, Fr[] c2) {
		
		Fr[] result = new Fr[c1.length + c2.length - 1];
		for (int i = 0; i < result.length; i++) 
			result[i] = new Fr(0); //fill with 0s to prevent NullPointerException later
		for (int i = 0; i < c1.length; i++) {
			for (int j = 0; j < c2.length; j++) {
				int degree = c1.length + c2.length - i - j - 2;
				Fr product = new Fr();
				Mcl.mul(product, c1[i], c2[j]);
				Mcl.add(result[result.length - degree - 1], result[result.length - degree - 1], product);	
			}
		}
		
		return result;
		
	}
	
	//adds two polynomials
	//similar functionality as multiply
	public static Fr[] add(Fr[] c1, Fr[] c2) {
		
		Fr[] result = new Fr[Math.max(c1.length, c2.length)];
		for (int i = 0; i < result.length; i++) 
			result[i] = new Fr(0); //fill with 0s to prevent NullPointerException later
		int lengthDiff = Math.abs(c1.length - c2.length);
		if (c1.length > c2.length) {
			for (int i = 0; i < c1.length; i++)
				result[i] = c1[i];
			for (int j = c2.length - 1; j >= 0; j--) {
				Mcl.add(result[j + lengthDiff], result[j + lengthDiff], c2[j]);
			}
		}
		else if (c2.length > c1.length) {
			for (int i = 0; i < c2.length; i++)
				result[i] = c2[i];
			for (int j = c1.length - 1; j >= 0; j--) {
				Mcl.add(result[j + lengthDiff], result[j + lengthDiff], c1[j]);
			}
		}
		else {
			for (int i = 0; i < c1.length; i++) {
				Mcl.add(result[i], c1[i], c2[i]);
			}
		}
		
		return result;
	}
	
	public static Fr[] negate(Fr[] polynomial) {
		Fr[] result = new Fr[polynomial.length];
		for (int i = 0; i < result.length; i++) {
			result[i] = new Fr();
			Mcl.mul(result[i], polynomial[i], new Fr(-1));
		}
		return result;
	}
	
	//returns the coefficients of the LaGrange polynomial
	//precondition: xValues.length == yValues.length
	//default: n = (# of data points - 1)
	//Runtime: O(n^3)
	public static Fr[] laGrange(Fr[] xValues, Fr[] yValues, int n) {
		
		Fr[] coefficients = new Fr[n + 1];
		Fr[] initialYValue = new Fr[1];
		initialYValue[0] = new Fr(1);
		Fr[] sum = new Fr[1];
		sum[0] = new Fr(0);
		for (int j = 0; j <= n; j++) {
			Fr[] yj = new Fr[1];
			yj[0] = yValues[j];
			Fr[] lx = lx(j, xValues, n);
			Fr[] product = multiply(yj, lx);
			sum = add(sum, product);
		}
		
		return sum;
		
	}
	
	
	//LaGrange helper function
	public static Fr[] lx(int j, Fr[] xValues, int n) {
		
		Fr[] numerator = {new Fr(1)};
		Fr denominator = new Fr(1);
		for (int m = 0; m <= n; m++) {
			if (m == j)
				continue;
			//else
			Fr secondElement = new Fr();
			Mcl.mul(secondElement, xValues[m], new Fr(-1));
			Fr[] next = {new Fr(1), secondElement};
			numerator = multiply(numerator, next);
			Fr difference = new Fr();
			Mcl.sub(difference, xValues[j],xValues[m]);
			Mcl.mul(denominator, denominator, difference);
		}
		
		for (int i = 0; i < numerator.length; i++) {
			Mcl.div(numerator[i], numerator[i], denominator);
		}
		
		return numerator;
		
	}
	
	//Lagrange Interpolation faster method
	//helper function to precompute the numerator and denominator
	private static Fr[] precomputedNumerator(Fr[] xVals, int n) {
		Fr[] numerator = {new Fr(1)};
		for (int m = 0; m <= n; m++) {
			Fr secondElement = new Fr();
			Mcl.mul(secondElement, xVals[m], new Fr(-1));
			Fr[] next = {new Fr(1), secondElement};
			numerator = multiply(numerator, next);
		}
		return numerator;
	}
	
	private static Fr[] lxFaster(int j, Fr[] xVals, Fr[] numerator, int n) {
		Fr[] polynomial = LagrangeInterpolationZp.syntheticDivideWithoutRemainder(numerator, xVals[j]);
		
		//extract the polynomial from the divided numerator
		Fr denominator = LagrangeInterpolationZp.computeFxHorner(polynomial, xVals[j]);

		for (int i = 0; i < polynomial.length; i++) {
			Mcl.div(polynomial[i], polynomial[i], denominator);
			
		}

		return polynomial;
	}
	
	//Time: O(n^2)
	public static Fr[] laGrangeFaster(Fr[] xValues, Fr[] yValues, int n) {
		
		Fr[] coefficients = new Fr[n + 1];
		Fr[] numerator = precomputedNumerator(xValues, n);
		Fr[] initialYValue = {new Fr(1)};
		Fr[] sum = {new Fr(0)};
		
		for (int j = 0; j <= n; j++) {
			Fr[] yj = {yValues[j]};
			Fr[] lx = lxFaster(j, xValues, numerator, n);
			Fr[] product = multiply(yj, lx);
			sum = add(sum, product);
		}
		
		return sum;
		
	}
	
	//computes f(x) given the coefficients of f(x) and the value of x
	//Time: O(n log(n))
	public static Fr computeFx(Fr[] coefficients, Fr x) {
		Fr sum = new Fr(0);
		Fr xExponentiated = new Fr(1);
		for (int degree = coefficients.length - 1; degree >= 0; degree--) {
			Fr product = new Fr();
			Mcl.mul(product, coefficients[degree], xExponentiated);
			Mcl.add(sum, sum, product);
			Mcl.mul(xExponentiated, xExponentiated, x);
		}
			
		return sum;
	}
	
	//computes f(x) given the coefficients of f(x) and the value of x -- fastest method
	//Time: O(n)
	public static Fr computeFxHorner(Fr[] coefficients, Fr x) {
			
			Fr result = new Fr(coefficients[0]);
			for (int i = 1; i < coefficients.length; i++) {
				Mcl.mul(result, result, x);
				Mcl.add(result, result, coefficients[i]);
			}
				
			return result;
		}
	
	//computes f(x) given the roots of a polynomial in form (x - a)(x - b)...
	public static Fr computeFxFromRoots(Fr[] roots, Fr x) {
		Fr product = new Fr(1);
		for (int i = 0; i < roots.length; i++) {
			Fr diff = new Fr();
			Mcl.sub(diff, x, roots[i]);
			Mcl.mul(product, product, diff);
		}
		return product;
	}
	
	//returns an Fr[] containing the resultant polynomial and the remainder
	//ex. if root = 1 then dividing by (x - 1)
	public static Fr[] syntheticDivide(Fr[] polynomial, Fr root) {
		Fr[] result = new Fr[polynomial.length];
		result[0] = new Fr(polynomial[0]);
		Fr previous = result[0];
		for (int i = 1; i < polynomial.length; i++) {
			Fr addend = new Fr();
			Mcl.mul(addend, root, previous);
			Mcl.add(addend, addend, polynomial[i]);
			result[i] = addend;
			previous = result[i];
		}
		return result;
	}
	
		//returns an Fr[] containing the resultant polynomial
		//ex. if root = 1 then dividing by (x - 1)
		public static Fr[] syntheticDivideWithoutRemainder(Fr[] polynomial, Fr root) {
			Fr[] result = new Fr[polynomial.length - 1];
			result[0] = new Fr(polynomial[0]);
			Fr previous = result[0];
			for (int i = 1; i < polynomial.length - 1; i++) {
				Fr addend = new Fr();
				Mcl.mul(addend, root, previous);
				Mcl.add(addend, addend, polynomial[i]);
				result[i] = addend;
				previous = result[i];
			}
			return result;
		}
	
	public static void main(String[] args) {
		File lib = new File("../../lib/libmcljava.dylib");
		System.load(lib.getAbsolutePath());
		Mcl.SystemInit(Mcl.BN254);
		//check how fast the new LaGrange evaluation algorithm is compared to the old one
		Fr[] xValues = new Fr[101];
		Fr[] yValues = new Fr[101];
		for (int i = 0; i < 101; i++) {
			xValues[i] = new Fr(i + 1);
			Fr c = new Fr();
			c.setByCSPRNG();
			yValues[i] = c;
		}
		long start = System.nanoTime();
		Fr[] laGrange = laGrangeFaster(xValues, yValues, 100);
		double secondsFaster =  ((double) (System.nanoTime() - start))/1E9;
		long start2 = System.nanoTime();
		Fr[] laGrange2 = laGrange(xValues, yValues, 100);
		double secondsFaster2 =  ((double) (System.nanoTime() - start2))/1E9;
		System.out.println("Old LaGrange: " + secondsFaster2);
		System.out.println("New LaGrange: " + secondsFaster);
		
	}
	
}
