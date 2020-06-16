package helperclasses;

import com.herumi.mcl.*;

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
	
	//returns the coefficients of the LaGrange polynomial
	//precondition: xValues.length == yValues.length
	//default: n = (# of data points - 1)
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
			//System.out.println(numerator[i] + " ");
		}
		
		return numerator;
		
	}
	
	
	//computes f(x) given the coefficients of f(x) and the value of x
	public static Fr computeFx(Fr[] coefficients, Fr x) {
		
		Fr sum = new Fr(0);
		int i = 0;
		
		for (int degree = coefficients.length - 1; degree >= 0; degree--) {
			Fr xExponentiated = Tools.power(x, degree);
			Fr product = new Fr();
			Mcl.mul(product, coefficients[i], xExponentiated);
			Mcl.add(sum, sum, product);
			i++;
		}
		
		return sum;
		
	}
	
	
	
	public static void main(String[] args) {
		
		System.load("/Users/arushchhatrapati/Documents/mcl/lib/libmcljava.dylib");
		Mcl.SystemInit(Mcl.BN254);
		Fr[] xValues = {new Fr(1), new Fr(2), new Fr(3), new Fr(4), new Fr(5), new Fr(6), new Fr(7)};
		Fr[] yValues = {new Fr(3), new Fr(5), new Fr(7), new Fr(9), new Fr(1), new Fr(1), new Fr(1)};
		Fr[] laGrange = laGrange(xValues, yValues, 6);
		for (Fr d: laGrange) {
			System.out.print(d + " ");
		}
	
		
	}
		
	
	
}
