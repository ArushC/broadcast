package helperclasses;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;

public class LagrangeInterpolation {
	
	//multiplies two polynomials
	//input: coefficients of 1st polynomial, coefficients of second polynomial
	//output: coefficients of resultant polynomial
	//ex. (x - 1) *  (x - 3) = x^2 - 4x + 3  --> multiply([1, -1], [1, -3]) = [1, -4, 3]
	public static BigDecimal[] multiply(BigDecimal[] c1, BigDecimal[] c2) {
		
		BigDecimal[] result = new BigDecimal[c1.length + c2.length - 1];
		for (int i = 0; i < result.length; i++) 
			result[i] = new BigDecimal(0); //fill with 0s to prevent NullPointerException later
		for (int i = 0; i < c1.length; i++) {
			for (int j = 0; j < c2.length; j++) {
				int degree = c1.length + c2.length - i - j - 2;
				result[result.length - degree - 1] = result[result.length - degree - 1].add(c1[i].multiply(c2[j]));	
			}
		}
		
		return result;
		
	}
	
	//adds two polynomials
	//similar functionality as multiply
	public static BigDecimal[] add(BigDecimal[] c1, BigDecimal[] c2) {
		
		BigDecimal[] result = new BigDecimal[Math.max(c1.length, c2.length)];
		for (int i = 0; i < result.length; i++) 
			result[i] = new BigDecimal(0); //fill with 0s to prevent NullPointerException later
		int lengthDiff = Math.abs(c1.length - c2.length);
		if (c1.length > c2.length) {
			for (int i = 0; i < c1.length; i++)
				result[i] = c1[i];
			for (int j = c2.length - 1; j >= 0; j--) {
				result[j + lengthDiff] = result[j + lengthDiff].add(c2[j]);
			}
		}
		else if (c2.length > c1.length) {
			for (int i = 0; i < c2.length; i++)
				result[i] = c2[i];
			for (int j = c1.length - 1; j >= 0; j--) {
				result[j + lengthDiff] = result[j + lengthDiff].add(c1[j]);
			}
		}
		else {
			for (int i = 0; i < c1.length; i++) {
				result[i] = c1[i].add(c2[i]);
			}
		}
		
		return result;
	}
	
	//returns the coefficients of the LaGrange polynomial
	//precondition: xValues.length == yValues.length
	//default: n = (# of data points - 1)
	public static BigDecimal[] Lx(BigDecimal[] xValues, BigDecimal[] yValues, int n) {
		
		BigDecimal[] coefficients = new BigDecimal[n + 1];
		BigDecimal[] initialYValue = new BigDecimal[1];
		initialYValue[0] = new BigDecimal(1);
		BigDecimal[] sum = new BigDecimal[1];
		sum[0] = new BigDecimal(0);
		for (int j = 0; j <= n; j++) {
			BigDecimal[] yj = new BigDecimal[1];
			yj[0] = yValues[j];
			BigDecimal[] lx = lx(j, xValues, n);
			BigDecimal[] product = multiply(yj, lx);
			sum = add(sum, product);
		}
		
		return sum;
		
	}
	
	//LaGrange helper function
	public static BigDecimal[] lx(int j, BigDecimal[] xValues, int n) {
		
		BigDecimal[] numerator = {new BigDecimal(1.0)};
		BigDecimal denominator = new BigDecimal(1.0);
		for (int m = 0; m <= n; m++) {
			if (m == j)
				continue;
			//else
			BigDecimal[] next = {new BigDecimal(1.0), xValues[m].multiply(new BigDecimal(-1))};
			numerator = multiply(numerator, next);
			denominator = denominator.multiply((xValues[j].subtract(xValues[m])));
		}
		
		for (int i = 0; i < numerator.length; i++) {
			numerator[i] = numerator[i].divide(denominator, MathContext.DECIMAL128);
			//System.out.println(numerator[i] + " ");
		}
		
		return numerator;
		
	}
	
	
	public static void main(String[] args) {
		
	//x: 1, 2, 3, 6
	//y: 1, 6, -2, -2
		
	BigDecimal[] xValues = {new BigDecimal(2), new BigDecimal(8), new BigDecimal(5), new BigDecimal(7), new BigDecimal(6), new BigDecimal(-6), new BigDecimal(3)};
	BigDecimal[] yValues = {new BigDecimal(3), new BigDecimal(5), new BigDecimal(7), new BigDecimal(9), new BigDecimal(1), new BigDecimal(1), new BigDecimal(1)};
	//BigDecimal[] zValues = multiply(xValues, yValues);
	//xValues[0] = xValues[0].add(new BigDecimal(4));
	/*BigDecimal[] add = multiply(xValues, yValues);
	for (BigDecimal b: add) {
		System.out.print(b + " ");
	}*/
	BigDecimal[] laGrange = Lx(xValues, yValues, 6);
	for (BigDecimal d: laGrange) {
		System.out.printf("%.4f ", d);
	}
	
		
	}
	
	
		
	
	
}
