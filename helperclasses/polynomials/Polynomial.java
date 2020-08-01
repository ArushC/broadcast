package helperclasses.polynomials;

//This polynomial class supports FFT multiplication
public class Polynomial {
	
	private Complex[] values;
	
	public Polynomial(Complex[] values) {
		this.values = values.clone();
	}
	
	public Complex[] getValues() {
		return this.values;
	}
	
	public Polynomial multiply(Polynomial other) {
		return new Polynomial(FFT.convolve(this.values, other.values));
	}
	
	public int[] getCoefficients() {
		int[] res = new int[values.length];
		for (int i = 0; i < values.length; i++) {
			res[i] = (int) Math.round(values[i].re());
		}
		return res;
	}
	
}
