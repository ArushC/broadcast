package helperclasses.structures;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ThreadLocalRandom;
import com.herumi.mcl.*;

public class VectorND {
	
	private ArrayList<Fr> coordinates = new ArrayList<Fr>();
	private int n;
	public final static int VECTOR_RANDOM = 0;
	public final static int VECTOR_ZERO = 1;
	public final static int VECTOR_ZERO_ONE = 2;
	
	public VectorND(int n) { //default random vector
		this(VECTOR_RANDOM, n);
	}
	
	public VectorND(VectorND vector, Fr... otherElements) {
		this.coordinates.addAll(vector.coordinates);
		this.coordinates.addAll(Arrays.asList(otherElements));
	}
	
	public VectorND(int type, int n) { //initialize vector with specified type
		
		if (type == VECTOR_RANDOM) {
			init(n);	
			this.n = n;
		}
		else if (type == VECTOR_ZERO) {
			for (int i = 0; i < n; i++)
				coordinates.add(new Fr(0));
		}
		else if (type == VECTOR_ZERO_ONE) { //random vector consisting of only zeroes and ones
			for (int i = 0; i < n; i++) {
				int p = ThreadLocalRandom.current().nextInt(0, 2);
				coordinates.add(p == 0 ? new Fr(0) : new Fr(1));
			}
		}
		
	}
	
	public VectorND(Fr... coords) { //accepts arbitrary number of Fr inputs
		for (Fr coord: coords)
			coordinates.add(coord);
	}
	
	public VectorND(ArrayList<Fr> coords) { //accepts the list of coordinates directly
		for (Fr coord: coords)
			coordinates.add(coord);
	}
	
	private void init(int n) {
		
		for (int i = 0; i < n; i++) {
			Fr c = new Fr();
			c.setByCSPRNG();
			coordinates.add(c);
		}
		this.n = n;
	}
	
	public void setByCSPRNG() {
		init(n); //to set randomly, just reinitialize
	}
	
	public VectorND multiply(Fr scalar) {
		
		ArrayList<Fr> resCoords = new ArrayList<Fr>();
		for (Fr coord: coordinates) {
			Fr product = new Fr();
			Mcl.mul(product, coord, scalar);
			resCoords.add(product);
		}
		
		return new VectorND(resCoords);
	}
	
	public VectorND add(VectorND v) {
		
		ArrayList<Fr> resCoords = new ArrayList<Fr>();
		for (int i = 0; i < this.coordinates.size(); i++) {
			Fr sum = new Fr();
			Mcl.add(sum, v.getCoords().get(i), this.coordinates.get(i));
			resCoords.add(sum);
		}
		
		return new VectorND(resCoords);
		
	}
	
	public Fr dotProduct(VectorND v) {
		
		Fr sum = new Fr(0);
		for (int i = 0; i < this.coordinates.size(); i++) {
			Fr product = new Fr();
			Mcl.mul(product, v.getCoords().get(i), this.getCoords().get(i));
			Mcl.add(sum, sum, product);
		}
		
		return sum;
	}
	
	public static ArrayList<G1> exponentiate(G1 g, VectorND v) {
		
		ArrayList<G1> res = new ArrayList<G1>();
		ArrayList<Fr> coords = v.getCoords();
		for (Fr coord: coords) {
			G1 exponentiated = new G1();
			Mcl.mul(exponentiated, g, coord);
			res.add(exponentiated);
		}
		
		return res;
		
	}
	
	//Overloaded method -- one for G1, one for G2
	public static ArrayList<G2> exponentiate(G2 g, VectorND v) {
		
		ArrayList<G2> res = new ArrayList<G2>();
		ArrayList<Fr> coords = v.getCoords();
		for (Fr coord: coords) {
			G2 exponentiated = new G2();
			Mcl.mul(exponentiated, g, coord);
			res.add(exponentiated);
		}
		
		return res;
		
	}
	
	//precondition: v1.size() = v2.size()
	public static GT vectorPairing(ArrayList<G1> v1, ArrayList<G2> v2) {
		GT result = new GT();
		Mcl.pairing(result, v1.get(0), v2.get(0));
		for (int i = 1; i < v1.size(); i++) {
			GT e = new GT();
			Mcl.pairing(e, v1.get(i), v2.get(i));
			Mcl.mul(result, result, e);
		}
		return result;
	}
	
	//getter method
	public ArrayList<Fr> getCoords() {
		return this.coordinates;
	}
	
	@Override
	public String toString() {
		String res = new String();
		for (int i = 0; i < coordinates.size(); i++) {
			res += coordinates.get(i).toString();
			if (!(i == coordinates.size() - 1))
				res +=  ", ";
		}
		
		return "(" + res + ")";
	}
	
}
