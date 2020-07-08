package schemes;
import java.util.ArrayList;
import java.util.Arrays;

import com.herumi.mcl.*;
//El Gamal Encryption With Shared Parameters
public class BaselineScheme {

	private static int N;
	
	//output public key and private key
	public static Object[] setup(int lambda, int N) {
		
		Mcl.SystemInit(lambda);
		BaselineScheme.N = N;
		
		//1. Random generator g
		G1 g = new G1();
		Mcl.hashAndMapToG1(g, "abc".getBytes());
		
		//instantiate public key, secret key
		ArrayList<G1> PK = new ArrayList<G1>();
		PK.add(g);
		ArrayList<Fr> SK = new ArrayList<Fr>();
		
		//2. random exponents x_1, x_2, ..., x_n
		//also compute h1, h2, ..., hn (hi = g^(xi))
		
		for (int i = 0; i < N; i++) {
			Fr xi = new Fr();
			xi.setByCSPRNG();
			SK.add(xi);
			G1 hi = new G1();
			Mcl.mul(hi, g, xi);
			PK.add(hi);
		}
		
		Object[] result = {PK, SK};
		return result;
	}
	
	//input: public key PK, symmetric key M, set S
	//encrypt to a set S, outputs the ciphertext CT
	public static G1[] enc(ArrayList<G1> PK, ArrayList<Integer> S, G1 M) {
		
		G1[] CT = new G1[N + 1];
		
		//random y in Z_q
		Fr y = new Fr();
		y.setByCSPRNG();
		
		//add g^y to CT
		G1 e1 = new G1();
		Mcl.mul(e1, PK.get(0), y);
		CT[0] = e1;
		
		//compute zi = (hi^y) * m for each i in S
		for (int i: S) {
			G1 hi = PK.get(i);
			G1 zi = new G1();
			Mcl.mul(zi, hi, y);
			Mcl.add(zi, zi, M);
			CT[i] = zi;
		}
		
		return CT;	
	}
	
	//input:  user i, the secret key for user i (xi), the ciphertext CT
	//output: the symmetric key M
	public static G1 dec(G1[] CT, Fr xi, int i) {
		G1 zi = CT[i];
		G1 c = CT[0];
		G1 res = new G1();
		Mcl.mul(res, c, xi);
		Mcl.sub(res, zi, res);	
		return res;	
	}
	
	public static void main(String[] args) {
		System.load("/Users/arushchhatrapati/Documents/mcl/lib/libmcljava.dylib");
		Object[] setup = setup(Mcl.BN254, 100);
		ArrayList<G1> PK = (ArrayList<G1>) setup[0];
		ArrayList<Fr> SK = (ArrayList<Fr>) setup[1];
		ArrayList<Integer> S = new ArrayList<Integer>();
		S.addAll(Arrays.asList(1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 99));
		G1 M = new G1();
		Mcl.hashAndMapToG1(M, "abc".getBytes());
		System.out.println("M = " + M);
		int i = 99;
		G1[] CT = enc(PK, S, M);
		Fr xi = SK.get(i - 1);
		G1 M1 = dec(CT, xi, i);
		System.out.println("M1 = " + M1);
	}
	
	
	
}
