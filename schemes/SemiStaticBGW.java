package schemes;

import com.herumi.mcl.*;
import java.io.*;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import helperclasses.Tools;

//The scheme: https://eprint.iacr.org/2008/268.pdf (3.1, page 8)
public class SemiStaticBGW {

	private static G2 gg;
	private static G1 g;
	private static int n;
	private static GT K;
	private static Fr alpha, t;
	//input: n = # of users
	//output: Object[] containing the public key PK and private key SK
	public static Object[] setup(int n) {
		
		//instead of creating a GroupGen function, the groups are generated here
		
		SemiStaticBGW.n = n; //save n so it can be used for later functions
		
		g = new G1();
		Mcl.hashAndMapToG1(g, "abc".getBytes());
		
		gg = new G2();
		Mcl.hashAndMapToG2(gg, "def".getBytes());
		
		//generate random alpha and t, precompute K
		alpha = new Fr();
		alpha.setByCSPRNG();
		
		t = new Fr();
		t.setByCSPRNG();
		
		precompute();
		
		//add g to public key
		ArrayList<Object> PK = new ArrayList<Object>();
		PK.add(g);
		
		//add pairing to public key
		GT e = new GT();
		Mcl.pairing(e, g, gg);
		Mcl.pow(e, e, alpha);
		PK.add(e);
		
		//generate random group elements h1, h2, ..., hn
		
		for (int i = 0; i < n; i++) {
			G1 h = new G1();
			Mcl.hashAndMapToG1(h, Tools.generateRandomBytes());
			PK.add(h);
		}
		
		//initialize SK
		G1 SK = new G1();
		Mcl.mul(SK, g, alpha);
		
		//add PK and SK to result (SK also contains PK)
		Object[] result = new Object[2];
		result[0] = PK;
		Object[] fullSK = new Object[2];
		fullSK[0] = SK;
		fullSK[1] = PK;
		result[1] = fullSK;
		return result;	
	}
	
	//input:  secret key SK and int i
	//output: ArrayList<Object> d, which contains all the individual secret keys
	public static ArrayList<Object> keyGen(int i, Object[] SK) {
		
		ArrayList<Object> di = new ArrayList<Object>();
		//extract from SK
		
		G1 MSK = new G1((G1) SK[0]);
		
		ArrayList<Object> PK = (ArrayList<Object>) SK[1];
		
		//generate random exponent ri
		Fr ri = new Fr();
		ri.setByCSPRNG();
		
		//calculate di0 and add to result
		G2 di0 = new G2();
		Mcl.mul(di0, gg, ri);
		Mcl.mul(di0, di0, new Fr(-1));
		di.add(di0);
		
		//calculate the other di's

		for (int j = 1; j <= n; j++) {
			if (j == i) {
				G1 dii = new G1();
				G1 hi = new G1((G1) PK.get(j + 1));
				Mcl.mul(hi, hi, ri);
				Mcl.add(dii, MSK, hi);
				di.add(dii);
			}
			else {
				G1 hj = new G1();
				Mcl.mul(hj, (G1) PK.get(j + 1), ri);
				di.add(hj);
			}
		
		}
		
		return di;
		
	}
	
	public static Object[] enc(ArrayList<Integer> S, ArrayList<Object> PK) {
		
		//calculate C1
		G2 C1 = new G2();
		Mcl.mul(C1, gg, t);
		
		//calculate C2
		G1 product = new G1((G1) PK.get(S.get(0) + 1));
		for (int k = 1; k < S.size(); k++) {
			int index = S.get(k) + 1;
			Mcl.add(product, product, (G1) PK.get(index));
		}
		G1 C2 = new G1();
		Mcl.mul(C2, product, t);
		
		//output (Hdr = (C1, C2), K)
		Object[] Hdr = new Object[2];
		Hdr[0] = C1;
		Hdr[1] = C2;
		Object[] result = new Object[2];
		result[0] = Hdr;
		result[1] = K;
		return result;
		
	}
	
	public static GT decrypt(ArrayList<Integer> S, int i, ArrayList<Object> di, Object[] Hdr, ArrayList<Object> PK) {
		
		//1. calculate the messy product in the first pairing
		G1 product = new G1((G1) di.get(i));
		
		for (int j: S) {
			if (j == i) 
				continue;
			G1 pork = new G1((G1) di.get(j));
			Mcl.add(product, product, pork);
		}
		
		//2. extract C1 and C2 from Hdr
		G2 C1 = (G2) Hdr[0];
		G1 C2 = (G1) Hdr[1];
		
		//calculate first pairing
		GT e1 = new GT();
		Mcl.pairing(e1, product, C1);

		//calculate second pairing
		GT e2 = new GT();
		Mcl.pairing(e2, C2, (G2) di.get(0));
		
		//result = e1 * e2
		GT result = new GT();
		Mcl.mul(result, e1, e2);
		return result;
		
	}
	
	
	//precomputes K = e(g, g)^(alpha * t)
	private static void precompute() {
		K = new GT();
		Mcl.pairing(K, g, gg);
		Fr exp = new Fr();
		Mcl.mul(exp, alpha, t);
		Mcl.pow(K, K, exp);
	}
	
	
	public static void main(String[] args) {
		
		System.load("/Users/arushchhatrapati/Documents/mcl/lib/libmcljava.dylib");
		Mcl.SystemInit(Mcl.BN254);
		Object[] setup = setup(100);
		ArrayList<Object> PK = (ArrayList<Object>) setup[0];
		Object[] SK = (Object[]) setup[1];
		int i = 75;
		ArrayList<Object> di = keyGen(i, SK);
		ArrayList<Integer> S = new ArrayList<Integer>();
		S.addAll(Arrays.asList(13, 75, 42, 14, 53, 12, 1));
		Object[] C = enc(S, PK);
		Object[] Hdr = (Object[]) C[0];
		GT K1 = decrypt(S, i, di, Hdr, PK);
		System.out.println("K = " + K);
		System.out.println("K1 = " + K1);
		
		
	}

	
	
}
