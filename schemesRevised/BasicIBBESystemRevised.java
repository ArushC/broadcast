package schemesRevised;
import helperclasses.polynomials.LagrangeInterpolationZp;
import helperclasses.structures.VectorND;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import com.herumi.mcl.*;

//The scheme: https://eprint.iacr.org/2019/038.pdf, 3.1 (page 11)
//No changes were made to the described scheme (it is written in the Type-3 Pairing setting)
public class BasicIBBESystemRevised {

	private static int m;
	
	//input: security parameter lambda, maximal subset size m
	//output: Object[] containing public key PP and master secret key MSK
	public static Object[] setup(int lambda, int m) {
		
		Mcl.SystemInit(lambda);
		BasicIBBESystemRevised.m = m; //save so it can be used for other functions
		
		//1. generate u1 and u2
		VectorND u1 = new VectorND(m + 1); //automatically set randomly by the class
		VectorND u2 = new VectorND(m + 1);
		
		//2. generate random exponents
		Fr alpha1 = new Fr();
		alpha1.setByCSPRNG();
		
		Fr alpha2 = new Fr();
		alpha2.setByCSPRNG();
		
		Fr w1 = new Fr();
		w1.setByCSPRNG();
		
		Fr w2 = new Fr();
		w2.setByCSPRNG();
		
		Fr b = new Fr();
		b.setByCSPRNG();
		
		Fr beta = new Fr();
		beta.setByCSPRNG();
		
		VectorND u = u2.multiply(b); //u = u1 + b(u2)
		u = u.add(u1);
		
		Fr w = new Fr(w2); //w = w1 + b(w2)
		Mcl.mul(w, w, b);
		Mcl.add(w, w1, w);
		
		Fr alpha = new Fr(alpha2); //alpha = alpha1 + b(alpha2)
		Mcl.mul(alpha, alpha, b);
		Mcl.add(alpha, alpha1, alpha);
		
		//generate g1 and g2
		G1 g1 = new G1();
		Mcl.hashAndMapToG1(g1, "abc".getBytes());
		
		G2 g2 = new G2();
		Mcl.hashAndMapToG2(g2, "def".getBytes());
		
		//calculate U1, W1, gT
		ArrayList<G1> U1 = VectorND.exponentiate(g1, u);
		
		G1 W1 = new G1();
		Mcl.mul(W1, g1, w);
		
		GT gT = new GT();
		Mcl.pairing(gT, g1, g2);
		Mcl.pow(gT, gT, alpha);
		
		//Add to public key PP
		ArrayList<Object> PP = new ArrayList<Object>();
		
		//e1 = g1, skip to e2
		G1 e2 = new G1();
		Mcl.mul(e2, g1, b);
		
		//e3, e4, e5, e6 are U1, W1, gT, g2, skip to e7
		ArrayList<G2> e7 = VectorND.exponentiate(g2, u1);
		ArrayList<G2> e8 = VectorND.exponentiate(g2, u2);
		
		G2 e9 = new G2();
		Mcl.mul(e9, g2, w1);
		
		G2 e10 = new G2();
		Mcl.mul(e10, g2, w2);
		
		G1 e11 = new G1();
		Mcl.mul(e11, g1, beta);
		Mcl.mul(e11, e11, alpha);
		
		G2 e12 = new G2();
		Mcl.mul(e12, g2, beta);
		Mcl.mul(e12, e12, alpha1);
		
		G2 e13 = new G2();
		Mcl.mul(e13, g2, beta);
		Mcl.mul(e13, e13, alpha2);
		
		G2 e14 = new G2();
		Fr exp14 = new Fr();
		Mcl.div(exp14, new Fr(1), beta);
		Mcl.mul(e14, g2, exp14);
		
		PP.addAll(Arrays.asList(g1, e2, U1, W1, gT, g2, e7, e8, e9, e10, e11, e12, e13, e14));
		
		//Add to MSK		
		G2 e1S = new G2();
		Mcl.mul(e1S, g2, alpha1);
		
		G2 e2S = new G2();
		Mcl.mul(e2S, g2, alpha2);
		
		Object[] MSK = {e1S, e2S};
		
		//return public and private key
		Object[] result = {PP, MSK};
		return result;
	}
	
	
	//input: public key PP, master secret key MSK, user ID
	public static Object[] keyGen(ArrayList<Object> PP, Object[] MSK, Fr ID) {
		
		//choose random r and tags
		Fr r = new Fr();
		r.setByCSPRNG();
		VectorND kTags = new VectorND(m);
		
		
		//calculate K1, K2, K3, ...
		G2 K1 = new G2();
		G2 K1Helper = new G2();
		Mcl.mul(K1Helper, (G2) PP.get(8), r);
		Mcl.add(K1, (G2) MSK[0], K1Helper);
		
		G2 K2 = new G2();
		G2 K2Helper = new G2();
		Mcl.mul(K2Helper, (G2) PP.get(9), r);
		Mcl.add(K2, (G2) MSK[1], K2Helper);
		
		G2 K3 = new G2();
		Mcl.mul(K3, (G2) PP.get(5), r);
		
		G2[] K4 = new G2[m];
		G2[] K5 = new G2[m];
		
		Fr expDen = new Fr(ID);
		
		for (int i = 1; i <= m; i++) {
			G2 K4i = new G2();
			G2 K5i = new G2();
			G2 K4iDen = new G2(((ArrayList<G2>) PP.get(6)).get(0));
			G2 K5iDen = new G2(((ArrayList<G2>) PP.get(7)).get(0));
			Mcl.mul(K4iDen, K4iDen, expDen);
			Mcl.mul(K5iDen, K5iDen, expDen);
			G2 K4NumHelper = new G2();
			G2 K5NumHelper = new G2();
			Mcl.mul(K4NumHelper, (G2) PP.get(8), kTags.getCoords().get(i - 1));
			Mcl.mul(K5NumHelper, (G2) PP.get(9), kTags.getCoords().get(i - 1));
			Mcl.add(K4i, K4NumHelper, ((ArrayList<G2>) PP.get(6)).get(i));
			Mcl.add(K5i, K5NumHelper, ((ArrayList<G2>) PP.get(7)).get(i));
			Mcl.sub(K4i, K4i, K4iDen);
			Mcl.sub(K5i, K5i, K5iDen);
			Mcl.mul(K4i, K4i, r);
			Mcl.mul(K5i, K5i, r);
			K4[i - 1] = K4i;
			K5[i -1] = K5i;
			Mcl.mul(expDen, expDen, ID);
		}
		
		Object[] SKID = {K1, K2, K3, K4, K5, kTags};
		return SKID;
	}
	
	
	public static Object[] encap(ArrayList<Object> PP, ArrayList<Fr> S) {
		
		//compute y
		ArrayList<Fr> y = computeY(S);
		
		//pick random s and ctag
		Fr s = new Fr();
		s.setByCSPRNG();
		
		Fr ctag = new Fr();
		ctag.setByCSPRNG();
		
		//compute C1, C2, C3
		
		G1 C1 = new G1();
		Mcl.mul(C1, (G1) PP.get(0), s);
		
		G1 C2 = new G1();
		Mcl.mul(C2, (G1) PP.get(1), s);
		
		G1 C3 = new G1();
		Mcl.mul(C3, (G1) PP.get(3), ctag);
		ArrayList<G1> U = (ArrayList<G1>) PP.get(2);
		for (int i = 0; i <= S.size(); i++) {
			G1 prod = new G1();
			Mcl.mul(prod, U.get(i), y.get(i));
			Mcl.add(C3, C3, prod);
		}
		
		Mcl.mul(C3, C3, s);
		
		GT K = new GT();
		Mcl.pow(K, (GT) PP.get(4), s);
		
		Object[] Hdr = {C1, C2, C3, ctag};
		Object[] result = {Hdr, K};
		
		return result;	
	}
	
	public static GT decap(ArrayList<Object> PP, ArrayList<Fr> S, Object[] Hdr, Object[] SKID) {
		
		//extract
		ArrayList<Fr> y = computeY(S);
		VectorND kTags = (VectorND) SKID[5];
		
		
		//compute ktag
		Fr ktag = new Fr(0);
		for (int i = 1; i <= m; i++) {
			Fr product = new Fr();
			Mcl.mul(product, kTags.getCoords().get(i - 1), y.get(i));
			Mcl.add(ktag, ktag, product);
		}
		
		//compute K4 and K5 products
		G2[] K4Values = (G2[]) SKID[3];
		G2[] K5Values = (G2[]) SKID[4];
		G2 K4Product = new G2();
		Mcl.mul(K4Product, K4Values[0], y.get(1)); //initialize
		G2 K5Product = new G2();
		Mcl.mul(K5Product, K5Values[0], y.get(1)); //initialize
		
		for (int i = 2; i <= m; i++) {
			G2 inner4 = new G2();
			Mcl.mul(inner4, K4Values[i - 1], y.get(i));
			G2 inner5 = new G2();
			Mcl.mul(inner5, K5Values[i - 1], y.get(i));
			Mcl.add(K4Product, inner4, K4Product);
			Mcl.add(K5Product, inner5, K5Product);
		}
		
		//Finally, compute the pairings
		//First extract everything that is needed
		G1 C1 = (G1) Hdr[0];
		G1 C2 = (G1) Hdr[1];
		G1 C3 = (G1) Hdr[2];
		Fr ctag = (Fr) Hdr[3];
		G2 K1 = (G2) SKID[0];
		G2 K2 = (G2) SKID[1];
		G2 K3 = (G2) SKID[2];
		
		GT e1 = new GT();
		Mcl.pairing(e1, C1, K4Product);
		
		GT e2 = new GT();
		Mcl.pairing(e2, C2, K5Product);
		
		GT e3 = new GT();
		Mcl.pairing(e3, C3, K3);
		Mcl.pow(e3, e3, new Fr(-1));
		
		GT A = new GT();
		Mcl.mul(A, e1, e2);
		Mcl.mul(A, A, e3);
		
		Fr exp = new Fr();
		Fr expDen = new Fr();
		Mcl.sub(expDen, ktag, ctag);
		Mcl.div(exp, new Fr(1), expDen);
		Mcl.pow(A, A, exp);
		
		GT e4 = new GT();
		Mcl.pairing(e4, C1, K1);
		
		GT e5 = new GT();
		Mcl.pairing(e5, C2, K2);
		
		GT K = new GT();
		Mcl.mul(K, e4, e5);
		Mcl.pow(A, A, new Fr(-1));
		Mcl.mul(K, K, A);
		
		return K;
	}
	
	//computes the vector y from the subset S
	private static ArrayList<Fr> computeY(ArrayList<Fr> S) {
		
		Fr[] PZ = {new Fr(1)};
		for (int i = 0; i < S.size(); i++) {
			Fr IDJ = new Fr(S.get(i));
			Mcl.mul(IDJ, IDJ, new Fr(-1));
			Fr[] next = {new Fr(1), IDJ};
			PZ = LagrangeInterpolationZp.multiply(PZ, next);
		}
		
		ArrayList<Fr> y = new ArrayList<Fr>(Arrays.asList(PZ));
		Collections.reverse(y);
		
		for (int i = y.size(); i <= m; i++)
			y.add(new Fr(0)); //pad with zeroes
		
		return y;
	}
	
	public static long[] printRuntimes(int N, int subsetSize, int lambda) {
		
		Mcl.SystemInit(lambda);
		long startSetup = System.nanoTime();
		Object[] setup = setup(Mcl.BN254, N);
		long elapsedSetup = System.nanoTime() - startSetup;
		double secondsSetup = ((double) elapsedSetup) / 1E9;
		Fr ID = new Fr();
		ID.setByCSPRNG();
		ArrayList<Object> PP = (ArrayList<Object>) setup[0];
		Object[] MSK = (Object[]) setup[1];
		
		
		//generate subset S
		ArrayList<Fr> S = new ArrayList<Fr>();
		for (int i = 0; i < subsetSize - 1; i++) {
			Fr d = new Fr();
			d.setByCSPRNG();
			S.add(d);
		}
		//random user ID to test decryption
		S.add(ID);
		
		long startKeygen = System.nanoTime();
		Object[] SKID = keyGen(PP, MSK, ID);
		long elapsedKeygen = System.nanoTime() - startKeygen;
		double secondsKeygen = ((double) elapsedKeygen) / 1E9;
		
		long startEnc = System.nanoTime();
		Object[] C = encap(PP, S);
		long elapsedEnc = System.nanoTime() - startEnc;
		double secondsEnc = ((double) elapsedEnc) / 1E9;
		
		Object[] Hdr = (Object[]) C[0];
		GT K = (GT) C[1];
		long startDec = System.nanoTime();
		GT K1 = decap(PP, S, Hdr, SKID);
		long elapsedDec = System.nanoTime() - startDec;
		double secondsDec = ((double) elapsedDec) / 1E9;
		
		String success = (K1.equals(K)) ? "SUCCESSFUL DECRYPTION" : "FAILED DECRYPTION";
		System.out.println(success + ": " + "n = " + N + ", subset size = " + subsetSize);
		System.out.println(); //padding
		System.out.println("setup took " + secondsSetup + " seconds");
		System.out.println("encrypt took " + secondsEnc + " seconds");
		System.out.println("keygen took " + secondsKeygen + " seconds");
		System.out.println("decrypt took " + secondsDec + " seconds");
		System.out.println(); //more padding
		
		long[] elapsedTimes = new long[4];
		elapsedTimes[0] = elapsedSetup;
		elapsedTimes[1] = elapsedEnc;
		elapsedTimes[2] = elapsedKeygen;
		elapsedTimes[3] = elapsedDec;
		
		return elapsedTimes;
		
	}
	
	//Test the runtimes of the algorithms
	public static void testRuntimes(int lambda, int percent) {
		for (int N = 100; N <= 1000000; N *= 10) {
			int subsetSize = (int) (0.01 * percent * N);
			printRuntimes(N, subsetSize, lambda);
		}
	}
	
	
	public static void main(String[] args) {
		
		File lib = new File("../../lib/libmcljava.dylib");
		System.load(lib.getAbsolutePath());
		printRuntimes(100000, 100, Mcl.BN254);
		//testRuntimes(Mcl.BN254, 10);
	}
	
}

