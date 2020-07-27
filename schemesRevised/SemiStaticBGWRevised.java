package schemesRevised;
import com.herumi.mcl.*;
import java.io.*;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ThreadLocalRandom;
import helperclasses.Tools;

//The scheme: https://eprint.iacr.org/2008/268.pdf (3.1, page 8)
//Changes made: optimized selection of G1 and G2, added both to the public key
//KeyGen takes the public key as a parameter (a slight modification to the described scheme)
//K is precomputed in the setup phase so that encryption does not require pairings
//Secret key contains g^(alpha) and gg^(alpha) 

public class SemiStaticBGWRevised {

	private static int n;
	private static GT K;
	
	//input: n = # of users
	//output: Object[] containing the public key PK and private key SK
	public static Object[] setup(int n) {
		
		//instead of creating a GroupGen function, the groups are generated here
		
		SemiStaticBGWRevised.n = n; //save n so it can be used for later functions
		
		G1 g = new G1();
		Mcl.hashAndMapToG1(g, "abc".getBytes());
		
		G2 gg = new G2();
		Mcl.hashAndMapToG2(gg, "def".getBytes());
		
		//generate random alpha and t, precompute K
		Fr alpha = new Fr();
		alpha.setByCSPRNG();
		
		precompute(g, gg, alpha);
		
		//add g to public key
		ArrayList<Object> PK = new ArrayList<Object>();
		PK.add(g);
		PK.add(gg);
		
		//add pairing to public key
		GT e = new GT();
		Mcl.pairing(e, g, gg);
		Mcl.pow(e, e, alpha);
		PK.add(e);
		
		//generate random group elements h1, h2, ..., hn
		byte[] bts = Tools.generateRandomBytes(n + 2);
		for (int i = 0; i < n; i++) {
			G1 h = new G1();
			byte[] randomBytes = Arrays.copyOfRange(bts, i, i + 3);
			Mcl.hashAndMapToG1(h, randomBytes);
			PK.add(h);
		}
		
		//initialize SK
		G1 sk1 = new G1();
		Mcl.mul(sk1, g, alpha);
		G2 sk2 = new G2();
		Mcl.mul(sk2, gg, alpha);
		Object[] SK = {sk1, sk2};
		
		//add PK and SK to result
		Object[] result = {PK, SK};
		return result;	
	}
	
	//input:  secret key SK and int i
	//output: ArrayList<Object> d, which contains all the individual secret keys
	public static ArrayList<Object> keyGen(int i, ArrayList<Object> PK, Object[] SK) {
		
		ArrayList<Object> di = new ArrayList<Object>();
		
		//generate random exponent ri
		Fr ri = new Fr();
		ri.setByCSPRNG();
		
		//calculate di0 and add to result
		G2 di0 = new G2();
		Mcl.mul(di0, (G2) PK.get(1), ri);
		Mcl.mul(di0, di0, new Fr(-1));
		di.add(di0);
		
		//calculate the other di's

		for (int j = 1; j <= n; j++) {
			if (j == i) {
				G1 dii = new G1();
				G1 hi = new G1((G1) PK.get(j + 2));
				Mcl.mul(hi, hi, ri);
				Mcl.add(dii, (G1) SK[0], hi);
				di.add(dii);
			}
			else {
				G1 hj = new G1();
				Mcl.mul(hj, (G1) PK.get(j + 2), ri);
				di.add(hj);
			}
		
		}
		
		return di;
	}
	
	public static Object[] enc(ArrayList<Integer> S, ArrayList<Object> PK) {
		
		Fr t = new Fr();
		t.setByCSPRNG();
		Mcl.pow(K, K, t);
		
		//calculate C1
		G2 C1 = new G2();
		Mcl.mul(C1, (G2) PK.get(1), t);
		
		//calculate C2
		G1 product = new G1((G1) PK.get(S.get(0) + 2));
		for (int k = 1; k < S.size(); k++) {
			int index = S.get(k) + 2;
			Mcl.add(product, product, (G1) PK.get(index));
		}
		
		G1 C2 = new G1();
		Mcl.mul(C2, product, t);
		
		//output (Hdr = (C1, C2), K)
		Object[] Hdr = {C1, C2};
		Object[] result = {Hdr, K};
		return result;
		
	}
	
	public static GT decrypt(ArrayList<Integer> S, int i, ArrayList<Object> di, Object[] Hdr, ArrayList<Object> PK) {
		
		//1. calculate the messy product in the first pairing
		G1 product = new G1((G1) di.get(i));
		
		for (int j: S) {
			if (j == i) 
				continue;
			G1 addend = new G1((G1) di.get(j));
			Mcl.add(product, product, addend);
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
	
	
	//precomputes K = e(g, g)^alpha
	private static void precompute(G1 g, G2 gg, Fr alpha) {
		K = new GT();
		Mcl.pairing(K, g, gg);
		Mcl.pow(K, K, alpha);
	}
	
		
private static long[] printRuntimes(int N, int subsetSize) {
		
		long startSetup = System.nanoTime();
		Object[] setup = setup(N);
		long elapsedSetup = System.nanoTime() - startSetup;
		double secondsSetup = ((double) elapsedSetup) / 1E9;

		ArrayList<Integer> S = new ArrayList<Integer>();
		
		//randomly generate numbers to put in the subset S (NO REPEATS)
		ArrayList<Integer> randomNums = new ArrayList<Integer>();
			for (int i = 0; i < N; i++) { 
				randomNums.add(i + 1);
			}
			for (int i = 0; i < subsetSize; i++) {
				int randomIndex = ThreadLocalRandom.current().nextInt(1, randomNums.size());
				int randomID = randomNums.get(randomIndex);
				S.add(randomID);
				randomNums.remove(randomIndex);
			}
			
		//Get random user ID u to test the decryption
		int i = S.get(ThreadLocalRandom.current().nextInt(0, S.size()));
		ArrayList<Object> PK = (ArrayList<Object>) setup[0];
		Object[] SK = (Object[]) setup[1];
		
		long startKeyGen = System.nanoTime();
		ArrayList<Object> di = keyGen(i, PK, SK);
		long elapsedKeyGen = System.nanoTime() - startKeyGen;
		double secondsKeyGen = ((double) elapsedKeyGen)/1E9;
		
		long startEnc = System.nanoTime();
		Object[] C = enc(S, PK);
		long elapsedEnc = System.nanoTime() - startEnc;
		double secondsEnc = ((double) elapsedEnc)/1E9;
		
		Object[] Hdr = (Object[]) C[0];
		long startDecrypt = System.nanoTime();
		GT K1 = decrypt(S, i, di, Hdr, PK);
		long elapsedDecrypt = System.nanoTime() - startDecrypt;
		double secondsDecrypt = ((double) elapsedDecrypt)/1E9;
		
		
		String result = (K.equals(K1)) ? "SUCCESSFUL DECRYPTION" : "FAILED DECRYPTION";
		System.out.println(result + ": " + "N = " + N + ", subset size = " + subsetSize);
		System.out.println(); //padding
		System.out.println("setup took " + secondsSetup + " seconds");
		System.out.println("encryption took " + secondsEnc + " seconds");
		System.out.println("key generation took " + secondsKeyGen + " seconds");
		System.out.println("decryption took " + secondsDecrypt + " seconds");
		System.out.println(); //more padding
		
		long[] elapsedTimes = new long[4];
		elapsedTimes[0] = elapsedSetup;
		elapsedTimes[1] = elapsedEnc;
		elapsedTimes[2] = elapsedKeyGen;
		elapsedTimes[3] = elapsedDecrypt;
		
		return elapsedTimes;
		
	}

	
	public static void testRuntimes(int percent) {
		for (int N = 10; N <= 1000000; N *= 10) {
			int subsetSize = (int) (0.01 * percent * N);
			printRuntimes(N, subsetSize);
		}
	}
	
	public static void main(String[] args) {
		File lib = new File("../../lib/libmcljava.dylib");
		System.load(lib.getAbsolutePath());
		Mcl.SystemInit(Mcl.BN254);
		testRuntimes(10);
	}

}
