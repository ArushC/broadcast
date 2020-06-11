import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.io.*;

import com.herumi.mcl.*;

//The scheme: https://eprint.iacr.org/2009/385.pdf (bottom of page 6)
public class WatersDualEncryption {
	
	private static G2 gg;
	private static Fr b; 
	
	//input: security parameter lambda (Mcl.BN254 or Mcl.BLS_381)
	//output: (public key PK, secret key MSK)
	public static Object[] setup(int lambda) { 
		
		Mcl.SystemInit(lambda); // curveType = Mcl.BN254 or Mcl.BLS12_381
		
		//Step 1: choose random generators
		
		//This is the g in G1 that is used in the private key
		G1 g = new G1();
		Mcl.hashAndMapToG1(g, "abc".getBytes());
		
		//This is the "real" g that is used for calculations in the public key
		gg = new G2();
		Mcl.hashAndMapToG2(gg, "def".getBytes());
		
		G1 v = new G1();
		Mcl.hashAndMapToG1(v, "ghi".getBytes());
		
		G1 v1 = new G1();
		Mcl.hashAndMapToG1(v1, "jkl".getBytes());

		G1 v2 = new G1();
		Mcl.hashAndMapToG1(v2, "mno".getBytes());
		
		G1 w = new G1();
		Mcl.hashAndMapToG1(w, "pqr".getBytes());
		
		G1 u = new G1();
		Mcl.hashAndMapToG1(u, "stu".getBytes());
		
		G1 h = new G1();
		Mcl.hashAndMapToG1(h, "vwx".getBytes());
		
		//Step 2: choose random exponents in Z_p
		
		Fr a1 = new Fr();
		a1.setByCSPRNG();
		
		Fr a2 = new Fr();
		a2.setByCSPRNG();
		
		b = new Fr();
		b.setByCSPRNG();
		
		Fr alpha = new Fr();
		alpha.setByCSPRNG();
		
		//Step 3: compute tau1 and tau2
		
		G1 tau1 = new G1(); 
		Mcl.mul(tau1, v1, a1);
		Mcl.add(tau1, v, tau1);
		
		G1 tau2 = new G1();
		Mcl.mul(tau2, v2, a2);
		Mcl.add(tau2, v, tau2);
		
		//Step 4: instantiate PK and add elements to it: g^b, g^a1, g^a2, ... etc. see the paper
		
		ArrayList<Object> PK = new ArrayList<Object>();
		
		G2 e0 = new G2();
		Mcl.mul(e0, gg, b);
		
		G2 e1 = new G2();
		Mcl.mul(e1, gg, a1);
		
		G2 e2 = new G2();
		Mcl.mul(e2, gg, a2);
		
		G2 e3 = new G2();
		Fr exp3 = new Fr();
		Mcl.mul(exp3, b, a1);
		Mcl.mul(e3, gg, exp3);
		
		G2 e4 = new G2();
		Fr exp4 = new Fr();
		Mcl.mul(exp4, b, a2);
		Mcl.mul(e4, gg, exp4);
		
		//e5 and e6 are just tau1 and tau2, so skip to e7
		G1 e7 = new G1();
		Mcl.mul(e7, tau1, b);
		
		G1 e8 = new G1();
		Mcl.mul(e8, tau2, b);
		
		//e9, e10, and e11 are w, u, and h, respectively, so skip to e12
		GT e12 = new GT();
		Fr exp12 = new Fr();
		Mcl.mul(exp12, exp3, alpha);
		Mcl.pairing(e12, g, gg); //e(g, g)^(alpha * a1 * b) = e(g, gg)^(alpha * a1 * b)
		Mcl.pow(e12, e12, exp12);
		
		PK.addAll(Arrays.asList(e0, e1, e2, e3, e4, tau1, tau2, e7, e8, w, u, h, e12));
		
		//Step 5: instantiate MSK and add elements to it: g, g^alpha, g^(alpha * a1), ... etc. (see paper)
		
		ArrayList<Object> MSK = new ArrayList<Object>();
		
		//e0S is just g, so skip to e1S
		G1 e1S = new G1();
		Mcl.mul(e1S, g, alpha);
		
		G1 e2S = new G1();
		Fr exp2S = new Fr();
		Mcl.mul(exp2S, alpha, a1);
		Mcl.mul(e2S, g, exp2S);
		
		//e3S, e4S, and e5S are v, v1, and v2, respectively
		MSK.addAll(Arrays.asList(g, e1S, e2S, v, v1, v2));
		MSK.addAll(PK); //add public parameters to MSK

		//return the public and private key
		Object[] result = new Object[2];
		result[0] = PK;
		result[1] = MSK;
		return result;
		
	}
	
	//input: public key PK, secret I, message M
	public static ArrayList<Object> encrypt(ArrayList<Object> PK, Fr I, GT M) {
		
		//Step 1: choose random values in Z_p
		
		Fr s1 = new Fr();
		s1.setByCSPRNG();
		
		Fr s2 = new Fr();
		s2.setByCSPRNG();
		
		Fr t = new Fr();
		t.setByCSPRNG();
		
		Fr tagC = new Fr();
		tagC.setByCSPRNG();
		
		//Step 2: calculate C0, C1, C2, ...
		
		GT C0 = new GT();
		Mcl.pow(C0, (GT) PK.get(12), s2);
		Mcl.mul(C0, M, C0);
		
		G2 C1 = new G2();
		Fr exp1 = new Fr();
		Mcl.add(exp1, s1, s2);
		Mcl.mul(C1, (G2) PK.get(0), exp1);
		
		G2 C2 = new G2();
		Mcl.mul(C2, (G2) PK.get(3), s1);
		
		G2 C3 = new G2();
		Mcl.mul(C3, (G2) PK.get(1), s1);
		
		G2 C4 = new G2();
		Mcl.mul(C4, (G2) PK.get(4), s2);
		
		G2 C5 = new G2();
		Mcl.mul(C5, (G2) PK.get(2), s2);
		
		G1 C6 = new G1();
		G1 part6 = new G1();
		Mcl.mul(part6, (G1) PK.get(6), s2);
		Mcl.mul(C6, (G1) PK.get(5), s1);
		Mcl.add(C6, C6, part6);
		
		G1 C7 = new G1();
		G1 part17 = new G1();
		Mcl.mul(part17, (G1) PK.get(7), s1);
		G1 part27 = new G1();
		Mcl.mul(part27, (G1) PK.get(8), s2);
		G1 part37 = new G1();
		Fr exp37 = new Fr();
		Mcl.mul(exp37, t, new Fr(-1));
		Mcl.mul(part37, (G1) PK.get(9), exp37);
		Mcl.add(C7, part17, part27);
		Mcl.add(C7, C7, part37);
		
		G1 E1 = new G1();
		G1 part1E1 = new G1();
		Mcl.mul(part1E1, (G1) PK.get(10), I);
		G1 part2E1 = new G1();
		Mcl.mul(part2E1, (G1) PK.get(9), tagC);
		Mcl.add(E1, part1E1, part2E1);
		Mcl.add(E1, E1, (G1) PK.get(11));
		Mcl.mul(E1, E1, t);
		
		G2 E2 = new G2(); //Has to be of type G2 to compute the pairing in decryption
		Mcl.mul(E2, gg, t);
		
		//Add all to the ciphertext
		ArrayList<Object> CT = new ArrayList<Object>();
		CT.addAll(Arrays.asList(C0, C1, C2, C3, C4, C5, C6, C7, E1, E2, tagC));
		return CT;
	}
	
	
	public static ArrayList<Object> keyGen(ArrayList<Object> MSK, Fr I) {
		
		//Step 1: choose random values in Z_p
		
		Fr r1 = new Fr();
		r1.setByCSPRNG();
		
		Fr r2 = new Fr();
		r2.setByCSPRNG();
		
		Fr z1 = new Fr();
		z1.setByCSPRNG();
		
		Fr z2 = new Fr();
		z2.setByCSPRNG();
		
		Fr tagK = new Fr();
		tagK.setByCSPRNG();
		
		Fr r = new Fr();
		Mcl.add(r, r1, r2);
		
		
		//Step 2: compute D1, D2, ...
		
		G1 g = (G1) MSK.get(0);
		
		G1 D1 = new G1();
		G1 helper1 = new G1();
		Mcl.mul(helper1, (G1) MSK.get(3), r);
		Mcl.add(D1, (G1) MSK.get(2), helper1);
		
		G1 D2 = new G1();
		G1 part12 = new G1();
		Mcl.mul(part12, (G1) MSK.get(1), new Fr(-1));
		G1 part22 = new G1();
		Mcl.mul(part22, (G1) MSK.get(4), r);
		G1 part32 = new G1();
		Mcl.mul(part32, g, z1);
		Mcl.add(D2, part12, part22);
		Mcl.add(D2, D2, part32);
		
		G1 D3 = new G1();
		Mcl.mul(D3, g, b);
		Mcl.mul(D3, D3, z1);
		Mcl.mul(D3, D3, new Fr(-1));
		
		G1 D4 = new G1();
		G1 part14 = new G1();
		Mcl.mul(part14, (G1) MSK.get(5), r); 
		G1 part24 = new G1();
		Mcl.mul(part24, g, z2);
		Mcl.add(D4, part14, part24);
		
		G1 D5 = new G1();
		Mcl.mul(D5, g, b);
		Mcl.mul(D5, D5, z2);
		Mcl.mul(D5, D5, new Fr(-1));
		
		G2 D6 = new G2((G2) MSK.get(6));
		Mcl.mul(D6, D6, r2);
		
		G2 D7 = new G2();
		Mcl.mul(D7, gg, r1); //Has to be of type G2 to compute the pairing in decryption
		
		G1 K = new G1();
		G1 part1K = new G1();
		Mcl.mul(part1K, (G1) MSK.get(16), I);
		G1 part2K = new G1();
		Mcl.mul(part2K, (G1) MSK.get(15), tagK);
		Mcl.add(K, part1K, part2K);
		Mcl.add(K, K, (G1) MSK.get(17));
		Mcl.mul(K, K, r1);
		
		//Finally, add all to the secret key
		ArrayList<Object> SK = new ArrayList<Object>();
		SK.addAll(Arrays.asList(D1, D2, D3, D4, D5, D6, D7, K, tagK));
		return SK;
		
	}
	
	//input: ciphertext CT, secret key SKI
	//output: message M (GT)
	public static GT decrypt(ArrayList<Object> CT, ArrayList<Object> SKI) {
		
		//Step 1: calculate A1 = e(C1, D1) * e(C2, D2) * ... * e(C5, D5)
		GT initialPairing = new GT();
		Mcl.pairing(initialPairing, (G1) SKI.get(0), (G2) CT.get(1));
		GT A1 = new GT(initialPairing);
		
		for (int i = 1; i < 5; i++) {
			GT e = new GT();
			Mcl.pairing(e, (G1) SKI.get(i), (G2) CT.get(i + 1));
			Mcl.mul(A1, A1, e);
		}
		
		//Step 2: calculate A2 = e(C6, D6) * e(C7, D7)
		GT p1 = new GT();
		Mcl.pairing(p1, (G1) CT.get(6), (G2) SKI.get(5));
		GT p2 = new GT();
		Mcl.pairing(p2, (G1) CT.get(7), (G2) SKI.get(6));
		GT A2 = new GT();
		Mcl.mul(A2, p1, p2);
		
		//Step 3: A3 = A1/A2
		GT A3 = new GT();
		Mcl.pow(A2, A2, new Fr(-1));
		Mcl.mul(A3, A1, A2);
		
		//Step 4: A4 = (e(E1, D7)/e(E2, K))^(1/(tagC - tagK))
		GT A4 = new GT();
		GT e1d7 = new GT();
		Mcl.pairing(e1d7, (G1) CT.get(8), (G2) SKI.get(6));
		GT e2K = new GT();
		Mcl.pairing(e2K, (G1) SKI.get(7), (G2) CT.get(9));
		Mcl.pow(e2K, e2K, new Fr(-1));
		Mcl.mul(A4, e1d7, e2K);
		Fr tagC = (Fr) CT.get(CT.size() - 1);
		Fr tagK = (Fr) SKI.get(SKI.size() - 1);
		Fr exp = new Fr();
		Mcl.sub(exp, tagC, tagK);
		Mcl.div(exp, new Fr(1), exp);
		Mcl.pow(A4, A4, exp);

		
		//FINALLY: M = (C0/(A3/A4))
		GT denominator = new GT();
		Mcl.pow(A4, A4, new Fr(-1));
		Mcl.mul(denominator, A3, A4);
		Mcl.pow(denominator, denominator, new Fr(-1));
		GT M = new GT();
		Mcl.mul(M, (GT) CT.get(0), denominator);
		
		return M;
		
	}
		
	
	//generates a random byte array -- helper function used to obtain random generators in G: NOT NECESSARY
	private static byte[] generateRandomBytes() {
		
		SecureRandom random = new SecureRandom();
		byte[] bytes = new byte[20];
		random.nextBytes(bytes);
		return bytes;
		
	}
	
	private static long[] printRuntimes(int lambda) {
		
		long startSetup = System.nanoTime();
		Object[] setup = setup(lambda);
		long elapsedSetup = System.nanoTime() - startSetup;
		double secondsSetup = ((double) elapsedSetup) / 1E6;
		
		ArrayList<Object> PK = (ArrayList<Object>) setup[0];
		ArrayList<Object> MSK = (ArrayList<Object>) setup[1];
		
		//generate a nondegenerate message M using pairings
		G1 g1 = new G1();
		Mcl.hashAndMapToG1(g1, "abc".getBytes());
		G2 g2 = new G2();
		Mcl.hashAndMapToG2(g2, "def".getBytes());
		GT M = new GT();
		Mcl.pairing(M, g1, g2);
		
		//generate a random value for the key I
		Fr I = new Fr();
		I.setByCSPRNG();
		
		
		//generate ciphertext and keys
		long startEncrypt = System.nanoTime();
		ArrayList<Object> CT = encrypt(PK, I, M);
		long elapsedEncrypt = System.nanoTime() - startEncrypt;
		double secondsEncrypt = ((double) elapsedEncrypt) / 1E6;
		long startKeyGen = System.nanoTime();
		ArrayList<Object> SKI = keyGen(MSK, I);
		long elapsedKeyGen = System.nanoTime() - startKeyGen;
		double secondsKeyGen = ((double) elapsedKeyGen) / 1E6;
		
		//finally, decrypt
		long startDecrypt = System.nanoTime();
		GT M1 = decrypt(CT, SKI);
		long elapsedDecrypt = System.nanoTime() - startDecrypt;
		double secondsDecrypt = ((double) elapsedDecrypt) / 1E6;
		
		String result = (M.equals(M1)) ? "SUCCESSFUL DECRYPTION" : "FAILED DECRYPTION";
		System.out.println(result);
		System.out.println(); //padding
		System.out.println("setup took " + secondsSetup + " milliseconds");
		System.out.println("encryption took " + secondsEncrypt + " milliseconds");
		System.out.println("key generation took " + secondsKeyGen + " milliseconds");
		System.out.println("decryption took " + secondsDecrypt + " milliseconds");
		System.out.println(); //more padding
		
		long[] elapsedTimes = new long[4];
		elapsedTimes[0] = elapsedSetup;
		elapsedTimes[1] = elapsedEncrypt;
		elapsedTimes[2] = elapsedKeyGen;
		elapsedTimes[3] = elapsedDecrypt;
		
		return elapsedTimes;
		
	}
	
	//n = # of times to test
	//lambda = Mcl.BLS12_381 or Mcl.BN254
	public static double[] testRuntimes(int lambda, int n) {
		
		System.out.println("Testing runtimes, n = " + n);
		
		
		long totalSetupTime = 0, totalEncryptionTime = 0, totalKeyGenTime = 0, totalDecryptionTime = 0;
		
		for (int i = 0; i < n; i++) {
			long[] elapsedTimes = printRuntimes(lambda);
			if (i == 0) { continue; } //for some reason, the setup time for the VERY FIRST run is always an outlier, do not include in average
			
			totalSetupTime += elapsedTimes[0];
			totalEncryptionTime += elapsedTimes[1];
			totalKeyGenTime += elapsedTimes[2];
			totalDecryptionTime += elapsedTimes[3];
		}
		
		double averageSetupTime = ((double) totalSetupTime) / (1E6 * n);
		double averageEncryptionTime = ((double) totalEncryptionTime) / (1E6 * n);
		double averageKeyGenTime = ((double) totalKeyGenTime) / (1E6 * n);
		double averageDecryptionTime = ((double) totalDecryptionTime) / (1E6 * n);
		
		System.out.println("Average setup time: " + averageSetupTime + " milliseconds");
		System.out.println("Average encryption time: " + averageEncryptionTime + " milliseconds");
		System.out.println("Average key generation time: " + averageKeyGenTime + " milliseconds");
		System.out.println("Average decryption time: " + averageDecryptionTime + " milliseconds");
		System.out.println();
		System.out.println("The average times do not include the data from the very first run because the first setup time always tends to be an outlier");
		
		//return a double array -- could run some statistical analyses on the sampling distributions of mean times for a sample size of n = 100
		double[] result = new double[4];
		result[0] = averageSetupTime;
		result[1] = averageEncryptionTime;
		result[2] = averageKeyGenTime;
		result[3] = averageDecryptionTime;
		
		return result;
		
	}
	
	public static void main(String[] args) {
		//change the file directory here
		System.load("/Users/arushchhatrapati/Documents/mcl/lib/libmcljava.dylib");
		testRuntimes(Mcl.BN254, 100);
		testRuntimes(Mcl.BLS12_381, 100); //this should be a lot slower
	}
			
}
