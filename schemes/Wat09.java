package schemes;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ThreadLocalRandom;
import java.io.*;
import helperclasses.miscellaneous.Tools;
import com.herumi.mcl.*;

//The scheme: https://eprint.iacr.org/2009/385.pdf (pages 28-29)
//Added g2 to the public key (g1 is in the private key), both g1^(b) and g2^(b) are included in the public key
//In paper: Explain how random elements u1, ..., uN were chosen
//mention that the KeyGen algorithm only generates the key for one user
public class Wat09 { 
	
	//input: security parameter lambda (Mcl.BN254 or Mcl.BLS_381)
	//output: (public key PK, secret key MSK)
	public static Object[] setup(int lambda, int n) { 
		
		Mcl.SystemInit(lambda); // curveType = Mcl.BN254 or Mcl.BLS12_381
		
		//Step 1: choose random generators
		
		G1 g = new G1();
		Mcl.hashAndMapToG1(g, "abc".getBytes());
		
		G2 gg = new G2();
		Mcl.hashAndMapToG2(gg, "def".getBytes());
		
		G1 v = new G1();
		Mcl.hashAndMapToG1(v, "ghi".getBytes());
		
		G1 v1 = new G1();
		Mcl.hashAndMapToG1(v1, "jkl".getBytes());

		G1 v2 = new G1();
		Mcl.hashAndMapToG1(v2, "mno".getBytes());
		
		G1 w = new G1();
		Mcl.hashAndMapToG1(w, "pqr".getBytes());
		
		//Generate random group elements u1, u2, ..., un, put them into an array
		byte[] bts = Tools.generateRandomBytes(n + 2);
		G1[] uElements = new G1[n];
		for (int i = 0; i < n; i++) {
			G1 u = new G1();
			byte[] randomBytes = Arrays.copyOfRange(bts, i, i + 3);
			Mcl.hashAndMapToG1(u, randomBytes);
			uElements[i] = u;
		}
		
		G1 h = new G1();
		Mcl.hashAndMapToG1(h, "vwx".getBytes());
		
		//Step 2: choose random exponents in Z_p
		
		Fr a1 = new Fr();
		a1.setByCSPRNG();
		
		Fr a2 = new Fr();
		a2.setByCSPRNG();
		
		Fr b = new Fr();
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
		
		//Element #0 in the public key is gg, so skip to e1
		G1 e1 = new G1();
		Mcl.mul(e1, g, b);
		
		G2 e2 = new G2();
		Mcl.mul(e2, gg, b);
		
		G2 e3 = new G2();
		Mcl.mul(e3, gg, a1);
		
		G2 e4 = new G2();
		Mcl.mul(e4, gg, a2);
		
		G2 e5 = new G2();
		Mcl.mul(e5, e3, b);
		
		G2 e6 = new G2();
		Mcl.mul(e6, e4, b);
		
		//e7 and e8 are just tau1 and tau2, so skip to e9
		G1 e9 = new G1();
		Mcl.mul(e9, tau1, b);
		
		G1 e10 = new G1();
		Mcl.mul(e10, tau2, b);
		
		//e11, e12, and e13 are w, uElements, and h, respectively, so skip to e14
		GT e14 = new GT();
		Mcl.pairing(e14, g, e5); //e(g, gg)^(alpha * a1 * b) = e(g, gg^(a1 * b))^alpha
		Mcl.pow(e14, e14, alpha);
		
		PK.addAll(Arrays.asList(gg, e1, e2, e3, e4, e5, e6, tau1, tau2, e9, e10, w, uElements, h, e14));
		
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
	public static ArrayList<Object> encrypt(ArrayList<Object> PK, ArrayList<Integer> S, GT M) {
		
		//Step 1: choose random values in Z_p
		
		Fr s1 = new Fr();
		s1.setByCSPRNG();
		
		Fr s2 = new Fr();
		s2.setByCSPRNG();
		
		Fr t = new Fr();
		t.setByCSPRNG();
		
		//Step 2: calculate C0, C1, C2, ...
		
		GT C0 = new GT();
		Mcl.pow(C0, (GT) PK.get(14), s2);
		Mcl.mul(C0, M, C0);
		
		G2 C1 = new G2();
		Fr exp1 = new Fr();
		Mcl.add(exp1, s1, s2);
		Mcl.mul(C1, (G2) PK.get(2), exp1);
		
		G2 C2 = new G2();
		Mcl.mul(C2, (G2) PK.get(5), s1);
		
		G2 C3 = new G2();
		Mcl.mul(C3, (G2) PK.get(3), s1);
		
		G2 C4 = new G2();
		Mcl.mul(C4, (G2) PK.get(6), s2);
		
		G2 C5 = new G2();
		Mcl.mul(C5, (G2) PK.get(4), s2);
		
		G1 C6 = new G1();
		G1 part6 = new G1();
		Mcl.mul(part6, (G1) PK.get(8), s2);
		Mcl.mul(C6, (G1) PK.get(7), s1);
		Mcl.add(C6, C6, part6);
		
		G1 C7 = new G1();
		G1 part17 = new G1();
		Mcl.mul(part17, (G1) PK.get(9), s1);
		G1 part27 = new G1();
		Mcl.mul(part27, (G1) PK.get(10), s2);
		G1 part37 = new G1();
		Fr exp37 = new Fr();
		Mcl.mul(exp37, t, new Fr(-1));
		Mcl.mul(part37, (G1) PK.get(11), exp37);
		Mcl.add(C7, part17, part27);
		Mcl.add(C7, C7, part37);
		
		
		G1[] uElements = (G1[]) PK.get(12);
		G1 E1 = new G1((G1) uElements[S.get(0) - 1]);
		for (int i = 1; i < S.size(); i++) {
			int j = S.get(i);
			Mcl.add(E1, E1, uElements[j - 1]);
		}
		Mcl.mul(E1, E1, t);
		
		G2 E2 = new G2(); //Has to be of type G2 to compute the pairing in decryption
		Mcl.mul(E2, (G2) PK.get(0), t);
		
		//Add all to the ciphertext
		ArrayList<Object> CT = new ArrayList<Object>();
		CT.addAll(Arrays.asList(C0, C1, C2, C3, C4, C5, C6, C7, E1, E2));
		return CT;
	}
	
	//generates a key for the kth user
	public static ArrayList<Object> keyGen(ArrayList<Object> MSK, int k) {
		
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
		
		G1 D3 = new G1((G1) MSK.get(7));
		Mcl.mul(D3, D3, z1);
		Mcl.mul(D3, D3, new Fr(-1));
		
		G1 D4 = new G1();
		G1 part14 = new G1();
		Mcl.mul(part14, (G1) MSK.get(5), r); 
		G1 part24 = new G1();
		Mcl.mul(part24, g, z2);
		Mcl.add(D4, part14, part24);
		
		G1 D5 = new G1((G1) MSK.get(7));
		Mcl.mul(D5, D5, z2);
		Mcl.mul(D5, D5, new Fr(-1));
		
		G2 D6 = new G2((G2) MSK.get(8));
		Mcl.mul(D6, D6, r2);
		
		G2 D7 = new G2();
		Mcl.mul(D7, (G2) MSK.get(6), r1); //Has to be of type G2 to compute the pairing in decryption
		
		G1 K = new G1();
		G1[] uElements = (G1[]) MSK.get(18);
		G1 part1K = new G1();
		Mcl.add(K, (G1) MSK.get(17), uElements[k - 1]);
		Mcl.mul(K, K, r1);
		G1[] keys = new G1[uElements.length];
		
		for (int i = 0; i < uElements.length; i++) {
			if (i == k - 1) {
				keys[i] = K;
			}
			else {
				G1 Ki = new G1();
				Mcl.mul(Ki, uElements[i], r1);
				keys[i] = Ki;
			}
		}
		
		//Finally, add all to the secret key
		ArrayList<Object> SK = new ArrayList<Object>();
		SK.addAll(Arrays.asList(D1, D2, D3, D4, D5, D6, D7, keys));
		return SK;
		
	}
	
	//input: ciphertext CT, secret key SKI
	//output: message M (GT)
	public static GT decrypt(ArrayList<Object> CT, ArrayList<Integer> S, ArrayList<Object> SKK) {
		
		//Step 1: calculate A1 = e(C1, D1) * e(C2, D2) * ... * e(C5, D5)
		GT initialPairing = new GT();
		Mcl.pairing(initialPairing, (G1) SKK.get(0), (G2) CT.get(1));
		GT A1 = new GT(initialPairing);
		
		for (int i = 1; i < 5; i++) {
			GT e = new GT();
			Mcl.pairing(e, (G1) SKK.get(i), (G2) CT.get(i + 1));
			Mcl.mul(A1, A1, e);
		}
		
		//Step 2: calculate A2 = e(C6, D6) * e(C7, D7)
		GT p1 = new GT();
		Mcl.pairing(p1, (G1) CT.get(6), (G2) SKK.get(5));
		GT p2 = new GT();
		Mcl.pairing(p2, (G1) CT.get(7), (G2) SKK.get(6));
		GT A2 = new GT();
		Mcl.mul(A2, p1, p2);
		
		//Step 3: A3 = A1/A2
		GT A3 = new GT();
		Mcl.pow(A2, A2, new Fr(-1));
		Mcl.mul(A3, A1, A2);
		
		//Step 4: Calculate A4
		GT A4 = new GT();
		GT e1d7 = new GT();
		Mcl.pairing(e1d7, (G1) CT.get(8), (G2) SKK.get(6));
		GT e2K = new GT();
		G1[] keys = (G1[]) SKK.get(7);
		G1 product = new G1((G1) keys[S.get(0) - 1]);
		for (int i = 1; i < S.size(); i++) {
			int j = S.get(i);
			Mcl.add(product, product, keys[j - 1]);
		}
		Mcl.pairing(e2K, product, (G2) CT.get(9));
		Mcl.pow(e1d7, e1d7, new Fr(-1));
		Mcl.mul(A4, e2K, e1d7);
		
		//FINALLY: M = (C0/(A3/A4))
		GT denominator = new GT();
		Mcl.pow(A4, A4, new Fr(-1));
		Mcl.mul(denominator, A3, A4);
		Mcl.pow(denominator, denominator, new Fr(-1));
		GT M = new GT();
		Mcl.mul(M, (GT) CT.get(0), denominator);
		
		return M;
		
	}
		
	
	private static long[] printRuntimes(int lambda, int n, int subsetSize) {
		
		long startSetup = System.nanoTime();
		Object[] setup = setup(lambda, n);
		long elapsedSetup = System.nanoTime() - startSetup;
		double secondsSetup = ((double) elapsedSetup) / 1E9;
		
		ArrayList<Object> PK = (ArrayList<Object>) setup[0];
		ArrayList<Object> MSK = (ArrayList<Object>) setup[1];
		
		//generate a nondegenerate message M using pairings
		G1 g1 = new G1();
		Mcl.hashAndMapToG1(g1, "abc".getBytes());
		G2 g2 = new G2();
		Mcl.hashAndMapToG2(g2, "def".getBytes());
		GT M = new GT();
		Mcl.pairing(M, g1, g2);
		
		//generate a random subset S
		ArrayList<Integer> S = new ArrayList<Integer>();
		
		//randomly generate numbers to put in the subset S (NO REPEATS)
		ArrayList<Integer> randomNums = new ArrayList<Integer>();
			for (int i = 0; i < n; i++) { 
				randomNums.add(i + 1);
			}
			for (int i = 0; i < subsetSize; i++) {
				int randomIndex = ThreadLocalRandom.current().nextInt(0, randomNums.size());
				int randomID = randomNums.get(randomIndex);
				S.add(randomID);
				randomNums.remove(randomIndex);
			}
				
		//Get random user ID i to test the decryption
		int i = S.get(ThreadLocalRandom.current().nextInt(0, S.size()));
		
		
		//generate ciphertext and keys
		long startEncrypt = System.nanoTime();
		ArrayList<Object> CT = encrypt(PK, S, M);
		long elapsedEncrypt = System.nanoTime() - startEncrypt;
		double secondsEncrypt = ((double) elapsedEncrypt) / 1E9;
		long startKeyGen = System.nanoTime();
		ArrayList<Object> SKK = keyGen(MSK, i);
		long elapsedKeyGen = System.nanoTime() - startKeyGen;
		double secondsKeyGen = ((double) elapsedKeyGen) / 1E9;
		
		//finally, decrypt
		long startDecrypt = System.nanoTime();
		GT M1 = decrypt(CT, S, SKK);
		long elapsedDecrypt = System.nanoTime() - startDecrypt;
		double secondsDecrypt = ((double) elapsedDecrypt) / 1E9;
		
		String result = (M.equals(M1)) ? "SUCCESSFUL DECRYPTION" : "FAILED DECRYPTION";
		System.out.println(result + ": n = " + n + ", subset size = " + subsetSize);
		System.out.println(); //padding
		System.out.println("setup took " + secondsSetup + " seconds");
		System.out.println("encryption took " + secondsEncrypt + " seconds");
		System.out.println("key generation took " + secondsKeyGen + " seconds");
		System.out.println("decryption took " + secondsDecrypt + " seconds");
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
	public static void testRuntimes(int lambda, int percent) {
		for (int N = 100; N <= 1000000; N *= 10) {
			int subsetSize = (int) (0.01 * percent * N);
			printRuntimes(lambda, N, subsetSize);
		}
	}
	
	public static void main(String[] args) {

		File lib = new File("../../lib/libmcljava.dylib");
		System.load(lib.getAbsolutePath());
		testRuntimes(Mcl.BN254, 10);
	}
	
		
}


