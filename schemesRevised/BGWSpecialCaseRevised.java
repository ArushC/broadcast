package schemesRevised;
import java.io.*;
import java.util.concurrent.ThreadLocalRandom;
import java.util.ArrayList;
import java.util.Arrays;
import com.herumi.mcl.*;
import helperclasses.Tools;

//This is the BGW scheme: https://eprint.iacr.org/2005/018.pdf (3.1)
//Changes made: optimized selection of G1 and G2, included gg^(alpha^i) for all i in the public key,
//              both g1 and g2 are included in the public key
//MENTION IN PAPER: precomputation of K
public class BGWSpecialCaseRevised {

	private static int n;
	private static GT K; //to be precomputed in the setup function
	
	//Input: n = # of receivers
	//Output: Object[]
	//first element:  list containing public key PK
	//second element: list containing private keys: d1, d2, ... dn 
	public static Object[] setup(int n) {
		
		BGWSpecialCaseRevised.n = n; //save n so it can be used for other functions
		
		//initialize random generator g in G
		G1 g = new G1();
		Mcl.hashAndMapToG1(g, "abc".getBytes());
		
		//get corresponding element in G2
		G2 gg = new G2();
		Mcl.hashAndMapToG2(gg, "def".getBytes());
		
		//random alpha and t in Z_p 
		Fr alpha = new Fr();
		alpha.setByCSPRNG(); 
		
		//precompute K
		precompute(g, gg, alpha);
		
		//Instantiate public key PK
		Object[] PK = new Object[3 * n + 3];
		PK[0] = g;
		PK[1] = gg;
		
		//from i = 1 to i = 2n (ALL), add g^(alpha^i)
		Fr exp = new Fr(alpha);
		for (int i = 1; i <= 2*n; i++) {
			G1 pub = new G1();
			Mcl.mul(pub, g, exp); // g_n = g^(exp)
			PK[i + 1] = pub;
			if (i <= n) {
				G2 pubG2 = new G2();
				Mcl.mul(pubG2, gg, exp);
				PK[2 * n + i + 1] = pubG2; //add g2^(alpha^i) to public key
			}
			Mcl.mul(exp, exp, alpha);
		}

		//random beta in Z_p
		Fr beta = new Fr();
		beta.setByCSPRNG();
		
		//compute v
		G1 v = new G1();
		Mcl.mul(v, g, beta); //v = g^beta
		PK[3 * n + 2] = v;
		
		//calculate private keys and put them into a list
		ArrayList<G1> privateKeys = new ArrayList<G1>();
		for (int i = 1; i <= n; i++) {
			G1 priv = new G1();
			Mcl.mul(priv, (G1) PK[i + 1], beta); // d_i = (g_i)^(beta)
			privateKeys.add(priv);
		}
		//return public key & private key
		Object[] result = new Object[2];
		result[0] = PK;
		result[1] = privateKeys;
		return result;
	}
	
	//precomputes K = e(gn+1, g)
	private static void precompute(G1 g, G2 gg, Fr alpha) {
		//calculate e(g, g_(n+1))
		K = new GT();
		G2 gNPlus1 = new G2();
		Fr exp = Tools.power(alpha, n+1);
		Mcl.mul(gNPlus1, gg, exp); //gn = g^(alpha^n), so gn+1 = g^(alpha^(n+1)) = g^(exp)
		Mcl.pairing(K, g, gNPlus1);
	}
	
	//Input: S = subset to which the message is brodcast, PK = public key
	//Output: Hdr = header (broadcast ciphertext), K = message encryption key
	public static Object[] encrypt(ArrayList<Integer> S, Object[] PK) {
			 		
		Fr t = new Fr();
		t.setByCSPRNG();
		Mcl.pow(K, K, t);
		
		//calculate C_0 (first element in Hdr)
		G2 c0 = new G2();
		Mcl.mul(c0, (G2) PK[1], t); //C_0 = gg^t
		
		//calculate C_1 (second element in Hdr)
		G1 v = (G1) PK[PK.length - 1];
		G1 product = new G1(v);
		
		for (int i = 0; i < S.size(); i++) {
			int j = S.get(i);
			Mcl.add(product, product, (G1) PK[n + 2 - j]); //product *= g_(n + 1 - j)
		}
		
		G1 c1 = new G1();
		Mcl.mul(c1, product, t);      //c1 = c1^(t)
		
		//return Hdr and K
		Object[] Hdr = new Object[2];
		Hdr[0] = c0;
		Hdr[1] = c1;
		Object[] result = new Object[2];
		result[0] = Hdr;
		result[1] = K;
		return result;
			
	}
	
	
	//Input: S = subset, i = user id, di = user private key, Hdr = header, PK = public key
	//Output: if i is in S, output message encryption key K. Use to decode C (brodcast body)
	public static GT decrypt(ArrayList<Integer> S, int i, G1 di, Object[] Hdr, Object[] PK) {
		
		//Calculate e(g_i, C1) = e(g, C1^(alpha^i))
		GT e1 = new GT();
		G1 c1 = (G1) Hdr[1];
		G2 gi = (G2) PK[2 * n + i + 1];
		Mcl.pairing(e1, c1, gi);
		
		//Calculate e(big messy expression, C0) --> see the paper, page 6: https://eprint.iacr.org/2005/018.pdf 
		//In this case, C0 is of type G2, big messy expression is of type G1
		
		G2 c0 = (G2) Hdr[0];
		
		//calculate big messy expression
		G1 product = new G1(di);
		
		for (int k = 0; k < S.size(); k++) {	
			int j = S.get(k);
			if (j == i) {
				continue;
			}
			Mcl.add(product, product, (G1) PK[n + 2 - j + i]); //product *= g_(n + 1 - j + i)
		}
		
		//finally, compute the pairing
		GT e2 = new GT();
		Mcl.pairing(e2, product, c0);
		
		GT K_R = new GT();
		//to compute e1 / e2, calculate e2^(-1) and output e1 * e2^(-1)
		Mcl.pow(e2, e2, new Fr(-1)); //CONFIRMED this works (after testing)
		Mcl.mul(K_R, e1, e2);
		
		return K_R;
		
	}
	
	//returns long[] containing elapsed time for setup, encrypt, and decrypt, respectively
	private static long[] printRuntimes(int n, int subsetSize) {
		
		//Get elapsed time for setup(n) 
		long startSetup = System.nanoTime();
		Object[] setup = setup(n);
		long elapsedSetup = System.nanoTime() - startSetup;
		double secondsSetup = ((double) elapsedSetup) / 1E9;
		//extract public/private key from setup
		Object[] PK = (Object[]) setup[0];
		ArrayList<G1> privateKeys = (ArrayList<G1>) setup[1];
		
		
		ArrayList<Integer> S = new ArrayList<Integer>();
		//randomly generate numbers to put in the subset S (NO REPEATS)
		ArrayList<Integer> randomNums = new ArrayList<Integer>();
		for (int i = 0; i < n; i++) { 
			randomNums.add(i + 1);
		}
		for (int i = 0; i < subsetSize; i++) {
			int randomIndex = ThreadLocalRandom.current().nextInt(1, randomNums.size());
			int randomID = randomNums.get(randomIndex);
			S.add(randomID);
			randomNums.remove(randomIndex);
		}
		
		//Get elapsed time for encrypt
		long startEncrypt = System.nanoTime();
		Object[] encrypted = encrypt(S, PK);
		long elapsedEncrypt = System.nanoTime() - startEncrypt;
		double secondsEncrypt = ((double) elapsedEncrypt) / 1E9;
		//Get random user ID i to test the decryption
		int i = S.get(ThreadLocalRandom.current().nextInt(0, S.size()));
		Object[] Hdr = (Object[]) encrypted[0];
		G1 di = privateKeys.get(i - 1);
		
		//Get elapsed time for decrypt
		long startDecrypt = System.nanoTime();
		GT K1 = decrypt(S, i, di, Hdr, PK);
		long elapsedDecrypt = System.nanoTime() - startDecrypt;
		double secondsDecrypt = ((double) elapsedDecrypt) / 1E9;
		
		//Finally, print out the results
		String success = (K1.equals(K)) ? "SUCCESSFUL DECRYPTION" : "FAILED DECRYPTION";
		System.out.println(success + ": " + "n = " + n + ", subset size = " + subsetSize);
		System.out.println(); //padding
		System.out.println("setup took " + secondsSetup + " seconds");
		System.out.println("encryption took " + secondsEncrypt + " seconds");
		System.out.println("decryption took " + secondsDecrypt + " seconds (i = " + i + ")");
		System.out.println(); //more padding
		
		long[] elapsedTimes = new long[3];
		elapsedTimes[0] = elapsedSetup;
		elapsedTimes[1] = elapsedEncrypt;
		elapsedTimes[2] = elapsedDecrypt;
		
		return elapsedTimes;
	}
	
	
	//Test the runtimes of the algorithms
	public static void testRuntimes(int percent) {
		for (int N = 10; N <= 1000000; N *= 10) {
			int subsetSize = (int) (0.01 * percent * N);
			printRuntimes(N, subsetSize);
		}
		
	}
	
	public static void main(String[] args) {

		File lib = new File("../../lib/libmcljava.dylib");
		System.load(lib.getAbsolutePath());
		Mcl.SystemInit(Mcl.BN254); // curveType = Mcl.BN254 or Mcl.BLS12_381
		testRuntimes(10);
		
	}

}
