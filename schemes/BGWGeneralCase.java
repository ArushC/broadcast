package schemes;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ThreadLocalRandom;
import helperclasses.Tools;
import com.herumi.mcl.*;

//How to use the library: https://github.com/herumi/mcl/blob/master/ffi/java/java.md
//This is the BGW scheme: https://eprint.iacr.org/2005/018.pdf
//Description of BLS: https://crypto.stanford.edu/~dabo/pubs/papers/BLSmultisig.html

//This implements a special case in which A = 1 and B = n
public class BGWGeneralScheme {

	
	//From https://books.google.com/books?id=kb7ZzFrJi48C&pg=PA106&lpg=PA106&dq=calculate+e(g,+g)+pairing&source=bl&ots=Taa2npE3kj&sig=ACfU3U1zNODtV4H_2pUbgscbcQsYljeujg&hl=en&sa=X&ved=2ahUKEwjD_q6hj-LpAhXLTN8KHS3eAU8Q6AEwCXoECAcQAQ#v=onepage&q=calculate%20e(g%2C%20g)%20pairing&f=falsev
	//e(g, g^(ab)) = e(g, g)^(ab) = e(g^a, g^b)
	//so instantiate g as a G1 group and make all the others G2 groups
	
	private static Fr alpha, t;
	private static G1 g;
	private static G2 gg;
	private static int n, A, B;
	private static GT K; //to be precomputed in the setup function
	
	//Input: n = # of receivers
	//Output: Object[]
	//first element:  list containing public key PK
	//second element: list containing private keys: d1, d2, ... dn 
	public static Object[] setup(int n) {
		
		BGWGeneralScheme.n = n; //save n so it can be used for other functions
		
		//instantiate B = floor(sqrt(n)),  A = ceil(n/B)
		B = (int) Math.floor(Math.pow((double) n, 0.5));
		A = (int) Math.ceil(((double) n) / B);

		//initialize random generator g in G
		g = new G1();
		Mcl.hashAndMapToG1(g, "abc".getBytes());
		
		//get corresponding element in G2
		gg = new G2();
		Mcl.hashAndMapToG2(gg, "def".getBytes());
		
		//random alpha and t in Z_p 
		alpha = new Fr();
		alpha.setByCSPRNG(); 
		t = new Fr();
		t.setByCSPRNG();
		
		//precompute K
		precompute();
		
		//Instantiate public key PK
		ArrayList<G2> PK = new ArrayList<G2>();
		PK.add(gg); //add group element gg (type G2)
		
		//from i = 1 to i = 2n (ALL)
		for (int i = 1; i <= 2*B; i++) {
			G2 pub = new G2();
			Fr exp = Tools.power(alpha, i); //exp = alpha^i
			Mcl.mul(pub, gg, exp); // g_n = g^(exp)
			PK.add(pub);
		}
		
		//generate random betas in Z_p and use to generate v1, v2, ..., vA
		ArrayList<Fr> randomBetas = new ArrayList<Fr>();
		for (int i = 1; i <= A; i++) {
			Fr beta = new Fr();
			beta.setByCSPRNG();
			G2 v = new G2();
			Mcl.mul(v, gg, beta);
			randomBetas.add(beta);
			PK.add(v);
		}
		
		
		//calculate private keys and put them into a list
		ArrayList<G2> privateKeys = new ArrayList<G2>();
		for (int i = 1; i <= n; i++) {
			G2 priv = new G2();
			int b = (i % B == 0) ? B : (i % B);
			int a = (int) Math.ceil(((double) i) / B); 
			Mcl.mul(priv, (G2) PK.get(b), randomBetas.get(a - 1)); // d_i = (g_i)^(beta)
			privateKeys.add(priv);
		}
		//return public key & private key
		Object[] result = new Object[2];
		result[0] = PK;
		result[1] = privateKeys;
		return result;
	}
	
	//precomputes K = e(gB+1, g)^t
	private static void precompute() {
		//calculate e(g, g_(B+1))
		GT e = new GT();
		G2 gBPlus1 = new G2();
		Fr exp = Tools.power(alpha, B+1);
		Mcl.mul(gBPlus1, gg, exp); //gn = g^(alpha^n), so gn+1 = g^(alpha^(n+1)) = g^(exp)
		Mcl.pairing(e, g, gBPlus1);
		K = new GT();
		Mcl.pow(K, e, t); // K = e(g, gn+1)^t
	}
	
	private static ArrayList<ArrayList<Integer>> computeSLSubsets(ArrayList<Integer> S, ArrayList<G2> PK) {
		
		ArrayList<ArrayList<Integer>> partSubsets = new ArrayList<ArrayList<Integer>>();
		
		for (int l = 1; l <= A; l++) {
			ArrayList<Integer> part = new ArrayList<Integer>();
			for (int i = 1; i <= B; i++) {
				part.add((l - 1) * B + i);  //part = {(l - 1)B + 1, (l - 1)B + 2, ..., lB} for l = 1, 2, ..., A
			}
			partSubsets.add(part);
		}
		
		ArrayList<ArrayList<Integer>> sHatLSubsets = new ArrayList<ArrayList<Integer>>(); //sHatL subsets = elements in S and each part subset
		for (ArrayList<Integer> part : partSubsets) {
			ArrayList<Integer> sHatL = new ArrayList<Integer>();
			for (Integer i: part) {
				if (S.contains(i))
					sHatL.add(i);
			}
			sHatLSubsets.add(sHatL);
		}
		
		//System.out.println(sHatLSubsets);
		
		int l = 1;
		ArrayList<ArrayList<Integer>> sLSubsets = new ArrayList<ArrayList<Integer>>();
		for (ArrayList<Integer> sHatL : sHatLSubsets) {
			ArrayList<Integer> sL = new ArrayList<Integer>();
			for (Integer x: sHatL) {
				sL.add(x - l * B + B); //check to make sure this is a subset of {1, ..., B}
			}
			l += 1;
			sLSubsets.add(sL);
		}
		
		return sLSubsets;
		
	}
	
	//Input: S = subset to which the message is brodcast, PK = public key
	//Output: Hdr = header (broadcast ciphertext), K = message encryption key
	public static Object[] encrypt(ArrayList<Integer> S, ArrayList<G2> PK) {
			 		
		ArrayList<ArrayList<Integer>> sLSubsets = computeSLSubsets(S, PK);
		ArrayList<Object> Hdr = new ArrayList<Object>();
		
		//calculate C_0 (first element in Hdr) and add to Hdr
		G1 c0 = new G1();
		Mcl.mul(c0, g, t); //C_0 = g^t
		Hdr.add(c0);
		
		//calculate rest of Hdr
		
		for (int i = 1; i <= A; i++) {
			G2 v = (G2) PK.get(2 * B + i);
			G2 product = new G2(v);
			for (int j : sLSubsets.get(i - 1)) {
				Mcl.add(product, product, (G2) PK.get(B + 1 - j)); //product *= g_(B + 1 - j)
			}
			Mcl.mul(product, product, t);
			Hdr.add(product);
		}
			
		//return Hdr and K
		Object[] result = new Object[2];
		result[0] = Hdr;
		result[1] = K;
		return result;
			
	}
	
	
	//Input: S = subset, i = user id, di = user private key, Hdr = header, PK = public key
	//Output: if i is in S, output message encryption key K. Use to decode C (brodcast body)
	public static GT decrypt(ArrayList<Integer> S, int i, G2 di, ArrayList<Object> Hdr, ArrayList<G2> PK) {
		
		//get a and b (i = (a - 1)b + B)
		int b = (i % B == 0) ? B : (i % B);
		int a = (int) Math.ceil(((double) i) / B);
		System.out.println("i = " + i + ": " + "a = " + a + ", b = " + b + ", B = " + B);
		
		//calculate the first pairing: e(g_b, C_a)
		GT e1 = new GT();
		G2 cA = (G2) Hdr.get(a);
		G1 gB = new G1();
		Fr exp = Tools.power(alpha, b);
		Mcl.mul(gB, g, exp);
		Mcl.pairing(e1, gB, cA);
		
		//calculate the second pairing: e(c0, big messy expression) --> see the paper, page 8: https://eprint.iacr.org/2005/018.pdf 
		G1 c0 = (G1) Hdr.get(0);
		ArrayList<Integer> sA = computeSLSubsets(S, PK).get(a - 1);
		
		G2 product = new G2(di);
		for (Integer j: sA) {	
			if (j == b) { //j != b
				continue;
			}
			Mcl.add(product, product, (G2) PK.get(B + 1 - j + b)); //product *= g_(B + 1 - j + b)
		}
		
		//finally, compute the pairing
		GT e2 = new GT();
		Mcl.pairing(e2, c0, product);
		
		GT K_R = new GT();
		//to compute e1 / e2, calculate e2^(-1) and output e1 * e2^(-1)
		Mcl.pow(e2, e2, new Fr(-1)); //CONFIRMED this works (after testing)
		Mcl.mul(K_R, e1, e2);
		
		return K_R;
		
	}
	
	//TESTING -------------------------------------------------------------------------------
	
	public static void assertEquals(String msg, String x, String y) {
		if (x.equals(y)) {
			System.out.println("OK : " + msg);
		} 
		else {
			System.out.println("NG : " + msg + ", x = " + x + ", y = " + y);
		}
	}
	
	public static Object[] testSetup(int n) {
		//TEST SETUP: returns the setup Object[]
		Object[] setup = setup(n);
		//EXTRACT EVERYTHING
		ArrayList<Object> PK = (ArrayList<Object>) setup[0];
		ArrayList<G2> privateKeys = (ArrayList<G2>) setup[1];
		for (int i = 1; i <= n; i++) {
			//Verify for each key: di = vi (v^(alpha^(i)))
			G2 di = privateKeys.get(i - 1); // (i - 1) because private keys starts from i = 1
			int b = (i % B == 0) ? B : (i % B);
			int a = (int) Math.ceil(((double) i) / B); 
			Fr exp = Tools.power(alpha, b);
			System.out.println("Exponent = " + exp.toString());
			G2 vi = new G2();
			Mcl.mul(vi, (G2) PK.get(2 * B + a), exp); //vi = v^(exp)
			assertEquals("i = " + i + " ", di.toString(), vi.toString());
		}
			
		return setup;
	}
	
	//must be called AFTER setup function has been called, otherwise instance variables will not be instantiated
	public static void testDecrypt(ArrayList<Integer> S, int i, G2 di, ArrayList<Object> Hdr, ArrayList<G2> PK) {
		GT K1 = decrypt(S, i, di, Hdr, PK);
		System.out.println("K = " + K.toString());
		System.out.println("K1 = " + K1.toString());
	}
	
	
	//Test the runtimes of the algorithms
	public static void testRuntimes() {
		
		//see how runtime changes with constant n and increasing subset size
		long totalSetupTime = 0;
		
		for (int i = 100; i <= 2000; i+=100) {
			long[] elapsedTimes = printRuntimes(10000, i);
			totalSetupTime += elapsedTimes[0];
		}
		
		double averageSetupTime = ((double) totalSetupTime) / (1E9 * 20);
		
		
		long totalEncryptionTime = 0;
		long totalDecryptionTime = 0;
		//see how runtime changes with increasing n, constant subset size = 100
		for (int i = 1000; i <= 20000; i += 1000) {
			long[] elapsedTimes =  printRuntimes(i, 100);
			totalEncryptionTime += elapsedTimes[1];
			totalDecryptionTime += elapsedTimes[2];
		}
		
		double averageEncryptionTime = ((double) totalEncryptionTime) / (1E9 * 20);
		double averageDecryptionTime = ((double) totalDecryptionTime) / (1E9 * 20);
		
		System.out.println("Average setup time, constant n = 10000: " + averageSetupTime + " seconds");
		System.out.println("Average encryption time, constant subset size = 100: " + averageEncryptionTime + " seconds");
		System.out.println("Average decryption time, constant subset size = 100: " + averageDecryptionTime + " seconds");
	
	}
	
		
	//returns long[] containing elapsed time for setup, encrypt, and decrypt, respectively
	private static long[] printRuntimes(int n, int subsetSize) {
		
		//Get elapsed time for setup(n) 
		long startSetup = System.nanoTime();
		Object[] setup = setup(n);
		long elapsedSetup = System.nanoTime() - startSetup;
		double secondsSetup = ((double) elapsedSetup) / 1E9;
		//extract public/private key from setup
		ArrayList<G2> PK = (ArrayList<G2>) setup[0];
		ArrayList<G2> privateKeys = (ArrayList<G2>) setup[1];
		
		
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
		ArrayList<Object> Hdr = (ArrayList<Object>) encrypted[0];
		G2 di = privateKeys.get(i - 1);
		
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
	
	public static void main(String[] args) {
		//change the file directory here
		System.load("/Users/arushchhatrapati/Documents/mcl/lib/libmcljava.dylib");
		Mcl.SystemInit(Mcl.BN254); // curveType = Mcl.BN254 or Mcl.BLS12_381
		testRuntimes();
		/*Object[] setup = testSetup(100);
		ArrayList<G2> PK = (ArrayList<G2>) setup[0];
		ArrayList<G2> privateKeys = (ArrayList<G2>) setup[1];
		int i = 40;
		G2 di = privateKeys.get(i - 1);
		ArrayList<Integer> S = new ArrayList<Integer>();
		S.addAll(Arrays.asList(1, 5, 31, 30, 40, 45, 22, 21, 20, 19, 18, 100));
		Object[] C = encrypt(S, PK);
		ArrayList<Object> Hdr = (ArrayList<Object>) C[0];
		testDecrypt(S, i, di, Hdr, PK);*/
		
		
		
		
	}

}


