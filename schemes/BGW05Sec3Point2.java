package schemes;
import java.io.*;
import java.util.*;
import java.util.concurrent.ThreadLocalRandom;
import com.herumi.mcl.*;

//This is the general case for the BGW scheme: https://eprint.iacr.org/2005/018.pdf (3.2)
//Changes made: optimized selection of G1 and G2, included g2^(alpha^i) for all i in the public key,
//both g1 and g2 are included in the public key, used HASHSETS to calculate the intersection of SHatL and SL
//and added a key generation function instead of computing the private keys for all N users in the setup function
//MAKE SURE TO MENTION IN PAPER: precomputation
public class BGW05Sec3Point2 {

	private static int n, A, B;
	private static GT K; //to be precomputed in the setup function
	
	//Input: n = # of receivers
	//Output: Object[]
	//first element:  public key PK
	//second element: MSK
	public static Object[] setup(int n) {
		
		BGW05Sec3Point2.n = n; //save n so it can be used for other functions
		
		//instantiate B = floor(sqrt(n)),  A = ceil(n/B)
		B = (int) Math.floor(Math.pow((double) n, 0.5));
		A = (int) Math.ceil(((double) n) / B);

		//initialize random generator g in G
		G1 g = new G1();
		Mcl.hashAndMapToG1(g, "abc".getBytes());
		
		//get corresponding element in G2
		G2 gg = new G2();
		Mcl.hashAndMapToG2(gg, "def".getBytes());
		
		//random alpha and t in Z_p 
		Fr alpha = new Fr();
		alpha.setByCSPRNG(); 
		
		//Instantiate public key PK
		Object[] PK = new Object[3 * B + A + 2];
		PK[0] = g;
		PK[1] = gg;
		
		//from i = 1 to i = 2B (except i = B + 1)
		Fr exp = new Fr(alpha);
		for (int i = 1; i <= 2*B; i++) {
			if (i == B + 1) {
				Mcl.mul(exp, exp, alpha);
				continue;
			}
			G1 pub = new G1();
			Mcl.mul(pub, g, exp); // g_n = g^(exp)
			PK[i + 1] = pub;
			if (i <= B) {
				G2 pubG2 = new G2();
				Mcl.mul(pubG2, gg, exp);
				PK[2 * B + i + 1] = pubG2; //add gg^(alpha^i) to public key for i = 1, 2, ..., B
			}
			Mcl.mul(exp, exp, alpha);
		}
		
		//precompute K
		precompute((G1) PK[B + 1], (G2) PK[2 * B + 2]);
		
		//generate random Gammas in Z_p and use to generate v1, v2, ..., vA
		ArrayList<Fr> randomGammas = new ArrayList<Fr>();
		for (int i = 1; i <= A; i++) {
			Fr gamma = new Fr();
			gamma.setByCSPRNG();
			G1 v = new G1();
			Mcl.mul(v, g, gamma);
			randomGammas.add(gamma);
			PK[3 * B + i + 1] = v;
		}
		
		
		//define a master secret key
		Object[] MSK = {randomGammas, PK};
		
		//return public key & private key
		Object[] result = new Object[2];
		result[0] = PK;
		result[1] = MSK;
		return result;
	}
	
	
	
	//generates the key for user i
	public static G1 keyGen(Object[] MSK, int i) {
		
		//extract from MSK
		ArrayList<Fr> randomGammas = (ArrayList<Fr>) MSK[0];
		Object[] PK = (Object[]) MSK[1];
		
		G1 priv = new G1();
		int b = (i % B == 0) ? B : (i % B); //i = (a-1)B + b
		int a = (int) Math.ceil(((double) i) / B); 
		Mcl.mul(priv, (G1) PK[b + 1], randomGammas.get(a - 1)); // d_i = (g_b)^(gamma)
		return priv;
	}
	
	
	//Input: S = subset to which the message is brodcast, PK = public key
	//Output: Hdr = header (broadcast ciphertext), K = message encryption key
	public static Object[] encrypt(HashSet<Integer> S, Object[] PK) {
			 		
		//generate the list of S_l subsets for l = 1, 2, ..., A
		ArrayList<ArrayList<Integer>> sLSubsets = new ArrayList<ArrayList<Integer>>();
		for (int l = 1; l <= A; l++)
			sLSubsets.add(getSLSubset(S, l));
		
		ArrayList<Object> Hdr = new ArrayList<Object>();
		
		Fr t = new Fr();
		t.setByCSPRNG();
		Mcl.pow(K, K, t);
		
		//calculate C_0 (first element in Hdr) and add to Hdr
		G2 c0 = new G2();
		Mcl.mul(c0, (G2) PK[1], t); //C_0 = g^t
		Hdr.add(c0);
		
		//calculate rest of Hdr
		
		for (int i = 1; i <= A; i++) {
			G1 v = (G1) PK[3 * B + i + 1];
			G1 product = new G1(v);
			for (int j : sLSubsets.get(i - 1)) {
				Mcl.add(product, product, (G1) PK[B + 2 - j]); //product *= g_(B + 1 - j)
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
	public static GT decrypt(HashSet<Integer> S, int i, G1 di, ArrayList<Object> Hdr, Object[] PK) {
		
		//get a and b (i = (a - 1)b + B)
		int b = (i % B == 0) ? B : (i % B);
		int a = (int) Math.ceil(((double) i) / B);
		
		//calculate the first pairing: e(g_b, C_a)
		GT e1 = new GT();
		G1 cA = (G1) Hdr.get(a);
		G2 gB = (G2) PK[2 * B + b + 1];
		Mcl.pairing(e1, cA, gB);
		
		//calculate the second pairing: e(c0, big messy expression) --> see the paper, page 8: https://eprint.iacr.org/2005/018.pdf 
		G2 c0 = (G2) Hdr.get(0);
		ArrayList<Integer> sA = getSLSubset(S, a);
		
		//calculate number of values in all subsets
		
		G1 product = new G1(di);
		for (Integer j: sA) {	
			if (j == b) { //j != b
				continue;
			}
			Mcl.add(product, product, (G1) PK[B + 2 - j + b]); //product *= g_(B + 1 - j + b)
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
	
	//precomputes K = e(g_(B+1), g) as e(g_B, gg_1)
	private static void precompute(G1 gB, G2 gg1) {
		//calculate e(g, g_(B+1))
		K = new GT();
		Mcl.pairing(K, gB, gg1);
	}

	private static ArrayList<Integer> getSLSubset(HashSet<Integer> S, int l) {
		
			HashSet<Integer> sHatL = new HashSet<Integer>();
			for (int i = 1; i <= B; i++) {
				sHatL.add((l - 1) * B + i);  //part = {(l - 1)B + 1, (l - 1)B + 2, ..., lB}
			}
				
			sHatL.retainAll(S); //calculate intersection of sHatL and S
			ArrayList<Integer> sL = new ArrayList<Integer>();
			for (Integer x: sHatL) {
				sL.add(x - l * B + B); //check to make sure this is a subset of {1, ..., B}
			}
		
			return sL;
	}

	//Test the runtimes of the algorithms
	public static void testRuntimes(int percent) {
		for (int N = 100; N <= 1000000; N *= 10) {
			int subsetSize = (int) (0.01 * percent * N);
			printRuntimes(N, subsetSize);
		}
	}
	
		
	//returns long[] containing elapsed time for setup, encrypt, and decrypt, respectively
	private static long[] printRuntimes(int n, int subsetSize) {
		
		//Get elapsed time for setup(n) 
		long startSetup = System.nanoTime();
		Object[] setup = setup(n);
		long elapsedSetup = System.nanoTime() - startSetup;
		double secondsSetup = ((double) elapsedSetup) / 1E9;
		//extract public key & MSK from setup
		Object[] PK = (Object[]) setup[0];
		Object[] MSK = (Object[]) setup[1];
		
		int randomID = 0;
		HashSet<Integer> S = new HashSet<Integer>();
		//randomly generate numbers to put in the subset S (NO REPEATS)
		ArrayList<Integer> randomNums = new ArrayList<Integer>();
			for (int i = 0; i < n; i++) { 
				randomNums.add(i + 1);
			}
			for (int i = 0; i < subsetSize; i++) {
				int randomIndex = ThreadLocalRandom.current().nextInt(1, randomNums.size());
				randomID = randomNums.get(randomIndex);
				S.add(randomID);
				randomNums.remove(randomIndex);
			}

		//Get elapsed time for encrypt
		long startEncrypt = System.nanoTime();
		Object[] encrypted = encrypt(S, PK);
		long elapsedEncrypt = System.nanoTime() - startEncrypt;
		double secondsEncrypt = ((double) elapsedEncrypt) / 1E9;
		
		//Get user ID i to test the decryption
		int i = randomID; //int i is the last ID that was added to the set S
		ArrayList<Object> Hdr = (ArrayList<Object>) encrypted[0];
		
		//generate key for user i
		long startKeyGen = System.nanoTime();
		G1 di = keyGen(MSK, i);
		long elapsedKeyGen = System.nanoTime() - startKeyGen;
		double secondsKeyGen = ((double) elapsedKeyGen)/1E9;
		
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
		System.out.println("key generation took " + secondsKeyGen + " seconds");
		System.out.println("decryption took " + secondsDecrypt + " seconds (i = " + i + ")");
		System.out.println(); //more padding
		
		long[] elapsedTimes = new long[4];
		elapsedTimes[0] = elapsedSetup;
		elapsedTimes[1] = elapsedEncrypt;
		elapsedTimes[2] = elapsedKeyGen;
		elapsedTimes[3] = elapsedDecrypt;
		
		return elapsedTimes;
	}
	
	public static void main(String[] args) {
		File lib = new File("../../lib/libmcljava.dylib");
		System.load(lib.getAbsolutePath());
		Mcl.SystemInit(Mcl.BN254);
		testRuntimes(10);
	}

}
