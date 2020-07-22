package schemes;

import java.io.*;
import java.util.concurrent.ThreadLocalRandom;
import java.util.ArrayList;
import java.util.Arrays;

import com.herumi.mcl.*;
import helperclasses.Tools;

//How to use the library: https://github.com/herumi/mcl/blob/master/ffi/java/java.md
//This is the BGW scheme: https://eprint.iacr.org/2005/018.pdf
//Description of BLS: https://crypto.stanford.edu/~dabo/pubs/papers/BLSmultisig.html
//This implements a special case in which A = 1 and B = n
public class BGWSpecialCase {

	//From https://books.google.com/books?id=kb7ZzFrJi48C&pg=PA106&lpg=PA106&dq=calculate+e(g,+g)+pairing&source=bl&ots=Taa2npE3kj&sig=ACfU3U1zNODtV4H_2pUbgscbcQsYljeujg&hl=en&sa=X&ved=2ahUKEwjD_q6hj-LpAhXLTN8KHS3eAU8Q6AEwCXoECAcQAQ#v=onepage&q=calculate%20e(g%2C%20g)%20pairing&f=falsev
	//e(g, g^(ab)) = e(g, g)^(ab) = e(g^a, g^b)
	//so instantiate g as a G1 group and make all the others G2 groups
	
	private static Fr alpha, t;
	private static String randomGenerator = "abc";
	private static G1 g;
	private static G2 gg;
	private static int n;
	private static GT K; //to be precomputed in the setup function
	
	//Input: n = # of receivers
	//Output: Object[]
	//first element:  list containing public key PK
	//second element: list containing private keys: d1, d2, ... dn 
	public static Object[] setup(int n) {
		
		BGWSpecialCase.n = n; //save n so it can be used for other functions
		
		//initialize random generator g in G
		g = new G1();
		Mcl.hashAndMapToG1(g, randomGenerator.getBytes());
		
		//get corresponding element in G2
		gg = new G2();
		Mcl.hashAndMapToG2(gg, randomGenerator.getBytes());
		
		//random alpha and t in Z_p 
		alpha = new Fr();
		alpha.setByCSPRNG(); 
		t = new Fr();
		t.setByCSPRNG();
		
		//precompute K
		precompute();
		
		//Instantiate public key PK
		ArrayList<G2> PK = new ArrayList<G2>();
		PK.add(gg); //add group element g (type G2)
		
		//from i = 1 to i = 2n (ALL)
		for (int i = 1; i <= 2*n; i++) {
			G2 pub = new G2();
			Fr exp = Tools.power(alpha, i); //exp = alpha^i
			Mcl.mul(pub, gg, exp); // g_n = g^(exp)
			PK.add(pub);
		}
		
		//random beta in Z_p
		Fr beta = new Fr();
		beta.setByCSPRNG();
		
		//compute v
		G2 v = new G2();
		Mcl.mul(v, gg, beta); //v = g^beta
		PK.add(v);
		
		//calculate private keys and put them into a list
		ArrayList<G2> privateKeys = new ArrayList<G2>();
		for (int i = 1; i <= n; i++) {
			G2 priv = new G2();
			Mcl.mul(priv, (G2) PK.get(i), beta); // d_i = (g_i)^(beta)
			privateKeys.add(priv);
		}
		//return public key & private key
		Object[] result = new Object[2];
		result[0] = PK;
		result[1] = privateKeys;
		return result;
	}
	
	//precomputes K = e(gn+1, g)^t
	private static void precompute() {
		//calculate e(g, g_(n+1))
		GT e = new GT();
		G2 gNPlus1 = new G2();
		Fr exp = Tools.power(alpha, n+1);
		Mcl.mul(gNPlus1, gg, exp); //gn = g^(alpha^n), so gn+1 = g^(alpha^(n+1)) = g^(exp)
		Mcl.pairing(e, g, gNPlus1);
		K = new GT();
		Mcl.pow(K, e, t); // K = e(g, gn+1)^t
	}
	
	//Input: S = subset to which the message is brodcast, PK = public key
	//Output: Hdr = header (broadcast ciphertext), K = message encryption key
	public static Object[] encrypt(ArrayList<Integer> S, ArrayList<G2> PK) {
			 		
		//calculate C_0 (first element in Hdr)
		G1 c0 = new G1();
		Mcl.mul(c0, g, t); //C_0 = g^t
		
		//calculate C_1 (second element in Hdr)
		G2 v = PK.get(PK.size() - 1);
		G2 product = new G2(v);
		
		for (int i = 0; i < S.size(); i++) {
			int j = S.get(i);
			Mcl.add(product, product, PK.get(n + 1 - j)); //product *= g_(n + 1 - j)
		}
		
		G2 c1 = new G2();
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
	public static GT decrypt(ArrayList<Integer> S, int i, G2 di, Object[] Hdr, ArrayList<Object> PK) {
		
		//Calculate e(g_i, C1) = e(g, C1^(alpha^i))
		GT e1 = new GT();
		G2 c1 = (G2) Hdr[1];
		Fr exp = Tools.power(alpha, i); //compute g^(alpha^i) for all i and include in PK
		G1 gi = new G1();               //g2 multiplication/exponentation is more expensive
		Mcl.mul(gi, g, exp); //gi = g^(alpha^i)
		Mcl.pairing(e1, gi, c1);
		
		//Calculate e(big messy expression, C0) --> see the paper, page 6: https://eprint.iacr.org/2005/018.pdf 
		//In this case, C0 is of type G1, big messy expression of type G2
		
		G1 c0 = (G1) Hdr[0];
		
		//calculate big messy expression
		G2 product = new G2(di);
		
		for (int k = 0; k < S.size(); k++) {	
			int j = S.get(k);
			if (j == i) {
				continue;
			}
			Mcl.add(product, product, PK.get(n + 1 - j + i)); //product *= g_(n + 1 - j + i)
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
	
//TESTING --------------------------------------------------------------------------------------------------------------------------
	
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
			ArrayList<G2> PK = (ArrayList<G2>) setup[0];
			ArrayList<G2> privateKeys = (ArrayList<G2>) setup[1];
			G2 v = PK.get(PK.size() - 1);
			for (int i = 1; i <= n; i++) {
				//Verify for each key: di = vi (v^(alpha^(i)))
				G2 di = privateKeys.get(i - 1); // (i - 1) because private keys starts from i = 1
				Fr exp = Tools.power(alpha, i);
				System.out.println("Exponent = " + exp.toString());
				G2 vi = new G2();
				Mcl.mul(vi, v, exp); //vi = v^(exp)
				assertEquals("i = " + i + " ", di.toString(), vi.toString());
			}
			
			return setup;
	}
	
	//must be called AFTER setup function has been called, otherwise instance variables will not be instantiated
	//Expected: K.toString() == K1.toString()
	public static void testDecrypt(ArrayList<Integer> S, int i, G2 di, Object[] Hdr, ArrayList<Object> PK) {
		GT K1 = decrypt(S, i, di, Hdr, PK);
		System.out.println("K = " + K.toString());
		System.out.println("K1 = " + K1.toString());
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
		Object[] Hdr = (Object[]) encrypted[0];
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
	
	public static void main(String[] args) {

		File lib = new File("../../lib/libmcljava.dylib");
		System.load(lib.getAbsolutePath());
		Mcl.SystemInit(Mcl.BN254); // curveType = Mcl.BN254 or Mcl.BLS12_381
		testRuntimes();
		/*Object[] setup = testSetup(100);
		ArrayList<Object> PK = (ArrayList<G2>) setup[0];
		ArrayList<G2> privateKeys = (ArrayList<G2>) setup[1];
		ArrayList<Integer> S = new ArrayList<Integer>();
		S.addAll(Arrays.asList(1, 6, 17, 25, 33, 49));
		int i = 6;
		G2 di = privateKeys.get(i - 1);
		Object[] encrypted = encrypt(S, PK);
		Object[] Hdr = (Object[]) encrypted[0];
		System.out.println();
		System.out.println("DECRYPT TEST: ");
		testDecrypt(S, i, di, Hdr, PK); //YES! IT WORKED! FINALLY (after so much debugging hahaha)*/
		
	}

}
