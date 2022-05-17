package schemes;
import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.ThreadLocalRandom;
import helperclasses.miscellaneous.Tools;
import com.herumi.mcl.*;

//Implementation of the baseline scheme over the group G1 used for the pairing-based schemes
public class BaselineOverPairings {

	private static int N;
	
	//output public key and private key
	public static Object[] setup(int lambda, int N) {
		
		Mcl.SystemInit(lambda);
		BaselineOverPairings.N = N;
		
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
		try {
			Mcl.sub(res, zi, res);
		}
		catch (NullPointerException e) {
			return null; //if cannot be decrypted, return null
		}
			
		return res;	
	}
	
	//returns long[] containing elapsed time for setup, encrypt, and decrypt, respectively
		private static long[] printRuntimes(int lambda, int N, int subsetSize) {
			
			//Get elapsed time for setup(n) 
			long startSetup = System.nanoTime();
			Object[] setup = setup(lambda, N);
			long elapsedSetup = System.nanoTime() - startSetup;
			double secondsSetup = ((double) elapsedSetup) / 1E9;
			//extract public/private key from setup
			ArrayList<G1> PK = (ArrayList<G1>) setup[0];
			ArrayList<Fr> privateKeys = (ArrayList<Fr>) setup[1];
			
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
			
			//Generate random message M
			G1 M = new G1();
			Mcl.hashAndMapToG1(M, Tools.generateRandomBytes(3));
			
			//Get elapsed time for encrypt
			long startEncrypt = System.nanoTime();
			G1[] CT = enc(PK, S, M);
			long elapsedEncrypt = System.nanoTime() - startEncrypt;
			double secondsEncrypt = ((double) elapsedEncrypt) / 1E9;
			
			//Get random user ID i to test the decryption
			int i = S.get(ThreadLocalRandom.current().nextInt(0, S.size()));
			Fr xi = privateKeys.get(i - 1);
			
			//Get elapsed time for decrypt
			long startDecrypt = System.nanoTime();
			G1 M1 = dec(CT, xi, i);
			long elapsedDecrypt = System.nanoTime() - startDecrypt;
			double secondsDecrypt = ((double) elapsedDecrypt) / 1E9;
			//Finally, print out the results
			
			String success = (M.equals(M1)) ? "SUCCESSFUL DECRYPTION" : "FAILED DECRYPTION";
			System.out.println(success + ": " + "n = " + N + ", subset size = " + subsetSize);
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

