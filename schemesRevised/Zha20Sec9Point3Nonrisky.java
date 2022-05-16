package schemesRevised;
import helperclasses.structures.VectorND;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.ThreadLocalRandom;
import com.herumi.mcl.Fr;
import com.herumi.mcl.GT;
import com.herumi.mcl.Mcl;

//nonrisky broadcast and trace scheme
//Zhandry's O(N^(1 - a), N^(1 - a), N^a) scheme -- setting a = 2/3 yields a scheme with public/private key of N^(1/3)
//The scheme: https://eprint.iacr.org/2020/954, section 9.3, post risk-mitigation compiler & user-expansion compiler
public class Zha20Sec9Point3Nonrisky {
	
	private static int nOverT, MT, T, N;
	
	public static Object[] setup(int N, int M, int T, int lambda) {
		
		//set up N/T instances of an (N/T) * (MT) broadcast multischeme
		int nOverT = (int) Math.ceil(((double) N)/T);
		
		//save to be used in later functions
		Zha20Sec9Point3Nonrisky.nOverT = nOverT;
		Zha20Sec9Point3Nonrisky.MT = M * T;
		Zha20Sec9Point3Nonrisky.T = T;
		Zha20Sec9Point3Nonrisky.N = N;
		
		int n = 2;
		int v = nOverT;
		int u = 1, t = 1;
		Object[] setup = Zha20Sec9Point3Risky.genMTB(u, v, t, n, lambda);
		ArrayList<Object> PK = (ArrayList<Object>) setup[0];
		Object[] MSK = (Object[]) setup[1];
		ArrayList<Object> secretKeys = new ArrayList<Object>();
		
		//initialize the tracing key
		ArrayList<Object> tk = new ArrayList<Object>();
		tk.add(MSK);
		
		//for each instance, choose random iStarJ in {1, 2, ..., N/T}
		for (int j = 1; j <= nOverT * MT; j++) {
			int iStarJ = ThreadLocalRandom.current().nextInt(1, nOverT + 1);
			tk.add(iStarJ);
		}
		
		Object[] result = {PK, tk};
		return result;

	}
		
	//setup authority generates the key for an individual user -- note that secret key has size O(N/T)
	public static ArrayList<Object> keyGen(int ID, ArrayList<Object> tk) {
			
			//extract from ID and tracing key
			Object[] MSK = (Object[]) tk.get(0);
			ArrayList<Object> skID = new ArrayList<Object>();
			int i = (ID % nOverT == 0) ? nOverT : (ID % nOverT);
			int j = (int) Math.ceil(((double) ID) / nOverT);
			
			//X = 1, 2, ..., N/T instances of the broadcast multischeme
			for (int X = 1; X <= nOverT; X++) {

				int newJ = j + MT * (X - 1); //calculate j and ID for the Xth instance
				int newID = (newJ - 1) * nOverT + i;
				
				int iStarJ = (int) tk.get(newJ);
				
				//calculate wJI value
				VectorND wJI;
				if (i < iStarJ) {
					wJI = new VectorND(new Fr(1), new Fr(1));
				}
				else if (i == iStarJ) {
					wJI = new VectorND(new Fr(1), new Fr(0));
				}
				else {
					wJI = new VectorND(new Fr(0), new Fr(0));
				}
				
				//finally, generate the key for the Xth instance
				Set<Integer> U = new HashSet<Integer>();
				U.add(newID);
				ArrayList<Object> ski = Zha20Sec9Point3Risky.extractMTB(MSK, U, wJI);
				skID.add(ski);

			}
			return skID;
	}
	
	public static Object[] enc(ArrayList<Object> PK, int k, Set<Integer> S) {
		
		ArrayList<Object> ciphertext = new ArrayList<Object>();
		
		//select a random instance of the multischeme
		int X = ThreadLocalRandom.current().nextInt(1, nOverT + 1);
		int startJ = k * T - (T - 1);
		startJ += MT * (X - 1);
		
		//compute the subsets of the S subset (partition based on the T instances)
		ArrayList<HashSet<Integer>> sSubsets = new ArrayList<HashSet<Integer>>();
		for (int x = 0; x < T; x++)
			sSubsets.add(new HashSet<Integer>()); //initialize the list containing the subsets
		for (int s: S) {
			int i = (s % nOverT == 0) ? nOverT : (s % nOverT);
			int j = (int) Math.ceil(((double) s)/nOverT);
			sSubsets.get(j - 1).add(i);
		}
		
		
		for (int l = startJ; l < startJ + T; l++) {
			
			Set<Integer> TjS = new HashSet<Integer>();
			for (int i: sSubsets.get((l - 1) % T)) {
				int ID = (l - 1) * nOverT + i;
				TjS.add(ID);
			}
			ciphertext.add(Zha20Sec9Point3Risky.encMTB(PK, TjS));
		}
			
		Object[] res = {ciphertext, X};
		return res;
		
	}
	
	//decrypt ciphertext c encrypted to instance X for a user with secret key skID 
	public static GT decrypt(int userID, ArrayList<Object> skID, ArrayList<Object> helperKey, Set<Integer> S, int X, Object[] c) {
		
		//extract j, i, and the secret key for this instance
		int jID = (int) Math.ceil(((double) userID) / nOverT);
		int iID = (userID % nOverT == 0) ? nOverT : (userID % nOverT);
		ArrayList<Object> ski = (ArrayList<Object>) skID.get(X - 1);
		
		//calculate the subset of the given instance
		HashSet<Integer> sSubset = new HashSet<Integer>();
		for (int i = 1; i <= nOverT; i++)
			sSubset.add((((jID - 1) * nOverT + i) % N) == 0 ? N : ((jID - 1) * nOverT + i) % N);
		
		sSubset.retainAll(S); //contains IDs of users in the given instance instead of i-values (faster)
		
		int jOfXInstance = jID + MT * (X - 1);
		int IDXInstance = (jOfXInstance - 1) * nOverT + iID;
		
		//calculate TjS
		Set<Integer> TjS = new HashSet<Integer>();
		for (int s: sSubset) {
			int i = (s % nOverT == 0) ? nOverT : (s % nOverT); //extract i-values when calculating TjS
			int ID = (jOfXInstance - 1) * nOverT + i;
			TjS.add(ID);
		}

		Set<Integer> U = new HashSet<Integer>();
		U.add(IDXInstance);	
		return Zha20Sec9Point3Risky.decMTB(ski, helperKey, TjS, U, c);
	}
	
	//To get the N^(1/3) scheme, set a = 2/3. In general, we get a scheme of size (N^(1 - a), N^(1 - a), N^a)
	public static void printRuntimes(int N, int M, double a, int subsetSize, int lambda) {
		
		int T = (int) Math.floor(Math.pow(N, a));
		long startSetup = System.nanoTime();
		Object[] setup = Zha20Sec9Point3Nonrisky.setup(N, M, T, Mcl.BN254);
		long elapsedSetup = System.nanoTime() - startSetup;
		double secondsSetup = ((double) elapsedSetup)/1E9;
		
		//random ID for the encryption & decryption
		int ID = ThreadLocalRandom.current().nextInt(1, N + 1);
		int k = (ID % N == 0) ? ID/N : ID/N + 1;
		int j = (int) Math.ceil(((double) ID) / nOverT);
		int i = (ID % nOverT == 0) ? nOverT : (ID % nOverT);
		
		Set<Integer> S = new HashSet<Integer>();
		S.add(ID);
		
		//randomly generate numbers to put in the subset S (NO REPEATS)
		ArrayList<Integer> randomNums = new ArrayList<Integer>();
		int randomID = 0;
		for (int z = 0; z < N; z++) {
			if (z + 1 == ID) continue; //already in the subset
			randomNums.add(z + 1);
		}
		for (int z = 1; z < subsetSize; z++) {
			int randomIndex = ThreadLocalRandom.current().nextInt(0, randomNums.size());
			randomID = randomNums.get(randomIndex);
			S.add(randomID);
			randomNums.remove(randomIndex);
		}
		
		//extract public/tracing key and helper key
		ArrayList<Object> PK = (ArrayList<Object>) setup[0];
		ArrayList<Object> tk = (ArrayList<Object>) setup[1];
		ArrayList<Object> helperKey = (ArrayList<Object>) PK.get(PK.size() - 1);
		
		long startEnc = System.nanoTime();
		Object[] C = Zha20Sec9Point3Nonrisky.enc(PK, k, S);
		long elapsedEnc = System.nanoTime() - startEnc;
		double secondsEnc = ((double) elapsedEnc)/1E9;
		
		//extract from ciphertext
		int X = (int) C[1];
		ArrayList<Object> ciphertext = (ArrayList<Object>) C[0];
		int ciphertextIndex = (j % T == 0) ? (T - 1) : (j % T - 1);
		Object[] enc = (Object[]) ciphertext.get(ciphertextIndex);
		GT encapsulatedKey = (GT) enc[1];
		Object[] c = (Object[]) enc[0];
		
		long startKeyGen = System.nanoTime();
		ArrayList<Object> skID = Zha20Sec9Point3Nonrisky.keyGen(ID, tk);
		long elapsedKeyGen = System.nanoTime() - startKeyGen;
		double secondsKeyGen = ((double) elapsedKeyGen)/1E9;
		
		//Finally, decrypt
		long startDec = System.nanoTime();
		GT encapsulatedKey1 = Zha20Sec9Point3Nonrisky.decrypt(ID, skID, helperKey, S, X, c);
		long elapsedDec = System.nanoTime() - startDec;
		double secondsDec = ((double) elapsedDec)/1E9;
		
		String success = (encapsulatedKey.equals(encapsulatedKey1)) ? "SUCCESSFUL DECRYPTION (" + (N * M) + " total users)": "FAILED DECRYPTION";
		System.out.println(success + ": " + "N = " + N + ", M = " + M + ", T = " + T + ", subset size = " + subsetSize);
		System.out.println(); //padding
		System.out.println("setup took " + secondsSetup + " seconds");
		System.out.println("encryption took " + secondsEnc + " seconds");
		System.out.println("key generation took " + secondsKeyGen + " seconds");
		System.out.println("decryption took " + secondsDec + " seconds (u = " + ID + ")");
		System.out.println(); //more padding
		
	}
	
	//test runtimes for M = 1 instances and a given value of a
	public static void testRuntimes(double a, int percent, int lambda) {
		for (int N = 100; N <= 1000000; N *= 10) {
			int subsetSize = (int) (0.01 * percent * N);
			printRuntimes(N, 1, a, subsetSize, lambda);
		}
	}
	
	public static void main(String[] args) {
		File lib = new File("../../lib/libmcljava.dylib");
		System.load(lib.getAbsolutePath());
		double a = 2.0/3;
		testRuntimes(a, 10, Mcl.BN254);
	}
	
}
