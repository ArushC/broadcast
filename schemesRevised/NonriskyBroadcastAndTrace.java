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

public class NonriskyBroadcastAndTrace {
	
	private static int nOverT, MT, T, N;
	
	public static Object[] setup(int N, int M, int T, int lambda) {
		
		//set up N/T instances of an (N/T) * (MT) broadcast multischeme
		int nOverT = N/T;
		
		//save to be used in later functions
		NonriskyBroadcastAndTrace.nOverT = nOverT;
		NonriskyBroadcastAndTrace.MT = M * T;
		NonriskyBroadcastAndTrace.T = T;
		NonriskyBroadcastAndTrace.N = N;
		
		int n = 2;
		int v = nOverT;
		int u = 1, t = 1;
		Object[] setup = RiskyMTBRevised.genMTB(u, v, t, n, lambda);
		ArrayList<Object> PK = (ArrayList<Object>) setup[0];
		Object[] MSK = (Object[]) setup[1];
		ArrayList<Object> secretKeys = new ArrayList<Object>();
		
		//initialize the tracing key
		ArrayList<Object> tk = new ArrayList<Object>();
		tk.add(MSK);
		
		//for each instance, choose random iStarJ in {1, 2, ..., N/T}
		for (int j = 1; j <= N * M; j++) {
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
				
				
				int newJ = j + MT * (X - 1);
				int newID = (newJ - 1) * nOverT + i;
				System.out.println("New ID for instance X = " + X + ": " + newID);
				
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
				ArrayList<Object> ski = RiskyMTBRevised.extractMTB(MSK, U, wJI);
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
			System.out.println("TjS, l = " + l + ": " + TjS);
			ciphertext.add(RiskyMTBRevised.encMTB(PK, TjS));
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
			sSubset.add(((jID - 1) * nOverT + i) % N == 0 ? N : ((jID - 1) * nOverT + i));
		
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
			
		System.out.println("TjS decryption, j = " + jOfXInstance + ": " + TjS);
		Set<Integer> U = new HashSet<Integer>();
		U.add(IDXInstance);	
		return RiskyMTBRevised.decMTB(ski, helperKey, TjS, U, c);
	}
	
	public static void main(String[] args) {
		File lib = new File("../../lib/libmcljava.dylib");
		System.load(lib.getAbsolutePath());
		int N = 27;
		int M = 3;
		int T = (int) Math.round(Math.pow(N, 2.0/3));
		Object[] setup = NonriskyBroadcastAndTrace.setup(N, M, T, Mcl.BN254);
		int ID = 26;
		int k = (ID == N) ? ID/N : ID/N + 1;
		int j = (int) Math.ceil(((double) ID) / nOverT); //FIX N OVER T
		int i = (ID % nOverT == 0) ? nOverT : (ID % nOverT);
		System.out.println("k = " + k + ", j = " + j + ", i = " + i);
		Set<Integer> S = new HashSet<Integer>();
		S.addAll(Arrays.asList(1, 3, 19, 26, 27));
		ArrayList<Object> PK = (ArrayList<Object>) setup[0];
		ArrayList<Object> tk = (ArrayList<Object>) setup[1];
		Object[] C = NonriskyBroadcastAndTrace.enc(PK, k, S);
		int X = (int) C[1];
		System.out.println("X = " + X);
		ArrayList<Object> ciphertext = (ArrayList<Object>) C[0];
		int ciphertextIndex = (j % T == 0) ? (T - 1) : (j % T - 1);
		Object[] enc = (Object[]) ciphertext.get(ciphertextIndex);
		GT encapsulatedKey = (GT) enc[1];
		Object[] c = (Object[]) enc[0];
		System.out.println("Encapsulated Key: " + encapsulatedKey);
		ArrayList<Object> skID = NonriskyBroadcastAndTrace.keyGen(ID, tk);
		ArrayList<Object> helperKey = (ArrayList<Object>) PK.get(PK.size() - 1);
		GT encapsulatedKey1 = NonriskyBroadcastAndTrace.decrypt(ID, skID, helperKey, S, X, c);
		System.out.println("Encapsulated Key 1: " + encapsulatedKey1);
		
	}
	
}
