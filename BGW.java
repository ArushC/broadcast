import java.io.*;
import java.math.BigInteger;
import java.util.ArrayList;
import java.math.*;

import com.herumi.mcl.*;

//How to use the library: https://github.com/herumi/mcl/blob/master/ffi/java/java.md
//This is the BGW scheme: https://eprint.iacr.org/2005/018.pdf
//Description of BLS: https://crypto.stanford.edu/~dabo/pubs/papers/BLSmultisig.html

public class BGW {

	private static Fr alpha;
	private static boolean allTestsPassed = false;
	
	//Input: n = # of receivers
	//Output: Object[]
	//first element:  list containing public key PK
	//second element: list containing private keys: d1, d2, ... dn 
	public static Object[] setup(int n) {
		
		//get random generator g in G
		G2 g = new G2();
		Mcl.hashAndMapToG2(g, "abc".getBytes());
		
		//random alpha in Z_p
		alpha = new Fr();
		alpha.setByCSPRNG(); 
		
		
		//Instantiate public key PK
		ArrayList<G2> PK = new ArrayList<G2>();
		BigInteger a_value = new BigInteger(alpha.toString());
		PK.add(g); //add group g
		
		//from i = 1 to i = n
		for (int i = 1; i < n; i++) {
			G2 pub = new G2();
			Fr exp = new Fr();
			Mcl.mul(exp, alpha, new Fr(i)); //exp = alpha^(i) 
			//System.out.println(exp);
			Mcl.mul(pub, g, exp); // g_n = g^(exp)
			//System.out.println(pub);
			PK.add(pub);
		}
		//from i = n to i = 2n (step = 2)
		for (int i = n; i <= 2*n; i+=2) {
			G2 pub = new G2();
			Fr exp = new Fr();
			Mcl.mul(exp, alpha, new Fr(i)); //exp = alpha^(i) 
			Mcl.mul(pub, g, exp);
			PK.add(pub);
		}
		
		//random beta in Z_p
		Fr beta = new Fr();
		beta.setByCSPRNG();
		//compute v
		G2 v = new G2();
		Mcl.mul(v, g, beta); //v = g^beta
		PK.add(v);
		
		//calculate private keys and put them into a list
		ArrayList<G2> privateKeys = new ArrayList<G2>();
		for (int i = 1; i <= n; i++) {
			G2 priv = new G2();
			Mcl.mul(priv, PK.get(i), beta); // d_i = (g_i)^(beta)
			privateKeys.add(priv);
		}
		//return public key & private key
		Object[] result = new Object[2];
		result[0] = PK;
		result[1] = privateKeys;
		return result;
	}
		
	//Input: S = subset to which the message is brodcast, PK = public key
	//Output: Hdr = header (broadcast ciphertext), K = message encryption key
	//public static int[] encrypt(ArrayList<Integer> S, int PK) {
			
	//}
	
	//Input: S = subset, i = user id, di = user private key, Hdr = header, PK = public key
	//Output: if i is in S, output message encryption key K. Use to decode C (brodcast body)
	//public static int[] decrypt(ArrayList<Integer> S, int i, int di, int Hdr, int PK) {
		
	//}
	
	//TESTING -------------------------------------------------------------------------------
	
	public static void assertEquals(String msg, String x, String y) {
		if (x.equals(y)) {
			System.out.println("OK : " + msg);
		} 
		else {
			System.out.println("NG : " + msg + ", x = " + x + ", y = " + y);
			allTestsPassed = false;
		}
	}
	
	public static void testSetup(int n) {
		//TEST SETUP
			Object[] setup = setup(n); //102 elements (public keys + g_0 + v)
			//EXTRACT EVERYTHING
			ArrayList<G2> PK = (ArrayList<G2>) setup[0];
			ArrayList<G2> privateKeys = (ArrayList<G2>) setup[1];
			G2 v = PK.get(PK.size() - 1);
			for (int i = 1; i <= n; i++) {
				//Verify for each key: di = vi (v^(alpha^(i)))
				G2 di = privateKeys.get(i - 1); // (i - 1) because private keys starts from i = 1
				Fr exp = new Fr();
				Mcl.mul(exp, alpha, new Fr(i)); //exp = (alpha^(i))
				System.out.println("Exponent = " + exp.toString());
				G2 vi = new G2();
				Mcl.mul(vi, v, exp); //vi = v^(exp)
				assertEquals("i = " + i + " ", di.toString(), vi.toString());
			}
	}
	
	
	
	public static void main(String[] args) {
		//change the file directory here
		System.load("/Users/arushchhatrapati/Documents/mcl/lib/libmcljava.dylib");
		Mcl.SystemInit(Mcl.BN254); // curveType = Mcl.BN254 or Mcl.BLS12_381
		testSetup(100); //SUCCESS!!
		
	}

}
