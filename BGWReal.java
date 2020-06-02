
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

import com.herumi.mcl.*;

//How to use the library: https://github.com/herumi/mcl/blob/master/ffi/java/java.md
//This is the BGW scheme: https://eprint.iacr.org/2005/018.pdf
//Description of BLS: https://crypto.stanford.edu/~dabo/pubs/papers/BLSmultisig.html

public class BGWReal {

	
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
		
		BGWReal.n = n; //save n so it can be used for other functions
		
		//initialize random generator g in G
		g = new G1();
		Mcl.hashAndMapToG1(g, randomGenerator.getBytes());
		
		//get corresponding element in G2
		gg = new G2();
		Mcl.hashAndMapToG2(gg, randomGenerator.getBytes());
		
		//random alpha in Z_p 
		alpha = new Fr();
		alpha.setByCSPRNG(); 
		t = new Fr();
		t.setByCSPRNG();
		
		//precompute K
		precompute();
		
		//Instantiate public key PK
		ArrayList<Object> PK = new ArrayList<Object>();
		PK.add(g); //add group element g (type G1)
		
		//from i = 1 to i = n
		for (int i = 1; i < n; i++) {
			G2 pub = new G2();
			Fr exp = new Fr();
			Mcl.mul(exp, alpha, new Fr(i)); //exp = alpha^(i) 
			//System.out.println(exp);
			Mcl.mul(pub, gg, exp); // g_n = g^(exp)
			//System.out.println(pub);
			PK.add(pub);
		}
		//from i = n to i = 2n (step = 2)
		for (int i = n; i <= 2*n; i+=2) {
			G2 pub = new G2();
			Fr exp = new Fr();
			Mcl.mul(exp, alpha, new Fr(i)); //exp = alpha^(i) 
			Mcl.mul(pub, gg, exp);
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
		Fr exp = new Fr();
		Mcl.mul(exp, alpha, new Fr(n + 1));
		Mcl.mul(gNPlus1, gg, exp); //gn = g^(alpha^n), so gn+1 = g^(alpha^(n+1)) = g^(exp)
		//will this work?
		Mcl.pairing(e, g, gNPlus1);
		K = new GT();
		Mcl.pow(K, e, t); // K = e(g, gn+1)^t
	}
	
	//Input: S = subset to which the message is brodcast, PK = public key
	//Output: Hdr = header (broadcast ciphertext), K = message encryption key
	public static Object[] encrypt(ArrayList<Integer> S, ArrayList<Object> PK) {
			 		
		//calculate C_0 (first element in Hdr)
		G1 c0 = new G1();
		Mcl.mul(c0, g, t); //C_0 = g^t
		
		//calculate C_1 (second element in Hdr)
		G2 v = (G2) PK.get(PK.size() - 1);
		
		int initialIndex = n + 1 - S.get(0); //S.get(0) = j
		G2 product = new G2((G2) PK.get(initialIndex));
		
		for (int i = 1; i < S.size(); i++) {
			int j = S.get(i);
			Mcl.add(product, product, (G2) PK.get(n + 1 - j)); //product *= g_(n + 1 - j)
		}
		
		G2 c1 = new G2();
		Mcl.add(c1, v, product); // c1 = v * product
		Mcl.mul(c1, c1, t);      //c1 = c1^(t)
		
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
		Fr exp = new Fr();
		Mcl.mul(exp, alpha, new Fr(i));
		Mcl.pairing(e1, g, c1);
		
		//Calculate e(big messy expression, C0) --> see the paper, page 6: https://eprint.iacr.org/2005/018.pdf 
		//In this case, C0 is of type G1, big messy expression of type G2
		
		G1 c0 = (G1) Hdr[0];
		
		//calculate big messy expression
		int initialIndex = n + 1 - S.get(0) + i; //S.get(0) = j
		G2 product = new G2((G2) PK.get(initialIndex));
		
		
		for (int k = 1; k < S.size(); k++) {
			if (k == i) {
				continue;
			}
			//else	
			int j = S.get(k);
			Mcl.add(product, product, (G2) PK.get(n + 1 - j + i)); //product *= g_(n + 1 - j + i)
		}
		
		G2 result = new G2();
		Mcl.add(result, di, product); //result = di * product
		
		//finally, compute the pairing
		GT e2 = new GT();
		Mcl.pairing(e2, c0, result);
		
		GT K = new GT();
		//to compute e1 / e2, calculate e2^(-1) and output e1 * e2^(-1)
		Mcl.pow(e2, e2, new Fr(-1)); //CONFIRMED this works (after testing)
		Mcl.mul(K, e1, e2);
		
		return K;
		
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
			
			return setup;
	}
	
	//must be called AFTER setup function has been called, otherwise instance variables will not be instantiated
	public static void testDecrypt(ArrayList<Integer> S, int i, G2 di, Object[] Hdr, ArrayList<Object> PK) {
		GT K1 = decrypt(S, i, di, Hdr, PK);
		System.out.println("K = " + K.toString());
		System.out.println("K1 = " + K1.toString());
	}

	
	
	public static void main(String[] args) {
		//change the file directory here
		System.load("/Users/arushchhatrapati/Documents/mcl/lib/libmcljava.dylib");
		Mcl.SystemInit(Mcl.BN254); // curveType = Mcl.BN254 or Mcl.BLS12_381
		Object[] setup = testSetup(100);
		ArrayList<Object> PK = (ArrayList<Object>) setup[0];
		ArrayList<G2> privateKeys = (ArrayList<G2>) setup[1];
		ArrayList<Integer> S = new ArrayList<Integer>();
		S.addAll(Arrays.asList(1, 6, 17, 25, 33));
		int i = 61;
		G2 di = privateKeys.get(i - 1);
		Object[] encrypted = encrypt(S, PK);
		Object[] Hdr = (Object[]) encrypted[0];
		testDecrypt(S, i, di, Hdr, PK); //Unfortunately, this doesn't work (K should be equal to K1)
										//Gotta do a lot of debugging now...
		
		
	}

}
