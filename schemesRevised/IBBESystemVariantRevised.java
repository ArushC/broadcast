package schemesRevised;
import helperclasses.polynomials.LagrangeInterpolationZp;
import java.io.File;
import java.util.ArrayList;
import helperclasses.miscellaneous.CustomPRF;
import helperclasses.miscellaneous.Tools;
import com.herumi.mcl.*;

//The scheme: https://eprint.iacr.org/2008/268.pdf (4.3.1, page 12)
//mention that K is precomputed in the setup function, key generation is for a single user
//Additional changes: fixed computation of g^P(alpha)
//Included g1 and g1^alpha in the public key (do not need g1^(alpha^2), g1^(alpha^3), ..., etc, only g1 and g1^alpha)
//Also note that this IBBE system is not meant to work when l = 1
public class IBBESystemVariantRevised {

	private static GT K;
	private static int lambda;
	private static G2 g2;
	
	//PRECONDITION: G is of order p >= n + l, l <= n
	//input: n = # of users, l = maximal size of a broadcast recipient group
	//output: Object[] containing the public key PK and private key SK
	public static Object[] setup(int n, int l, int lambda) {
		
		Mcl.SystemInit(lambda);
		IBBESystemVariantRevised.lambda = lambda;
		//instead of creating a GroupGen function, the groups are generated here
		
		//generate random g1, g2
		G1 g1 = new G1();
		Mcl.hashAndMapToG1(g1, "abc".getBytes());
		g2 = new G2();
		Mcl.hashAndMapToG2(g2, "def".getBytes());
		
		//generate random exponents alpha, beta, gamma, t
		Fr alpha = new Fr();
		alpha.setByCSPRNG();
		
		Fr beta = new Fr();
		beta.setByCSPRNG();
		
		Fr gamma = new Fr();
		gamma.setByCSPRNG();
		
		//compute gHat1 = g1^(beta) and gHat2 = g2^(beta)
		G1 gHat1 = new G1();
		Mcl.mul(gHat1, g1, beta);

		G2 gHat2 = new G2();
		Mcl.mul(gHat2, g2, beta);
		
		//precompute K
		precompute(g1, gHat2, alpha, gamma, l);
		
		//add n, l, g1^(gamma), g1^(gamma * alpha) to PK
		ArrayList<Object> PK = new ArrayList<Object>();
		PK.add(n);
		PK.add(l);
		G1 element = new G1();
		Mcl.mul(element, g1, gamma); //add g1^(gamma)
		PK.add(element);
		G1 element2 = new G1();
		Mcl.mul(element2, g1, gamma);
		Mcl.mul(element2, element2, alpha); //add g1^(gamma * alpha)
		PK.add(element2);
		
		PK.add(g1);
		
		G1 element4 = new G1(); //add g1^alpha
		Mcl.mul(element4, g1, alpha);

		//second part of PK is going to be a set containing the set elements: {[gHat1^(alpha^j), gHat2^(alpha^k)}
		//for all j in [0, l] and all k in [0, l - 1]
		ArrayList<Object> setElements = new ArrayList<Object>();
		
		Fr exp = new Fr(1);
		for (int j = 0; j <= l; j++) {
				Object[] elmnt = (j > l - 2) ? new Object[1] : new Object[2];
				G1 e1 = new G1();
				Mcl.mul(e1, gHat1, exp);
				elmnt[0] = e1;
				if (j <= l - 2) {
					G2 e2 = new G2();
					Mcl.mul(e2, gHat2, exp);
					elmnt[1] = e2;
				}
				setElements.add(elmnt);	
				Mcl.mul(exp, exp, alpha);
		}
		
		PK.add(element4);
		PK.add(setElements);
		
		//Generate random key kappa for a pseudorandom function psi
		Fr kappa = new Fr();
		kappa.setByCSPRNG();
		
		//add alpha, gamma, kappa to SK
		Fr[] SK = new Fr[3];
		SK[0] = alpha;
		SK[1] = gamma;
		SK[2] = kappa;
		
		//return (PK, SK)
		Object[] result = new Object[2];
		result[0] = PK;
		result[1] = SK;
		return result;		
	}
	
	//input:  secret key SK and user ID
	//output: ArrayList<Object> di, the user's secret key
	public static Object[] keyGen(Fr ID, Fr[] SK) {
		
		//Extract from SK
		Fr alpha = SK[0];
		Fr gamma = SK[1];
		Fr kappa = SK[2];
		
		//1. compute ri = phi(kappa, i) where phi is a PRF
		CustomPRF phi = new CustomPRF(kappa, ID);
		Fr ri = phi.compute(lambda);
		
		//2. compute hi = g2^((gamma - ri)/(alpha - i))
		G2 hi = new G2();
		Fr exp = new Fr();
		Fr expNum = new Fr();
		Fr expDen = new Fr();
		Mcl.sub(expNum, gamma, ri);
		Mcl.sub(expDen, alpha, ID);
		Mcl.div(exp, expNum, expDen);
		Mcl.mul(hi, g2, exp);
		
		//return (ri, hi)
		Object[] di = {ri, hi};
		return di;
		
	}
	
	
	public static Object[] enc(ArrayList<Fr> S, ArrayList<Object> PK) {
		
		//extract l and n, compute k
		int n = (int) PK.get(0);
		int l = (int) PK.get(1);
		
		ArrayList<Object> setElements = (ArrayList<Object>) PK.get(6);
		Object[] setElement = (Object[]) setElements.get(0);
		G1 gHat1 = (G1) setElement[0];
		
		//compute Px
		Fr[] Px = computePx(n, l, new Fr(-1), S);
		
		Fr t = new Fr();
		t.setByCSPRNG();
		Mcl.pow(K, K, t);
		
		//compute C1, C2, C3, C4
		G1 C1 = new G1(gHat1);
		Mcl.mul(C1, C1, Px[Px.length - 1]);
		//calculate gHat1^(P(alpha))
		int index = 1;
		for (int i = Px.length - 1; i > 0; i--) {
			G1 addend = new G1();
			Mcl.mul(addend, (G1) (((Object[]) setElements.get(index))[0]) , Px[i - 1]);
			Mcl.add(C1, C1, addend);
			index += 1;
		}
		
		Mcl.mul(C1, C1, t);
		
		G1 C2 = new G1();
		Mcl.mul(C2, (G1) PK.get(2), t);
		
		G1 C3 = new G1();
		Mcl.mul(C3, (G1) PK.get(4), t);
		
		GT C4 = new GT();
		Mcl.pairing(C4, (G1) PK.get(5), (G2) (((Object[]) setElements.get(l - 2))[1]));
		Mcl.pow(C4, C4, t);
		
		//add to Hdr and result
		Object[] Hdr = {C1, C2, C3, C4};
		Object[] result = {Hdr, K};
		return result;
		
	}
	
	
	//returns the key K1 of type GT
	public static GT decrypt(ArrayList<Fr> S, Fr ID, Object[] di, Object[] Hdr, ArrayList<Object> PK) {
		
		//1. extract from PK and compute P(x)
		int n = (int) PK.get(0);
		int l = (int) PK.get(1);
		Fr[] Px = computePx(n, l, ID, S); //P(x)/(x - i)
		//2. extract from Hdr and di
		G1 C1 = (G1) Hdr[0];
		G1 C2 = (G1) Hdr[1];
		G1 C3 = (G1) Hdr[2];
		GT C4 = (GT) Hdr[3];
		Fr ri = (Fr) di[0];
		G2 hi = (G2) di[1];
		
		//3. compute gHat2^(Pi(alpha))
		ArrayList<Object> setElements = (ArrayList<Object>) PK.get(6);
		Fr[] negated = LagrangeInterpolationZp.negate(Px); //note: need to ignore item at position 0 in this array because this is an l - 2 degree polyn.
		Object[] setElement = (Object[]) setElements.get(0);
		G2 gHat2 = (G2) setElement[1]; //if l = 1, then gHat2 does not need to be extracted
		G2 fin = new G2(gHat2);
		Mcl.mul(fin, fin, negated[negated.length - 1]);
		int index = 1;
		for (int j = negated.length - 1; j > 1; j--) {
			G2 addend = new G2();
			Mcl.mul(addend, (G2) (((Object[]) setElements.get(index))[1]) , negated[j - 1]);
			Mcl.add(fin, fin, addend);
			index += 1;
		}
		
		Fr ei = new Fr();
		Mcl.mul(ei, ri, new Fr(-1));
		
		//compute the pairings
		GT e1 = new GT();
		Mcl.pairing(e1, C1, hi);
		
		GT e2 = new GT();
		G1 firstE2 = new G1();
		Mcl.mul(firstE2, C3, ei);
		Mcl.add(firstE2, C2, firstE2);
		Mcl.pairing(e2, firstE2, fin);
		
		GT denominator = new GT();
		Mcl.pow(denominator, C4, ei);
		Mcl.pow(denominator, denominator, new Fr(-1));
		
		//multiply all together to get the result
		GT result = new GT();
		Mcl.mul(result, e1, e2);
		Mcl.mul(result, result, denominator);
		
		return result;
		
	}

	//HELPER FUNCTIONS -----------------------------------------------------------------------------------------------------------------------------
	
	//computes Px as defined in the tagEncrypt function
	private static Fr[] computePx(int n, int l, Fr i, ArrayList<Fr> S) {
		
		Fr NEGATIVEONE = new Fr(-1);
		int k = S.size();
		ArrayList<Fr> ijValues = new ArrayList<Fr>();
		for (int j = k + 1; j <= l; j++) {
			ijValues.add(new Fr(n + j));
		}
		
		ArrayList<Fr> allIJValues = new ArrayList<Fr>();
		allIJValues.addAll(S);
		allIJValues.addAll(ijValues);
		
		//compute P(x) (or Pi(x), if i is in allIJValues)
		Fr[] Px = {new Fr(1)};
		for (int j = 1; j <= l; j++) {
			Fr item = allIJValues.get(j - 1);
			if (item.equals(i))
				continue;
			//else
			Fr secondElement = new Fr();
			Mcl.mul(secondElement, item, NEGATIVEONE);
			Fr[] next = {new Fr(1), secondElement};
			Px = LagrangeInterpolationZp.multiply(Px, next);
		}
		
		return Px;	
	}
	
	
	//precomputes K = e(gHat2, g1)^(alpha * gamma^(l - 1) * t)
	private static void precompute(G1 g1, G2 gHat2, Fr alpha, Fr gamma, int l) {

		K = new GT();
		Mcl.pairing(K, g1, gHat2);
		Fr exp = new Fr();
		Fr alphaExp = Tools.power(alpha, l - 1);
		Mcl.mul(exp, gamma, alphaExp);
		Mcl.pow(K, K, exp);
		
	}
	
	private static long[] printRuntimes(int n, int l, int lambda) {
		
		long startSetup = System.nanoTime();
		Object[] setup = setup(n, l, lambda);
		long elapsedSetup = System.nanoTime() - startSetup;
		double secondsSetup = ((double) elapsedSetup) / 1E9;
		
		ArrayList<Object> PK = (ArrayList<Object>) setup[0];
		Fr[] MSK = (Fr[]) setup[1];
		
		//random user ID to test decryption
		Fr ID = new Fr();
		ID.setByCSPRNG();
				
		//generate subset S
		ArrayList<Fr> S = new ArrayList<Fr>();
		for (int i = 0; i < l - 1; i++) {
			Fr d = new Fr();
			d.setByCSPRNG();
			S.add(d);
		}
				
		S.add(ID);
		
		//generate ciphertext and keys
		long startEncrypt = System.nanoTime();
		Object[] C = enc(S, PK);
		long elapsedEncrypt = System.nanoTime() - startEncrypt;
		double secondsEncrypt = ((double) elapsedEncrypt) / 1E9;
		long startKeyGen = System.nanoTime();
		Object[] di = keyGen(ID, MSK);
		long elapsedKeyGen = System.nanoTime() - startKeyGen;
		double secondsKeyGen = ((double) elapsedKeyGen) / 1E9;
		
		//finally, decrypt
		Object[] Hdr = (Object[]) C[0];
		GT K1 = (GT) C[1];
		long startDecrypt = System.nanoTime();
		GT K2 = decrypt(S, ID, di, Hdr, PK);
		long elapsedDecrypt = System.nanoTime() - startDecrypt;
		double secondsDecrypt = ((double) elapsedDecrypt) / 1E9;
		
		String result = (K1.equals(K2)) ? "SUCCESSFUL DECRYPTION" : "FAILED DECRYPTION";
		System.out.println(result + ": n = " + n + ", l = " + l);
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
	

	public static void testRuntimes(int lambda, int percent) {
		for (int N = 100; N <= 1000000; N *= 10) {
			int subsetSize = (int) (0.01 * percent * N);
			printRuntimes(N, subsetSize, lambda);
		}
	}	
	
	public static void main(String[] args) {
		File lib = new File("../../lib/libmcljava.dylib");
		System.load(lib.getAbsolutePath());
		testRuntimes(Mcl.BN254, 10);
	}
	
}
