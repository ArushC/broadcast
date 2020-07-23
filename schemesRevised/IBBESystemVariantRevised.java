package schemesRevised;
import com.herumi.mcl.*;

import helperclasses.LagrangeInterpolationZp;
import helperclasses.CustomPRF;
import helperclasses.Tools;

import java.io.*;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ThreadLocalRandom;

//The scheme: https://eprint.iacr.org/2008/268.pdf (4.3.1, page 12)
//mention that K is precomputed in the setup function
//Additional changes: fixed computation of g^P(alpha)
//Added e(g1, gHat2)^(alpha^(l - 1)) to the public key
public class IBBESystemVariantRevised {

	private static CustomPRF phi;
	private static Fr gamma, t;
	private static GT K;
	private static Fr alpha;
	private static int lambda;
	private static G1 g1;
	private static G2 g2;
	
	//PRECONDITION: G is of order p >= n + l, l <= n
	//input: n = # of users, l = maximal size of a broadcast recipient group
	//output: Object[] containing the public key PK and private key SK
	public static Object[] setup(int n, int l, int lambda) {
		
		Mcl.SystemInit(lambda);
		IBBESystemVariantRevised.lambda = lambda;
		//instead of creating a GroupGen function, the groups are generated here
		
		//generate random g1, g2
		g1 = new G1();
		Mcl.hashAndMapToG1(g1, "abc".getBytes());
		g2 = new G2();
		Mcl.hashAndMapToG2(g2, "def".getBytes());
		
		//generate random exponents alpha, beta, gamma, t
		alpha = new Fr();
		alpha.setByCSPRNG();
		
		Fr beta = new Fr();
		beta.setByCSPRNG();
		
		gamma = new Fr();
		gamma.setByCSPRNG();
		
		t = new Fr();
		t.setByCSPRNG();
		
		
		//compute gHat1 = g1^(beta) and gHat2 = g2^(beta)
		G1 gHat1 = new G1();
		Mcl.mul(gHat1, g1, beta);

		G2 gHat2 = new G2();
		Mcl.mul(gHat2, g2, beta);
		
		//precompute K
		precompute(g1, gHat2, alpha, l, t);
		
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
		
		//second part of PK is going to be a set containing the set elements: {[gHat1^(alpha^j), gHat2^(alpha^k)}
		//for all j in [0, l] and all k in [0, l - 2]
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
		
		//add e(g1, gHat2)^(alpha^(l - 1)) to PK
		GT element3 = new GT();
		Mcl.pairing(element3, g1, gHat2);
		Fr exp3 = Tools.power(alpha, l - 1);
		Mcl.pow(element3, element3, exp3);
		PK.add(element3);
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
	
	//input:  secret key SK and int i
	//output: ArrayList<Object> d, which contains all the individual secret keys
	public static Object[] keyGen(int i, Fr[] SK) {
		
		//Extract from SK
		Fr alpha = SK[0];
		Fr gamma = SK[1];
		Fr kappa = SK[2];
		
		//1. compute ri = phi(kappa, i) where phi is a PRF
		phi = new CustomPRF(kappa, i);
		Fr ri = phi.compute(lambda);
		
		//2. compute hi = g2^((gamma - ri)/(alpha - i))
		G2 hi = new G2();
		Fr exp = new Fr();
		Fr expNum = new Fr();
		Fr expDen = new Fr();
		Mcl.sub(expNum, gamma, ri);
		Mcl.sub(expDen, alpha, new Fr(i));
		Mcl.div(exp, expNum, expDen);
		Mcl.mul(hi, g2, exp);
		
		//return (ri, hi)
		Object[] di = new Object[2];
		di[0] = ri;
		di[1] = hi;
		return di;
		
	}
	
	
	public static Object[] enc(ArrayList<Integer> S, ArrayList<Object> PK) {
		
		//extract l and n, compute k
		int n = (int) PK.get(0);
		int l = (int) PK.get(1);
		
		ArrayList<Object> setElements = (ArrayList<Object>) PK.get(5);
		Object[] setElement = (Object[]) setElements.get(0);
		G1 gHat1 = (G1) setElement[0];
		
		//compute Px
		Fr[] Px = computePx(n, l, -1, S);
		
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
		Mcl.mul(C3, g1, t);
		
		GT C4 = new GT();
		Mcl.pow(C4, (GT) PK.get(4), t);
		
		//add to Hdr and result
		Object[] Hdr = {C1, C2, C3, C4};
		Object[] result = {Hdr, K};
		return result;
		
	}
	
	
	//returns the key K1 of type GT
	public static GT decrypt(ArrayList<Integer> S, int i, Object[] di, Object[] Hdr, ArrayList<Object> PK) {
		
		//1. extract from PK and compute P(x)
		int n = (int) PK.get(0);
		int l = (int) PK.get(1);
		Fr[] Px = computePx(n, l, i, S); //P(x)/(x - i)
		//2. extract from Hdr and di
		G1 C1 = (G1) Hdr[0];
		G1 C2 = (G1) Hdr[1];
		G1 C3 = (G1) Hdr[2];
		GT C4 = (GT) Hdr[3];
		Fr ri = (Fr) di[0];
		G2 hi = (G2) di[1];
		
		//3. compute gHat2^(Pi(alpha))
		ArrayList<Object> setElements = (ArrayList<Object>) PK.get(5);
		Fr[] negated = LagrangeInterpolationZp.negate(Px); //note: need to ignore item at position 0 in this array because this is an l - 2 degree polyn.
		Object[] setElement = (Object[]) setElements.get(0);
		G2 gHat2 = (l > 1) ? (G2) setElement[1] : new G2(); //if l = 1, then gHat2 does not need to be extracted
		G2 fin = new G2(gHat2);
		Mcl.mul(fin, fin, negated[negated.length - 1]);
		//calculate gHat1^(P(alpha))
		int index = 1;
		for (int j = negated.length - 1; j > 1; j--) {
			G2 addend = new G2();
			Mcl.mul(addend, (G2) (((Object[]) setElements.get(index))[1]) , negated[j - 1]);
			Mcl.add(fin, fin, addend);
			index += 1;
		}	
		
		if (l == 1) //special case
			Mcl.mul(fin, fin, new Fr(0));
		
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
	private static Fr[] computePx(int n, int l, int i, ArrayList<Integer> S) {
		
		int k = S.size();
		ArrayList<Integer> ijValues = new ArrayList<Integer>();
		for (int j = k + 1; j <= l; j++) {
			ijValues.add(n +  j);
		}
		
		ArrayList<Integer> allIJValues = new ArrayList<Integer>();
		allIJValues.addAll(S);
		allIJValues.addAll(ijValues);
		
		//compute P(x) (or Pi(x), if i is in allIJValues)
		Fr[] Px = {new Fr(1)};
		for (int j = 1; j <= l; j++) {
			Fr item = new Fr(allIJValues.get(j - 1));
			if (item.equals(new Fr(i)))
				continue;
			//else
			Fr secondElement = new Fr();
			Mcl.mul(secondElement, item, new Fr(-1));
			Fr[] next = {new Fr(1), secondElement};
			Px = LagrangeInterpolationZp.multiply(Px, next);
		}
		
		return Px;	
	}
	
	
	//precomputes K = e(gHat2, g1)^(alpha * gamma^(l - 1) * t)
	private static void precompute(G1 g1, G2 gHat2, Fr alpha, int l, Fr t) {

		K = new GT();
		Mcl.pairing(K, g1, gHat2);
		Fr exp = new Fr();
		Mcl.mul(exp, gamma, t);
		Fr alphaExp = Tools.power(alpha, l - 1);
		Mcl.mul(exp, exp, alphaExp);
		Mcl.pow(K, K, exp);
		
	}
	
	private static long[] printRuntimes(int n, int l, int lambda) {
		
		long startSetup = System.nanoTime();
		Object[] setup = setup(n, l, lambda);
		long elapsedSetup = System.nanoTime() - startSetup;
		double secondsSetup = ((double) elapsedSetup) / 1E9;
		
		ArrayList<Object> PK = (ArrayList<Object>) setup[0];
		Fr[] MSK = (Fr[]) setup[1];
		
		ArrayList<Integer> S = new ArrayList<Integer>();
		//randomly generate numbers to put in the subset S (NO REPEATS)
		ArrayList<Integer> randomNums = new ArrayList<Integer>();
		for (int i = 0; i < n; i++) { 
			randomNums.add(i + 1);
		}
		for (int i = 0; i < l; i++) {
			int randomIndex = ThreadLocalRandom.current().nextInt(0, randomNums.size());
			int randomID = randomNums.get(randomIndex);
			S.add(randomID);
			randomNums.remove(randomIndex);
		}
		
		//generate ciphertext and keys
		int i = S.get(ThreadLocalRandom.current().nextInt(0, S.size()));
		long startEncrypt = System.nanoTime();
		Object[] C = enc(S, PK);
		long elapsedEncrypt = System.nanoTime() - startEncrypt;
		double secondsEncrypt = ((double) elapsedEncrypt) / 1E9;
		long startKeyGen = System.nanoTime();
		Object[] di = keyGen(i, MSK);
		long elapsedKeyGen = System.nanoTime() - startKeyGen;
		double secondsKeyGen = ((double) elapsedKeyGen) / 1E9;
		
		//finally, decrypt
		Object[] Hdr = (Object[]) C[0];
		GT K1 = (GT) C[1];
		long startDecrypt = System.nanoTime();
		GT K2 = decrypt(S, i, di, Hdr, PK);
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
	for (int N = 10; N <= 1000000; N *= 10) {
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
