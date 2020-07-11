package schemes;
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
public class IBBESystemVariant {

	private static CustomPRF phi;
	private static Fr alpha, gamma, t;
	private static GT K;
	private static int lambda;
	private static G2 g2, gHat2;
	
	//PRECONDITION: G is of order p >= n + l, l <= n
	//input: n = # of users, l = maximal size of a broadcast recipient group
	//output: Object[] containing the public key PK and private key SK
	public static Object[] setup(int n, int l, int lambda) {
		
		Mcl.SystemInit(lambda);
		IBBESystemVariant.lambda = lambda;
		//instead of creating a GroupGen function, the groups are generated here
		
		//generate random g1, g2
		G1 g1 = new G1();
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

		gHat2 = new G2();
		Mcl.mul(gHat2, g2, beta);
		
		//precompute K
		precompute(g1, gHat2, l, t);
		
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
		
		//second part of PK is going to be a set containing the g1 set elements: {[g1^(alpha^j), gHat1^(alpha^j) for all 0 <= j <= l}
		ArrayList<Object> g1Set = new ArrayList<Object>();
	
		for (int j = 0; j <= l; j++) {
				Object[] elmnt1 = new Object[2];
				G1 e1 = new G1();
				G1 e2 = new G1();
				Mcl.mul(e1, g1, Tools.power(alpha, j));
				Mcl.mul(e2, gHat1, Tools.power(alpha, j));
				elmnt1[0] = e1;
				elmnt1[1] = e2;
				g1Set.add(elmnt1);				
		}

		PK.add(g1Set);
		
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
		int k = S.size();
		
		//extract g1, gHat1, gHat2
		ArrayList<Object> g1Parts = (ArrayList<Object>) PK.get(4);
		Object[] setElement = (Object[]) g1Parts.get(0);
		G1 g1 = (G1) setElement[0];
		G1 gHat1 = (G1) setElement[1];
		
		//compute Px1
		Fr[] Px = computePx(n, l, k, S);
		
		//compute C1, C2, C3, C4
		G1 C1 = new G1();
		Fr exp1 = LagrangeInterpolationZp.computeFx(Px, alpha); //P(alpha)
		Mcl.mul(exp1, exp1, t);
		Mcl.mul(C1, gHat1, exp1);
		
		G1 C2 = new G1();
		Mcl.mul(C2, (G1) PK.get(2), t);
		
		G1 C3 = new G1();
		Mcl.mul(C3, g1, t);
		
		GT C4 = new GT();
		Fr exp4 = new Fr(t);
		Mcl.mul(exp4, exp4, Tools.power(alpha, l - 1));
		Mcl.pairing(C4, g1, gHat2);
		Mcl.pow(C4, C4, exp4);
		
		//add to Hdr
		Object[] Hdr = new Object[4];
		Hdr[0] = C1;
		Hdr[1] = C2;
		Hdr[2] = C3;
		Hdr[3] = C4;
		
		//add to and return result
		Object[] result = new Object[3];
		result[0] = Hdr;
		result[1] = K;
		return result;
		
	}
	
	
	//returns the key K1 of type GT
	public static GT decrypt(ArrayList<Integer> S, int i, Object[] di, Object[] Hdr, ArrayList<Object> PK) {
		
		//1. extract from PK and compute P(x)
		int n = (int) PK.get(0);
		int l = (int) PK.get(1);
		int k = S.size();
		Fr[] Px = computePx(n, l, k, S);
		
		//2. extract from Hdr and di
		G1 C1 = (G1) Hdr[0];
		G1 C2 = (G1) Hdr[1];
		G1 C3 = (G1) Hdr[2];
		GT C4 = (GT) Hdr[3];
		Fr ri = (Fr) di[0];
		G2 hi = (G2) di[1];
		
		//3. compute Pi(alpha) and ei
		Fr Pialpha = new Fr(Tools.power(alpha, l - 1));
		Fr numeratorPi = LagrangeInterpolationZp.computeFx(Px, alpha);
		Fr denominatorPi = new Fr();
		Mcl.sub(denominatorPi, alpha, new Fr(i));
		Fr fracPi = new Fr();
		Mcl.div(fracPi, numeratorPi, denominatorPi);
		Mcl.sub(Pialpha, Pialpha, fracPi);
		
		Fr ei = new Fr();
		Mcl.mul(ei, ri, new Fr(-1));
		
		//compute the pairings
		GT e1 = new GT();
		Mcl.pairing(e1, C1, hi);
		
		GT e2 = new GT();
		G1 firstE2 = new G1();
		Mcl.mul(firstE2, C3, ei);
		Mcl.add(firstE2, C2, firstE2);
		G2 secondE2 = new G2();
		Mcl.mul(secondE2, gHat2, Pialpha);
		Mcl.pairing(e2, firstE2, secondE2);
		
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
	private static Fr[] computePx(int n, int l, int k, ArrayList<Integer> S) {
		
		
		ArrayList<Integer> ijValues = new ArrayList<Integer>();
		
		for (int j = k + 1; j <= l; j++) {
			ijValues.add(n +  j);
		}
		
		ArrayList<Integer> allIJValues = new ArrayList<Integer>();
		allIJValues.addAll(S);
		allIJValues.addAll(ijValues);
		
		//compute P(x)
		Fr[] Px = {new Fr(1)};
		for (int j = 1; j <= l; j++) {
			Fr secondElement = new Fr();
			Mcl.mul(secondElement, new Fr(allIJValues.get(j - 1)), new Fr(-1));
			Fr[] next = {new Fr(1), secondElement};
			Px = LagrangeInterpolationZp.multiply(Px, next);
		}
		
		return Px;
		
	}
	
	
	//precomputes K = e(gHat2, g1)^(alpha * gamma^(l - 1) * t)
	private static void precompute(G1 g1, G2 gHat2, int l, Fr t) {
		
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
	

public static double[] testRuntimes(int lambda) {
		
		long totalSetupTime = 0, totalEncryptionTime = 0, totalKeyGenTime = 0, totalDecryptionTime = 0;
		
		for (int n = 5000; n < 100000; n+= 5000) {
			long[] elapsedTimes = printRuntimes(n, 1000, lambda); //constant l = 1000; changing n does not affect the other times
			totalSetupTime += elapsedTimes[0];
			totalEncryptionTime += elapsedTimes[1];
			totalKeyGenTime += elapsedTimes[2];
			totalDecryptionTime += elapsedTimes[3];
		}
		
		int count = 0;
		for (int l = 1000; l <= 10000; l += 1000) {
			printRuntimes(100000, l, lambda); //value of n really doesn't matter, as long as it's constant
			count++;
		}
		double averageSetupTime = ((double) totalSetupTime) / (1E9 * count);
		double averageEncryptionTime = ((double) totalEncryptionTime) / (1E9 * count);
		double averageKeyGenTime = ((double) totalKeyGenTime) / (1E9 * count);
		double averageDecryptionTime = ((double) totalDecryptionTime) / (1E9 * count);
		
		System.out.println("Average setup time, constant l = 1000: " + averageSetupTime + " seconds");
		System.out.println("Average encryption time, constant l = 1000: " + averageEncryptionTime + " seconds");
		System.out.println("Average key generation time, constant l = 1000: " + averageKeyGenTime + " seconds");
		System.out.println("Average decryption time, constant l = 1000: " + averageDecryptionTime + " seconds");
		
		//return a double array -- could run some statistical analyses on the sampling distributions of mean times for a sample size of n = 100
		double[] result = new double[4];
		result[0] = averageSetupTime;
		result[1] = averageEncryptionTime;
		result[2] = averageKeyGenTime;
		result[3] = averageDecryptionTime;
		
		return result;
		
}
	
	
	
	public static void main(String[] args) {
		System.load("/Users/arushchhatrapati/Documents/mcl/lib/libmcljava.dylib");
		testRuntimes(Mcl.BN254);	
	}
	
	
	
}
