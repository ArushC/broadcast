package schemesRevised;
import com.herumi.mcl.*;
import helperclasses.polynomials.LagrangeInterpolationZp;
import helperclasses.miscellaneous.*;
import java.io.*;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ThreadLocalRandom;

//The scheme: https://eprint.iacr.org/2008/268.pdf (4.1, page 10)
//Note one part that was very difficult to figure out because it is not explained in the paper:
//The only way C4 can be computed without accessing the private key is as e((g1^alpha)^F(alpha), (gHat2)^(alpha^(l-2)))^t
//Also note that this IBBE system does not work when l = 1
public class IBBESystemRevised {

	private static G2 g2;
	private static GT K;
	private static int lambda;
	
	//PRECONDITION: G is of order p >= n + l, l <= n
	//input: n = # of users, l = maximal size of a broadcast recipient group
	//output: Object[] containing the public key PK and private key SK
	public static Object[] setup(int n, int l, int lambda) {
		
		Mcl.SystemInit(lambda);
		IBBESystemRevised.lambda = lambda;
		
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
		
		//second part of PK is going to be the following
		//first element of the object array is a set containing the g1 set elements: {[g1^(alpha^j), gHat1^(alpha^j) for all 0 <= j <= l}
		//second element of the object array is a set containing the g2 set elements {[g2^(alpha^j), gHat2^(alpha^j) for all 0 <= j <= l - 2}
		ArrayList<Object> g1Set = new ArrayList<Object>();
		ArrayList<Object> g2Set = new ArrayList<Object>();
		
		Fr exp = new Fr(1);
		for (int j = 0; j <= l; j++) {
				G1 e1 = new G1();
				G1 e2 = new G1();
				Mcl.mul(e1, g1, exp);
				Mcl.mul(e2, gHat1, exp);
				Object[] elmnt1 = {e1, e2};
				g1Set.add(elmnt1);
				
				if (j <= l - 2) {
					G2 e3 = new G2();
					G2 e4 = new G2();
					Mcl.mul(e3, g2, exp);
					Mcl.mul(e4, gHat2, exp);
					Object[] elmnt2 = {e3, e4};
					g2Set.add(elmnt2);
				}	
				Mcl.mul(exp, exp, alpha);
		}
		Object[] pkSecondPart = {g1Set, g2Set};
		PK.add(pkSecondPart);
		
		//Generate random key kappa for a pseudorandom function psi
		Fr kappa = new Fr();
		kappa.setByCSPRNG();
		
		//add alpha, gamma, kappa to SK
		Fr[] SK = {alpha, gamma, kappa};
		
		//return (PK, SK)
		Object[] result = {PK, SK};
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
		CustomPRF phi = new CustomPRF(kappa, i);
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
		Object[] di = {ri, hi};
		return di;
		
	}
	
	//outputs a random l - 1 degree polynomial
	public static Fr[] tagGen(ArrayList<Integer> S, ArrayList<Object> PK) {
		
		//extract l and n, compute k
		int n = (int) PK.get(0);
		int l = (int) PK.get(1);
		int k = S.size();
		
		ArrayList<Fr> xVals = new ArrayList<Fr>();
		ArrayList<Fr> yVals = new ArrayList<Fr>();
		
		for (int j = k + 1; j <= l; j++) {
			xVals.add(new Fr(n + j));
			yVals.add(new Fr(1)); //F(n + j) = 1 for j in [k + 1, l]
		}
		
		//this is the "random" part: the x values are the i-values in {1, ..., k}, F(i) = random Fr
		for (int i = 1; i <= k; i++) {
			Fr fi = new Fr();
			fi.setByCSPRNG();
			xVals.add(new Fr(i));
			yVals.add(new Fr(fi));
		}
		
		//calculate the lagrange interpolation based on these points
		Fr[] xValues = xVals.toArray(new Fr[0]);
		Fr[] yValues = yVals.toArray(new Fr[0]);
		Fr[] tau = LagrangeInterpolationZp.laGrange(xValues, yValues, l - 1);
		
		return tau;
		
	}
	
	public static Object[] tagEncrypt(Fr[] tau, ArrayList<Integer> S, ArrayList<Object> PK) {
		
		//extract l and n, compute k
		int n = (int) PK.get(0);
		int l = (int) PK.get(1);
		
		//extract g1, gHat1, gHat2
		Object[] pkSecondPart = (Object[]) PK.get(4);
		ArrayList<Object> g1Parts = (ArrayList<Object>) pkSecondPart[0];
		ArrayList<Object> g2Parts = (ArrayList<Object>) pkSecondPart[1];
		Object[] setElement = (Object[]) g1Parts.get(0);
		G1 g1 = (G1) setElement[0];
		G1 gHat1 = (G1) setElement[1];
		G1 g1ToTheAlpha = (G1) ((Object[]) g1Parts.get(1))[0];
		G2 gHat2 = (G2) ((Object[]) g2Parts.get(0))[1];
		
		//compute Px
		Fr[] Px = computePx(n, l, -1, S);
		
		Fr t = new Fr();
		t.setByCSPRNG();
		Mcl.pow(K, K, t);
		
		//compute C1, C2, C3, C4
		G1 C1 = new G1(gHat1);
		Mcl.mul(C1, C1, Px[Px.length - 1]);
		//compute gHat1^P(alpha)
		int index = 1;
		for (int i = Px.length - 1; i > 0; i--) {
			G1 addend = new G1();
			Mcl.mul(addend, (G1) (((Object[]) g1Parts.get(index))[1]) , Px[i - 1]);
			Mcl.add(C1, C1, addend);
			index += 1;
		}
		Mcl.mul(C1, C1, t);
		
		G1 C2 = new G1();
		Mcl.mul(C2, (G1) PK.get(2), t);
		
		G1 C3 = new G1(g1);
		Mcl.mul(C3, C3, tau[tau.length - 1]);
		//compute g1^(F(alpha))
		index = 1;
		for (int i = tau.length - 1; i > 0; i--) {
			G1 addend = new G1();
			Mcl.mul(addend, (G1) (((Object[]) g1Parts.get(index))[0]) , tau[i - 1]);
			Mcl.add(C3, C3, addend);
			index += 1;
		}
		Mcl.mul(C3, C3, t);
		
		GT C4 = new GT();
		//1. calculate (gHat2)^(alpha^(l-2))
		G2 g2ToTheAlphaPower = (G2) (((Object[]) g2Parts.get(g2Parts.size() - 1)))[1];
		//2. calculate (g1^alpha)^F(alpha)
		G1 g1ToThePower = new G1(g1ToTheAlpha);
		Mcl.mul(g1ToThePower, g1ToThePower, tau[tau.length - 1]);
		index = 2;
		for (int i = tau.length - 1; i > 0; i--) {
			G1 addend = new G1();
			Mcl.mul(addend, (G1) (((Object[]) g1Parts.get(index))[0]) , tau[i - 1]);
			Mcl.add(g1ToThePower, g1ToThePower, addend);
			index++;
		}
		
		Mcl.pairing(C4, g1ToThePower, g2ToTheAlphaPower);
		Mcl.pow(C4, C4, t);
		
		//add to Hdr
		Object[] Hdr = {C1, C2, C3, C4};
		
		//add to and return result
		Object[] result = {tau, Hdr, K};
		return result;
		
	}
	
	//returns an array containing (tau, Hdr, K)
	public static Object[] enc(ArrayList<Integer> S, ArrayList<Object> PK) {
		Fr[] tau = tagGen(S, PK);
		return tagEncrypt(tau, S, PK);
	}
	
	//returns the key K1 of type GT
	public static GT decrypt(ArrayList<Integer> S, int i, Object[] di, Fr[] tau, Object[] Hdr, ArrayList<Object> PK) {
		
		//1. extract from PK and compute P(x)
		Object[] pkSecondPart = (Object[]) PK.get(4);
		ArrayList<Object> g2Parts = (ArrayList<Object>) pkSecondPart[1];
		int n = (int) PK.get(0);
		int l = (int) PK.get(1);
		Fr[] Px = computePx(n, l, i, S);
		
		//2. extract from Hdr and di
		G1 C1 = (G1) Hdr[0];
		G1 C2 = (G1) Hdr[1];
		G1 C3 = (G1) Hdr[2];
		GT C4 = (GT) Hdr[3];
		Fr ri = (Fr) di[0];
		G2 hi = (G2) di[1];
		
		//3. compute g2^(Fi(alpha))
		Fr Fi = LagrangeInterpolationZp.computeFxHorner(tau, new Fr(i));
		Fr[] subtracted = tau.clone();
		Fr constant = new Fr();
		Mcl.sub(constant, subtracted[subtracted.length - 1], Fi);
		subtracted[subtracted.length - 1] = constant;
		Fr[] FiPoly = LagrangeInterpolationZp.syntheticDivide(subtracted, new Fr(i));
		G2 fiFin = new G2(g2);
		Mcl.mul(fiFin, fiFin, FiPoly[FiPoly.length - 2]); //ignore last element because it is the remainder
		int index = 1;
		for (int j = FiPoly.length - 2; j > 0; j--) {
			G2 addend = new G2();
			Mcl.mul(addend, (G2) (((Object[]) g2Parts.get(index))[0]) , FiPoly[j - 1]);
			Mcl.add(fiFin, fiFin, addend);
			index += 1;
		}
		
		//4. compute gHat2^(Pi(alpha))
		Fr[] negated = LagrangeInterpolationZp.negate(Px); //note: need to ignore item at position 0 in this array because this is an l - 2 degree polyn.
		G2 gHat2 = (G2) (((Object[]) ((ArrayList<Object>) pkSecondPart[1]).get(0))[1]); //if l = 1, then gHat2 does not need to be extracted
		G2 piFin = new G2(gHat2);
		Mcl.mul(piFin, piFin, negated[negated.length - 1]);
		index = 1;
		for (int j = negated.length - 1; j > 1; j--) {
			G2 addend = new G2();
			Mcl.mul(addend, (G2) (((Object[]) g2Parts.get(index))[1]) , negated[j - 1]);
			Mcl.add(piFin, piFin, addend);
			index += 1;
		}
		
		Fr ei = new Fr();
		Mcl.div(ei, ri, Fi);
		Mcl.mul(ei, ei, new Fr(-1));
		
		//compute the pairings
		GT e1 = new GT();
		G2 secondE1 = new G2();
		Mcl.mul(secondE1, fiFin, ei);
		Mcl.add(secondE1, hi, secondE1);
		Mcl.pairing(e1, C1, secondE1);
		
		GT e2 = new GT();
		G1 firstE2 = new G1();
		Mcl.mul(firstE2, C3, ei);
		Mcl.add(firstE2, C2, firstE2);
		Mcl.pairing(e2, firstE2, piFin);
		
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
	
	
	//precomputes K = e(gHat2, g1)^(gamma * alpha^(l - 1))
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
		
		ArrayList<Integer> S = new ArrayList<Integer>();
		//randomly generate numbers to put in the subset S (NO REPEATS)
		ArrayList<Integer> randomNums = new ArrayList<Integer>();
		for (int i = 0; i < n; i++) { 
			randomNums.add(i + 1);
		}
		for (int i = 0; i < l; i++) {
			int randomIndex = ThreadLocalRandom.current().nextInt(1, randomNums.size());
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
		Fr[] tau = (Fr[]) C[0];
		Object[] Hdr = (Object[]) C[1];
		GT K1 = (GT) C[2];
		long startDecrypt = System.nanoTime();
		GT K2 = decrypt(S, i, di, tau, Hdr, PK);
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
