package schemesRevised;
import helperclasses.polynomials.LagrangeInterpolationZp;

import helperclasses.structures.NDMatrix;
import helperclasses.structures.VectorND;

import java.io.File;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.ThreadLocalRandom;
import com.herumi.mcl.*;


//The scheme: https://eprint.iacr.org/2020/954, section 9.3, page 58
//"we use a single helper key for all users and include the helper key in the public parameters
//to obtain a risky broadcast multi-scheme of size (N, 1, 1)"
//CHANGES MADE: include R^(-1) in the master secret key so it does not have to be recomputed during secret key extraction
//NOTE: The scheme described in Zhandry's paper does not work. The pairings do not cancel out during decryption
//We modified the construction slightly by doing the following:
//1. in the public key, second element: g1^((beta * gamma, 0, 0) * R^(-1)) --> g1^((beta, 0, 0) * R^(-1))
//2. in the vector exponent of hTheta: ((beta - tTheta)/(beta * gamma)) --> (beta - tTheta)/gamma
//3. in the vector exponent of c1: (alpha * beta * gamma) --> (alpha * gamma)
public class RiskyMTBRevised {
	
	private static int u, t, n, v;
	private static G1 g1;
	private static G2 g2;
	
	//input: v = number of users, n = vector dimension, u = t = 1
	//the larger t is, the less likely it is that decryption for a random user will be possible
	public static Object[] genMTB(int u, int v, int t, int n, int lambda) {
		
		Mcl.SystemInit(lambda);
		
		//save so it can be used in later functions
		RiskyMTBRevised.u = u;
		RiskyMTBRevised.t = t;
		RiskyMTBRevised.n = n;
		RiskyMTBRevised.v = v;
		
		//choose random beta, gamma
		Fr beta = new Fr();
		beta.setByCSPRNG();
		
		Fr gamma = new Fr();
		gamma.setByCSPRNG();
		//Choose random g1, g2
		g1 = new G1();
		Mcl.hashAndMapToG1(g1, "abc".getBytes());
		
		g2 = new G2();
		Mcl.hashAndMapToG2(g2, "def".getBytes());
		
		Object[] MSK = new Object[5]; //instantiate master secret key
		
		NDMatrix R = new NDMatrix(NDMatrix.MATRIX_RANDOM, n + 2, n + 2);
		
		//compute and add R inverse
		NDMatrix RInverted = new NDMatrix(R);
		RInverted.invert();
		MSK[3] = RInverted;
		
		//transpose and add R
		R = R.transpose();
		
		MSK[0] = beta;
		MSK[1] = gamma;
		MSK[2] = R;
		
		//add all tau_theta to the master secret key
		//note that theta = {1, 2, ..., v} because we set the bound on the number of helper keys equal to v
		Fr tauTheta = new Fr();
		tauTheta.setByCSPRNG();
		MSK[4] = tauTheta;
		
		//calculate the public key
		ArrayList<Object> PK = new ArrayList<Object>();
		
		GT e1 = new GT();
		Mcl.pairing(e1, g1, g2);
		Mcl.pow(e1, e1, beta);
		
		VectorND exp2NT = new VectorND(VectorND.VECTOR_ZERO, n + 2);
		exp2NT.getCoords().set(0, gamma);
		VectorND exp2 = NDMatrix.transform(RInverted, exp2NT);
		ArrayList<G2> e2 = VectorND.exponentiate(g2, exp2);
		
		PK.add(e1); //add first two elements to public key
		PK.add(e2);
		
		Fr gammaToTheI = new Fr(1);
		for (int i = 0; i <= v; i++) {
			VectorND exp3NT = new VectorND(VectorND.VECTOR_ZERO, n + 2);
			exp3NT.getCoords().set(0, gammaToTheI);
			VectorND exp3 = NDMatrix.transform(R, exp3NT);
			ArrayList<G1> pkJ = VectorND.exponentiate(g1, exp3);
			PK.add(pkJ);
			Mcl.mul(gammaToTheI, gammaToTheI, gamma);
		}
		
		//add the common helper key to the public key
		//compute the helper key
		ArrayList<Object> hkTheta = new ArrayList<Object>();
		VectorND expHVNT = new VectorND(VectorND.VECTOR_ZERO, n + 2);
		Fr e1ExpHVNT = new Fr();
		Mcl.sub(e1ExpHVNT, beta, tauTheta);
		Mcl.div(e1ExpHVNT, e1ExpHVNT, gamma);
		expHVNT.getCoords().set(0, e1ExpHVNT);
		VectorND expHV = NDMatrix.transform(R, expHVNT);
		ArrayList<G1> hTheta = VectorND.exponentiate(g1, expHV);	
		hkTheta.add(hTheta);
				
		//compute all the g2 vectors
		Fr gammaToThePowerHelper = new Fr(1);
		for (int j = 0; j <= v; j++) {
				//nontransformed
			VectorND expHVENT = new VectorND(VectorND.VECTOR_ZERO, n + 2);
			Fr e1ExpHVENT = new Fr();
			Mcl.mul(e1ExpHVENT, tauTheta, gammaToThePowerHelper);
			expHVENT.getCoords().set(0, e1ExpHVENT);
			//transformed
			VectorND expHVE = NDMatrix.transform(R, expHVENT);
			ArrayList<G1> elementHelper = VectorND.exponentiate(g1, expHVE);
			hkTheta.add(elementHelper);
			Mcl.mul(gammaToThePowerHelper, gammaToThePowerHelper, gamma);
		}
		
		PK.add(hkTheta);
		
		Object[] res = {PK, MSK};
		return res;
	}
	
	public static ArrayList<Object> extractMTB(Object[] MSK, Set<Integer> U, VectorND x) {
		
		
		Fr beta = (Fr) MSK[0];  //extract from MSK
		Fr gamma = (Fr) MSK[1]; 
		NDMatrix R = (NDMatrix) MSK[2];
		NDMatrix RInverted = (NDMatrix) MSK[3];
		Fr tTheta = (Fr) MSK[4];
		
		//calculate the product in the numerator
		Fr product = new Fr(1);
		Fr ONE = new Fr(1);
		for (int s: U) {
			Fr c = new Fr();
			Mcl.div(c, gamma, new Fr(s));	
			Mcl.sub(c, ONE, c);
			Mcl.mul(product, product, c);
		}

		ArrayList<Object> sk = new ArrayList<Object>();
		
		//compute all of the secret keys
		Fr gammaToThePower = new Fr(1);
		for (int l = 0; l <= u - t; l++) {
			
			Fr etaL = new Fr();
			etaL.setByCSPRNG();
			
			Fr e1ExpNT = new Fr();
			Mcl.mul(e1ExpNT, tTheta, gammaToThePower);
			Mcl.div(e1ExpNT, e1ExpNT, product);
			
			NDMatrix Dl = new NDMatrix(NDMatrix.MATRIX_DIAGONAL, n, n);
			VectorND vectorInExponent = NDMatrix.transform(Dl, x);
			
			//create ArrayList of Fr's to put the elements in the vector in order
			ArrayList<Fr> expNTVals = new ArrayList<Fr>();
			expNTVals.add(e1ExpNT);
			expNTVals.addAll(vectorInExponent.getCoords());
			expNTVals.add(etaL);
			
			VectorND expNT = new VectorND(expNTVals);
			VectorND exp = NDMatrix.transform(RInverted, expNT);
			ArrayList<G2> element = VectorND.exponentiate(g2, exp);
			sk.add(element);
			Mcl.mul(gammaToThePower, gammaToThePower, gamma);
		}
		
		return sk;
	}
	
	//encrypt with public key
	public static Object[] encMTB(ArrayList<Object> PK, Set<Integer> S) {
		
		//random alpha
		Fr alpha = new Fr();
		alpha.setByCSPRNG();
		
		ArrayList<G2> c1 = new ArrayList<G2>((ArrayList<G2>) PK.get(1));
		
		for (int l = 0; l < c1.size(); l++) {
			G2 elC1 = new G2();
			Mcl.mul(elC1, c1.get(l), alpha);
			c1.set(l, elC1);
		}
		
		
		GT encapsulatedKey = new GT((GT) PK.get(0));
		Mcl.pow(encapsulatedKey, encapsulatedKey, alpha);
		
		//calculate c2
		//first compute coefficients of polynomial in numerator
		Fr[] coefficients = {new Fr(1)};
		for (int s: S) {
			Fr[] next = {new Fr(-1), new Fr(1)};
			Mcl.div(next[0], next[0], new Fr(s));
			coefficients = LagrangeInterpolationZp.multiply(coefficients, next);
		}
		
		//then use the coefficients to calculate c2
		ArrayList<G1> c2 = new ArrayList<G1>((ArrayList<G1>) PK.get(2));
		int k = coefficients.length;
		for (int z = 0; z < c2.size(); z++) {
			G1 helper = new G1();
			Mcl.mul(helper, c2.get(z), coefficients[k - 1]);
			c2.set(z, helper);
		}
		
		for (int i = k - 2; i >= 0; i--) {
			//if (coefficients.length == 1) break;
			ArrayList<G1> element = new ArrayList<G1>((ArrayList<G1>) PK.get(k - i + 1));
			for (int y = 0; y < c2.size(); y++) {
				G1 p = new G1();
				Mcl.mul(p, element.get(y), coefficients[i]);
				element.set(y, p);
			}
			
			for (int j = 0; j < element.size(); j++) {
				G1 b = new G1();
				Mcl.add(b, c2.get(j), element.get(j));
				c2.set(j, b);
			}
			
		}
		
		for (int l = 0; l < c2.size(); l++) {
			G1 b = new G1();
			Mcl.mul(b, c2.get(l), alpha);
			c2.set(l, b);
		}
		
		Object[] c = {c1, c2};
		Object[] res = {c, encapsulatedKey};
		return res;
	}
	
	public static GT decMTB(ArrayList<Object> sk, ArrayList<Object> helperKey, Set<Integer> S, Set<Integer> U, Object[] c) {
		
		//extract from the ciphertext and secret key
		ArrayList<G1> hTheta = (ArrayList<G1>) helperKey.get(0);
		ArrayList<G2> c1 = (ArrayList<G2>) c[0];
		ArrayList<G1> c2 = (ArrayList<G1>) c[1];
		
		Set<Integer> SDiff = new HashSet<Integer>(S);
		Set<Integer> UDiff = new HashSet<Integer>(U);
		
		SDiff.removeAll(U); //SDiff = S \ U
		UDiff.removeAll(S); //UDiff = U \ S
		
		//compute the polynomial Q
		Fr[] Q = {new Fr(1)};
		for (int j: SDiff) {
			Fr[] next = {new Fr(-1), new Fr(1)};
			Mcl.div(next[0], next[0], new Fr(j));
			Q = LagrangeInterpolationZp.multiply(Q, next);
		}
		
		//compute the polynomial R
		Fr[] R = {new Fr(1)};
		for (int i: UDiff) {
			Fr[] next = {new Fr(-1), new Fr(1)};
			Mcl.div(next[0], next[0], new Fr(i));
			R = LagrangeInterpolationZp.multiply(R, next);
		}
		
		//compute the polynomial P
		Fr[] PNumerator = LagrangeInterpolationZp.negate(Q);
		Fr[] P = (PNumerator.length == 1) ? new Fr[1] : new Fr[PNumerator.length - 1];
		
		for (int i = 0; i < PNumerator.length - 1; i++)
			P[i] = PNumerator[i]; //divides perfectly, so the coefficients of P are just the  
								  //coefficients of PNumerator, excluding the last coefficient
		
		if (PNumerator.length == 1)
			P[0] = new Fr(0); //special case, P is a zero polynomial
		
		ArrayList<G1> JThetaP = new ArrayList<G1>((ArrayList<G1>) helperKey.get(1));
		int p = P.length;
		for (int z = 0; z < JThetaP.size(); z++) {
			G1 helperJTP = new G1();
			Mcl.mul(helperJTP, JThetaP.get(z), P[p - 1]);
			JThetaP.set(z, helperJTP);
		}
		
		for (int i = p - 2; i >= 0; i--) {
			ArrayList<G1> element = new ArrayList<G1>((ArrayList<G1>) helperKey.get(p - i));
			for (int y = 0; y < element.size(); y++) {
				G1 a = new G1();
				Mcl.mul(a, element.get(y), P[i]);
				element.set(y, a);
			}
			//add the element by the vector JThetaP
			for (int j = 0; j < element.size(); j++) {
				G1 b = new G1();
				Mcl.add(b, JThetaP.get(j), element.get(j));
				JThetaP.set(j, b);
			}
		}
		
		//compute HThetaR
		ArrayList<G2> HThetaR = new ArrayList<G2>((ArrayList<G2>) sk.get(0));
		int r = R.length;
		for (int z = 0; z < HThetaR.size(); z++)  {
			G2 helperHTR = new G2();
			Mcl.mul(helperHTR, HThetaR.get(z), R[r - 1]);
			HThetaR.set(z, helperHTR);
		}

		for (int i = r - 2; i >= 0; i--) {
			if (U.size() == 1) break;
			ArrayList<G2> element = new ArrayList<G2>((ArrayList<G2>) sk.get(r - i - 1));
			for (int y = 0; y < element.size(); y++) {
				G2 a = new G2();
				Mcl.mul(a, element.get(y), R[i]);
				element.set(y, a);
			}
			
			//add the element by the vector JThetaP
			for (int j = 0; j < element.size(); j++) {
				G2 b = new G2();
				Mcl.add(b, HThetaR.get(j), element.get(j));
				HThetaR.set(j, b);
			}
		}
		
		GT e1 = VectorND.vectorPairing(hTheta, c1);
		GT e2 = VectorND.vectorPairing(JThetaP, c1);
		GT e3 = VectorND.vectorPairing(c2, HThetaR);
		Mcl.mul(e1, e1, e2);
		Mcl.mul(e1, e1, e3);
		return e1;
	}
	
	//note that this scheme sets up M instances of a broadcast encryption system with N users, for a total of N * M users
	//encryption happens to the jth instance. In each instance, only the users in that instance with i-values in the subset S can decrypt
	public static void printRuntimes(int N, int M, int subsetSize, int lambda) {
		
		//random ID for the encryption & decryption
		int ID = ThreadLocalRandom.current().nextInt(1, N + 1);
		int j = (int) Math.ceil(((double) ID) / N);
		int i = (ID % N == 0) ? N : (ID % N);
		
		Set<Integer> S = new HashSet<Integer>();
		S.add(i);
		
		//randomly generate numbers to put in the subset S (NO REPEATS)
		ArrayList<Integer> randomNums = new ArrayList<Integer>();
		int randomID = 0;
		for (int z = 0; z < N; z++) {
			if (z + 1 == i) continue; //already in the subset
			randomNums.add(z + 1);
		}
		for (int z = 1; z < subsetSize; z++) {
			int randomIndex = ThreadLocalRandom.current().nextInt(0, randomNums.size());
			randomID = randomNums.get(randomIndex);
			S.add(randomID);
			randomNums.remove(randomIndex);
		}
		
		long startSetup = System.nanoTime();
		Object[] setup = RiskyBroadcastMultischeme.setup(N, M, Mcl.BN254);
		long elapsedSetup = System.nanoTime() - startSetup;
		double secondsSetup = ((double) elapsedSetup)/1E9;
		
		//extract public/tracing key
		ArrayList<Object> PK = (ArrayList<Object>) setup[0];
		ArrayList<Object> tk = (ArrayList<Object>) setup[1];
		
		long startEncrypt = System.nanoTime();
		Object[] enc = RiskyBroadcastMultischeme.enc(PK, j, S);
		long elapsedEncrypt = System.nanoTime() - startEncrypt;
		double secondsEncrypt = ((double) elapsedEncrypt)/1E9;
		
		//extract from ciphertext
		GT encapsulatedKey = (GT) enc[1];
		Object[] c = (Object[]) enc[0];
		
		long startKeyGen = System.nanoTime();
		ArrayList<Object> skID = RiskyBroadcastMultischeme.keyGen(ID,  tk);
		long elapsedKeyGen = System.nanoTime() - startKeyGen;
		double secondsKeyGen = ((double) elapsedKeyGen)/1E9;
		ArrayList<Object> helperKey = (ArrayList<Object>) PK.get(PK.size() - 1);
		
		//finally, decrypt
		long startDec = System.nanoTime();
		GT encapsulatedKey1 = RiskyBroadcastMultischeme.decrypt(skID, helperKey, S, ID, c);
		long elapsedDec = System.nanoTime() - startDec;
		double secondsDec = ((double) elapsedDec)/1E9;
		
		//Finally, print out the results
		String success = (encapsulatedKey.equals(encapsulatedKey1)) ? "SUCCESSFUL DECRYPTION (" + (N * M) + " total users)": "FAILED DECRYPTION";
		System.out.println(success + ": " + "N = " + N + ", M = " + M + ", subset size = " + subsetSize);
		System.out.println(); //padding
		System.out.println("setup took " + secondsSetup + " seconds");
		System.out.println("encryption took " + secondsEncrypt + " seconds");
		System.out.println("key generation took " + secondsKeyGen + " seconds");
		System.out.println("decryption took " + secondsDec + " seconds (u = " + ID + ")");
		System.out.println(); //more padding
		
		
	}
	
	public static void testRuntimes(int percent, int lambda) {
		for (int N = 100; N <= 1000000; N *= 10) {
			//int n = (int) Math.ceil(Math.pow(N, 0.5));
			//int m = (int) Math.ceil(((double) N)/n);
			int subsetSize = (int) Math.round(0.01 * percent * N);
			printRuntimes(N, 1, subsetSize, lambda);
		}
		
	}
	
	//TEST THE MTB SCHEME
	/*public static void main(String[] args) {
		File lib = new File("../../lib/libmcljava.dylib");
		System.load(lib.getAbsolutePath());
		Object[] gen = genMTB(10, 1000, 1, 10, Mcl.BN254);
		ArrayList<Object> PK = (ArrayList<Object>) gen[0];
		Object[] MSK = (Object[]) gen[1];
			
		Set<Integer> S = new HashSet<Integer>();
		S.addAll(Arrays.asList(1, 4, 6, 7, 10));
		Set<Integer> U = new HashSet<Integer>();
		U.addAll(Arrays.asList(1));
		VectorND x = new VectorND(VectorND.VECTOR_ZERO, 10);
		ArrayList<Object> sk = extractMTB(MSK, U, x);
		ArrayList<Object> helperKey = (ArrayList<Object>) PK.get(PK.size() - 1);
		Object[] enc = encMTB(PK, S);
		Object[] C = (Object[]) enc[0];
		GT encapsulatedKey = (GT) enc[1];
		GT encapsulatedKey1 = decMTB(sk, helperKey, S, U, C);
		System.out.println("Encapsulated Key Actual: " + encapsulatedKey);
		System.out.println("Decrypted Encapsulated Key: " + encapsulatedKey1);
	}*/
	
	public static void main(String[] args) {
		File lib = new File("../../lib/libmcljava.dylib");
		System.load(lib.getAbsolutePath());
		testRuntimes(10, Mcl.BN254);
	}
	
}

//this is the risky broadcast & trace scheme given in section 9.3 of https://eprint.iacr.org/2020/954
//note: page 36 is helpful
class RiskyBroadcastMultischeme {
	
	private static int N, M;
	
	//N = # of users in broadcast system, M = # of instances of the system, lambda = security parameter
	public static Object[] setup(int N, int M, int lambda) {
		
		//save to be used in later functions
		RiskyBroadcastMultischeme.N = N;
		RiskyBroadcastMultischeme.M = M;
		
		int n = 2;
		int v = N;
		int u = 1, t = 1;
		Object[] setup = RiskyMTBRevised.genMTB(u, v, t, n, lambda);
		ArrayList<Object> PK = (ArrayList<Object>) setup[0];
		Object[] MSK = (Object[]) setup[1];
		ArrayList<Object> secretKeys = new ArrayList<Object>();
		
		//initialize the tracing key
		ArrayList<Object> tk = new ArrayList<Object>();
		tk.add(MSK);
		
		//for each j in M, choose random iStarJ in {1, 2, ..., N}
		for (int j = 1; j <= M; j++) {
			int iStarJ = ThreadLocalRandom.current().nextInt(1, N + 1);
			tk.add(iStarJ);
		}
		
		Object[] result = {PK, tk};
		return result;

	}
	
	//setup authority generates the key for an individual user
	public static ArrayList<Object> keyGen(int ID, ArrayList<Object> tk) {
		
		//extract from ID and tracing key
		int j = (int) Math.ceil(((double) ID) / N);
		int i = (ID % N == 0) ? N : (ID % N);
		Object[] MSK = (Object[]) tk.get(0);
		int iStarJ = (int) tk.get(j);
		
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
		
		//finally, generate the key
		Set<Integer> U = new HashSet<Integer>();
		U.add(ID);
		ArrayList<Object> skID = RiskyMTBRevised.extractMTB(MSK, U, wJI);
		return skID;
	}
	
	public static Object[] enc(ArrayList<Object> PK, int j, Set<Integer> S) {
		
		Set<Integer> TjS = new HashSet<Integer>();
		
		for (int i: S) {
			int ID = (j - 1) * N + i;
			TjS.add(ID);
		}
		
		return RiskyMTBRevised.encMTB(PK, TjS);
		
	}
	
	//decrypt for user u with secret key skID
	public static GT decrypt(ArrayList<Object> skID, ArrayList<Object> helperKey, Set<Integer> S, int userID, Object[] c) {
		
		int j = (int) Math.ceil(((double) userID) / N);
		
		Set<Integer> TjS = new HashSet<Integer>();
		
		for (int i: S) {
			int ID = (j - 1) * N + i;
			TjS.add(ID);
		}
		
		Set<Integer> U = new HashSet<Integer>();
		U.add(userID);
		
		return RiskyMTBRevised.decMTB(skID, helperKey, TjS, U, c);
	}
	
}
