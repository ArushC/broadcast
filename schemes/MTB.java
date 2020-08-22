package schemes;
import helperclasses.polynomials.LagrangeInterpolationZp;
import helperclasses.structures.NDMatrix;
import helperclasses.structures.VectorND;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import com.herumi.mcl.*;

//
//The scheme: (include link to Zhandry's paper here once it is published), section 9.2
//CHANGES MADE: include R^(-1) in the master secret key so it does not have to be recomputed during secret key extraction
//NOTE: The scheme described in Zhandry's paper does not work. The pairings do not cancel out during decryption
//We modified the construction slightly by doing the following:
//1. in the public key, second element: g1^((beta * gamma, 0, 0) * R^(-1)) --> g1^((beta, 0, 0) * R^(-1))
//2. in the vector exponent of hTheta: ((beta - tTheta)/(beta * gamma)) --> (beta - tTheta)/gamma
//3. in the vector exponent of c1: (alpha * beta * gamma) --> (alpha * gamma)
public class MTB {
	
	private static ArrayList<Integer> chi = new ArrayList<Integer>(); //helper keys
	private static int u, t, n, v;
	private static G1 g1;
	private static G2 g2;
	
	//input: v = number of users, n = vector dimension, u = t = 1
	//the larger t is, the less likely it is that decryption for a random user will be possible
	public static Object[] genMTB(int u, int v, int t, int n, int lambda) {
		
		Mcl.SystemInit(lambda);
		
		//save so it can be used in later functions
		MTB.u = u;
		MTB.t = t;
		MTB.n = n;
		MTB.v = v;
		chi.addAll(Arrays.asList(1, 2, 3, 4, 5, 6, 7, 8, 9, 10));
		
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
		
		Object[] MSK = new Object[chi.size() + 4]; //instantiate master secret key
		
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
		
		//add tau_theta for each theta in chi to the master secret key
		for (int theta: chi) {
			Fr tauTheta = new Fr();
			tauTheta.setByCSPRNG();
			MSK[theta + 3] = tauTheta; //assuming first is tau_1, ..., tau_n
		}
		
		//calculate the public key
		ArrayList<Object> PK = new ArrayList<Object>();
		
		GT e1 = new GT();
		Mcl.pairing(e1, g1, g2);
		Mcl.pow(e1, e1, beta);
		
		VectorND exp2NT = new VectorND(VectorND.VECTOR_ZERO, n + 2);
		exp2NT.getCoords().set(0, gamma);
		VectorND exp2 = NDMatrix.transform(RInverted, exp2NT);
		ArrayList<G1> e2 = VectorND.exponentiate(g1, exp2);
		
		PK.add(e1); //add first two elements to public key
		PK.add(e2);
		
		Fr gammaToTheI = new Fr(1);
		for (int i = 0; i <= v; i++) {
			VectorND exp3NT = new VectorND(VectorND.VECTOR_ZERO, n + 2);
			exp3NT.getCoords().set(0, gammaToTheI);
			VectorND exp3 = NDMatrix.transform(R, exp3NT);
			ArrayList<G2> pkJ = VectorND.exponentiate(g2, exp3);
			PK.add(pkJ);
			Mcl.mul(gammaToTheI, gammaToTheI, gamma);
		}
		
		Object[] res = {PK, MSK};
		return res;
	}
	
	public static Object[] extractMTB(Object[] MSK, Set<Integer> U, VectorND x) {
		
		
		Fr beta = (Fr) MSK[0];  //extract from MSK
		Fr gamma = (Fr) MSK[1]; 
		NDMatrix R = (NDMatrix) MSK[2];
		NDMatrix RInverted = (NDMatrix) MSK[3];
		
		//the function f will iterate over ints in U, first int that could be a valid index will be used, otherwise use 0
		//in the same loop, compute the product in the denominator of the exponent
		int index = 0;
		Fr product = new Fr(1);
		Fr ONE = new Fr(1);
		boolean firstFound = false;
		for (int s: U) {
			
			if (s < chi.size() && !firstFound) {
				index = s;
				firstFound = true;
			}
			
			Fr c = new Fr();
			Mcl.div(c, gamma, new Fr(s));	
			Mcl.sub(c, ONE, c);
			Mcl.mul(product, product, c);
		}
		
		int theta = chi.get(index);
		Fr tTheta = (Fr) MSK[theta + 3]; //get tTheta
		
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
			ArrayList<G1> element = VectorND.exponentiate(g1, exp);
			sk.add(element);
			Mcl.mul(gammaToThePower, gammaToThePower, gamma);
		}
		
		//compute the helper key
		ArrayList<Object> hkTheta = new ArrayList<Object>();
		VectorND expHVNT = new VectorND(VectorND.VECTOR_ZERO, n + 2);
		Fr e1ExpHVNT = new Fr();
		Mcl.sub(e1ExpHVNT, beta, tTheta);
		Mcl.div(e1ExpHVNT, e1ExpHVNT, gamma);
		expHVNT.getCoords().set(0, e1ExpHVNT);
		VectorND expHV = NDMatrix.transform(R, expHVNT);
		ArrayList<G2> hTheta = VectorND.exponentiate(g2, expHV);	
		hkTheta.add(hTheta);
		
		//compute all the g2 vectors
		Fr gammaToThePowerHelper = new Fr(1);
		for (int j = 0; j <= v; j++) {
			//nontransformed
			VectorND expHVENT = new VectorND(VectorND.VECTOR_ZERO, n + 2);
			Fr e1ExpHVENT = new Fr();
			Mcl.mul(e1ExpHVENT, tTheta, gammaToThePowerHelper);
			expHVENT.getCoords().set(0, e1ExpHVENT);
			//transformed
			VectorND expHVE = NDMatrix.transform(R, expHVENT);
			ArrayList<G2> elementHelper = VectorND.exponentiate(g2, expHVE);
			hkTheta.add(elementHelper);
			Mcl.mul(gammaToThePowerHelper, gammaToThePowerHelper, gamma);
		}
		
		Object[] SK = {sk, hkTheta};
		return SK;
	}
	
	//encrypt with MSK
	public static Object[] encMTB(Object[] MSK, Set<Integer> S, VectorND y) {
		
		Fr beta = (Fr) MSK[0];  //extract from MSK
		Fr gamma = (Fr) MSK[1]; 
		NDMatrix R = (NDMatrix) MSK[2];
		NDMatrix RInverted = (NDMatrix) MSK[3];
		
		//random alpha
		Fr alpha = new Fr();
		alpha.setByCSPRNG();
		
		Fr product = new Fr(1);
		Fr ONE = new Fr(1);
		for (int s: S) { //compute product in numerator
			Fr c = new Fr();
			Mcl.div(c, gamma, new Fr(s));	
			Mcl.sub(c, ONE, c);
			Mcl.mul(product, product, c);
		}
		
		//compute C1
		VectorND expC1NT = new VectorND(VectorND.VECTOR_ZERO, n + 2);
		Fr e1ExpC1NT = new Fr();
		Mcl.mul(e1ExpC1NT, alpha, gamma);
		expC1NT.getCoords().set(0, e1ExpC1NT);
		VectorND expC1 = NDMatrix.transform(RInverted, expC1NT);
		ArrayList<G1> c1 = VectorND.exponentiate(g1, expC1);
		
		//compute C2
		NDMatrix E = new NDMatrix(NDMatrix.MATRIX_DIAGONAL, n, n);
		VectorND vectorInExponent = NDMatrix.transform(E, y);
		Fr e1ExpC2 = new Fr();
		Mcl.mul(e1ExpC2, product, alpha);
		ArrayList<Fr> expC2NTElements = new ArrayList<Fr>();
		expC2NTElements.add(e1ExpC2);
		expC2NTElements.addAll(vectorInExponent.getCoords());
		expC2NTElements.add(new Fr(0));
		VectorND expC2NT = new VectorND(expC2NTElements);
		VectorND expC2 = NDMatrix.transform(R, expC2NT);
		ArrayList<G2> c2 = VectorND.exponentiate(g2, expC2);
		
		//compute encapsulated key
		GT encapsulatedKey = new GT();
		Mcl.pairing(encapsulatedKey, g1, g2);
		Mcl.pow(encapsulatedKey, encapsulatedKey, beta);
		Mcl.pow(encapsulatedKey, encapsulatedKey, alpha);
		
		Object[] c = {c1, c2};
		Object[] res = {c, encapsulatedKey};
		return res;
	}
	
	//encrypt with public key
	public static Object[] encMTB(ArrayList<Object> PK, Set<Integer> S) {
		
		//random alpha
		Fr alpha = new Fr();
		alpha.setByCSPRNG();
		
		ArrayList<G1> c1 = new ArrayList<G1>((ArrayList<G1>) PK.get(1));
		G1 e1C1 = new G1();
		Mcl.mul(e1C1, c1.get(0), alpha);
		c1.set(0, e1C1);
		
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
		ArrayList<G2> c2 = new ArrayList<G2>((ArrayList<G2>) PK.get(2));
		int k = coefficients.length;
		for (int z = 0; z < c2.size(); z++) {
			G2 helper = new G2();
			Mcl.mul(helper, c2.get(z), coefficients[k - 1]);
			c2.set(z, helper);
		}
		
		for (int i = k - 2; i >= 0; i--) {
			ArrayList<G2> element = new ArrayList<G2>((ArrayList<G2>) PK.get(k - i + 1));
			for (int y  = 0; y < c2.size(); y++) {
				G2 p = new G2();
				Mcl.mul(p, element.get(y), coefficients[i]);
				element.set(y, p);
			}
			for (int j = 0; j < element.size(); j++) {
				G2 b = new G2();
				Mcl.add(b, c2.get(j), element.get(j));
				c2.set(j, b);
			}
			
		}
		
		for (int l = 0; l < c2.size(); l++) {
			G2 b = new G2();
			Mcl.mul(b, c2.get(l), alpha);
			c2.set(l, b);
		}
		
		Object[] c = {c1, c2};
		Object[] res = {c, encapsulatedKey};
		return res;
	}
	
	public static GT decMTB(Object[] SK, Object[] MSK, Set<Integer> S, Set<Integer> U, Object[] c) {
		
		Fr gamma = (Fr) MSK[1]; 
		NDMatrix RInverted = (NDMatrix) MSK[3];
		
		//extract from the ciphertext and secret key
		ArrayList<Object> sk = (ArrayList<Object>) SK[0];
		ArrayList<Object> hkTheta = (ArrayList<Object>) SK[1];
		ArrayList<G2> hTheta = (ArrayList<G2>) hkTheta.get(0);
		ArrayList<G1> c1 = (ArrayList<G1>) c[0];
		ArrayList<G2> c2 = (ArrayList<G2>) c[1];
		
		
		Set<Integer> SDiff = new HashSet<Integer>(S);
		Set<Integer> UDiff = new HashSet<Integer>(U);
		Set<Integer> intersect = new HashSet<Integer>(U);
		
		SDiff.removeAll(U); //SDiff = S \ U
		UDiff.removeAll(S); //UDiff = U \ S
		intersect.retainAll(S);
		
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
		Fr lastCoefficient = new Fr();
		Mcl.add(lastCoefficient, PNumerator[PNumerator.length - 1], new Fr(1));
		PNumerator[PNumerator.length - 1] = lastCoefficient;
		Fr[] P = LagrangeInterpolationZp.syntheticDivideWithoutRemainder(PNumerator, new Fr(0));
		ArrayList<G2> JThetaP = new ArrayList<G2>((ArrayList<G2>) hkTheta.get(1));
		int p = P.length;
		for (int z = 0; z < JThetaP.size(); z++) {
			G2 helperJTP = new G2();
			Mcl.mul(helperJTP, JThetaP.get(z), P[p - 1]);
			JThetaP.set(z, helperJTP);
		}
		
		for (int i = p - 2; i >= 0; i--) {
			ArrayList<G2> element = new ArrayList<G2>((ArrayList<G2>) hkTheta.get(p - i));
			for (int y = 0; y < element.size(); y++) {
				G2 a = new G2();
				Mcl.mul(a, element.get(y), P[i]);
				element.set(y, a);
			}
			//add the element by the vector JThetaP
			for (int j = 0; j < element.size(); j++) {
				G2 b = new G2();
				Mcl.add(b, JThetaP.get(j), element.get(j));
				JThetaP.set(j, b);
			}
		}
		
		//compute HThetaR
		ArrayList<G1> HThetaR = new ArrayList<G1>((ArrayList<G1>) sk.get(0));
		int r = R.length;
		for (int z = 0; z < HThetaR.size(); z++)  {
			G1 helperHTR = new G1();
			Mcl.mul(helperHTR, HThetaR.get(z), R[r - 1]);
			HThetaR.set(z, helperHTR);
		}

		for (int i = r - 2; i >= 0; i--) {
			ArrayList<G1> element = new ArrayList<G1>((ArrayList<G1>) sk.get(r - i - 1));
			for (int y = 0; y < element.size(); y++) {
				G1 a = new G1();
				Mcl.mul(a, element.get(y), R[i]);
				element.set(y, a);
			}
			
			//add the element by the vector JThetaP
			for (int j = 0; j < element.size(); j++) {
				G1 b = new G1();
				Mcl.add(b, HThetaR.get(j), element.get(j));
				HThetaR.set(j, b);
			}
		}
		
		GT e1 = VectorND.vectorPairing(c1, hTheta);
		GT e2 = VectorND.vectorPairing(c1, JThetaP);
		GT e3 = VectorND.vectorPairing(HThetaR, c2);
		Mcl.mul(e1, e1, e2);
		Mcl.mul(e1, e1, e3);
		return e1;
	}
		
	
	public static void main(String[] args) {
		File lib = new File("../../lib/libmcljava.dylib");
		System.load(lib.getAbsolutePath());
		Object[] gen = genMTB(10, 100, 1, 10, Mcl.BN254);
		ArrayList<Object> PK = (ArrayList<Object>) gen[0];
		Object[] MSK = (Object[]) gen[1];
		
		Set<Integer> S = new HashSet<Integer>();
		S.addAll(Arrays.asList(1, 4, 6, 7, 10));
		Set<Integer> U = new HashSet<Integer>();
		U.addAll(Arrays.asList(1, 2, 3, 4));
		
		VectorND y = new VectorND(VectorND.VECTOR_ZERO_ONE, 10);
		//calculate attribute x (such that x dot y = 0) -- note both x & y are binary
		Fr ZERO = new Fr(0);
		ArrayList<Fr> xVals = new ArrayList<Fr>();
		for (int i = 0; i < y.getCoords().size(); i++)
			xVals.add(y.getCoords().get(i).equals(ZERO) ? new Fr(1) : new Fr(0));
		VectorND x = new VectorND(xVals);
		Object[] SK = extractMTB(MSK, U, x);
		Object[] enc = encMTB(MSK, S, y);
		Object[] C = (Object[]) enc[0];
		GT encapsulatedKey = (GT) enc[1];
		GT encapsulatedKey1 = decMTB(SK, MSK, S, U, C);
		System.out.println("Encapsulated Key Actual: " + encapsulatedKey);
		System.out.println("Decrypted Encapsulated Key: " + encapsulatedKey1);

	}
	
}


