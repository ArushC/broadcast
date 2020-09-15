package schemesRevised;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.concurrent.ThreadLocalRandom;
import java.io.*;
import miscellaneous.Tools;
import helperclasses.structures.Vector2D;
import com.herumi.mcl.*;

//The scheme: https://eprint.iacr.org/2009/532.pdf (5.2, page 10)
//Changes made: added a master secret key to the scheme which contains the public parameters, in addition to:
//r1, r2, ..., r_m, c1, c2, ..., c_m, alpha1, alpha2, ..., alpha_m
//also added a key generation function to the ABBE scheme
//finally, changed the elements in g1 and g2 to obtain the most efficient scheme

public class AugBERevised {
	
	private static int m;
	
	private static int m;
	
	//ouput: public key PK, private keys SK
	public static Object[] setupABBE(int N, int lambda) {
		
		Mcl.SystemInit(lambda);
		
		//1. Create random generators and calculate m
		G1 g1 = new G1();
		Mcl.hashAndMapToG1(g1, "abc".getBytes());
	
		G2 g2 = new G2();
		Mcl.hashAndMapToG2(g2, "def".getBytes());
		
		int m = (int) Math.ceil(Math.pow(N, 0.5));
		AugBERevisedV2.m = m; //save this value so it can be used in other functions

		//2. Choose random exponents and store in arrays
		Fr[] cExponents = new Fr[m];
		Fr[] alphaExponents = new Fr[m];
		Fr[] rExponents = new Fr[m];
		
		for (int i = 0; i < m; i++) {
			Fr c = new Fr();
			Fr alpha = new Fr();
			Fr r = new Fr();
			r.setByCSPRNG();
			c.setByCSPRNG();
			alpha.setByCSPRNG();
			cExponents[i] = c;
			alphaExponents[i] = alpha;
			rExponents[i] = r;
		}
		
		
		//3. Generate the public key
		//The public key is going to be an Object[] containing three elements
		//First two elements: g1, g2
		//Third element is a (4 * m) matrix: {{E1, ..., Em}, {G1, ..., Gm}, {H1, ..., Hm}, {u1, ..., um}}
		
		Object[] PK = new Object[3];
		PK[0] = g1;
		PK[1] = g2;
		
		//Third element: a 4 * m matrix
		Object[][] pkThirdPart = new Object[4][m];
		
		byte[] b = Tools.generateRandomBytes(m + 2);
		
		//Add everything in this one gigantic for loop (repeat this single-column addition m times)
		for (int i = 0; i < m; i++) {
			
			G1 Ei = new G1(); //CHECK THIS if nothing else is the issue
			Mcl.mul(Ei, g1, rExponents[i]);
			
			GT Gi = new GT();
			Mcl.pairing(Gi, g1, g2);
			Mcl.pow(Gi, Gi, alphaExponents[i]);
		
			
			G2 Hi = new G2();
			Mcl.mul(Hi, g2, cExponents[i]);
			
			G1 ui = new G1();
			byte[] randomBytes = Arrays.copyOfRange(b, i, i + 3);
			Mcl.hashAndMapToG1(ui, randomBytes);
			
			//Finally, add all the elements to their places in the matrix
			pkThirdPart[0][i] = Ei;
			pkThirdPart[1][i] = Gi;
			pkThirdPart[2][i] = Hi;
			pkThirdPart[3][i] = ui;
		}
		
		PK[2] = pkThirdPart;
		
		//define a master secret key to contain rExponents, alphaExponents, cExponents, in addition to the public parameters
		Object[] MSK = {PK, rExponents, alphaExponents, cExponents};
		
		Object[] result = {PK, MSK};
		return result;
	}
	
	//generates the key for user u
	private static Object[] keyGenABBE(int u, Object[] MSK) {
			
		//extract from MSK
		Object[] PK = (Object[]) MSK[0];
		G1 g1 = (G1) PK[0];
		G2 g2 = (G2) PK[1];
		Object[] uValues = (Object[]) (((Object[][]) PK[2])[3]);
		Fr[] rExponents = (Fr[]) MSK[1];
		Fr[] alphaExponents = (Fr[]) MSK[2];
		Fr[] cExponents = (Fr[]) MSK[3];
			
		//1. extract x and y from u (u = (x-1)m + y), generate random delta(x, y)
		int y = (u % m == 0) ? m : (u % m);
		int x = (int) Math.ceil(((double) u) / m);
			
		Fr deltaXY = new Fr();
		deltaXY.setByCSPRNG();
			
		//2. instantiate and add the first two elements to the arraylist
		Object[] SK = new Object[m+2];
			
		G1 e1 = new G1();
		G1 part21 = new G1();
		G1 part31 = new G1();
		Mcl.mul(part21, g1, rExponents[x - 1]);
		Mcl.mul(part21, part21, cExponents[y - 1]);
		Mcl.mul(part31, (G1) uValues[y - 1], deltaXY);
		Mcl.mul(e1, g1, alphaExponents[x - 1]);
		Mcl.add(e1, e1, part21);
		Mcl.add(e1, e1, part31);
			
		G2 e2 = new G2();
		Mcl.mul(e2, g2, deltaXY);
			
		SK[0] = e1;
		SK[1] = e2;

		//3. add u1^(deltaXY), ..., (u_(y-1))^(deltaXY), (u_(y+1))^(deltaXY), ... (u_m)^(deltaXY)  //CURRENTLY DO NOT KNOW WHETHER TO ADD (u_(y))^(deltaXY)
		for (int i = 1; i <= m; i++) {
			if (i == y)
				continue;
			//else
			G1 element = new G1();
			Mcl.mul(element, (G1) uValues[i - 1], deltaXY);
			SK[i + 1] = new G1(element);
		}
			
		return SK;	
	}
	
	//input: subset S, public key PK, subset of users S, encrypt to position u = (i, j), message M in GT
	public static Object[][] encryptABBE(HashSet<Integer> S, Object[] PK, int u, GT M) {
		
		//extract data
		int j = (u % m == 0) ? m : (u % m);
		int i = (int) Math.ceil(((double) u) / m);
		
		//1. generate random exponents & vectors
		Fr t = new Fr();
		t.setByCSPRNG();
		
		Fr eta = new Fr();
		eta.setByCSPRNG();
		
		Fr[] sExponents = new Fr[m];
		Vector2D[] wVectors = new Vector2D[m];
		
		for (int k = 0; k < m; k++) {
			Fr s = new Fr();
			Vector2D w = new Vector2D();
			s.setByCSPRNG();
			w.setByCSPRNG();
			sExponents[k] = s;
			wVectors[k] = w;
		}
		
		//more vectors
		Vector2D v1 = new Vector2D();
		v1.setByCSPRNG();
		
		
		Vector2D vc = new Vector2D();
		vc.setByCSPRNG();
		
		//get random vector v2 such that v1 dot v2 = 0
		Fr randomScalar = new Fr();
		randomScalar.setByCSPRNG();
		
		Fr v2X = new Fr();
		Mcl.mul(v2X, v1.getY(), new Fr(-1));
		Vector2D v2 = new Vector2D(v2X, v1.getX()); //if v1 = (x, y), v2 = (-y, x) * randomScalar
		v2 = v2.multiply(randomScalar);
		
		//calculate vPrimeC = vc + vcr * v2
		Fr vcr = new Fr();
		vcr.setByCSPRNG();
		Vector2D product = v2.multiply(vcr);
		Vector2D vPrimeC = vc.add(product);
		
		//2. get Sx (all the y-values in S) -- helper function below
		HashSet<Integer> Sx = getSx(S);
		
		//3. add to the ciphertext using the helper functions
		Object[][] C = new Object[2][m];
		for (int k = 0; k < m; k++) //x components [0] - [m-1]
			C[0][k] = getXCiphertextComponents(k+1, i, vc, v1, t, PK, Sx, eta, sExponents, M);
		for (int k = 0; k < m; k++) //y components [m] - [2m - 1]
			C[1][k] = getYCiphertextComponents(k+1, j, PK, vc, vPrimeC, eta, t, wVectors);
		
		return C;
	}
	
	
	
	public static GT decryptABBE(Object[][] C, HashSet<Integer> S, Object[] SK, int u)     {
		
		//extract data
		int y = (u % m == 0) ? m : (u % m);
		int x = (int) Math.ceil(((double) u) / m);
		
		//1. compute KPrimeXY
		G1 product = new G1((G1) SK[0]);
		HashSet<Integer> Sx = getSx(S);
		for (int k: Sx) {
			if (k == y)
				continue;	
			//else
			Mcl.add(product, product, (G1) SK[k+1]);
		}
		G1 KPrimeXY = new G1(product);
		
		//2. extract everything that is needed
		Object[] xCiphertext = (Object[]) C[0][x - 1];
		Object[] yCiphertext = (Object[]) C[1][y-1];
		G1[] Rx = (G1[]) xCiphertext[0];
		G2 Ax = (G2) xCiphertext[1];
		G1 Tx = (G1) xCiphertext[2];
		G1[] RSquigglex = (G1[]) xCiphertext[3];
		GT Bx = (GT) xCiphertext[4];
		G2[] Cy = (G2[]) yCiphertext[0];
		G2[] CSquiggley = (G2[]) yCiphertext[1];
		G2 dDoublePrimeXY = (G2) SK[1];
		
		//3. compute the pairings
		GT numerator = new GT();
		GT e1 = computeVectorPairing(Rx, Cy);
		Mcl.mul(numerator, Bx, e1);
		
		GT denominator = new GT();
		GT e2 = computeVectorPairing(RSquigglex, CSquiggley);
		GT e3 = new GT();
		Mcl.pairing(e3, Tx, dDoublePrimeXY);
		Mcl.pow(e3, e3, new Fr(-1));
		GT e4 = new GT();
		Mcl.pairing(e4, KPrimeXY, Ax);
		Mcl.mul(denominator, e2, e3);
		Mcl.mul(denominator, denominator, e4);
		
		GT result = new GT();
		Mcl.pow(denominator, denominator, new Fr(-1));
		Mcl.mul(result, numerator, denominator);
		
		return result;
	}
	
	
	
	private static GT computeVectorPairing(G1[] g1Vals, G2[] g2Vals) {
		
		G1 g1 = g1Vals[0];
		G2 g2 = g2Vals[0]; 
		GT e1 = new GT();
		Mcl.pairing(e1, g1, g2);
		
		G1 g3 = g1Vals[1];
		G2 g4 = g2Vals[1];
		GT e2 = new GT();
		Mcl.pairing(e2, g3, g4);
		
		GT result = new GT();
		Mcl.mul(result, e1, e2);
		return result;
	}
	
	private static HashSet<Integer> getSx(HashSet<Integer> S) {
		//note the k-values returned must be UNIQUE
		HashSet<Integer> Sx = new HashSet<Integer>();
		HashSet<Integer> uniqueKValues = new HashSet<Integer>();
		for (int u: S)  {
			int y = (u % m == 0) ? m : (u % m);
			if (!(uniqueKValues.contains(y))) {
				Sx.add(y);
				uniqueKValues.add(y);
			}
				
		}
		return Sx;
	}
	
	//SEE PAGE 13 (VERY TOP): https://eprint.iacr.org/2006/298.pdf
	private static Object[] getXCiphertextComponents(int x, int i, Vector2D vc, Vector2D v1, Fr t, Object[] PK, HashSet<Integer> Sx, Fr eta, Fr[] sExponents, GT M) {
		
		G1[] Rx = new G1[2]; //G1 vector
		G2 Ax = new G2();
		G1 Tx = new G1();
		G1[] RSquigglex = new G1[2]; //another G1 vector
		GT Bx = new GT();
		
		//extract from the public key
		G1 g1 = (G1) PK[0];
		G2 g2 = (G2) PK[1];
		Object[][] pkThirdPart = (Object[][]) PK[2];
		
		//Compute the product for Tx
		G1 product = new G1();
		Mcl.mul(product, product, new Fr(0));
		for (int k: Sx) {
			G1 uk = new G1((G1) pkThirdPart[3][k-1]);
			Mcl.add(product, product, uk);
		}
		
		
		if (x < i) {
			
			//random vector zx and exponents ax and bx
			Vector2D zx = new Vector2D();
			zx.setByCSPRNG();
			
			Fr ax = new Fr();
			ax.setByCSPRNG();
			
			Fr bx = new Fr();
			bx.setByCSPRNG();
			
			G1 Rxx = new G1();
			Mcl.mul(Rxx, g1, zx.getX());
			G1 Rxy = new G1();
			Mcl.mul(Rxy, g1, zx.getY());
			Rx[0] = Rxx;
			Rx[1] = Rxy;
			
			Mcl.mul(Ax, g2, ax); 
			
			Mcl.mul(Tx, product, ax);
			
			Vector2D expRSquiggleX = zx.multiply(eta);
			G1 RSquiggleXx = new G1();
			Mcl.mul(RSquiggleXx, g1, expRSquiggleX.getX());
			G1 RSquiggleXy = new G1();
			Mcl.mul(RSquiggleXy, g1, expRSquiggleX.getY());
			RSquigglex[0] = RSquiggleXx;
			RSquigglex[1] = RSquiggleXy;
			
			Mcl.pow(Bx, (GT) pkThirdPart[1][x - 1], bx);

		}
		
		else if (x == i) {
			
			//random vector vi
			Vector2D vi = new Vector2D();
			vi.setByCSPRNG();
			
			Vector2D expRx = vi.multiply(sExponents[x - 1]);
			G1 Rxx = new G1((G1) pkThirdPart[0][x - 1]);
			Mcl.mul(Rxx, Rxx, expRx.getX());
			G1 Rxy = new G1((G1) pkThirdPart[0][x - 1]);
			Mcl.mul(Rxy, Rxy, expRx.getY());
			Rx[0] = Rxx;
			Rx[1] = Rxy;
			
			Fr expAi = vi.dotProduct(vc);
			Mcl.mul(expAi, expAi, t);
			Mcl.mul(expAi, expAi, sExponents[x - 1]);
			Mcl.mul(Ax, g2, expAi);
			
			Mcl.mul(Tx, product, expAi);
			
			Vector2D expRSquiggleX = vi.multiply(eta);
			expRSquiggleX = expRSquiggleX.multiply(sExponents[x - 1]);
			G1 RSquiggleXx = new G1((G1) pkThirdPart[0][x - 1]);
			Mcl.mul(RSquiggleXx, RSquiggleXx, expRSquiggleX.getX());
			G1 RSquiggleXy = new G1((G1) pkThirdPart[0][x - 1]);
			Mcl.mul(RSquiggleXy, RSquiggleXy, expRSquiggleX.getY());
			RSquigglex[0] = RSquiggleXx;
			RSquigglex[1] = RSquiggleXy;
			
			Mcl.pow(Bx, (GT) pkThirdPart[1][x - 1], expAi);
			Mcl.mul(Bx, M, Bx);	
			
		}
		
		else { // (x > i)
			
			//pick random vectors
			Fr vPrimeX = new Fr();
			vPrimeX.setByCSPRNG();
			
			Vector2D vx = v1.multiply(vPrimeX);
			
			Vector2D expRx = vx.multiply(sExponents[x - 1]);
			G1 Rxx = new G1((G1) pkThirdPart[0][x - 1]);
			Mcl.mul(Rxx, Rxx, expRx.getX());
			G1 Rxy = new G1((G1) pkThirdPart[0][x - 1]);
			Mcl.mul(Rxy, Rxy, expRx.getY());
			Rx[0] = Rxx;
			Rx[1] = Rxy;
			
			Fr expAi = vx.dotProduct(vc);
			Mcl.mul(expAi, expAi, t);
			Mcl.mul(expAi, expAi, sExponents[x - 1]);
			Mcl.mul(Ax, g2, expAi);
			
			Mcl.mul(Tx, product, expAi);
			
			Vector2D expRSquiggleX = vx.multiply(eta);
			expRSquiggleX = expRSquiggleX.multiply(sExponents[x - 1]);
			G1 RSquiggleXx = new G1((G1) pkThirdPart[0][x - 1]);
			Mcl.mul(RSquiggleXx, RSquiggleXx, expRSquiggleX.getX());
			G1 RSquiggleXy = new G1((G1) pkThirdPart[0][x - 1]);
			Mcl.mul(RSquiggleXy, RSquiggleXy, expRSquiggleX.getY());
			RSquigglex[0] = RSquiggleXx;
			RSquigglex[1] = RSquiggleXy;
			
			Mcl.pow(Bx, (GT) pkThirdPart[1][x - 1], expAi);
			Mcl.mul(Bx, M, Bx);	
			
		}
		
		Object[] ciphertextX = {Rx, Ax, Tx, RSquigglex, Bx};		
		return ciphertextX;	
	}
	
	private static Object[] getYCiphertextComponents(int y, int j, Object[] PK, Vector2D vc, Vector2D vPrimeC, Fr eta, Fr t, Vector2D[] wVectors) {
		
		G2[] Cy = new G2[2]; //both of these are technically 2D G2 vectors
		G2[] CSquiggley = new G2[2];
		
		//extract the public key
		G2 g2 = (G2) PK[1];
		Object[][] pkThirdPart = (Object[][]) PK[2];
		
		//compute
		Vector2D v = (y < j) ? vPrimeC : vc;
		
		Vector2D expCy1 = v.multiply(t);
		Vector2D expCy2 = wVectors[y - 1].multiply(eta);
		G2 Cyx = new G2((G2) pkThirdPart[2][y - 1]);
		G2 helperCyx = new G2();
		G2 Cyy = new G2((G2) pkThirdPart[2][y - 1]);
		G2 helperCyy = new G2(); 
		Mcl.mul(helperCyx, g2, expCy2.getX());
		Mcl.mul(Cyx, Cyx, expCy1.getX());
		Mcl.add(Cyx, Cyx, helperCyx); //The "dot" operation in the paper in this case is defined as addition?
		Mcl.mul(helperCyy, g2, expCy2.getY()); 
		Mcl.mul(Cyy, Cyy, expCy1.getY());
		Mcl.add(Cyy, Cyy, helperCyy);
		Cy[0] = Cyx;
		Cy[1] = Cyy;
		
		Vector2D expCSquiggley = wVectors[y - 1];
		G2 CSquiggleyx = new G2();
		G2 CSquiggleyy = new G2();
		Mcl.mul(CSquiggleyx, g2, expCSquiggley.getX());
		Mcl.mul(CSquiggleyy, g2, expCSquiggley.getY());
		CSquiggley[0] = CSquiggleyx;
		CSquiggley[1] = CSquiggleyy;
		
		Object[] ciphertextY = {Cy, CSquiggley};
		return ciphertextY;
	}
	
	
	//returns long[] containing elapsed time for setup, encrypt, and decrypt, respectively
		private static long[] printRuntimes(int N, int subsetSize, int lambda) {
			
			//Get elapsed time for setup(n) 
			long startSetup = System.nanoTime();
			Object[] setup = setupABBE(N, lambda);
			long elapsedSetup = System.nanoTime() - startSetup;
			double secondsSetup = ((double) elapsedSetup) / 1E9;
			
			//generate random message M in GT
			GT M = new GT();
			G1 g1 = new G1();
			Mcl.hashAndMapToG1(g1, Tools.generateRandomBytes(3));
			G2 g2 = new G2();
			Mcl.hashAndMapToG2(g2, Tools.generateRandomBytes(3));
			Mcl.pairing(M, g1, g2);
			
			//extract PK & MSK from setup and initialize S
			Object[] PK = (Object[]) setup[0];
			Object[] MSK = (Object[]) setup[1];
			HashSet<Integer> S = new HashSet<Integer>();
			
			//randomly generate numbers to put in the subset S (NO REPEATS)
			ArrayList<Integer> randomNums = new ArrayList<Integer>();
			int randomID = 0;
				for (int i = 0; i < N; i++) { 
					randomNums.add(i + 1);
				}
				for (int i = 0; i < subsetSize; i++) {
					int randomIndex = ThreadLocalRandom.current().nextInt(0, randomNums.size());
					randomID = randomNums.get(randomIndex);
					S.add(randomID);
					randomNums.remove(randomIndex);
				}
				
			long startEncrypt = System.nanoTime();
			Object[][] C = encryptABBE(S, PK, 1, M);
			long elapsedEncrypt = System.nanoTime() - startEncrypt;
			double secondsEncrypt = ((double) elapsedEncrypt) / 1E9;
					
			//Get random user ID u to test the decryption
			int u = randomID;
			
			//Get elapsed time for keygen
			long startKeyGen = System.nanoTime();
			Object[] SK = keyGenABBE(u, MSK);
			long elapsedKeyGen = System.nanoTime() - startKeyGen;
			double secondsKeyGen = ((double) elapsedKeyGen) / 1E9;
			
			//get elabsed time for decrypt
			long startDecrypt = System.nanoTime();
			GT M1 = decryptABBE(C, S, SK, u);
			long elapsedDecrypt = System.nanoTime() - startDecrypt;
			double secondsDecrypt = ((double) elapsedDecrypt) / 1E9;
			
			//Finally, print out the results
			String success = (M1.equals(M)) ? "SUCCESSFUL DECRYPTION" : "FAILED DECRYPTION";
			System.out.println(success + ": " + "N = " + N + ", subset size = " + subsetSize);
			System.out.println(); //padding
			System.out.println("setup took " + secondsSetup + " seconds");
			System.out.println("encryption took " + secondsEncrypt + " seconds");
			System.out.println("key generation took " + secondsKeyGen + " seconds");
			System.out.println("decryption took " + secondsDecrypt + " seconds (u = " + u + ")");
			System.out.println(); //more padding
			
			long[] elapsedTimes = new long[3];
			elapsedTimes[0] = elapsedSetup;
			elapsedTimes[1] = elapsedEncrypt;
			elapsedTimes[2] = elapsedDecrypt;
			
			return elapsedTimes;
		}
	
		
	//Test the runtimes of the algorithms
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
