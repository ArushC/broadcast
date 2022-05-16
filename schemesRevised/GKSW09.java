package schemesRevised;
import helperclasses.structures.Vector2D;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.concurrent.ThreadLocalRandom;
import helperclasses.miscellaneous.Tools;
import com.herumi.mcl.*;

//This is a PLBE scheme
//The original Type-II pairing construction details are outlined on pages 7-9 of  
//https://eprint.iacr.org/eprint-bin/getfile.pl?entry=2009/532&version=20091104:184423&file=532.pdf
//Modifications made to convert this to a Type-III pairing setting were described in section 3.8.2 of our paper
//Changes made: added a master secret key to the scheme which contains the public paramaters, in addition to:
//r1, r2, ..., r_m, c1, c2, ..., c_m, alpha1, alpha2, ..., alpha_m
//also added a key generation function to the PLBE scheme
public class GKSW09 {
	
	private static int m;
	
	//ouput: public key PK, private key for each user in the system
	public static Object[] setupPLBE(int N, int lambda) {
		
		Mcl.SystemInit(lambda);
		
		//1. Create random generators and calculate m
		G1 g1 = new G1();
		Mcl.hashAndMapToG1(g1, "abc".getBytes());
	
		G2 g2 = new G2();
		Mcl.hashAndMapToG2(g2, "def".getBytes());
		
		int m = (int) Math.ceil(Math.pow(N, 0.5));
		GKSW09.m = m; //save this value so it can be used in other functions

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
		//Third element is a (3 * m) matrix: {{E1, ..., Em}, {G1, ..., Gm}, {H1, ..., Hm}}
		
		Object[] PK = new Object[3];
		PK[0] = g1;
		PK[1] = g2;
		
		//Third element: a 3 * m matrix
		Object[][] pkThirdPart = new Object[3][m];
		
		GT g1g2Pairing = new GT();
		Mcl.pairing(g1g2Pairing, g1, g2);
		
		//Add everything in this one gigantic for loop (repeat this single-column addition m times)
		for (int i = 0; i < m; i++) {
			
			G1 Ei = new G1();
			Mcl.mul(Ei, g1, rExponents[i]);
			
			GT Gi = new GT(g1g2Pairing);
			Mcl.pow(Gi, Gi, alphaExponents[i]);
		
			G2 Hi = new G2();
			Mcl.mul(Hi, g2, cExponents[i]);
			
			//Finally, add all the elements to their places in the matrix
			pkThirdPart[0][i] = Ei;
			pkThirdPart[1][i] = Gi;
			pkThirdPart[2][i] = Hi;
		}
		
		PK[2] = pkThirdPart;
		
		Object[] MSK = {PK, rExponents, cExponents, alphaExponents};
		Object[] result = {PK, MSK};
		return result;
	}
	
	//generates the key for user u
	public static G2 keyGenPLBE(int u, Object[] MSK) {
		
		//extract from MSK
		Object[] PK = (Object[]) MSK[0];
		G1 g1 = (G1) PK[0];
		G2 g2 = (G2) PK[1];
		Fr[] rExponents = (Fr[]) MSK[1];
		Fr[] cExponents = (Fr[]) MSK[2];
		Fr[] alphaExponents = (Fr[]) MSK[3];
		
		int y = (u % m == 0) ? m : (u % m);
		int x = (int) Math.ceil(((double) u) / m);
		G2 kP2 = new G2();
		Mcl.mul(kP2, g2, rExponents[x - 1]);
		Mcl.mul(kP2, kP2, cExponents[y - 1]);
		G2 SK = new G2();
		Mcl.mul(SK, g2, alphaExponents[x - 1]);
		Mcl.add(SK, SK, kP2);
		return SK;
	}
	

	//input: Object[] PK, message M
	//output: _____________________
	public static Object[][] encryptPLBE(Object[] PK, GT M) {
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
		
		//3. add to the ciphertext using the helper functions
		Object[][] C = new Object[2][m];
		for (int x = 1; x <= m; x++) //x components [0] - [m-1]
			C[0][x - 1] = getXCiphertextComponents(x, vc, v1, t, PK, eta, sExponents, M);
		for (int y = 1; y <= m; y++) //y components [m] - [2m - 1]
			C[1][y - 1] = getYCiphertextComponents(y, PK, vc, eta, t, wVectors);
		
		return C;
	}

	//encrypts to recipients with a u-value > u
	public static Object[][] trEncryptPLBE(Object[] PK, int u, GT M) {
			
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
			
		//3. add to the ciphertext using the helper functions
		Object[][] C = new Object[2][m];
		for (int k = 0; k < m; k++) //x components [0] - [m-1]
			C[0][k] = getTraceableXCiphertextComponents(k+1, i, vc, v1, t, PK, eta, sExponents, M);
		for (int k = 0; k < m; k++) //y components [m] - [2m - 1]
			C[1][k] = getTraceableYCiphertextComponents(k+1, j, PK, vc, vPrimeC, eta, t, wVectors);
			
		return C;
	}
	
	public static GT decryptPLBE(Object[][] C, G2 SK, int u)   {
		
		//extract data
		int y = (u % m == 0) ? m : (u % m);
		int x = (int) Math.ceil(((double) u) / m);
		
		//2. extract everything that is needed
		Object[] xCiphertext = (Object[]) C[0][x - 1];
		Object[] yCiphertext = (Object[]) C[1][y-1];
		G1[] Rx = (G1[]) xCiphertext[0];
		G1[] RSquigglex = (G1[]) xCiphertext[1];
		G1 Ax = (G1) xCiphertext[2];
		GT Bx = (GT) xCiphertext[3];
		G2[] Cy = (G2[]) yCiphertext[0];
		G2[] CSquiggley = (G2[]) yCiphertext[1];
		
		//3. compute the pairings
		GT numerator = new GT();
		GT e1 = computeVectorPairing(Rx, Cy);
		Mcl.mul(numerator, Bx, e1);
		
		GT denominator = new GT();
		GT e2 = computeVectorPairing(RSquigglex, CSquiggley);
		GT e3 = new GT();
		Mcl.pairing(e3, Ax, SK);
		Mcl.mul(denominator, e2, e3);
		
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

	//almost exactly the same as from the AugBE system -- just remove Tx
	private static Object[] getTraceableXCiphertextComponents(int x, int i, Vector2D vc, Vector2D v1, Fr t, Object[] PK, Fr eta, Fr[] sExponents, GT M) {
			
			G1[] Rx = new G1[2]; //G1 vector
			G1 Ax = new G1();
			G1[] RSquigglex = new G1[2]; //another G1 vector
			GT Bx = new GT();
			
			//extract from the public key
			G1 g1 = (G1) PK[0];
			G2 g2 = (G2) PK[1];
			Object[][] pkThirdPart = (Object[][]) PK[2];
	
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
				
				Mcl.mul(Ax, g1, ax); 
				
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
				Mcl.mul(Ax, g1, expAi);
				
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
				Mcl.mul(Ax, g1, expAi);
				
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
			
			Object[] ciphertextX = {Rx, RSquigglex, Ax, Bx};		
			return ciphertextX;	
		}

	//exactly the same as from the augmented broadcast system
	private static Object[] getTraceableYCiphertextComponents(int y, int j, Object[] PK, Vector2D vc, Vector2D vPrimeC, Fr eta, Fr t, Vector2D[] wVectors) {
	
		G2[] Cy = new G2[2]; //both of these are technically 2D G2 vectors
		G2[] CSquiggley = new G2[2];
		
		//extract the public key
		G1 g1 = (G1) PK[0];
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
		Mcl.add(Cyx, Cyx, helperCyx);
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
	
	//non-tracing encryption
	private static Object[] getXCiphertextComponents(int x, Vector2D vc, Vector2D v1, Fr t, Object[] PK, Fr eta, Fr[] sExponents, GT M) {
		
		G1[] Rx = new G1[2]; //G1 vector
		G1 Ax = new G1();
		G1[] RSquigglex = new G1[2]; //another G1 vector
		GT Bx = new GT();
		
		//extract from the public key
		G1 g1 = (G1) PK[0];
		G2 g2 = (G2) PK[1];
		Object[][] pkThirdPart = (Object[][]) PK[2];
		
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
		Mcl.mul(Ax, g1, expAi);	
			
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
				
		Object[] ciphertextX = {Rx, RSquigglex, Ax, Bx};		
		return ciphertextX;		
	}
	
	//non-tracing encryption
	private static Object[] getYCiphertextComponents(int y, Object[] PK, Vector2D vc, Fr eta, Fr t, Vector2D[] wVectors) {
		
		G2[] Cy = new G2[2]; //both of these are technically 2D G2 vectors
		G2[] CSquiggley = new G2[2];
		
		//extract the public key
		G2 g2 = (G2) PK[1];
		Object[][] pkThirdPart = (Object[][]) PK[2];
		
		Vector2D expCy1 = vc.multiply(t);
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
		private static long[] printRuntimes(int N, int lambda) {
			
			//Get elapsed time for setup(n) 
			long startSetup = System.nanoTime();
			Object[] setup = setupPLBE(N, lambda);
			long elapsedSetup = System.nanoTime() - startSetup;
			double secondsSetup = ((double) elapsedSetup) / 1E9;
			
			//generate random message M in GT
			GT M = new GT();
			G1 g1 = new G1();
			Mcl.hashAndMapToG1(g1, Tools.generateRandomBytes(3));
			G2 g2 = new G2();
			Mcl.hashAndMapToG2(g2, Tools.generateRandomBytes(3));
			Mcl.pairing(M, g1, g2);
			
			//extract PK & MSK from setup
			Object[] PK = (Object[]) setup[0];
			Object[] MSK = (Object[]) setup[1];
			
				
			long startEncrypt = System.nanoTime();
			Object[][] C = encryptPLBE(PK, M);
			long elapsedEncrypt = System.nanoTime() - startEncrypt;
			double secondsEncrypt = ((double) elapsedEncrypt) / 1E9;
					
			//Get random user ID u to test the decryption
			int u = ThreadLocalRandom.current().nextInt(1, N + 1);
			
			//Get elapsed time for keygen
			long startKeyGen = System.nanoTime();
			G2 SK = keyGenPLBE(u, MSK);
			long elapsedKeyGen = System.nanoTime() - startKeyGen;
			double secondsKeyGen = ((double) elapsedKeyGen) / 1E9;
			
			//get elabsed time for decrypt
			long startDecrypt = System.nanoTime();
			GT M1 = decryptPLBE(C, SK, u);
			long elapsedDecrypt = System.nanoTime() - startDecrypt;
			double secondsDecrypt = ((double) elapsedDecrypt) / 1E9;
			
			//Finally, print out the results
			String success = (M1.equals(M)) ? "SUCCESSFUL DECRYPTION" : "FAILED DECRYPTION";
			System.out.println(success + ": " + "N = " + N);
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
	public static void testRuntimes(int lambda) {
		for (int N = 100; N <= 1000000; N *= 10) {
			printRuntimes(N, lambda);
		}	
	}
	
	public static void main(String[] args) {
		File lib = new File("../../lib/libmcljava.dylib");
		System.load(lib.getAbsolutePath());
		testRuntimes(Mcl.BN254);
	}
}
