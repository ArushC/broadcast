import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.io.*;
import com.herumi.mcl.*;



//The scheme: https://eprint.iacr.org/2006/298.pdf (pages 11-13)
public class BonehWatersABE {
	
	private static G2 gg, ggp, ggq;
	private static ArrayList<G1> uList;
	

	//ouput: public key PK, private keys SK
	public static ArrayList<Object> setupABE(int N, int lambda) {
		
		Mcl.SystemInit(lambda);
		
		//1. Create random generators in G1
		G1 gp = new G1();
		Mcl.hashAndMapToG1(gp, "abc".getBytes());
	
		G1 hp = new G1();
		Mcl.hashAndMapToG1(hp, "def".getBytes());
		
		G1 gq = new G1();
		Mcl.hashAndMapToG1(gq, "ghi".getBytes());
		
		G1 hq = new G1();
		Mcl.hashAndMapToG1(hq, "jkl".getBytes());
		
		G1 g = new G1();
		Mcl.add(g, gp, gq); //g = gp * gq (G1)
		
		G1 h = new G1();
		Mcl.add(h, hp, hq);
		
		//2. create a corresponding set of random generators in G2
		//note hhq, hhp, and hh are not needed because they will not be used
		ggq = new G2(); 
		Mcl.hashAndMapToG2(ggq, "mno".getBytes());
		
		ggp = new G2();
		Mcl.hashAndMapToG2(ggp, "pqr".getBytes());
		
		gg = new G2();
		Mcl.add(gg, ggp, ggq); //gg = ggp * ggq (G2)
		
		
		//2. Choose random elements u
		int m = (int) Math.ceil(Math.pow(N, 0.5));
		ArrayList<G1> uqList = new ArrayList<G1>();
		uList = new ArrayList<G1>(); //u = up * uq
		
		for (int i = 0; i < m; i++) {
			//Instantiate and multiply
			G1 upi = new G1();
			Mcl.hashAndMapToG1(upi, generateRandomBytes());
			G1 uqi = new G1();
			Mcl.hashAndMapToG1(uqi, generateRandomBytes());
			G1 ui = new G1();
			Mcl.add(ui, upi, uqi);
			//Add all to their lists (up list not needed)
			uqList.add(uqi);
			uList.add(ui);
			
		}
		
		//3. Choose random exponents
		ArrayList<Fr> cExponents = new ArrayList<Fr>();
		ArrayList<Fr> alphaExponents = new ArrayList<Fr>();
		ArrayList<Fr> rExponents = new ArrayList<Fr>();
		
		for (int i = 0; i < m; i++) {
			Fr c = new Fr();
			Fr alpha = new Fr();
			Fr r = new Fr();
			r.setByCSPRNG();
			c.setByCSPRNG();
			alpha.setByCSPRNG();
			cExponents.add(c);
			alphaExponents.add(alpha);
			rExponents.add(r);
		}
		
		Fr delta = new Fr();
		delta.setByCSPRNG();
		
		Fr beta = new Fr();
		beta.setByCSPRNG();
		
		Fr gamma = new Fr();
		gamma.setByCSPRNG();
		
		
		//4. Generate the public key
		//The public key is going to be an Object[] containing two elements
		//First element is a single-dimensional array containing g, h, VSquiggle, and V
		//Second element is a (9 * m) matrix containing the rest of the elements
		
		Object[] PK = new Object[2];
		
		//First element: NOTE THAT Eq is the last element in this list, even though it is on the second line
		Object[] pkFirstPart = new Object[5];
		pkFirstPart[0] = g;
		pkFirstPart[1] = h;
		G2 VSquiggle = new G2();
		G2 helper = new G2();
		Mcl.mul(helper, gg, delta);
		Mcl.mul(VSquiggle, ggp, gamma);
		Mcl.add(VSquiggle, helper, VSquiggle);
		pkFirstPart[2] = VSquiggle;
		G1 V = new G1();
		Mcl.mul(V, h, delta);
		pkFirstPart[3] = V;
		G2 Eq = new G2();
		Mcl.mul(Eq, ggq, beta);
		pkFirstPart[4] = Eq;
		
		//Second element: a 9 * m matrix
		Object[][] pkSecondPart = new Object[9][m];
		
		//Add everything in this one gigantic for loop (repeat this single-column addition m times)
		for (int i = 0; i < m; i++) {
			
			G2 Ei = new G2(); //CHECK THIS if nothing else is the issue
			Mcl.mul(Ei, gg, rExponents.get(i));
			
			G2 Eqi = new G2();
			Mcl.mul(Eqi, ggq, beta);
			Mcl.mul(Eqi, Eqi, rExponents.get(i));
			
			G1 Fi = new G1();
			Mcl.mul(Fi, h, rExponents.get(i));
			
			G1 Fqi = new G1();
			Mcl.mul(Fqi, hq, beta);
			Mcl.mul(Fqi, Fqi, rExponents.get(i));
			
			GT Gi = new GT();
			Mcl.pairing(Gi, g, gg);
			Mcl.pow(Gi, Gi, alphaExponents.get(i));
			
			GT Gqi = new GT();
			Mcl.pairing(Gqi, gq, ggq);
			Mcl.pow(Gqi, Gqi, beta);
			Mcl.pow(Gqi, Gqi, alphaExponents.get(i));
			
			G1 Hi = new G1();
			Mcl.mul(Hi, g, cExponents.get(i));
			
			G1 Ui = new G1(uList.get(i));
			
			G1 Uqi = new G1();
			Mcl.mul(Uqi, uqList.get(i), beta);
			
			//Finally, add all the elements to their places in the matrix
			pkSecondPart[0][i] = Ei;
			pkSecondPart[1][i] = Eqi;
			pkSecondPart[2][i] = Fi;
			pkSecondPart[3][i] = Fqi;
			pkSecondPart[4][i] = Gi;
			pkSecondPart[5][i] = Gqi;
			pkSecondPart[6][i] = Hi;
			pkSecondPart[7][i] = Ui;
			pkSecondPart[8][i] = Uqi;
		}
		
		PK[0] = pkFirstPart;
		PK[1] = pkSecondPart;
		
		ArrayList<Object> result = new ArrayList<Object>();
		result.add(PK);
		
		//Now add all the secret keys to result: SK_u, where u = 1, 2, ..., N (This will take O(n) time)
		for (int i = 0; i < N; i++) {
			Object[] SK = getSK(i + 1, m, g, uList, alphaExponents, rExponents, cExponents);
			result.add(SK);
		}
		
		return result;
	}
	
	private static Object[] getSK(int u, int m, G1 g, ArrayList<G1> uList, ArrayList<Fr> alphaExponents, ArrayList<Fr> rExponents, ArrayList<Fr> cExponents) {
		
		//1. extract x and y from u (u = (x-1)m + y), generate random delta(x, y)
		int y = (u % m == 0) ? m : (u % m);
		int x = (int) Math.ceil(((double) u) / m);
		
		Fr sigmaXY = new Fr();
		sigmaXY.setByCSPRNG();
		
		//2. instantiate and add the first two elements to the arraylist
		Object[] SK = new Object[m+2];
		
		G1 e1 = new G1();
		G1 part21 = new G1();
		G1 part31 = new G1();
		Mcl.mul(part21, g, rExponents.get(x - 1));
		Mcl.mul(part21, part21, cExponents.get(y - 1));
		Mcl.mul(part31, uList.get(y - 1), sigmaXY);
		Mcl.mul(e1, g, alphaExponents.get(x - 1));
		Mcl.add(e1, e1, part21);
		Mcl.add(e1, e1, part31);
		
		G2 e2 = new G2();
		Mcl.mul(e2, gg, sigmaXY);
		
		SK[0] = e1;
		SK[1] = e2;
		
		/*System.out.println("x = " + x + ", y = " + y);
		System.out.println(m);*/
		
		//3. add u1^(deltaXY), ..., (u_(y-1))^(deltaXY), (u_(y+1))^(deltaXY), ... (u_m)^(deltaXY)  //CURRENTLY DO NOT KNOW WHETHER TO ADD (u_(y))^(deltaXY)
		
		for (int i = 1; i <= m; i++) {
			if (i == y)
				continue;
			//else
			G1 element = new G1();
			Mcl.mul(element, uList.get(i - 1), sigmaXY);
			SK[i + 1] = element;
		}
		
		return SK;
		
	}
	
	
	//input: subset S, public key PK, int[][] containing rows with x-value i [0] and y-value j [1], message M in GT
	public static ArrayList<Object> encryptABE(int[][] S, Object[] PK, int[] values, GT M) {
		
		//extract data
		int m = ((Object[][]) PK[1])[0].length;
		int i = values[0], j = values[1];
		
		//1. generate random exponents
		Fr t = new Fr();
		t.setByCSPRNG();
		
		Fr kappa = new Fr();
		kappa.setByCSPRNG();
		
		ArrayList<Fr> wExponents = new ArrayList<Fr>();
		ArrayList<Fr> sExponents = new ArrayList<Fr>();
		
		for (int k = 0; k < m; k++) { //w1, ..., w_m  and s1, ..., sm
			Fr w = new Fr();
			Fr s = new Fr();
			w.setByCSPRNG();
			s.setByCSPRNG();
			wExponents.add(w);
			sExponents.add(s);
		}
		
		ArrayList<Fr> bExponents = new ArrayList<Fr>();
		
		for (int k = 1; k <= j - 1; k++) { //b1, ..., b_(j-1)
			Fr b = new Fr();
			b.setByCSPRNG();
			bExponents.add(b);
		}
		
		Fr[][] vExponents = new Fr[i-1][3]; //vExponents are defined as an (i - 1) * 3 array: see page 12 of https://eprint.iacr.org/2006/298.pdf (very strange)
		for (int k = 0; k < vExponents.length; k++) {
			for (int l = 0; l < vExponents[0].length; l++) {
				Fr v = new Fr();
				v.setByCSPRNG();
				vExponents[k][l] = v;
			}
		}
		
		//2. get Sx (all the y-values in S) -- helper function below
		ArrayList<Integer> Sx = getSx(S);
		
		//3. add to the ciphertext using more helper functions
		ArrayList<Object> C = new ArrayList<Object>();
		for (int k = 0; k < m; k++) //x components [0] - [m-1]
			C.add(getXCiphertextComponents(k+1, i, PK, Sx, kappa, t, vExponents, sExponents, M));
		for (int k = 0; k < m; k++) //y components [m] - [2m - 1]
			C.add(getYCiphertextComponents(k+1, j, PK, kappa, t, wExponents, bExponents));
		
		return C;
	}
	
	public static GT decryptABE(int[][] S, int[] userCoordinates, Object[] SK, ArrayList<Object> C, Object[] PK)     {
		
		//1. compute KPrimeXY
		G1 product = new G1((G1) SK[0]);
		int x = userCoordinates[0], y = userCoordinates[1];
		ArrayList<Integer> Sx = getSx(S);
		//System.out.println("Sx = " + Sx + ", y = " + y);
		for (int k: Sx) {
			if (k == y) {
				continue;
			}
			//else
			Mcl.add(product, product, (G1) SK[k+1]);
		}
		
		G1 KPrimeXY = new G1(product);
		
		//2. extract everything that is needed
		int m = ((Object[][]) PK[1])[0].length;
		Object[] xCiphertext = (Object[]) C.get(x-1);
		//System.out.println("m = " + m);
		Object[] yCiphertext = (Object[]) C.get(m + y - 1);
		G2 Rx = (G2) xCiphertext[0];
		G1 RSquigglex = (G1) xCiphertext[1];
		G1 Tx = (G1) xCiphertext[2];
		G2 Ax = (G2) xCiphertext[3];
		GT Bx = (GT) xCiphertext[4];
		G1 Cy = (G1) yCiphertext[0];
		G2 CSquiggley = (G2) yCiphertext[1];
		G2 dDoublePrimeXY = (G2) SK[1];
		
		//3. compute the pairings
		GT numerator = new GT();
		GT e1 = new GT();
		Mcl.pairing(e1, KPrimeXY, Ax);
		GT e2 = new GT();
		Mcl.pairing(e2, RSquigglex, CSquiggley);
		Mcl.mul(numerator, e1, e2);
		
		GT denominator = new GT();
		GT e3 = new GT();
		Mcl.pairing(e3, Cy, Rx);
		GT e4 = new GT();
		Mcl.pairing(e4, Tx, dDoublePrimeXY);
		Mcl.mul(denominator, e3, e4);
		
		GT overallDenominator = new GT();
		Mcl.pow(denominator, denominator, new Fr(-1));
		Mcl.mul(overallDenominator, numerator, denominator);
		
		GT result = new GT();
		Mcl.pow(overallDenominator, overallDenominator, new Fr(-1));
		Mcl.mul(result, Bx, overallDenominator);
		
		return result;
	}
	
	
	private static ArrayList<Integer> getSx(int[][] S) {
		ArrayList<Integer> Sx = new ArrayList<Integer>();
		for (int k = 0; k < S.length; k++) 
			Sx.add(S[k][1]);
		return Sx;
	}
	
	//SEE PAGE 13 (VERY TOP): https://eprint.iacr.org/2006/298.pdf
	private static Object[] getXCiphertextComponents(int x, int i, Object[] PK, ArrayList<Integer> Sx, Fr kappa, Fr t, Fr[][] vExponents, ArrayList<Fr> sExponents, GT M) {
		
		G2 Rx = new G2();
		G1 Tx = new G1();
		G1 RSquigglex = new G1();
		GT Bx = new GT();
		G2 Ax = new G2();
		//extract from the public key
		Object[] pkFirstPart = (Object[]) PK[0];
		Object[][] pkSecondPart = (Object[][]) PK[1];
		
		if (x > i) {
			
			Mcl.mul(Rx, (G2) pkSecondPart[1][x-1] , sExponents.get(x - 1));
			
			Mcl.mul(RSquigglex, (G1) pkSecondPart[3][x-1], kappa);
			Mcl.mul(RSquigglex, RSquigglex, sExponents.get(x - 1));
			
			Mcl.mul(Ax, (G2) pkFirstPart[4], t);
			Mcl.mul(Ax, Ax, sExponents.get(x - 1));
			
			Mcl.pow(Bx, (GT) pkSecondPart[5][x-1], t);
			Mcl.pow(Bx, Bx, sExponents.get(x - 1));
			Mcl.mul(Bx, M, Bx);
			
			G1 product = (G1) pkSecondPart[8][Sx.get(0)-1];
			for (int l = 1; l < Sx.size(); l++) {
				int k = Sx.get(l);
				G1 Uqk = (G1) pkSecondPart[8][k-1];
				Mcl.add(product, product, Uqk);
			}
			
			Mcl.mul(Tx, product, t);
			Mcl.mul(Tx, Tx, sExponents.get(x - 1));
			
		}
		
		else if (x == i) {
			
			Mcl.mul(Rx, (G2) pkSecondPart[0][x-1], sExponents.get(x - 1));
			
			Mcl.mul(RSquigglex, (G1) pkSecondPart[2][x-1], kappa);
			Mcl.mul(RSquigglex, RSquigglex, sExponents.get(x - 1));
			
			Mcl.mul(Ax, gg, t);
			Mcl.mul(Ax, Ax, sExponents.get(x - 1));
			
			Mcl.pow(Bx, (GT) pkSecondPart[4][x-1], t);
			Mcl.pow(Bx, Bx, sExponents.get(x - 1));
			Mcl.mul(Bx, M, Bx);
			
			G1 product = (G1) pkSecondPart[7][Sx.get(0)-1];
			for (int l = 1; l < Sx.size(); l++) {
				int k = Sx.get(l);
				G1 Uk = (G1) pkSecondPart[7][k-1];
				Mcl.add(product, product, Uk);
			}
			
			Mcl.mul(Tx, product, t);
			Mcl.mul(Tx, Tx, sExponents.get(x - 1));
			
		}
		
		else { // (x < i)
			
			Mcl.mul(Rx, gg, vExponents[x-1][0]);
			
			Mcl.mul(RSquigglex, (G1) pkFirstPart[1], kappa);
			Mcl.mul(RSquigglex, RSquigglex, vExponents[x-1][0]);
			
			Mcl.mul(Ax, gg, vExponents[x-1][1]);
			
			Mcl.pairing(Bx, (G1) pkFirstPart[0], gg);
			Mcl.pow(Bx, Bx, vExponents[x-1][2]);
			
			G1 product = (G1) pkSecondPart[7][Sx.get(0)-1];
			for (int l = 1; l < Sx.size(); l++) {
				int k = Sx.get(l);
				G1 Uk = (G1) pkSecondPart[7][k-1];
				Mcl.add(product, product, Uk);
			}
			
			Mcl.mul(Tx, product, vExponents[x-1][1]);
			
		}
		
		Object[] ciphertextX = new Object[5];
		ciphertextX[0] = Rx;
		ciphertextX[1] = RSquigglex;
		ciphertextX[2] = Tx;
		ciphertextX[3] = Ax;
		ciphertextX[4] = Bx;
		
		return ciphertextX;	
	}
	
	private static Object[] getYCiphertextComponents(int y, int j, Object[] PK, Fr kappa, Fr t, ArrayList<Fr> wExponents, ArrayList<Fr> bExponents) {
		
		G1 Cy = new G1();
		G2 CSquiggley = new G2();
		//extract the public key
		Object[] pkFirstPart = (Object[]) PK[0];
		Object[][] pkSecondPart = (Object[][]) PK[1];
		
		//calculate the first part of Cy and CSquiggle y
		G1 helperCy = new G1();
		Mcl.mul(helperCy, (G1) pkFirstPart[1], kappa);
		Mcl.mul(helperCy, helperCy, wExponents.get(y - 1));
		Mcl.mul(Cy, (G1) pkSecondPart[6][y-1], t);
		Mcl.add(Cy, Cy, helperCy);
		Mcl.mul(CSquiggley, gg, wExponents.get(y - 1));
		
		if (y < j) {
			
			G1 anotherHelperCy = new G1();
			Mcl.mul(anotherHelperCy, (G1) pkFirstPart[3], kappa);
			Mcl.mul(anotherHelperCy, anotherHelperCy, bExponents.get(y - 1));
			Mcl.add(Cy, Cy, anotherHelperCy);
			
			G2 helperCSquiggley = new G2();
			Mcl.mul(helperCSquiggley, (G2) pkFirstPart[2], bExponents.get(y - 1));
			Mcl.add(CSquiggley, CSquiggley, helperCSquiggley);
		}
		
		Object[] ciphertextY = new Object[2];
		ciphertextY[0] = Cy;
		ciphertextY[1] = CSquiggley;
		return ciphertextY;
	}
	
	
	//generates a random byte array -- helper function used to obtain random generators in G 
	//RUNTIME: 0.021 seconds 
	private static byte[] generateRandomBytes() {
		
		byte[] bytes = new byte[20];
		new SecureRandom().nextBytes(bytes);
		return bytes;
			
		}
	
	
	public static void main(String[] args) {
		
		System.load("/Users/arushchhatrapati/Documents/mcl/lib/libmcljava.dylib");
		//Mcl.SystemInit(Mcl.BN254);
		ArrayList<Object> setup = setupABE(100, Mcl.BN254);
		Object[] PK = (Object[]) setup.get(0);
		int m = ((Object[][]) PK[1])[0].length;
		//System.out.println("m = " + m);
		int[] values = {3, 3};
		int[][] S = {{1, 2}, {4, 5}, {6, 8}, {9, 10}, {4, 3}};
		
		GT M = new GT();
		G1 g1 = new G1();
		Mcl.hashAndMapToG1(g1, generateRandomBytes());
		G2 g2 = new G2();
		Mcl.hashAndMapToG2(g2, generateRandomBytes());
		Mcl.pairing(M, g1, g2);
		System.out.println("M = " + M);
		ArrayList<Object> C = encryptABE(S, PK, values, M); //NO ERRORS THROWN -- that's a good sign!
		//ok now let's decrypt it
		int[] userCoordinates = {4, 5};
		int u = (userCoordinates[0] - 1) * m + userCoordinates[1];
		Object[] SK = (Object[]) setup.get(u);
		GT M1 = decryptABE(S, userCoordinates, SK, C, PK);
		System.out.println("M1 = " + M1); //This doesn't work... but WHY? I think everything is written correctly (checked twice).
		//M should be EQUAL to M1
		//Is there an underlying logic error in the code? Or a typo in the paper?
		
		//System.out.println(SKu[y + 1]);
		/*for (Object o: SKu)
			System.out.println(o);
		System.out.println();
		System.out.println();
		
		for (G1 uu: uList) {
			G1 ff = new G1();
			Mcl.mul(ff, uu, test);
			System.out.println(ff);
		}*/
		
		//ALL DIFFERENT OUTPUTS
		
		/*System.out.println();
		for (int u = 1; u <= 50; u++) {
			Object[] SKu = (Object[]) setup.get(u);
			int y = (u % m == 0) ? m : (u % m);
			int x = (int) Math.ceil(((double) u) / m);
			int[] userCoordinates = {x, y};
			GT M1 = decryptABE(S, userCoordinates, SKu, C, PK);
			System.out.println("u = " + u + ", x = " + x + ", y = " + y);
			System.out.println("M1 = " + M1);
			System.out.println();
		}*/
		
		
		
		
		
		//CONFIDENT ABOUT: setup, getSK, getXCiphertextComponents, getYCiphertextComponents, 
		//NOT CONFIDENT ABOUT: NONE?!
		
	}
	
	
}
