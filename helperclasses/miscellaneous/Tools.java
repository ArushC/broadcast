package helperclasses.miscellaneous;

import java.math.BigInteger;
import java.security.SecureRandom;

import com.herumi.mcl.Fr;
import com.herumi.mcl.Mcl;

public class Tools {
	
	private static final BigInteger CURVE_ORDER_BN254 = new BigInteger("16798108731015832284940804142231733909759579603404752749028378864165570215949");
	
	public static Fr power(Fr base, int exponent) {
		
		    Fr baseCopy = new Fr(base); //don't change the value of base
		    Fr res = new Fr(1);     // Initialize result 
		    while (exponent > 0) 
		    { 
		        // If y is odd, multiply x with result 
		        if (exponent % 2 == 1) 
		            Mcl.mul(res, res, baseCopy);
		        // y must be even now 
		        exponent = exponent >> 1; // y = y/2 
		        Mcl.mul(baseCopy, baseCopy, baseCopy);  // Change x to x^2 
		    } 
		    return res;
	}
	
	//calculate modular multiplicative inverse over BN254 using the BigInteger class
	public static Fr reciprocal(Fr fr) {
		BigInteger y = new BigInteger(fr.toString());
		BigInteger t = y.modInverse(CURVE_ORDER_BN254);
		return new Fr(t.toString());		
	}
	
	
	//generates a random byte array -- helper function used to obtain random generators in G 
	public static byte[] generateRandomBytes(int numberOfBytes) {
					
		byte[] bytes = new byte[numberOfBytes];
		new SecureRandom().nextBytes(bytes);
		return bytes;
						
		}	
	
	
}
