package helperclasses;

import java.security.SecureRandom;
import com.herumi.mcl.Fr;
import com.herumi.mcl.Mcl;

public class Tools {

	
	public static Fr power(Fr base, int exponent) {
		
		    Fr res = new Fr(1);     // Initialize result 
		    while (exponent > 0) 
		    { 
		        // If y is odd, multiply x with result 
		        if (exponent % 2 == 1) 
		            Mcl.mul(res, res, base);
		        // y must be even now 
		        exponent = exponent>>1; // y = y/2 
		        Mcl.mul(base, base, base);  // Change x to x^2 
		    } 
		    return res;
	}
	
	
	
	//generates a random byte array -- helper function used to obtain random generators in G 
	public static byte[] generateRandomBytes(int numberOfBytes) {
					
		byte[] bytes = new byte[numberOfBytes];
		new SecureRandom().nextBytes(bytes);
		return bytes;
						
		}
	
	
	
}
