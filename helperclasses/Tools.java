package helperclasses;

import java.security.SecureRandom;
import com.herumi.mcl.Fr;
import com.herumi.mcl.Mcl;

public class Tools {

	
	public static Fr power(Fr base, int exponent) {
		Fr res = new Fr(1);
		for (int i = 0; i < exponent; i++) {
			Mcl.mul(res, res, base); //res = res * base (do this exponent # of times)
		}
		return res;
	}
	
	
	
	//generates a random byte array -- helper function used to obtain random generators in G 
	public static byte[] generateRandomBytes() {
					
		byte[] bytes = new byte[20];
		new SecureRandom().nextBytes(bytes);
		return bytes;
						
		}
	
	
	
}
