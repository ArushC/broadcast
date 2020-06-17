package helperclasses;

import com.herumi.mcl.Fr;
import com.herumi.mcl.G1;
import com.herumi.mcl.Mcl;

//CustomPRF class for the IBBE construction
public class CustomPRF {
	
	private String key;
	private int input;
	
	
	//input: 1. key - an Fr that has been set randomly using setByCSPRNG()
	//       2. input in [1, n]
	public CustomPRF(Fr key, int input) {
		
		Mcl.add(key, key, new Fr(input));
		this.key = key.toString();
		this.input = input;
		
	}
	
	
	public void setKey(Fr newKey) {
		Mcl.add(newKey, newKey, new Fr(input));
		this.key = newKey.toString();
	}
	
	
	public Fr compute() {
		G1 g1 = new G1();
		Mcl.hashAndMapToG1(g1, key.getBytes()); 
		String g1String = g1.toString();
		int endIndex = g1String.indexOf(' ', 2); 
		Fr result = new Fr(g1String.substring(2, endIndex)); //result is defined as the first Fr that appears in the G1 string
		return result;
	}
	
	
}
