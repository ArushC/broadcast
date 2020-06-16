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
		
		this.key = new Fr(key).toString();
		this.input = input;
		
	}
	
	
	public void setKey(Fr newKey) {
		this.key = new Fr(newKey).toString();
	}
	
	
	public Fr compute() {
		String repeatedKey = new String(new char[input]).replace("\0", key);
		G1 g1 = new G1();
		Mcl.hashAndMapToG1(g1, repeatedKey.getBytes()); 
		String g1String = g1.toString();
		int endIndex = g1String.indexOf(' ', 2); 
		Fr result = new Fr(g1String.substring(2, endIndex)); //result is defined as the first Fr that appears in the G1 string
		return result;
	}
	
	
}

