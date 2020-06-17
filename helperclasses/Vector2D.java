package helperclasses;
import com.herumi.mcl.*;

public class Vector2D {

	private Fr x;
	private Fr y;
	
	public Vector2D() {
		//default is (1, 1)
		this.x = new Fr(1);
		this.y = new Fr(1);
	}
	
	public Vector2D(Fr x, Fr y) {	
		
		this.x = new Fr(x);
		this.y = new Fr(y);
		
	}
	
	
	//randomly set the x & y coordinates of this vector
	public void setByCSPRNG() {
		
		Fr x = new Fr();
		x.setByCSPRNG();
		Fr y = new Fr();
		y.setByCSPRNG();
		this.x = new Fr(x);
		this.y = new Fr(y);
		
	}
	
	//scalar multiplication
	public Vector2D multiply(Fr scalar) {
		
		Fr newX = new Fr();
		Fr newY = new Fr();
		Mcl.mul(newX, scalar, this.x);
		Mcl.mul(newY, scalar, this.y);
		return new Vector2D(newX, newY);
		
	}
	
	public Fr dotProduct(Vector2D v) {
		
		Fr result = new Fr();
		Fr vx = v.getX();
		Fr vy = v.getY();
		Fr xProduct = new Fr();
		Fr yProduct = new Fr();
		Mcl.mul(xProduct, this.x, vx);
		Mcl.mul(yProduct, this.y, vy);
		Mcl.add(result, xProduct, yProduct);
		return result;
		
	}
	
	@Override
	public String toString() {
		return "( " + this.x.toString() + ", " + this.y.toString() + ")";
	}
	
	//getter methods
	public Fr getX() {
		return this.x;
	}
	
	public Fr getY() {
		return this.y;
	}
	
}
