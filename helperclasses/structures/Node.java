package helperclasses.structures;

public class Node{
	private int datum;
	
	public Node(int datum){
		this.datum = datum;
	}
	
	public boolean equals(Node other) {
		return other.datum == this.datum;
	}
	
	public int getDatum() {
		return this.datum;
	}
	
	@Override
	public String toString() {
		return Integer.toString(datum);
	}
}
