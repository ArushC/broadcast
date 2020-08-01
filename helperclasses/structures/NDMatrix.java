package helperclasses.structures;
import java.util.ArrayList;

import com.herumi.mcl.*;

public class NDMatrix {
	
	public static final int INIT_ROW_VECTORS = 0;
	public final static int INIT_COLUMN_VECTORS = 1;
	
	private int rows, columns;
	private Fr[][] elements;
	private VectorND[] rowVectors, columnVectors;
	
	//initalize matrix from row vectors/column vectors
	public NDMatrix(VectorND[] vectors, int mode) {
		
		if (mode == INIT_ROW_VECTORS) {
				this.rows = vectors.length;
				this.columns = vectors[0].getCoords().size();
				elements = new Fr[rows][columns];
				rowVectors = vectors.clone();
				columnVectors = new VectorND[this.columns];
				
				//initialize column vectors
				for (int i = 0; i < this.columns; i++) {
					ArrayList<Fr> columnCoords = new ArrayList<Fr>();
					for (int j = 0; j < this.rows; j++) {
						Fr element = vectors[j].getCoords().get(i);
						elements[j][i] = element;
						columnCoords.add(element);
					}
					columnVectors[i] = new VectorND(columnCoords);
				}
		}
		
		else if (mode == INIT_COLUMN_VECTORS) {
			this.rows = vectors[0].getCoords().size();
			this.columns = vectors.length;
			elements = new Fr[rows][columns];
			columnVectors = vectors.clone();
			rowVectors = new VectorND[this.rows];
			
			//initialize row vectors
			for (int i = 0; i < this.rows; i++) {
				ArrayList<Fr> rowCoords = new ArrayList<Fr>();
				for (int j = 0; j < this.columns; j++) {
					Fr element = vectors[j].getCoords().get(i);
					rowCoords.add(element);
					elements[i][j] = element;
				}
				rowVectors[i] = new VectorND(rowCoords);
			}
			
		}
		
		else {
			throw new IllegalArgumentException("INVALID MODE: mode should be either NDMatrix.INIT_ROW_VECTORS or NDMatrix.INIT_COLUMN_VECTORS");
		}	
		
	}
	
	//initialize matrix directly from elements
	public NDMatrix(Fr[][] elements) {
		
		this.elements = elements.clone();
		this.rows = elements.length;
		this.columns = elements[0].length;
		columnVectors = new VectorND[this.columns];
		rowVectors = new VectorND[this.rows];
		initVectors();
		
	}
	
	//initialize a matrix with specified dimensions
	public NDMatrix(int rows, int columns) {
		this.rows = rows;
		this.columns = columns;
		rowVectors = new VectorND[this.rows];
		columnVectors = new VectorND[this.columns];
		elements = new Fr[rows][columns];
		//set the matrix to default values (0, 0, 0, ...)
		for (int i = 0; i < elements.length; i++) {
			for (int j = 0; j < elements[0].length; j++) {
				elements[i][j] = new Fr(0);
			}
		}
		initVectors(); //initialize the vectors based on the defined elements
	}
	
	//reinitialize the elements of the matrix randomly
	public void setByCSPRNG() {
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				Fr element = new Fr();
				element.setByCSPRNG();
				elements[i][j] = element;
			}
		}
		initVectors(); //reinitialize vectors
	}
	
	//Helper function to initialize the rows/column vectors if the matrix is defined by its elements
	private void initVectors() {
		
		//initialize row vectors
		for (int i = 0; i < elements.length; i++) {
			ArrayList<Fr> rowCoords = new ArrayList<Fr>();
			for (int j = 0; j < elements[0].length; j++) {
				rowCoords.add(elements[i][j]);
			}
			rowVectors[i] = new VectorND(rowCoords);
		}
		
		//initalize column vectors
		for (int i = 0; i < elements[0].length; i++) {
			ArrayList<Fr> columnCoords = new ArrayList<Fr>();
			for (int j = 0; j < elements.length; j++) {
				columnCoords.add(elements[j][i]);
			}
			columnVectors[i] = new VectorND(columnCoords);
		}
		
	}
	
	//transpose this matrix (rows = columns)
	public NDMatrix transpose() {
		VectorND[] newRowVectors = columnVectors.clone();
		return new NDMatrix(newRowVectors, NDMatrix.INIT_ROW_VECTORS);
	}
	
	//setter methods
	public void set(int x, int y, Fr element) {
		elements[x][y] = element;
		rowVectors[x].getCoords().set(y, element);
		columnVectors[y].getCoords().set(x, element);
	}
	
	public void setRow(int x, VectorND v) {
		VectorND newVector = new VectorND(v.getCoords().size());
		//change elements
		for (int i = 0; i < this.columns; i++) {
			Fr newElement = v.getCoords().get(i);
			elements[x][i] = newElement;
			newVector.getCoords().set(i, newElement);
		}
		//make sure to change the column vectors too!
		for (int i = 0; i < columnVectors.length; i++) {
			columnVectors[i].getCoords().set(x, v.getCoords().get(i));
		}
			
		rowVectors[x] = newVector;
	}
	
	public void setColumn(int y, VectorND v) {
		VectorND newVector = new VectorND(v.getCoords().size());
		//change elements
		for (int i = 0; i < this.rows; i++) {
			Fr newElement = v.getCoords().get(i);
			elements[i][y] = newElement;
			newVector.getCoords().set(i, newElement);
		}
		
		//make sure to change the row vectors too!
		for (int i = 0; i < rowVectors.length; i++) {
			rowVectors[i].getCoords().set(y, v.getCoords().get(i));
		}
		columnVectors[y] = newVector;
	}
	
	//multiply two matrices
	public NDMatrix multiply(NDMatrix matrix) {
		
		if (!(this.columns == matrix.rows))
			throw new IllegalArgumentException("Error: # of columns in first matrix not equal to # of rows in second matrix");
		
		int newRows = this.rows, newColumns = matrix.columns;
		Fr[][] newElements = new Fr[newRows][newColumns];
		
		for (int i = 0; i < newRows; i++) {
			VectorND v1 = this.getRowVectors()[i];
			for (int j = 0; j < newColumns; j++) {
				VectorND v2 = matrix.getColumnVectors()[j];
				newElements[i][j] = v1.dotProduct(v2);
			}
		}
		
		return new NDMatrix(newElements);
	}
	
	@Override
	public String toString() {
		String result = new String();
		for (VectorND row: rowVectors)  {
			result += row.toString() + "\n";
		}
		return result;	
	}
	
	//Getter methods
	public VectorND[] getRowVectors() {
		return this.rowVectors;
	}
	
	public VectorND[] getColumnVectors() {
		return this.columnVectors;
	}
	
}
