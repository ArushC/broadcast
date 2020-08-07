package helperclasses.structures;
import java.io.File;
import java.util.ArrayList;
import miscellaneous.Tools;
import com.herumi.mcl.Fr;
import com.herumi.mcl.Mcl;

public class NDMatrix {
	
	public static final int INIT_ROW_VECTORS = 0;
	public final static int INIT_COLUMN_VECTORS = 1;
	public static final int MATRIX_ZERO = 2;
	public static final int MATRIX_DIAGONAL = 3;
	public static final int MATRIX_RANDOM = 4;
	public static final int MATRIX_IDENTITY = 5;
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
	public NDMatrix(int type, int rows, int columns) {
		this.rows = rows;
		this.columns = columns;
		rowVectors = new VectorND[this.rows];
		columnVectors = new VectorND[this.columns];
		elements = new Fr[rows][columns];
		// a zero matrix
		if (type == NDMatrix.MATRIX_ZERO) {
			for (int i = 0; i < elements.length; i++) {
				for (int j = 0; j < elements[0].length; j++) {
					elements[i][j] = new Fr(0);
				}
			}
			initVectors();
		}
		else if (type == NDMatrix.MATRIX_RANDOM) //set randomly if it is specified
			this.setByCSPRNG();
		
		else if (type == NDMatrix.MATRIX_IDENTITY) { //otherwise instantiate as an identity matrix
			if (!(rows == columns)) throw new IllegalArgumentException("ERROR: Identity matrix must be a square matrix");
			for (int i = 0; i < elements.length; i++) {
				for (int j = 0; j < elements[0].length; j++) {
					elements[i][j] = (i == j) ? new Fr(1) : new Fr(0);
				}
			}
			initVectors();
		}
		else if (type == NDMatrix.MATRIX_DIAGONAL) { //or even a diagonal matrix...
			for (int i = 0; i < elements.length; i++) {
				for (int j = 0; j < elements[0].length; j++) {
					if (i == j) {
						Fr element = new Fr();
						element.setByCSPRNG();
						elements[i][j] = element;
					}
					else
						elements[i][j] = new Fr(0);
				}
			}
			initVectors();
		}
		else throw new IllegalArgumentException("INVALID TYPE: valid types include NDMatrix.MATRIX_DEFAULT, NDMatrix.MATRIX_RANDOM, etc.");
	}
	
	//instantiate a copy of another matrix
	public NDMatrix(NDMatrix matrix) {
		Fr[][] elements = new Fr[matrix.rows][matrix.columns];
		for(int i = 0; i < matrix.rows; i++)
		    elements[i] = matrix.getElements()[i].clone();
		this.elements = elements;
		this.rows = matrix.rows;
		this.columns = matrix.columns;
		rowVectors = new VectorND[this.rows];
		columnVectors = new VectorND[this.columns];
		initVectors();
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
	
	//based on the algorithm given on this website: https://www.codeproject.com/Articles/405128/Matrix-Operations-in-Java
	//very bad runtime -- O(n!). The better one is given in determinantAndRef()
	public static Fr laPlaceDeterminant(NDMatrix matrix) {
	    if (matrix.rows != matrix.columns)
	        throw new IllegalArgumentException("ERROR: cannot calculate determinant of a nonsquare matrix");
	    if (matrix.rows == 1) {
		return matrix.getElements()[0][0];
	    } 
	    else if (matrix.rows == 2) {
	    	Fr cPos = new Fr();
	    	Mcl.mul(cPos, matrix.getElements()[0][0],matrix.getElements()[1][1]);
	    	Fr cNeg = new Fr();
	    	Mcl.sub(cNeg, matrix.getElements()[0][1], matrix.getElements()[1][0]);
	        Mcl.sub(cPos, cPos, cNeg);
	        return cPos;
	    }
	    Fr sum = new Fr(0);
	    for (int i = 0; i < matrix.columns; i++) {
	    	Fr addend = new Fr();
	    	Mcl.mul(addend, (i % 2 == 0) ? new Fr(1) : new Fr(-1) , matrix.getElements()[0][i]);
	    	Mcl.mul(addend, addend, laPlaceDeterminant(matrix.minorMatrix(0, i)));
	    	Mcl.add(sum, sum, addend);
	    }
	    return sum;
	}
	
	// function to reduce matrix to reduced row echelon form (O(n^3))
	//Modified from the code given on this website: https://www.nayuki.io/res/gauss-jordan-elimination-over-any-field/Matrix.java
	public void reducedRowEchelonForm() {
		
		// Compute row echelon form (REF)
		int numPivots = 0;
		Fr ZERO = new Fr(0);
		for (int j = 0; j < columns && numPivots < rows; j++) {  // For each column
			// Find a pivot row for this column
			int pivotRow = numPivots;
			while (pivotRow < rows && elements[pivotRow][j].equals(ZERO))
				pivotRow++;
			if (pivotRow == rows)
				continue;  // Cannot eliminate on this column
			swapRows(numPivots, pivotRow);
			pivotRow = numPivots;
			numPivots++;
			
			// Simplify the pivot row
			multiplyRow(pivotRow, Tools.reciprocal(elements[pivotRow][j]));
			
			// Eliminate rows below
			for (int i = pivotRow + 1; i < rows; i++) {
				Fr negated = new Fr();
				Mcl.mul(negated, elements[i][j], new Fr(-1));
				addRows(pivotRow, i, negated);
			}
				
		}
		
		// Compute reduced row echelon form (RREF)
		for (int i = numPivots - 1; i >= 0; i--) {
			// Find pivot
			int pivotCol = 0;
			while (pivotCol < columns && elements[i][pivotCol].equals(ZERO))
				pivotCol++;
			if (pivotCol == columns)
				continue;  // Skip this all-zero row
			
			// Eliminate rows above
			for (int j = i - 1; j >= 0; j--) {
				Fr negated = new Fr();
				Mcl.mul(negated, elements[j][pivotCol], new Fr(-1));
				addRows(i, j, negated);
			}
				
		}
	}
	
	//O(n^3)) algorithm to find determinant of a square matrix
	//Modified from the code given on this website: https://www.nayuki.io/res/gauss-jordan-elimination-over-any-field/Matrix.java
	public Fr determinantAndRef() {
		if (rows != columns)
			throw new IllegalStateException("Matrix dimensions are not square");
		
		Fr ZERO = new Fr(0);
		Fr det = new Fr(1);
		// Compute row echelon form (REF)
		int numPivots = 0;
		for (int j = 0; j < columns; j++) {  // For each column
			// Find a pivot row for this column
			int pivotRow = numPivots;
			while (pivotRow < rows && elements[pivotRow][j].equals(ZERO))
				pivotRow++;
			
			if (pivotRow < rows) {
				// This column has a nonzero pivot
				if (numPivots != pivotRow) {
					swapRows(numPivots, pivotRow);
					Mcl.mul(det, det, new Fr(-1));
				}
				pivotRow = numPivots;
				numPivots++;
				
				// Simplify the pivot row
				Fr temp = elements[pivotRow][j];
				Fr rec = Tools.reciprocal(temp);
				multiplyRow(pivotRow, rec);
				Mcl.mul(det, temp, det);
				
				// Eliminate rows below
				for (int i = pivotRow + 1; i < rows; i++)  {
					Fr negated = new Fr();
					Mcl.mul(negated, elements[i][j], new Fr(-1));
					addRows(pivotRow, i, negated);
				}
			}
			
			// Update determinant
			Mcl.mul(det, elements[j][j], det);
		}
		return det;
	}
	
	public void invert() {
		
		Fr ZERO = new Fr(0);
		Fr ONE = new Fr(1);
		
		if (rows != columns)
			throw new IllegalStateException("Matrix dimensions are not square");
		
		// Build augmented matrix: [this | identity]
		NDMatrix temp = new NDMatrix(NDMatrix.MATRIX_ZERO, rows, columns * 2);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				temp.set(i, j, elements[i][j]);
				temp.set(i, j + columns, i == j ? new Fr(1) : new Fr(0));
			}
		}
		
		// Do the main calculation
		temp.reducedRowEchelonForm();
		
		// Check that the left half is the identity matrix
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				if (!temp.getElements()[i][j].equals(i == j ? ONE : ZERO))
					throw new IllegalStateException("Matrix is not invertible");
			}
		}
		
		// Extract inverse matrix from: [identity | inverse]
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++)
				set(i, j, temp.getElements()[i][j + columns]);
		}
	}
	
	//multiplies an entire row of the matrix by the factor
	public void multiplyRow(int row, Fr factor) {
		if (row < 0 || row >= rowVectors.length)
			throw new IndexOutOfBoundsException("Row index out of bounds");
		
		VectorND v = this.rowVectors[row].multiply(factor);
		this.setRow(row, v);	
	}
	
	//destRow += srcRow * factor
	public void addRows(int srcRow, int destRow, Fr factor) {
		if (srcRow < 0 || srcRow >= rowVectors.length || destRow < 0 || destRow >= rowVectors.length)
			throw new IndexOutOfBoundsException("Row index out of bounds");
		VectorND v = this.rowVectors[srcRow].multiply(factor);
		v = v.add(this.rowVectors[destRow]);
		this.setRow(destRow, v);
	}
	
	
	//swap two rows -- used for rref and determinant calculation
	public void swapRows(int row0, int row1) {
		if (row0 < 0 || row0 >= rowVectors.length || row1 < 0 || row1 >= rowVectors.length)
			throw new IndexOutOfBoundsException("Row index out of bounds");
		VectorND temp = rowVectors[row0];
		rowVectors[row0] = rowVectors[row1];
		rowVectors[row1] = temp;
	}
	
	public NDMatrix minorMatrix(int i, int j)  {
		int c = 0, d = 0;
		Fr[][] resVals = new Fr[rows - 1][columns - 1];
		for (int a = 0; a < rows; a++) {
			if (a != i) {
				d = 0;
				for (int b = 0; b < rows; b++) {
					if (b != j) {
						resVals[c][d] = new Fr(elements[a][b]);
						d ++;
					}	
				}
				c++;
			}
		}		
		return new NDMatrix(resVals);
				
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
			throw new IllegalArgumentException("ERROR: # of columns in first matrix not equal to # of rows in second matrix");
		
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
	
	
	//Preconditions: matrix must be a square matrix and # of rows = # of elements in vector
	public static VectorND transform(NDMatrix matrix, VectorND vector) {
		if (!(matrix.rows == matrix.columns))
			throw new IllegalArgumentException("Invalid parameters for vector transformation: matrix must be square");
		else if  (!(vector.getCoords().size() == matrix.columns))
			throw new IllegalArgumentException("Invalid parameters for vector transformation: # of columns in matrix != vector size");
		
		ArrayList<Fr> res = new ArrayList<Fr>();
		for (VectorND v: matrix.getRowVectors()) {
			res.add(v.dotProduct(vector));
		}
		
		return new VectorND(res);
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
	
	public Fr[][] getElements() {
		return this.elements;
	}
	
}
