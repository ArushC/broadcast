package helperclasses.polynomials;
import java.util.ArrayList;
import java.util.Arrays;

//I DID NOT WRITE MOST OF THIS. I modified the code that I found here: 
Reference: https://www.cs.princeton.edu/~wayne/kleinberg-tardos/pdf/05DivideAndConquerII.pdf

/******************************************************************************
 *  Compilation:  javac FFT.java
 *  Execution:    java FFT n
 *  Dependencies: Complex.java
 *
 *  Compute the FFT and inverse FFT of a length n complex sequence
 *  using the radix 2 Cooley-Tukey algorithm.

 *  Bare bones implementation that runs in O(n log n) time and O(n)
 *  space. Our goal is to optimize the clarity of the code, rather
 *  than performance.
 *
 *  This implementation uses the primitive root of unity w = e^(-2 pi i / n).
 *  Some resources use w = e^(2 pi i / n).
 *
 *  
 *
 *  Limitations
 *  -----------
 *   -  assumes n is a power of 2
 *
 *   -  not the most memory efficient algorithm (because it uses
 *      an object type for representing complex numbers and because
 *      it re-allocates memory for the subarray, instead of doing
 *      in-place or reusing a single temporary array)
 *  
 *  For an in-place radix 2 Cooley-Tukey FFT, see
 *  https://introcs.cs.princeton.edu/java/97data/InplaceFFT.java.html
 *
 ******************************************************************************/

public class FFT {

    // compute the FFT of x[], assuming its length n is a power of 2
    public static Complex[] fft(Complex[] x) {
        int n = x.length;

        // base case
        if (n == 1) return new Complex[] { x[0] };

        // radix 2 Cooley-Tukey FFT
        if (n % 2 != 0) {
            throw new IllegalArgumentException("n is not a power of 2");
        }

        // compute FFT of even terms
        Complex[] even = new Complex[n/2];
        for (int k = 0; k < n/2; k++) {
            even[k] = x[2*k];
        }
        Complex[] evenFFT = fft(even);

        // compute FFT of odd terms
        Complex[] odd  = even;  // reuse the array (to avoid n log n space)
        for (int k = 0; k < n/2; k++) {
            odd[k] = x[2*k + 1];
        }
        Complex[] oddFFT = fft(odd);

        // combine
        Complex[] y = new Complex[n];
        for (int k = 0; k < n/2; k++) {
            double kth = -2 * k * Math.PI / n;
            Complex wk = new Complex(Math.cos(kth), Math.sin(kth));
            y[k]       = evenFFT[k].plus (wk.times(oddFFT[k]));
            y[k + n/2] = evenFFT[k].minus(wk.times(oddFFT[k]));
        }
        return y;
    }


    // compute the inverse FFT of x[], assuming its length n is a power of 2
    public static Complex[] ifft(Complex[] x) {
        int n = x.length;
        Complex[] y = new Complex[n];

        // take conjugate
        for (int i = 0; i < n; i++) {
            y[i] = x[i].conjugate();
        }

        // compute forward FFT
        y = fft(y);

        // take conjugate again
        for (int i = 0; i < n; i++) {
            y[i] = y[i].conjugate();
        }

        // divide by n
        for (int i = 0; i < n; i++) {
            y[i] = y[i].scale(1.0 / n);
        }

        return y;

    }

    // compute the circular convolution of x and y
    public static Complex[] cconvolve(Complex[] x, Complex[] y) {

        // should probably pad x and y with 0s so that they have same length
        // and are powers of 2
        if (x.length != y.length) {
            throw new IllegalArgumentException("Dimensions don't agree");
        }

        int n = x.length;

        // compute FFT of each sequence
        Complex[] a = fft(x);
        Complex[] b = fft(y);

        // point-wise multiply
        Complex[] c = new Complex[n];
        for (int i = 0; i < n; i++) {
            c[i] = a[i].times(b[i]);
        }

        // compute inverse FFT
        return ifft(c);
    }
    
    //multiply two polynomials using a linear convolution -- the function provided in the source code did not work, so I wrote my own
    public static Complex[] convolve(Complex[] poly1, Complex[] poly2) {
    	
    	int max = Math.max(poly1.length, poly2.length);
    	while (!(isPowerOfTwo(max)))
    		max += 1;
    	
    	poly1 = zfill(poly1, max); //pad both to the appropriate number
    	poly2 = zfill(poly2, max); 
    	
    	//compute circular convolution
    	return cconvolve(poly1, poly2);
    		
    }
    
    /*public static Complex[] multiplyPoly(Polynomial... polynomials) {
    	int N = polynomials.length, M = polynomials[0].length;
    	if (N == 1) return polynomials[0];
    	if (N == 2) return multiplyPoly(polynomials[0], polynomials[1]);
    	else {
    		int n = (int) (N/2), nOp = N - n;
    		Complex[] firstHalf = new Complex[][M];
    	}
    }*/
    
    /*public static Complex[] multiplyPolyRecursive(Complex[][] f0, Complex[][] f1) {
    	int N1 = f0[0].length, N2 = f1[0].length;
    	if (f0.length == 1 && f1.length == 1) return multiplyPoly(f0[0], f1[1]);
    	if (f0.length == 1 && f1.length == 0) return f0[0];
    	if (f1.length == 1 && f0.length == 0) return f1[0];
    	else {
    		int n0 = f0.length/2, n1 = f1.length/2, nOp0 = f0.length - n1, nOp1 =  f1.length - n1;
    		Complex[][] n0Part = new Complex[n0][N1], n1Part = new Complex[n1][N1], nOp0Part = new Complex[nOp0][N2], nOp1Part = new Complex[nOp1][N2];
    		//copy all the elements to their arrays -- it is a pain, but this is the only way to do it
    		for (int i = 0; i < n0; i++)
    			n0Part[i] = f0[i];
    		for (int i = 0; i < n1; i++)
    			n1Part[i] = f1[i];
    		for (int i = 0; i < nOp0; i++) {
    			nOp0Part[i] = f0[i + n1];
    		}
    		for (int i = 0; i < nOp1; i++) {
    			nOp1Part[i] = f1[i + n1];
    		}
    		Complex[] r1 = multiplyPolyRecursive(n0Part, nOp0Part);
    		Complex[] r2 = multiplyPolyRecursive(n1Part, nOp1Part);
    	}
    }

    
    public static Complex[] multiplyPolyRecursive(Complex[][] polynomials, int low, int high) { 
    	int N = polynomials.length;
    	if (N == 1) return polynomials[0];
    	else if (N == 2) return multiplyPoly(polynomials[0], polynomials[1]);
    	else  {
    		
    	}
    }*/
    
    //check if some integer is a power of two
    public static boolean isPowerOfTwo(int n) 
    { 
        if(n==0) 
        return false; 
      
    return (int)(Math.ceil((Math.log(n) / Math.log(2)))) ==  
           (int)(Math.floor(((Math.log(n) / Math.log(2))))); 
    } 
    
    //pad an array of complexes
    public static Complex[] zfill(Complex[] x, int n) {
    	Complex[] res = new Complex[n];
    	for (int i = 0; i < x.length; i++)
    		res[i] = x[i];
    	for (int j = x.length; j < n; j++)
    		res[j] = new Complex(0, 0);
    	return res;
    }

    // compute the DFT of x[] via brute force (n^2 time)
    public static Complex[] dft(Complex[] x) {
        int n = x.length;
        Complex ZERO = new Complex(0, 0);
        Complex[] y = new Complex[n];
        for (int k = 0; k < n; k++) {
            y[k] = ZERO;
            for (int j = 0; j < n; j++) {
                int power = (k * j) % n;
                double kth = -2 * power *  Math.PI / n;
                Complex wkj = new Complex(Math.cos(kth), Math.sin(kth));
                y[k] = y[k].plus(x[j].times(wkj));
            }
        }
        return y;
    }

    // display an array of Complex numbers to standard output
    public static void show(Complex[] x, String title) {
        System.out.println(title);
        System.out.println("-------------------");
        for (int i = 0; i < x.length; i++) {
            System.out.println(x[i]);
        }
        System.out.println();
    }


   /***************************************************************************
    *  Test client and sample execution
    *
    *  % java FFT 4
    *  x
    *  -------------------
    *  -0.03480425839330703
    *  0.07910192950176387
    *  0.7233322451735928
    *  0.1659819820667019
    *
    *  y = fft(x)
    *  -------------------
    *  0.9336118983487516
    *  -0.7581365035668999 + 0.08688005256493803i
    *  0.44344407521182005
    *  -0.7581365035668999 - 0.08688005256493803i
    *
    *  z = ifft(y)
    *  -------------------
    *  -0.03480425839330703
    *  0.07910192950176387 + 2.6599344570851287E-18i
    *  0.7233322451735928
    *  0.1659819820667019 - 2.6599344570851287E-18i
    *
    *  c = cconvolve(x, x)
    *  -------------------
    *  0.5506798633981853
    *  0.23461407150576394 - 4.033186818023279E-18i
    *  -0.016542951108772352
    *  0.10288019294318276 + 4.033186818023279E-18i
    *
    *  d = convolve(x, x)
    *  -------------------
    *  0.001211336402308083 - 3.122502256758253E-17i
    *  -0.005506167987577068 - 5.058885073636224E-17i
    *  -0.044092969479563274 + 2.1934338938072244E-18i
    *  0.10288019294318276 - 3.6147323062478115E-17i
    *  0.5494685269958772 + 3.122502256758253E-17i
    *  0.240120239493341 + 4.655566391833896E-17i
    *  0.02755001837079092 - 2.1934338938072244E-18i
    *  4.01805098805014E-17i
    *
    ***************************************************************************/

    public static void main(String[] args) { 
        int n = 32;
        Complex[] x = new Complex[n];

        // original data
        for (int i = 0; i < n; i++) {
            x[i] = new Complex(i, 0);
        }
        show(x, "x");

        // FFT of original data
        Complex[] y = fft(x);
        show(y, "y = fft(x)");

        // FFT of original data
        Complex[] y2 = dft(x);
        show(y2, "y2 = dft(x)");

        // take inverse FFT
        Complex[] z = ifft(y);
        show(z, "z = ifft(y)");

        // circular convolution of x with itself
        Complex[] c = cconvolve(x, x);
        show(c, "c = cconvolve(x, x)");

        // linear convolution of x with itself
        Complex[] d = convolve(x, x);
        show(d, "d = convolve(x, x)");
        
        Complex[] f1 = {new Complex(4, 0), new Complex(5, 0), new Complex(3, 0), new Complex(1, 0), new Complex(7, 0)};
        Complex[] f2 = {new Complex(1, 0), new Complex(3, 0), new Complex(5, 0)};
        Complex[] res = convolve(f1, f2);
        show(res, "multiply f1 by f2");
        }

}
