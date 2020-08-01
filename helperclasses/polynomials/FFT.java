package helperclasses.polynomials;

import java.util.ArrayList;
import java.util.Arrays;
//THIS CODE IS BASED ON:
//Reference: https://www.cs.princeton.edu/~wayne/kleinberg-tardos/pdf/05DivideAndConquerII.pdf
//I had to modify the convolve algorithm because the way it was written, it only worked when x.length == y.length
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
 *  Reference: https://www.cs.princeton.edu/~wayne/kleinberg-tardos/pdf/05DivideAndConquerII.pdf
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
    
 // compute the linear convolution of x and y -- modified to work when x and y are not the same size
    public static Complex[] convolve(Complex[] x, Complex[] y) {
        Complex ZERO = new Complex(0, 0);
        int equalLimit = x.length + y.length - 1;
        int padLimit = equalLimit;
        while (!(isPowerOfTwo(padLimit))) {
        	padLimit += 1;
        }
        Complex[] a = new Complex[padLimit], b = new Complex[padLimit];
        //copy elements
        for (int i = 0; i < x.length; i++)
        	a[i] = x[i];
        for (int i = 0; i < y.length; i++)
        	b[i] = y[i];
        //pad with zeroes
        for (int i = x.length; i < padLimit; i++)
        	a[i] = ZERO;
        for (int i = y.length; i < padLimit; i++)
        	b[i] = ZERO;
        
        return cconvolve(a, b);
    }
    
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
        int n = 4;
        Complex[] x = new Complex[n];
        
        // original data
        for (int i = 0; i < n; i++) {
            x[i] = new Complex(i + 1, 0);
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
    
        }

}
