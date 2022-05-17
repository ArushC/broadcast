# Broadcast Encryption Implementations

This is a library containing implementations of broadcast encryption systems using the MCL Java pairings library. The source code for the implemented schemes can be found under the "schemes" directory. The runtimes were tested on a 2014 Macbook Air with a 1.4 GHz Dual-Core Intel Core i5 processor and 4GB RAM, and can be found under the "runtimes" directory. 

If you want to run the source code and test runtimes on your own system, follow these steps:

1. Download the MCL library from https://github.com/herumi/mcl. For the remainder of these steps, we will assume that the library was installed in a folder titled "mcl".
2. Build the library. Instructions for how to build the library are at the bottom of the page linked above, in the README doc.
3. Build the Java bindings for the library. Instructions for how to do this are at the bottom of the following page: https://github.com/herumi/mcl/tree/master/ffi/java.

IMPORTANT NOTE: by the end of steps 1-3, you should have successfully generated a .lib file titled "libmcljava.dylib". The lib file should be found at "mcl/lib". Without this library file, none of the schemes can be run.

WARNING: building this library might be difficult, depending on the architecture that you are running. I had to hardcode more than 12 modifications for the binary files to compile properly on a Macbook Pro running the M1 chip. If you are having trouble setting up the library, please don't hesitate to reach out.

4. Download this repository, and move the folders "schemes" and "helperclasses" to the directory "mcl/ffi/java".
5. Navigate to "mcl/ffi/java" in a terminal, and run the following commands:
   * javac -d . helperclasses/miscellaneous/*.java
   * javac -d . helperclasses/structures/*.java
   * javac -d . schemes/*.java
    
   Note that if you modify the source code in any of the files, you will have to rerun the above commands to see the changes reflected in the compiled code.
   
6. To test runtimes for one of the schemes, run the command "java schemes.**NAMEOFCLASS**", where **NAMEOFCLASS** is the name of any one of the files in the "schemes" directory, excluding the '.java' file extension.  

An analysis of observed runtimes can be found in the full version of the research paper based on this implementation: https://eprint.iacr.org/2021/1526.
