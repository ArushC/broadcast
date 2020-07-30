
// Run using g++ --std=c++11 baseline.cpp -I/usr/local/Cellar/openssl@1.1/1.1.1g/include/ -L/usr/local/Cellar/openssl@1.1/1.1.1g/lib/ -lssl -lcrypto
//Note filepaths may vary depending on the system
//Make sure you have the latest version of openssl installed
//Please excuse the messy code -- I am not a regular C++ coder
//  Created by Arush Chhatrapati on 7/15/20.
//  Copyright © 2020 Arush Chhatrapati. All rights reserved.

#include <openssl/ec.h>
#include <openssl/err.h>
#include <openssl/obj_mac.h>
#include <openssl/objects.h>
#include <openssl/rand.h>
#include <openssl/bn.h>
#include <openssl/opensslconf.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <ctime>

using namespace std;
using namespace std::chrono;

//Use these if testing for N = 1000000 and comment out lines #149, 150, and 171
//BIGNUM *private_key[1000000];
//EC_POINT *public_key[1000001];
//EC_POINT *ciphertext[1000000];

//function to get a random point (generator) on a curve (NOT NECESSARY)
EC_POINT* get_random_generator(EC_GROUP *group, BIGNUM *order, BN_CTX *ctx) {
    
    BIGNUM *x = BN_new(), *y = BN_new();
    bool FindAnotherPoint;
    EC_POINT *g = EC_POINT_new(group);
    FindAnotherPoint = true;
    do {
          
        // generate random point
        BN_rand_range(x, order);
        EC_POINT_set_compressed_coordinates_GFp(group, g, x, 1, ctx);
        EC_POINT_get_affine_coordinates_GFp(group, g, x, y, ctx);
            // make sure point is on curve and not zero

        if(BN_is_zero(x) || BN_is_zero(y)) {
            FindAnotherPoint = true;
            continue;
        }

        if(EC_POINT_is_on_curve(group, g, ctx)) {
            FindAnotherPoint = false;
        }
    } while(FindAnotherPoint);
       
    //free x & y from memory
    BN_free(x);
    BN_free(y);
    
    return g;
}

//input: security parameter lambda, int N = number of users, and
//all of the fields that will be populated/modified by this function
//Note public key size = N + 1, private key size = N
void setup(int lambda, int N, EC_GROUP *(&groupInp), BIGNUM *(&orderInp), BN_CTX *(&ctxInp), BIGNUM* private_key[], EC_POINT* public_key[]) {
    
    //initialize group and group element g
    EC_GROUP *group = EC_GROUP_new_by_curve_name(lambda);
    BIGNUM *order = BN_new(); //order
    BN_CTX *ctx = BN_CTX_new(); //cofactor
    EC_GROUP_get_order(group, order, ctx);
    
    //random generator g
    public_key[0] = get_random_generator(group, order, ctx);
   
    //calculate x_1, x_2, ..., x_N and h1, h2, ..., h_N
    for (int i = 0; i < N; i++)
    {
        private_key[i] = BN_new();
        BN_rand_range(private_key[i], order); //set x_i randomly
        
        public_key[i+1] = EC_POINT_dup(public_key[0], group);
        EC_POINT_mul(group, public_key[i + 1], private_key[i], NULL, NULL, ctx); //h_i = g^(x_i)
    
    }
    
    //change the references that were passed into the function
    groupInp = group;
    orderInp = order;
    ctxInp = ctx;
}

//input: reference to subset of users S, message m, public_key, ciphertext (to be populated), additional parameters that were defined in the setup
//output: ciphertext
void enc(int* S, int subsetSize, EC_POINT* m, EC_POINT* public_key[], EC_POINT *ciphertext[], int N, EC_GROUP *group, BIGNUM *order, BN_CTX *ctx) {
    
    //pick a random y in Zq
    BIGNUM *y = BN_new();
    BN_rand_range(y, order);
    
    ciphertext[0] = EC_POINT_dup(public_key[0], group);
    EC_POINT_mul(group, ciphertext[0], y, NULL, NULL, ctx); //c = g^y

    for (int j = 0; j < subsetSize; j++) {
        int i = S[j];
        EC_POINT *zi = EC_POINT_dup(public_key[i], group); //z_i = h_i
        EC_POINT_mul(group, zi, y, NULL, NULL, ctx);       //z_i = z_i^y
        ciphertext[i] = EC_POINT_new(group);
        EC_POINT_add(group, ciphertext[i], zi, m, ctx); //c[i] = z_i * m
    }
    
    //free y from memory
    BN_free(y);
}


//input: ciphertext, private_key, user i, additional parameters defined in setup
void dec(EC_POINT **ciphertext, BIGNUM *private_key[], EC_POINT *result, int i, EC_GROUP *group, BIGNUM *order, BN_CTX *ctx) {
    
    
    //calculate denominator = 1/(c^x_i)
    EC_POINT *denominator = EC_POINT_dup(ciphertext[0], group);
    EC_POINT_mul(group, denominator, private_key[i - 1], NULL, NULL, ctx);
    EC_POINT_invert(group, denominator, ctx);
    
    //calculate numerator = z_i
    EC_POINT *numerator = EC_POINT_dup(ciphertext[i], group);
    
    //finally, calculate result
    EC_POINT_add(group, result, numerator, denominator, ctx);
    
    //free numerator and denominator from memory
    EC_POINT_free(numerator);
    EC_POINT_free(denominator);
}

//start clock and end clock
#define BEGIN {start = high_resolution_clock::now();}
#define END(msg) {stop = high_resolution_clock::now();\
     duration = duration_cast<microseconds>(stop - start); \
       cout << msg << " " << ((double)duration.count()/1000000) << endl; }

void printRuntimes(int lambda, int N, int subsetSize) {
    
    cout <<"N = " << N << ", subset size = " << subsetSize << endl << endl;
    srand ( time(NULL) );
    EC_GROUP *group;
    BIGNUM *order;
    BN_CTX *ctx;
    BIGNUM *private_key[N];
    EC_POINT *public_key[N + 1];
    auto start = high_resolution_clock::now();
    setup(lambda, N, group, order, ctx, private_key, public_key);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout <<"setup took: " << ((double)duration.count()/1000000) << endl;
    
    vector<int> randomNums;
    for (int i = 0; i < N; i++) {
        randomNums.push_back(i + 1);
    }
    
    int sArr[subsetSize];
    int* S = sArr;
    
    for (int i = 0; i < subsetSize; i++) {
        int randIndex = (std::rand() % randomNums.size());
        S[i] = randomNums.at(randIndex);
        randomNums.erase(randomNums.begin() + randIndex);
    }
     //subset
    EC_POINT *ciphertext[N];
    EC_POINT *m = get_random_generator(group, order, ctx);

    BEGIN
    enc(S, subsetSize, m, public_key, ciphertext, N, group, order, ctx);
    END("encrypt took: ")

    int rd = (std::rand() % subsetSize);
    int i = S[rd];
    
    EC_POINT *result = EC_POINT_new(group);
    BEGIN
    dec(ciphertext, private_key, result, i, group, order, ctx);
    END("decrypt took: ")
    
    cout <<endl; //padding
}

void testRuntimes(int lambda, int percent) {
    for (int N = 100; N <= 100000; N *= 10) {
        int subsetSize = (int) (percent * 0.01 * N);
        printRuntimes(lambda, N, subsetSize);
    }
}

int main()
{
 //    {"P-192", NID_X9_62_prime192v1},
 //    {"P-224", NID_secp224r1},
 //    {"P-256", NID_X9_62_prime256v1},
 //    {"P-384", NID_secp384r1},
 //    {"P-521", NID_secp521r1}
    cout <<"All times are in seconds" <<endl << endl;
    int lambda = NID_X9_62_prime256v1;
    testRuntimes(lambda, 10); //only tests from N = 100 to N = 100000
    //printRuntimes(lambda, 1000000, 100000); //Uncomment to test N = 1000000
                                             //and follow instructions on line #24
}
