#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <string.h>
#include <stdlib.h>

#define ARRAY_SIZE(array) (sizeof((array))/sizeof((array[0])))

/*
Rémi Tournebize in collaboration with Yves Vigouroux
remi.tournebize at ird dot fr
start: 8 june 2015
*/

/*
Only 3 possible genotypes using a BEAGLE format:
NB. Assumes only diploid individuals!
0 = maj/maj
1 = maj/min
2 = min/min
*/

/*
IMPROVEMENTS
- nPop > 2 & nPop ==1
- only biallelic sites
- minInd
- propose to reduce the SFS if high rate of missing data
- speed up the algorithm by removing geno configures if only zeros at min/min across all inds
- include Hardy Weinberg maximization
*/


using namespace std;

int main()
{

    int nPop = 2;
    int popDiploidSizes[2] = { 4, 4 };


// Compute number of diploid individuals
    int nAllDiploInd = 0;
    for ( int i = 0; i < ARRAY_SIZE(popDiploidSizes); i++ ) { nAllDiploInd += popDiploidSizes[i]; }

// Initialize the inferior bound of columns in each population
    int popBounds[ARRAY_SIZE(popDiploidSizes)];
    popBounds[0] = 0;
    int init = 0;
    for ( int i = 1; i < ARRAY_SIZE(popDiploidSizes); i++ ) {
        popBounds[i] = init + popDiploidSizes[i-1] - 1;
        init += popDiploidSizes[i-1];
    }

// Computes a matrix containing all possible genotype combinations
    int rows = 1; for ( int i = 1; i <= nAllDiploInd; i++ ) { rows *= 3; }
    int cols = nAllDiploInd;

    int genotype = 0;
    int comb[rows][cols];
    int increment = 0;
    int lapse = 0;
    for ( int j = cols-1; j >= 0 ; j-- ) {
        lapse = 0;
        genotype = 0;
        for ( int i = 0; i < rows; i++ ) {
            comb[i][j] = genotype;
            lapse++;
            if ( lapse >= increment ) {
                genotype++;
                lapse = 0;
            }
            if ( genotype > 2 ) {
                genotype = 0;
            }
        }
        if ( increment == 0 ) { increment = 1; }
        increment = increment * 3;
    }

// Compute the population-specific allelic count
// CHECK IT OUT!!!
int sum[rows][nPop];
for ( int i = 0; i < rows; i++ ) { for ( int j = 0; j < nPop; j++ ) { sum[i][j] = 0; } } // Initialize at 0
int p = 0;
for ( int i = 0; i < rows; i++ ) {
    for ( int j = 0; j < cols ; j++ ) {
        p = 0;
        if ( j > popBounds[p+1] ) {
            p++;
        }
        sum[i][p] += comb[i][j];
    }
}

// Initialize the SFS
float sfs[popDiploidSizes[0]*2+1][popDiploidSizes[1]*2+1];
for ( int i = 0; i < popDiploidSizes[0]*2+1; i++ ) {
    for ( int j = 0; j < popDiploidSizes[1]*2+1; j++ ) {
        sfs[i][j] = 0;
    }
}

// Open input file (BEAGLE format)
std::ifstream infile("C:/Users/Windows/Desktop/THESE IRD/BAM/0.ANGSD/ownANGSD/GLb.beagle");
std::string line;
std::string partial;

int index = 0;
double prod;
int geno, pos, s;
double lar[3+nAllDiploInd*3+10]; // we add 10 more elements "by security"


// Iterate through the BEAGLE file to build the 2d-SFS
while (std::getline(infile, line) && index <= 100)
{

    // read input file line-by-line
    std::istringstream iss(line);
    std::string token;
    int w = 0;
    while(std::getline(iss, token, '\t')) { // assumes a tab delimitation
        double temp = ::atof(token.c_str());
        lar[w] = temp;
        w++;
    }

    if ( index != 0 ) { // do not read the first line

        for ( int cb = 0; cb < rows; cb++ ) {

            // is major?
            s = sum[cb][0]+sum[cb][1];

            // compute the overall GL
            for ( int j = 0; j < cols; j++ ) {
                geno = comb[cb][j];
                pos = j * 3 + geno + 3; // the last element is to scale regarding the beagle formatted line (3 tokens in the beginning)
                if ( j == 0 ) {
                    prod = lar[pos];
                } else {
                    prod *= lar[pos];
                }
            }

            // add while folding
            if ( s > nAllDiploInd ) {
                sfs[nAllDiploInd-sum[cb][0]][nAllDiploInd-sum[cb][1]] += prod;
            } else {
                sfs[sum[cb][0]][sum[cb][1]] += prod;
            }
        }

    }
    index++;
}

    infile.close();

    // Output the 2d-SFS
    ofstream osfs;
    osfs.open("C:/Users/Windows/Desktop/THESE IRD/BAM/0.ANGSD/ownANGSD/foldSFS/2dsfs.txt");
            for (int k = 0; k < 9; k++)
            {
                for (int j = 0; j < 9; j++) { osfs << sfs[k][j] << " "; }
                osfs << "\n";
            }
    osfs.close();


    cout << "Done!\n";
    return 0;
}

