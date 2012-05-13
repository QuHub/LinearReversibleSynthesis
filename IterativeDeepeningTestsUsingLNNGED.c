// Iterative Deepening Tests Using LNNGED Linear Reversible Circuit Synthesis 
// Copyright 2011, 2012 Ben Schaeffer
// Permission to copy this file is granted under the terms of the 
// GNU Lesser General Public License. See COPYING.LESSER.txt for details.
// Date: November 21, 2011
// Version: 0.3
//	
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
////
// Abstract: This program compares nearest neighbor linear reversible circuit
// synthesis methods with long range methods. The long range approaches
// were described in the article "Optimal Synthesis of Linear Reversal
// Circuits" by Patel, Markov, Hayes. The results represent average total
// nearest neighbor CNOT gates per synthesis of heavily randomized circuits.
//
// Can synthesize up to 64 wire circuits.
//
// Gate Encoding is based on the control wire and is described below:
// Format is CNOT(control, target)
// CNOT(0, 1) is encoded as 0, CNOT(1, 2) is encoded as 1, etc.
// CNOT(1, 0) is encoded as -1, CNOT(2, 1) is encoded as -2, etc.
//
// The synthesis output is an array of encoded gates and is 
// INVALID-terminated. Applying the output gate sequence to the 
// problem (i.e. input) circuit will change the circuit to the
// identity matrix.
//
// To solve the transposed circuit using LNN synthesis methods first
// convert the problem circuit with "TransposeCircuit", call the
// desired synthesis function, and then process the output by using
// "TransposeCNotList". After following these states the output gate
// order will be the same as the non-transposed approach, i.e. 
// applying the output gate sequence to the problem (i.e. input)
// circuit will change the circuit to the identity matrix.

//todo: clean up old
//Gaussian elimination code and experiment with 1 flag optimization,
//eliminate warnings, match names to paper, comment on reverse order

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#define FALSE (0)
#define TRUE (1)
#define TESTS_TO_RUN (100)
//This macro performs a nearest neighbor CNOT operation on circuit "a"
#define ApplyCNot(a,c) if (c < 0) a[-c-1] ^= a[-c]; else a[c+1] ^= a[c]
//#define ApplyCNot(circuit, gate) if (gate < 0) circuit[-gate - 1] ^= circuit[-gate]; else circuit[gate + 1] ^= circuit[gate]
#define INVALID (128)	//end of gate sequence marker
#define NSTART (8)
#define NEND (64)
#define NINCREMENT (8)
#define INVALID (128)
#define NMAX (64)
#define CNOTLISTSIZE (4*NMAX*NMAX)

typedef int bool;			
static uint64_t identity[NMAX];

void Display(int N, uint64_t * circuit, uint64_t cost, int gates, int depth, int cnot);
void TransposeCircuit(int N, uint64_t * source, uint64_t * destination); 
void TransposeCNotList(int * source, int * destination); //Transposes a INVALID-terminated gate list
void Randomize(int N, uint64_t * circuit); //performs 2*N*N operations on solved matrix
void CopyCircuit(int N, uint64_t * source, uint64_t * destination);
void Initialize(void);
void DisplayAlgorithm1Progress(int N, uint64_t * circuit, int gatecount, int controlrow,
		int targetrow);

// The following two functions modify the input circuit. Counts nearest neighbor gates.
int Algorithm_1_by_Patel_Markov_Hayes(int N, uint64_t * circuit, bool displayprogress); //returns gate count
int Long_Range_Gaussian_CNOT_Synthesis(int N, uint64_t * circuit, bool displayprogress); //returns gate count

// All subsequent functions make copies of the variable inputcircuit
// Linear Nearest Neighbor Gaussian Elimination using the "first achieve upper triangle matrix" approach
int LNNGE_FAUTM(int N, uint64_t * inputcircuit, int * cnotlist); //returns gate count
// Linear Nearest Neighbor Gaussian Elimination using the "first achieve lower triangle matrix" approach
int LNNGE_FALTM(int N, uint64_t * inputcircuit, int * cnotlist); //returns gate count
int LNNGED_FAUTM(int N, uint64_t * inputcircuit, int depth, int * cnotlist); //returns gate count
int LNNGED_FALTM(int N, uint64_t * inputcircuit, int depth, int * cnotlist); //returns gate count


int main(void) 
{
	uint64_t inputcircuit[NMAX], circuit[NMAX], transposedcircuit[NMAX];
	int i, d, lowestcandidate, temp, cnotlist[CNOTLISTSIZE];
	double nnDepthGaussianElimination_4runs_counts[5];
	bool displaysynthesis = FALSE;
	int N = NMAX;//Ranges from 4 to 64
	
	Initialize();
	N = 8;
	for (d = 0; d < 5; d++)
		nnDepthGaussianElimination_4runs_counts[d] = 0;
	for (i = 0;i < TESTS_TO_RUN; i++){
		Randomize(N, circuit);

		//Prepare transposed matrix for future function calls
		TransposeCircuit(N, circuit, transposedcircuit);

		for (d = 0; d < 5; d++){
			//Nearest Neighbor Gaussian Elimination with depth 1 subcolumn
			//searching during phase I, best of 4 approaches
			lowestcandidate = LNNGED_FAUTM(N, circuit, d, cnotlist);
			temp = LNNGED_FALTM(N, circuit, d, cnotlist);
			if (lowestcandidate > temp)
				lowestcandidate = temp;
			temp = LNNGED_FAUTM(N, transposedcircuit, d, cnotlist);
			if (lowestcandidate > temp)
				lowestcandidate = temp;
			temp = LNNGED_FALTM(N, transposedcircuit, d, cnotlist);
			if (lowestcandidate > temp)
				lowestcandidate = temp;
			nnDepthGaussianElimination_4runs_counts[d] += lowestcandidate;
			}		
		}
		printf("%d wire tests depth 0 through 4:\n%Lf, %Lf, %Lf, %Lf, %Lf\n", N, 
			nnDepthGaussianElimination_4runs_counts[0]/((double)TESTS_TO_RUN),
			nnDepthGaussianElimination_4runs_counts[1]/((double)TESTS_TO_RUN),
			nnDepthGaussianElimination_4runs_counts[2]/((double)TESTS_TO_RUN),
			nnDepthGaussianElimination_4runs_counts[3]/((double)TESTS_TO_RUN),
			nnDepthGaussianElimination_4runs_counts[4]/((double)TESTS_TO_RUN)
			);
	N = 16;	
	for (d = 0; d < 5; d++)
		nnDepthGaussianElimination_4runs_counts[d] = 0;
	for (i = 0;i < TESTS_TO_RUN; i++){
		Randomize(N, circuit);

		//Prepare transposed matrix for future function calls
		TransposeCircuit(N, circuit, transposedcircuit);

		for (d = 0; d < 5; d++){
			//Nearest Neighbor Gaussian Elimination with depth 1 subcolumn
			//searching during phase I, best of 4 approaches
			lowestcandidate = LNNGED_FAUTM(N, circuit, d, cnotlist);
			temp = LNNGED_FALTM(N, circuit, d, cnotlist);
			if (lowestcandidate > temp)
				lowestcandidate = temp;
			temp = LNNGED_FAUTM(N, transposedcircuit, d, cnotlist);
			if (lowestcandidate > temp)
				lowestcandidate = temp;
			temp = LNNGED_FALTM(N, transposedcircuit, d, cnotlist);
			if (lowestcandidate > temp)
				lowestcandidate = temp;
			nnDepthGaussianElimination_4runs_counts[d] += lowestcandidate;

			//printf("%d,", lowestcandidate);
			}		
			//printf("\n");
		}
		printf("%d wire tests depth 0 through 4:\n%Lf, %Lf, %Lf, %Lf, %Lf\n", N, 
			nnDepthGaussianElimination_4runs_counts[0]/((double)TESTS_TO_RUN),
			nnDepthGaussianElimination_4runs_counts[1]/((double)TESTS_TO_RUN),
			nnDepthGaussianElimination_4runs_counts[2]/((double)TESTS_TO_RUN),
			nnDepthGaussianElimination_4runs_counts[3]/((double)TESTS_TO_RUN),
			nnDepthGaussianElimination_4runs_counts[4]/((double)TESTS_TO_RUN)
			);
	
	return 0;
}

void TransposeCircuit(int N, uint64_t * source, uint64_t * destination){
	uint64_t destinationflag = 1, sourceflag;
	for (int i=0; i < N; i++)
		destination[i] = 0;
	for (int destinationcolumn = 0; destinationcolumn < N; destinationcolumn++)
	{
		sourceflag = 1;
		for (int sourcecolumn = 0; sourcecolumn < N; sourcecolumn++)
		{
			if (source[destinationcolumn] & sourceflag)
				destination[sourcecolumn] |= destinationflag;
			sourceflag<<=1;
		}
		destinationflag<<=1;
	}

}

void TransposeCNotList(int * source, int * destination){
	int lower = 0, higher = 0;
	while(source[higher] != INVALID)
		higher++;
	destination[higher] = INVALID;
	while(source[lower] != INVALID)
		destination[--higher] = -(source[lower++]+1);
}
 
void Initialize(void){
	uint64_t one = 1;
	for (int i = 0; i < NMAX; i++){
		identity[i] = one<<i;
	}
}

void Randomize(int N, uint64_t * circuit){ //performs 2*N ^ 2 operations on solved matrix
	int count = 2*N*N, x1, x2;

	for (int i = 0; i < N; i++)
		circuit[i] = identity[i];
	for (; count > 0; count--)
	{
		x1 = rand()%N;// get a number between 0 and N
		x2 = (x1 + rand()%(N-1) + 1)%N;// get a different number between 0 and N
		if(rand()%2)
		{//randomly use a cnot
			circuit[x1] ^= circuit[x2];
		}
		else
		{//randomly swap wires
			circuit[x1] ^= circuit[x2];
			circuit[x2] ^= circuit[x1];
			circuit[x1] ^= circuit[x2];
		}
	}
}

// "Algorithm 1" by Patel, Markov, Hayes, uses long-range gates
int Algorithm_1_by_Patel_Markov_Hayes(int N, uint64_t * circuit, bool displayprogress) //returns gate count
{
	uint64_t rowmask, one = 1;//handles section grouping
	int m, col, maxcol, row, targetrow, sub_row_pattern;
	int count, diagonal_one;
	count = 0;
	if (displayprogress)
		DisplayAlgorithm1Progress(N, circuit, 0, 0, 0);
	//first calculate m and row mask
	if (N < 32)
	{
	  m = 2;
	  rowmask = 0xFFFFFFFFFFFFFFFC;
	}
	else
	{
	  m = 3;
	  rowmask = 0xFFFFFFFFFFFFFFF8;
	}

	col = 0;
	while (col < N - 1) {
		maxcol = col + m; //variable to mark width of section
		if (maxcol > N - 1)
			maxcol = N - 1;
		//Step A
		for (row = col; row < N - 1; row++){
			sub_row_pattern = circuit[row] & ~rowmask;
			if (sub_row_pattern){//only search nonzero sub rows
				for (targetrow = row + 1; targetrow < N; targetrow++){
					if (sub_row_pattern == (circuit[targetrow] & ~rowmask)){
						circuit[targetrow] ^= circuit[row];
						if (targetrow == row + 1)
							count++;
						else
							count += ((targetrow - row) << 2) - 4;// cost of nearest neighbor conversion
						if (displayprogress)
							DisplayAlgorithm1Progress(N, circuit, count, row, targetrow);
					}
				}
			}
		}
		for (; col < maxcol; col++){//Step B
			if (circuit[col] & one)
				diagonal_one = TRUE;
			else
				diagonal_one = FALSE;
			for (targetrow = col + 1; targetrow < N; targetrow++){
				if (circuit[targetrow] & one){
					if (!diagonal_one) {
						diagonal_one = TRUE;
						circuit[col] ^= circuit[targetrow];
						if (targetrow == col + 1)
							count++;
						else
							count += ((targetrow - col) << 2) - 4;// cost of nearest neighbor conversion
						if (displayprogress)
							DisplayAlgorithm1Progress(N, circuit, count, targetrow, col);
					}
					//Step C
					circuit[targetrow] ^= circuit[col];
					if (targetrow == col + 1)					
						count++;
					else
						count += ((targetrow - col) << 2) - 4;// cost of nearest neighbor conversion
					if (displayprogress)
						DisplayAlgorithm1Progress(N, circuit, count, col, targetrow);
				}
			}
			//shift test flag for next iteration
			one <<= 1;
		}
		rowmask <<= m;
	}

	return count;
}

void DisplayAlgorithm1Progress(int N, uint64_t * circuit, int gatecount, int controlrow,
		int targetrow){
	uint64_t one = 1;
	int i = 0, j;
	char A='A', chr;

	if (!gatecount)
		printf ("Initial Circuit\n");
	else
		printf ("After CNOT(%d -> %d), Nearest Neighbor Gate Count = %d\n", 
				controlrow, targetrow, gatecount);
	
	for (i = 0; i < N; ++i) {
		for(j=0; j < N; j++){
			if (circuit[i] & (one<<j))
				chr = (char)j + A;
			else
				chr = ' ';
			printf("%c ", chr);//Output appropriate variable
		}
		printf("\n");//end of row
	}
	printf("\n");//end with an extra blank line
}

void CopyCircuit(int N, uint64_t * source, uint64_t * destination){
	for (int i = 0; i < N; i++)
		destination[i] = source[i];
}

int Long_Range_Gaussian_CNOT_Synthesis(int N, uint64_t * circuit, bool displayprogress){ //returns gate count
	uint64_t one = 1;//handles section grouping
	int col, targetrow;
	int count, diagonal_one;
	count = 0;
	if (displayprogress)
		DisplayAlgorithm1Progress(N, circuit, 0, 0, 0);
	//first calculate m and row mask

	col = 0;
	while (col < N - 1) {
		if (circuit[col] & one)
			diagonal_one = TRUE;
		else
			diagonal_one = FALSE;
		for (targetrow = col + 1; targetrow < N; targetrow++){
			if (circuit[targetrow] & one){
				if (!diagonal_one) {
					diagonal_one = TRUE;
					circuit[col] ^= circuit[targetrow];
					if (targetrow == col + 1)
						count++;
					else
						count += ((targetrow - col) << 2) - 4;// cost of nearest neighbor conversion
					if (displayprogress)
						DisplayAlgorithm1Progress(N, circuit, count, targetrow, col);
				}
				circuit[targetrow] ^= circuit[col];
				if (targetrow == col + 1)					
					count++;
				else
					count += ((targetrow - col) << 2) - 4;// cost of nearest neighbor conversion
				if (displayprogress)
					DisplayAlgorithm1Progress(N, circuit, count, col, targetrow);
			}
		}
		col++;
		one <<= 1;
	}

	while (col > 0) {
		for (targetrow = col - 1; targetrow >= 0; targetrow--){
			if (circuit[targetrow] & one){
				circuit[targetrow] ^= circuit[col];
				if (targetrow == col - 1)					
					count++;
				else
					count += ((-targetrow + col) << 2) - 4;// cost of nearest neighbor conversion
				if (displayprogress)
					DisplayAlgorithm1Progress(N, circuit, count, col, targetrow);
			}
		}
		col--;
		one >>= 1;
	}

	return count;
}

int LNNGE_FAUTM(int N, uint64_t * inputcircuit, int *cnotlist)
{   //returns gate count
    int totalgates = 0, cnot = INVALID, column = 0, row;
    uint64_t circuit[N];

    for (int i = 0; i < N; i++) //use copy of circuit   
        circuit[i] = inputcircuit[i];
    for (; column < N - 1; column++)    //first phase of Gaussian Elimination
        for (row = N - 1; row > column; row--)
            if (circuit[row] & ((uint64_t)1 << column))
            {
                if (circuit[row - 1] & ((uint64_t)1 << column))
                {
                    cnotlist[totalgates] = (row - 1);   //CNOT down gate                
                    totalgates++;
                    circuit[row] ^= circuit[row - 1];
                }
                else
                {
                    cnotlist[totalgates] = -(row);  //CNOT up gate
                    totalgates++;
                    cnotlist[totalgates] = (row - 1);   //CNOT down gate
                    totalgates++;
                    circuit[row - 1] ^= circuit[row];
                    circuit[row] ^= circuit[row - 1];
                }
            }
    for (; column > 0; column--)
    {   //second phase of Gaussian Elimination
        for (row = 0; row < column && !(circuit[row] & ((uint64_t)1 << column)); row++)
            ;   //search for top instance of variable associated with the column
        if (row != column)
        {   //First extend "1"s up
            for (int rowhelper = column; rowhelper - 1 > row; rowhelper--)
            {
                if (!(circuit[rowhelper - 1] & ((uint64_t)1 << column)))
                {
                    cnotlist[totalgates] = -(rowhelper);    //CNOT up gate
                    totalgates++;
                    circuit[rowhelper - 1] ^= circuit[rowhelper];
                }
            }
            for (; ++row <= column;)
            {   //Next eliminate continuing down
                cnotlist[totalgates] = -(row);  //CNOT up gate
                totalgates++;
                circuit[row - 1] ^= circuit[row];
            }
        }
    }
    cnotlist[totalgates] = INVALID;
    return totalgates;
}

int LNNGE_FALTM(int N, uint64_t * inputcircuit, int * cnotlist){//returns gate count
	int totalgates = 0;
	int lastpenalty;
	int cnot = INVALID;
	
	//use copy of circuit
	uint64_t circuit[N];
	for (int i = 0; i < N; i++)
		circuit[i] = inputcircuit[i];

	//lastpenalty = CalculateCostHeuristic(circuit);
	//Display(circuit, lastpenalty, totalgates, 0, cnot);
	if (1){
	
	// first look for a lower triangle CNOT gate, then upper triangle
		int column = N-1, row;
		for (;column > 0; column--)	
			for (row = 0; row < column; row++)
				if(circuit[row]&((uint64_t)1<<column)){
					//now check if for another instance of this variable
					//on the row above, necessitating a CNOT down gate
					if(circuit[row+1]&((uint64_t)1<<column)){
						cnotlist[totalgates] = -(row+1);					
						totalgates++;
						circuit[row] ^= circuit[row+1];						
					}
					else {
						//CNOT down gate needed, then up
						cnotlist[totalgates] = row;
						totalgates++;
						cnotlist[totalgates] = -(row+1);					
						totalgates++;
						circuit[row+1] ^= circuit[row];
						circuit[row] ^= circuit[row+1];
					}
				//	Display(circuit, lastpenalty, totalgates, 0, row);
				} 
		//now examine upper triangle starting from column N-1
		//due to the fact that the lower triangle is solved only CNOT
		//up are permitted, otherwise the gates might interfere with the
		//solved entries
		for (;column < N-1; column++)	{
			for (row = N-1; row > column && !(circuit[row]&((uint64_t)1<<column)); row--)
				;
			if (row!= column)
			{
				//First grow down
				for (int rowhelper = column; rowhelper + 1 < row; rowhelper++) {

					if (!(circuit[rowhelper + 1]&((uint64_t)1<<column))){
						cnotlist[totalgates] = rowhelper;
						totalgates++;
						circuit[rowhelper + 1] ^= circuit[rowhelper];
					}
				}
				//Next grow up
				for (;--row >= column;) {
						cnotlist[totalgates] = row;
						totalgates++;
						circuit[row+1] ^= circuit[row];
				//Display(circuit, lastpenalty, totalgates, 0, row);
				}
				
			}
			
		}

		//lastpenalty = CalculateCostHeuristic(circuit);
		//Display(circuit, lastpenalty, totalgates, 0, cnot);
	}
	cnotlist[totalgates] = INVALID;
	return totalgates;
}


int LNNGED_FAUTM(int N, uint64_t * inputcircuit, int depth, int * cnotlist){//returns gate count
	if (depth == 0) 
		return LNNGE_FAUTM(N, inputcircuit, cnotlist);
	int totalgates = 0;
	int lastpenalty;
	int cnot = INVALID;
	//use copy of circuit
	uint64_t circuit[N];
	for (int i = 0; i < N; i++)
		circuit[i] = inputcircuit[i];

	//lastpenalty = CalculateCostHeuristic(circuit);
	//Display(circuit, lastpenalty, totalgates, 0, cnot);
	if (1){

	// first look for a lower triangle CNOT gate, then upper triangle
		int column = 0, row;
		for (;column < N-1; column++)	
			for (row = N-1; row > column; row--)
				if(circuit[row]&((uint64_t)1<<column)){
					//now check if for another instance of this variable
					//on the row above, necessitating a CNOT down gate
					if(circuit[row-1]&((uint64_t)1<<column)){
						cnotlist[totalgates] = (row-1);
						//printf(" just now added %d, %d\n", cnotlist[totalgates], totalgates);
						totalgates++;
						circuit[row] ^= circuit[row-1];						
					}
					else {
						//Check for above instance
						int rowabove = row - 2, instancefound = FALSE;
						while (!instancefound && rowabove >= column) 
							if (circuit[rowabove]&((uint64_t)1<<column)) 
								instancefound = TRUE;							
							else
								rowabove--;
						if (!instancefound) {
							cnotlist[totalgates] = -(row);
						//printf(" just now added %d, %d\n", cnotlist[totalgates], totalgates);
							totalgates++;
							cnotlist[totalgates] = (row-1);
						//printf(" just now added %d, %d\n", cnotlist[totalgates], totalgates);
							totalgates++;
							circuit[row-1] ^= circuit[row];
							circuit[row] ^= circuit[row-1];
						}
						else {//choose best heuristic and adjust row
							int minimumheuristic, mincnotdown = 0, cnotdown = 0, rown, cd, temp;
							//first set minimumheuristic to all CNot up
							for (rown = row; rown - 1> rowabove; rown--)
								circuit[rown-1] ^= circuit[rown];
							//p= cnotlist+totalgates;
							lastpenalty = minimumheuristic = LNNGED_FAUTM(N, circuit, depth - 1, cnotlist+totalgates); 
							//restore circuit
							for (rown++; rown <= row; rown++)
								circuit[rown-1] ^= circuit[rown];

							//compare against rest
							for (cnotdown = 1; cnotdown < row - rowabove; cnotdown++){
								// compute up and down gates, find cost, and restore
								for (cd = 0; cd < cnotdown; cd++)
									circuit[rowabove + 1 + cd] ^= circuit[rowabove + cd];
								for (rown = row; rown - 1> rowabove+cnotdown; rown--)
									circuit[rown-1] ^= circuit[rown];
								temp = LNNGED_FAUTM(N, circuit, depth - 1, cnotlist+totalgates);
								if (temp < minimumheuristic){
									lastpenalty = minimumheuristic = temp;
									mincnotdown = cnotdown;
								}
								//restore circuit
								for (rown++; rown <= row; rown++)
									circuit[rown-1] ^= circuit[rown];
								for (cd--; cd >= 0; cd--)
									circuit[rowabove + 1 + cd] ^= circuit[rowabove + cd];
							}

							//choose best
							cnotdown = mincnotdown;
							for (cd = 0; cd < cnotdown; cd++) {
								cnotlist[totalgates] = (rowabove + cd);
						//printf(" just now added %d, %d\n", cnotlist[totalgates], totalgates);
								totalgates++;
								circuit[rowabove + 1 + cd] ^= circuit[rowabove + cd];
							}
							for (rown = row; rown - 1> rowabove+cnotdown; rown--){
								cnotlist[totalgates] = -(rown);
						//printf(" just now added %d, %d\n", cnotlist[totalgates], totalgates);
								totalgates++;
								circuit[rown-1] ^= circuit[rown];
							}

							row++;//adjustment so row calculation starts over
						}
					}
					//Display(circuit, lastpenalty, totalgates, 0, row);
				} 
		//now examine upper triangle starting from column N-1
		//due to the fact that the lower triangle is solved only CNOT
		//up are permitted, otherwise the gates might interfere with the
		//solved entries
 		for (;column > 0; column--)	{
			for (row = 0; row < column && !(circuit[row]&((uint64_t)1<<column)); row++)
				;
			if (row!= column)
			{
				//First grow up
				for (int rowhelper = column; rowhelper - 1 > row; rowhelper--) {

					if (!(circuit[rowhelper - 1]&((uint64_t)1<<column))){
						cnotlist[totalgates] = -(rowhelper);
						//printf(" just now added %d, %d\n", cnotlist[totalgates], totalgates);
						totalgates++;
						circuit[rowhelper - 1] ^= circuit[rowhelper];
					}
				}
				//Next grow down
				for (;++row <= column;) {
						cnotlist[totalgates] = -(row);
						//printf(" just now added %d, %d\n", cnotlist[totalgates], totalgates);
						totalgates++;
						circuit[row-1] ^= circuit[row];
				//Display(circuit, lastpenalty, totalgates, 0, row);
				}
				
			}
		}

		//lastpenalty = CalculateCostHeuristic(circuit);
		//Display(circuit, lastpenalty, totalgates, 0, cnot);
	}
	cnotlist[totalgates] = INVALID;
	return totalgates;
}

int LNNGED_FALTM(int N, uint64_t * inputcircuit, int depth, int * cnotlist){//returns gate count
	if (depth == 0) return LNNGE_FALTM(N, inputcircuit, cnotlist);
	int totalgates = 0;
	int lastpenalty;
	int cnot = INVALID;
	//use copy of circuit
	uint64_t circuit[N];
	for (int i = 0; i < N; i++)
		circuit[i] = inputcircuit[i];

	//lastpenalty = CalculateCostHeuristic(circuit);
	//Display(circuit, lastpenalty, totalgates, 0, cnot);
	if (1){
	
	// first look for a lower triangle CNOT gate, then upper triangle
		int column = N-1, row;
		for (;column > 0; column--)	
			for (row = 0; row < column; row++)
				if(circuit[row]&((uint64_t)1<<column)){
					//now check if for another instance of this variable
					//on the row above, necessitating a CNOT down gate
					if(circuit[row+1]&((uint64_t)1<<column)){
						cnotlist[totalgates] = -(row+1);
						totalgates++;
						circuit[row] ^= circuit[row+1];						
					}
					else {
						//Check for above instance
						int rowabove = row + 2, instancefound = FALSE;
						while (!instancefound && rowabove <= column) 
							if (circuit[rowabove]&((uint64_t)1<<column)) 
								instancefound = TRUE;							
							else
								rowabove++;
						if (!instancefound) {
							cnotlist[totalgates] = row;
							totalgates++;
							cnotlist[totalgates] = -(row+1);
							totalgates++;
							circuit[row+1] ^= circuit[row];
							circuit[row] ^= circuit[row+1];
						}
						else {//choose best heuristic and adjust row
							int minimumheuristic, mincnotdown = 0, cnotdown = 0, rown, cd, temp;
							//first set minimumheuristic to all CNot up
							for (rown = row; rown + 1< rowabove; rown++)
								circuit[rown+1] ^= circuit[rown];
							lastpenalty = minimumheuristic = LNNGED_FALTM(N, circuit, depth - 1, cnotlist+totalgates);
							//restore circuit
							for (rown--; rown >= row; rown--)
								circuit[rown+1] ^= circuit[rown];

							//compare against rest
							for (cnotdown = 1; cnotdown < rowabove - row; cnotdown++){
								// compute up and down gates, find cost, and restore		
								for (cd = 0; cd < cnotdown; cd++)
									circuit[rowabove - 1 - cd] ^= circuit[rowabove - cd];
								for (rown = row; rown+1<rowabove-cnotdown; rown++)
									circuit[rown+1] ^= circuit[rown];
								temp = LNNGED_FALTM(N, circuit, depth - 1, cnotlist+totalgates);
								if (temp < minimumheuristic){
									lastpenalty = minimumheuristic = temp;
									mincnotdown = cnotdown;
								}
								//restore circuit
								for (rown--; rown >= row; rown--)
									circuit[rown+1] ^= circuit[rown];
								for (cd--; cd >= 0; cd--)
									circuit[rowabove - 1 - cd] ^= circuit[rowabove - cd];
							}

							//choose best
							cnotdown = mincnotdown;
							for (cd = 0; cd < cnotdown; cd++) {
								cnotlist[totalgates] = -(rowabove - cd);
								totalgates++;
								circuit[rowabove - 1 - cd] ^= circuit[rowabove - cd];
							}
							for (rown = row; rown + 1< rowabove-cnotdown; rown++){
								cnotlist[totalgates] = rown;
								totalgates++;
								circuit[rown+1] ^= circuit[rown];
							}

							row--;//adjustment so row calculation starts over
						}
					}
					//Display(circuit, lastpenalty, totalgates, 0, row);
				} 
		//now examine upper triangle starting from column N-1
		//due to the fact that the lower triangle is solved only CNOT
		//up are permitted, otherwise the gates might interfere with the
		//solved entries
		for (;column < N-1; column++)	{
			for (row = N-1; row > column && !(circuit[row]&((uint64_t)1<<column)); row--)
				;
			if (row!= column)
			{
				//First grow down
				for (int rowhelper = column; rowhelper + 1 < row; rowhelper++) {

					if (!(circuit[rowhelper + 1]&((uint64_t)1<<column))){
						cnotlist[totalgates] = rowhelper;
						totalgates++;
						circuit[rowhelper + 1] ^= circuit[rowhelper];
					}
				}
				//Next grow up
				for (;--row >= column;) {
						cnotlist[totalgates] = row;
						totalgates++;
						circuit[row+1] ^= circuit[row];
				//Display(circuit, lastpenalty, totalgates, 0, row);
				}
				
			}
		}

		//lastpenalty = CalculateCostHeuristic(circuit);
		//Display(circuit, lastpenalty, totalgates, 0, cnot);
	}
	cnotlist[totalgates] = INVALID;
	return totalgates;
}



void Display(int N, uint64_t * circuit, uint64_t cost, int gates, int depth, int cnot) {
	 
	uint64_t one = 1;
	int i=0, j, cnotcontrol, cnottarget;
	char A='A', chr;
	 //return;
	if (cnot != INVALID)
	{
		if (cnot>=0)
		{
			cnotcontrol=cnot;
			cnottarget=cnot+1;
		}
		else
		{
			cnotcontrol=-cnot;
			cnottarget=-cnot-1;
		}
		printf ("After CNOT(%d -> %d): ",cnotcontrol,cnottarget);
	}
	printf("Cost %lld, Total Gates %d, Depth %d\n", cost, gates, depth);
	for (;i<N;++i) {
		for(j=0;j<depth;j++)
			printf("    ");//indentation based on depth
		for(j=0;j<N;j++){
			if (circuit[i]&(one<<j))
				chr = (char)j+A;
			else
				chr = ' ';

			printf("%c ",chr);//Output appropriate variable
		}
		printf("\n");//end of row
	}
	printf("\n");//end with extra blank line
	;
}

