// Comparisons of LNN Linear Reversible Circuit Synthesis Methods
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
//
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
void ReverseandTransposeCNOTList(int * source, int * destination); //Transposes a INVALID-terminated gate list
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
	int i, lowestcandidate, temp, cnotlist[CNOTLISTSIZE], gates;
	double algorithm1gatecounts;
	double lrGaussianEliminationcounts;
	double nnGaussianElimination_4runs_counts;
	double nnDepthGaussianElimination_4runs_counts;
	bool displaysynthesis = FALSE;
	int N = NMAX;//Valid between 4 and 64
	
	Initialize();
	//Output in ".csv" format, compatible with spreadsheet programs
	printf("Comparisons of LNN Linear Reversible Circuit Synthesis"
		"Methods (Average Adjacent CNOT Gate Counts):\n, Gaussian Elimination,"
		" Algorithm 1, LNNGE, LNNGED\n");
	printf("GE,Alg. 1,LNNGE,LNNGEx2,LNNGEx4,LNNGED1,LNNGED1x2,LNNGED1x4\n");
	for (N = NSTART; N <= NEND; N += NINCREMENT){
		algorithm1gatecounts = 0;
		lrGaussianEliminationcounts = 0;
		nnGaussianElimination_4runs_counts = 0;
		nnDepthGaussianElimination_4runs_counts = 0;
		printf("Random linear reversible circuit synthesis of %d wires\n", N);
		for (i = 0;i < TESTS_TO_RUN; i++){
			
			Randomize(N, circuit);

			CopyCircuit(N, circuit, inputcircuit);
			gates = Long_Range_Gaussian_CNOT_Synthesis(N, inputcircuit, displaysynthesis);
			lrGaussianEliminationcounts += gates;
			printf("%d,", gates);

			CopyCircuit(N, circuit, inputcircuit);
			gates = Algorithm_1_by_Patel_Markov_Hayes(N, inputcircuit, displaysynthesis);
			TransposeCircuit(N, inputcircuit, transposedcircuit);
			gates += Algorithm_1_by_Patel_Markov_Hayes(N, transposedcircuit, displaysynthesis);
			printf("%d,", gates);
			algorithm1gatecounts += gates;

			//Prepare transposed matrix for future function calls
			TransposeCircuit(N, circuit, transposedcircuit);

			//Nearest Neighbor Gaussian Elimination, best of 4 approaches
			lowestcandidate = LNNGE_FAUTM(N, circuit, cnotlist);
			printf("%d,", lowestcandidate);
			temp = LNNGE_FALTM(N, circuit, cnotlist);
			if (lowestcandidate > temp)
				lowestcandidate = temp;
			printf("%d,", lowestcandidate);
			
			temp = LNNGE_FAUTM(N, transposedcircuit, cnotlist);
			if (lowestcandidate > temp)
				lowestcandidate = temp;
			temp = LNNGE_FALTM(N, transposedcircuit, cnotlist);
			if (lowestcandidate > temp)
				lowestcandidate = temp;
			printf("%d,", lowestcandidate);
			
			nnGaussianElimination_4runs_counts += lowestcandidate;

			//Nearest Neighbor Gaussian Elimination with depth 1 subcolumn
			//searching during phase I, best of 4 approaches
			lowestcandidate = LNNGED_FAUTM(N, circuit, 1, cnotlist);
			printf("%d,", lowestcandidate);
			temp = LNNGED_FALTM(N, circuit, 1, cnotlist);
			if (lowestcandidate > temp)
				lowestcandidate = temp;
			printf("%d,", lowestcandidate);
			temp = LNNGED_FAUTM(N, transposedcircuit, 1, cnotlist);
			if (lowestcandidate > temp)
				lowestcandidate = temp;
			temp = LNNGED_FALTM(N, transposedcircuit, 1, cnotlist);
			if (lowestcandidate > temp)
				lowestcandidate = temp;
			nnDepthGaussianElimination_4runs_counts += lowestcandidate;
			printf("%d\n", lowestcandidate);
		}
		printf("\n\n\n\n\n");
	}
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

void ReverseandTransposeCNOTList(int * source, int * destination){
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
    int totalgates = 0, column = 0, row;
    uint64_t circuit[N], flag;

    for (int i = 0; i < N; i++) //use copy of circuit   
        circuit[i] = inputcircuit[i];
    for (; column < N - 1; column++)    //first phase of Gaussian Elimination
	{
		flag = (uint64_t)1 << column;
        for (row = N - 1; row > column; row--)
            if (circuit[row] & flag)
            {
                if (circuit[row - 1] & flag)
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
	}

    for (; column > 0; column--)
    {   //second phase of Gaussian Elimination
		flag = (uint64_t)1 << column;
        for (row = 0; row < column && !(circuit[row] & flag); row++)
            ;   //search for top instance of variable associated with the column
        if (row != column)
        {   //First extend "1"'s up
            for (int rowhelper = column; rowhelper - 1 > row; rowhelper--)
            {
                if (!(circuit[rowhelper - 1] & flag))
                {
                    cnotlist[totalgates] = -(rowhelper);    //CNOT up gate
                    totalgates++;
                    circuit[rowhelper - 1] ^= circuit[rowhelper];
                }
            }
            for (; ++row <= column;)
            {   //Next eliminate "1"'s 
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
	int totalgates = 0, column, row;
	uint64_t circuit[N], flag;
	
	for (int i = 0; i < N; i++) //use copy of circuit
		circuit[i] = inputcircuit[i];
	for (column = N-1; column > 0; column--) //first phase of Gaussian Elimination
	{
		flag = (uint64_t)1 << column;
		for (row = 0; row < column; row++)
			if(circuit[row] & flag)
			{
				if(!(circuit[row + 1] & flag))
				{
					cnotlist[totalgates] = row; //CNOT down gate
					totalgates++;
					circuit[row + 1] ^= circuit[row];
				}
				cnotlist[totalgates] = -(row + 1); //CNOT up gate
				totalgates++;
				circuit[row] ^= circuit[row + 1];
			//	Display(circuit, lastpenalty, totalgates, 0, row);
			} 
	}
	for (;column < N-1; column++)	
	{   //second phase of Gaussian Elimination
		flag = (uint64_t)1 << column;
		for (row = N-1; row > column && !(circuit[row] & flag); row--)
            ;   //search for lowest instance of variable associated with the column
		if (row != column)
        {   //First extend "1"'s down
			for (int rowhelper = column; rowhelper + 1 < row; rowhelper++)
			{
				if (!(circuit[rowhelper + 1] & flag))
				{
					cnotlist[totalgates] = rowhelper; //CNOT down gate
					totalgates++;
					circuit[rowhelper + 1] ^= circuit[rowhelper];
				}
			}
			for (;--row >= column;) 
            {   //Next eliminate "1"'s
				cnotlist[totalgates] = row; //CNOT down gate
				totalgates++;
				circuit[row + 1] ^= circuit[row];
			}			
		}			 
	}

	cnotlist[totalgates] = INVALID;
	return totalgates;
}

int LNNGED_FAUTM(int N, uint64_t * inputcircuit, int depth, int * cnotlist){//returns gate count
    int totalgates = 0, column = 0, row;
    uint64_t circuit[N], flag;
	if (depth == 0) 
		return LNNGE_FAUTM(N, inputcircuit, cnotlist);

    for (int i = 0; i < N; i++) //use copy of circuit   
        circuit[i] = inputcircuit[i];
	for (;column < N-1; column++)	
	{	//first phase of Gaussian Elimination
		flag = (uint64_t)1 << column;
		for (row = N-1; row > column; row--)
			if(circuit[row] & flag)
			{
				//now check if for another instance of this variable
				//on the row above, necessitating a CNOT down gate
				if(circuit[row - 1] & flag)
				{
					cnotlist[totalgates] = (row - 1); //CNOT down gate
					totalgates++;
					circuit[row] ^= circuit[row - 1];						
				}
				else 
				{
					//Check for higher instance of variable in the same column
					//i.e. row[0] represents a wire that physically is higher
					//than row[1]
					int rowabove = row - 2, instancefound = FALSE;
					while (!instancefound && rowabove >= column) 
						if (circuit[rowabove] & flag) 
							instancefound = TRUE;							
						else
							rowabove--;
					if (!instancefound) 
					{
						cnotlist[totalgates] = -(row); //CNOT up gate
						totalgates++;
						cnotlist[totalgates] = (row - 1); //CNOT down gate
						totalgates++;
						circuit[row - 1] ^= circuit[row];
						circuit[row] ^= circuit[row - 1];
					}
					else
					{//choose best heuristic and adjust row
						int minimumheuristic, mincnotdown = 0, cnotdown = 0, rown, cd, temp;
						//first set minimumheuristic to all CNOT up
						for (rown = row; rown - 1> rowabove; rown--)
							circuit[rown - 1] ^= circuit[rown];
						minimumheuristic = LNNGED_FAUTM(N, circuit, depth - 1, cnotlist+totalgates); 

						//compare against rest
						for (cnotdown = 1; cnotdown < row - rowabove; cnotdown++)
						{	// compute deltas, find cost, and ultimately restore
							circuit[rowabove + cnotdown] ^= circuit[rowabove + cnotdown + 1];
							circuit[rowabove + cnotdown] ^= circuit[rowabove + cnotdown - 1];
							temp = LNNGED_FAUTM(N, circuit, depth - 1, cnotlist+totalgates);
							if (temp < minimumheuristic){
								minimumheuristic = temp;
								mincnotdown = cnotdown;
							}
						}//restore circuit
						for (rown = row; rown - 1 > rowabove; rown--)
							circuit[rown - 1] ^= circuit[rown - 2];

						//choose best
						cnotdown = mincnotdown;
						for (cd = 0; cd < cnotdown; cd++)
						{
							cnotlist[totalgates] = (rowabove + cd); //CNOT down gate
							totalgates++;
							circuit[rowabove + 1 + cd] ^= circuit[rowabove + cd];
						}
						for (rown = row; rown - 1> rowabove+cnotdown; rown--)
						{
							cnotlist[totalgates] = -(rown); //CNOT up gate
							totalgates++;
							circuit[rown - 1] ^= circuit[rown];
						}

						row++; //adjustment so row calculation starts over
					}
				}
			} 
	}
    for (; column > 0; column--)
    {   //second phase of Gaussian Elimination
		flag = (uint64_t)1 << column;
        for (row = 0; row < column && !(circuit[row] & flag); row++)
            ;   //search for top instance of variable associated with the column
        if (row != column)
        {   //First extend "1"'s up
            for (int rowhelper = column; rowhelper - 1 > row; rowhelper--)
            {
                if (!(circuit[rowhelper - 1] & flag))
                {
                    cnotlist[totalgates] = -(rowhelper);    //CNOT up gate
                    totalgates++;
                    circuit[rowhelper - 1] ^= circuit[rowhelper];
                }
            }
            for (; ++row <= column;)
            {   //Next eliminate "1"'s 
                cnotlist[totalgates] = -(row);  //CNOT up gate
                totalgates++;
                circuit[row - 1] ^= circuit[row];
            }
        }
    }
    cnotlist[totalgates] = INVALID;
    return totalgates;
}

int LNNGED_FALTM(int N, uint64_t * inputcircuit, int depth, int * cnotlist){//returns gate count
    int totalgates = 0, column, row;
    uint64_t circuit[N], flag;
	if (depth == 0)  
		return LNNGE_FALTM(N, inputcircuit, cnotlist);

    for (int i = 0; i < N; i++) //use copy of circuit   
        circuit[i] = inputcircuit[i];
	for (column = N-1; column > 0; column--)	
	{	//first phase of Gaussian Elimination
		flag = (uint64_t)1 << column;
		for (row = 0; row < column; row++)
			if(circuit[row] & flag)
			{
				//now check if for another instance of this variable
				//on the row above, necessitating a CNOT up gate
				if(circuit[row+1] & flag)
				{
					cnotlist[totalgates] = -(row+1); //CNOT up gate
					totalgates++;
					circuit[row] ^= circuit[row+1];						
				}
				else 
				{
					//Check for lower instance of variable in the same column
					int rowbelow = row + 2, instancefound = FALSE;
					while (!instancefound && rowbelow <= column) 
						if (circuit[rowbelow] & flag) 
							instancefound = TRUE;							
						else
							rowbelow++;
					if (!instancefound)
					{
						cnotlist[totalgates] = row; //CNOT down gate
						totalgates++;
						cnotlist[totalgates] = -(row+1); //CNOT up gate
						totalgates++;
						circuit[row+1] ^= circuit[row];
						circuit[row] ^= circuit[row+1];
					}
					else
					{//choose best heuristic and adjust row
						int minimumheuristic, mincnotup = 0, cnotup = 0, rown, cu, temp;
						//first set minimumheuristic to all CNOT down
						for (rown = row; rown + 1 < rowbelow; rown++)
							circuit[rown + 1] ^= circuit[rown];
						minimumheuristic = LNNGED_FALTM(N, circuit, depth - 1, cnotlist+totalgates);
						
						//compare against rest
						for (cnotup = 1; cnotup < rowbelow - row; cnotup++)
						{
							// compute deltas, find cost, and eventually restore circuit
							circuit[rowbelow - cnotup] ^= circuit[rowbelow - cnotup + 1];
							circuit[rowbelow - cnotup] ^= circuit[rowbelow - cnotup - 1];
							temp = LNNGED_FALTM(N, circuit, depth - 1, cnotlist+totalgates);
							if (temp < minimumheuristic)
							{
								minimumheuristic = temp;
								mincnotup = cnotup;
							}
						}

						//restore circuit
						for (rown = row; rown + 1 < rowbelow; rown++)
							circuit[rown + 1] ^= circuit[rown + 2];
						
						//choose best
						cnotup = mincnotup;
						for (cu = 0; cu < cnotup; cu++) 
						{
							cnotlist[totalgates] = -(rowbelow - cu); //CNOT up gate
							totalgates++;
							circuit[rowbelow - 1 - cu] ^= circuit[rowbelow - cu];
						}
						for (rown = row; rown + 1 < rowbelow - cnotup; rown++)
						{
							cnotlist[totalgates] = rown; //CNOT down gate
							totalgates++;
							circuit[rown+1] ^= circuit[rown];
						}

						row--;//adjustment so row calculation starts over
					}
				}
			} 
	}
	for (;column < N-1; column++)	
	{   //second phase of Gaussian Elimination
		flag = (uint64_t)1 << column;
		for (row = N-1; row > column && !(circuit[row] & flag); row--)
            ;   //search for lowest instance of variable associated with the column
		if (row != column)
        {   //First extend "1"'s down
			for (int rowhelper = column; rowhelper + 1 < row; rowhelper++)
			{
				if (!(circuit[rowhelper + 1] & flag))
				{
					cnotlist[totalgates] = rowhelper; //CNOT down gate
					totalgates++;
					circuit[rowhelper + 1] ^= circuit[rowhelper];
				}
			}
			for (;--row >= column;) 
            {   //Next eliminate "1"'s
				cnotlist[totalgates] = row; //CNOT down gate
				totalgates++;
				circuit[row + 1] ^= circuit[row];
			}			
		}			 
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

