// Comparisons of LNN Linear Reversible Circuit Synthesis Methods
// Copyright 2011, 2012 Ben Schaeffer
// Permission to copy this file is granted under the terms of the 
// GNU Lesser General Public License. See COPYING.LESSER.txt for details.
// Date: October 24, 2012
// Version: 0.4
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
#define ApplyCNOT(a, c) if (c < 0) a[-c-1] ^= a[-c]; else a[c+1] ^= a[c]
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
// Linear Nearest Neighbor Gaussian Elimination using the "upper triangle matrix" approach
int LNNGE_UTM(int N, uint64_t * inputcircuit, int * cnotlist); //returns gate count
// Linear Nearest Neighbor Gaussian Elimination using the "lower triangle matrix" approach
int LNNGE_LTM(int N, uint64_t * inputcircuit, int * cnotlist); //returns gate count
// Linear Nearest Neighbor Gaussian Elimination with Depth using the "upper triangle matrix" approach
int LNNGED_UTM(int N, uint64_t * inputcircuit, int * cnotlist, int depth); //returns gate count
// Linear Nearest Neighbor Gaussian Elimination with Depth using the "lower triangle matrix" approach
int LNNGED_LTM(int N, uint64_t * inputcircuit, int * cnotlist, int depth); //returns gate count
// Linear Nearest Neighbor Alternating Elimination, solve for upper diagonal first
int LNNAE_U(int N, uint64_t * inputcircuit, int * cnotlist); //returns gate count
// Linear Nearest Neighbor Alternating Elimination, solve for lower diagonal first
int LNNAE_L(int N, uint64_t * inputcircuit, int * cnotlist); //returns gate count
// Linear Nearest Neighbor Alternating Elimination with Depth, solve for upper diagonal first
int LNNAED_U(int N, uint64_t * inputcircuit, int * cnotlist, int depth); //returns gate count
// Linear Nearest Neighbor Alternating Elimination with Depth, solve for lower diagonal first
int LNNAED_L(int N, uint64_t * inputcircuit, int * cnotlist, int depth); //returns gate count

int main(void) 
{
	uint64_t inputcircuit[NMAX], circuit[NMAX], transposedcircuit[NMAX],
	inversecircuit[NMAX], inversetransposedcircuit[NMAX];
	int i, lowestcandidate, temp, cnotlist[CNOTLISTSIZE], gates;
	double algorithm1total;
	double GaussianEliminationtotal;
	double LNNGEtotal;
	double LNNGEBestOf8total;
	double LNNAEtotal;
	double LNNAEBestOf8total;
	double LNNGEDN16BestOf8total[5];
	double LNNAEDN16BestOf8total[5];
	bool displaysynthesis = FALSE;
	int N = NMAX;//Valid between 4 and 64
	
	Initialize();
	//Output in ".csv" format, compatible with spreadsheet programs
	printf("Comparisons of LNN Linear Reversible Circuit Synthesis"
	"Methods (Average Adjacent CNOT Gate Counts):\n, Gaussian Elimination,"
	" Algorithm 1, LNNGE, LNNAE, LNNGE Best of 8, LNNAED Best of 8\n");
	for (N = NSTART; N <= NEND; N += NINCREMENT){
		algorithm1total = 0;
		GaussianEliminationtotal = 0;
		LNNGEtotal = 0;
		LNNGEBestOf8total = 0;
		LNNAEtotal = 0;
		LNNAEBestOf8total = 0;
		//printf("Random linear reversible circuit synthesis of %d wires\n", N);
		for (i = 0;i < TESTS_TO_RUN; i++){
			Randomize(N, circuit);
			
			CopyCircuit(N, circuit, inputcircuit);
			gates = Long_Range_Gaussian_CNOT_Synthesis(N, inputcircuit, displaysynthesis);
			GaussianEliminationtotal += gates;
			//printf("%d,", gates);

			CopyCircuit(N, circuit, inputcircuit);
			gates = Algorithm_1_by_Patel_Markov_Hayes(N, inputcircuit, displaysynthesis);
			TransposeCircuit(N, inputcircuit, transposedcircuit);
			gates += Algorithm_1_by_Patel_Markov_Hayes(N, transposedcircuit, displaysynthesis);
			//printf("%d,", gates);
			algorithm1total += gates;

			//Prepare transposed matrix for future function calls
			TransposeCircuit(N, circuit, transposedcircuit);

			//LNNGE best of 8 approaches
			lowestcandidate = LNNGE_UTM(N, circuit, cnotlist);
			LNNGEtotal += lowestcandidate;
			//printf("%d,", lowestcandidate);
			
			//Use result to compute matrix inverse
			CopyCircuit(N, identity, inversecircuit);
			for(temp = 0; cnotlist[temp] != INVALID; temp++)
			ApplyCNOT(inversecircuit, cnotlist[temp]);
			TransposeCircuit(N, inversecircuit, inversetransposedcircuit);

			temp = LNNGE_LTM(N, circuit, cnotlist);
			if (lowestcandidate > temp)
			lowestcandidate = temp;
			//printf("%d,", lowestcandidate);
			temp = LNNGE_UTM(N, transposedcircuit, cnotlist);
			if (lowestcandidate > temp)
			lowestcandidate = temp;
			temp = LNNGE_LTM(N, transposedcircuit, cnotlist);
			if (lowestcandidate > temp)
			lowestcandidate = temp;
			//printf("%d,", lowestcandidate);			
			temp = LNNGE_UTM(N, inversecircuit, cnotlist);
			if (lowestcandidate > temp)
			lowestcandidate = temp;
			//printf("%d,", lowestcandidate);
			temp = LNNGE_LTM(N, inversecircuit, cnotlist);
			if (lowestcandidate > temp)
			lowestcandidate = temp;
			//printf("%d,", lowestcandidate);
			temp = LNNGE_UTM(N, inversetransposedcircuit, cnotlist);
			if (lowestcandidate > temp)
			lowestcandidate = temp;
			temp = LNNGE_LTM(N, inversetransposedcircuit, cnotlist);
			if (lowestcandidate > temp)
			lowestcandidate = temp;
			//printf("%d,", lowestcandidate);
			LNNGEBestOf8total += lowestcandidate;

			//LNNAE best of 8 approaches
			lowestcandidate = LNNAE_U(N, circuit, cnotlist);
			LNNAEtotal += lowestcandidate;
			//printf("%d,", lowestcandidate);
			temp = LNNAE_L(N, circuit, cnotlist);
			if (lowestcandidate > temp)
			lowestcandidate = temp;
			//printf("%d,", lowestcandidate);
			temp = LNNAE_U(N, transposedcircuit, cnotlist);
			if (lowestcandidate > temp)
			lowestcandidate = temp;
			temp = LNNAE_L(N, transposedcircuit, cnotlist);
			if (lowestcandidate > temp)
			lowestcandidate = temp;
			//printf("%d,", lowestcandidate);
			temp = LNNAE_U(N, inversecircuit, cnotlist);
			if (lowestcandidate > temp)
			lowestcandidate = temp;
			//printf("%d,", lowestcandidate);
			temp = LNNAE_L(N, inversecircuit, cnotlist);
			if (lowestcandidate > temp)
			lowestcandidate = temp;
			//printf("%d,", lowestcandidate);
			temp = LNNAE_U(N, inversetransposedcircuit, cnotlist);
			if (lowestcandidate > temp)
			lowestcandidate = temp;
			temp = LNNAE_L(N, inversetransposedcircuit, cnotlist);
			if (lowestcandidate > temp)
			lowestcandidate = temp;
			//printf("%d,", lowestcandidate);
			LNNAEBestOf8total += lowestcandidate;
			//printf("%f,", LNNAEBestOf8total);
			
			//Iterative deepening test for LNNGED and LNNAED
			if (N == 16)
			{
				for (int d = 0; d < 5; d++) 
				{
					//LNNGED best of 8 approaches
					lowestcandidate = LNNGED_UTM(N, circuit, cnotlist, d);
					//printf("%d,", lowestcandidate);
					temp = LNNGED_LTM(N, circuit, cnotlist, d);
					if (lowestcandidate > temp)
					lowestcandidate = temp;
					//printf("%d,", lowestcandidate);
					temp = LNNGED_UTM(N, transposedcircuit, cnotlist, d);
					if (lowestcandidate > temp)
					lowestcandidate = temp;
					temp = LNNGED_LTM(N, transposedcircuit, cnotlist, d);
					if (lowestcandidate > temp)
					lowestcandidate = temp;
					//printf("%d,", lowestcandidate);			
					temp = LNNGED_UTM(N, inversecircuit, cnotlist, d);
					if (lowestcandidate > temp)
					lowestcandidate = temp;
					//printf("%d,", lowestcandidate);
					temp = LNNGED_LTM(N, inversecircuit, cnotlist, d);
					if (lowestcandidate > temp)
					lowestcandidate = temp;
					//printf("%d,", lowestcandidate);
					temp = LNNGED_UTM(N, inversetransposedcircuit, cnotlist, d);
					if (lowestcandidate > temp)
					lowestcandidate = temp;
					temp = LNNGED_LTM(N, inversetransposedcircuit, cnotlist, d);
					if (lowestcandidate > temp)
					lowestcandidate = temp;
					//printf("%d,", lowestcandidate);
					LNNGEDN16BestOf8total[d] += lowestcandidate;

					//LNNAED best of 8 approaches
					lowestcandidate = LNNAED_U(N, circuit, cnotlist, d);
					//printf("%d,", lowestcandidate);
					temp = LNNAED_L(N, circuit, cnotlist, d);
					if (lowestcandidate > temp)
					lowestcandidate = temp;
					//printf("%d,", lowestcandidate);
					temp = LNNAED_U(N, transposedcircuit, cnotlist, d);
					if (lowestcandidate > temp)
					lowestcandidate = temp;
					temp = LNNAED_L(N, transposedcircuit, cnotlist, d);
					if (lowestcandidate > temp)
					lowestcandidate = temp;
					//printf("%d,", lowestcandidate);
					temp = LNNAED_U(N, inversecircuit, cnotlist, d);
					if (lowestcandidate > temp)
					lowestcandidate = temp;
					//printf("%d,", lowestcandidate);
					temp = LNNAED_L(N, inversecircuit, cnotlist, d);
					if (lowestcandidate > temp)
					lowestcandidate = temp;
					//printf("%d,", lowestcandidate);
					temp = LNNAED_U(N, inversetransposedcircuit, cnotlist, d);
					if (lowestcandidate > temp)
					lowestcandidate = temp;
					temp = LNNAED_L(N, inversetransposedcircuit, cnotlist, d);
					if (lowestcandidate > temp)
					lowestcandidate = temp;
					//printf("%d,", lowestcandidate);
					LNNAEDN16BestOf8total[d] += lowestcandidate;
				}
			}
		}
		printf("%d, %f, %f, %f, %f, %f, %f\n", N, 
		GaussianEliminationtotal/TESTS_TO_RUN,
		algorithm1total/TESTS_TO_RUN,
		LNNGEtotal/TESTS_TO_RUN,
		LNNAEtotal/TESTS_TO_RUN,
		LNNGEBestOf8total/TESTS_TO_RUN,
		LNNAEBestOf8total/TESTS_TO_RUN);
	}
	printf("Iterative deepening comparison of LNNGED and LNNAED for n=16\n");
	printf("Depth, Average CNOT gate count");
	for (int d = 0; d < 5; d++) 
	printf("\n%f, %f", LNNGEDN16BestOf8total[d]/TESTS_TO_RUN,
	LNNAEDN16BestOf8total[d]/TESTS_TO_RUN);
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

int LNNGE_UTM(int N, uint64_t * inputcircuit, int *cnotlist)
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

int LNNGE_LTM(int N, uint64_t * inputcircuit, int * cnotlist){//returns gate count
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

int LNNGED_UTM(int N, uint64_t * inputcircuit, int * cnotlist, int depth){//returns gate count
	int totalgates = 0, column = 0, row;
	uint64_t circuit[N], flag;
	if (depth == 0) 
	return LNNGE_UTM(N, inputcircuit, cnotlist);

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
					minimumheuristic = LNNGED_UTM(N, circuit, cnotlist+totalgates, depth - 1); 

					//compare against rest
					for (cnotdown = 1; cnotdown < row - rowabove; cnotdown++)
					{	// compute deltas, find cost, and ultimately restore
						circuit[rowabove + cnotdown] ^= circuit[rowabove + cnotdown + 1];
						circuit[rowabove + cnotdown] ^= circuit[rowabove + cnotdown - 1];
						temp = LNNGED_UTM(N, circuit, cnotlist+totalgates, depth - 1);
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

int LNNGED_LTM(int N, uint64_t * inputcircuit, int * cnotlist, int depth){//returns gate count
	int totalgates = 0, column, row;
	uint64_t circuit[N], flag;
	if (depth == 0)  
	return LNNGE_LTM(N, inputcircuit, cnotlist);

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
					minimumheuristic = LNNGED_LTM(N, circuit, cnotlist+totalgates, depth - 1);
					
					//compare against rest
					for (cnotup = 1; cnotup < rowbelow - row; cnotup++)
					{
						// compute deltas, find cost, and eventually restore circuit
						circuit[rowbelow - cnotup] ^= circuit[rowbelow - cnotup + 1];
						circuit[rowbelow - cnotup] ^= circuit[rowbelow - cnotup - 1];
						temp = LNNGED_LTM(N, circuit, cnotlist+totalgates, depth - 1);
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

// Linear Nearest Neighbor Alternating Elimination, solve for upper diagonal first
int LNNAE_U(int N, uint64_t * inputcircuit, int * cnotlist) //returns gate count
{
	int totalgates = 0, transposedtotalgates = 0, column, row, transposedcnotlist[CNOTLISTSIZE];
	uint64_t circuit[N], transposedcircuit[N], flag;

	for (int i = 0; i < N; i++) //use copy of circuit   
	circuit[i] = inputcircuit[i];
	for (flag = 1, column = 0; column < N - 1; column++, flag <<= 1)    //first phase of Gaussian Elimination
	{
		for (row = N - 1; row > column; row--)
		{
			if (circuit[row] & flag)
			{
				if (!(circuit[row - 1] & flag))
				{
					cnotlist[totalgates] = -(row);  //CNOT up gate
					totalgates++;
					circuit[row - 1] ^= circuit[row];
				}
				cnotlist[totalgates] = (row - 1);   //CNOT down gate
				totalgates++;
				circuit[row] ^= circuit[row - 1];
			}
		}
		//1. Transpose circuit
		//2. Use forward substitution and backwards elimination on column
		//3. Transpose back
		TransposeCircuit(N, circuit, transposedcircuit); 
		for (row = N - 1; row > column; row--)
		{
			if (transposedcircuit[row] & flag)
			{
				if (!(transposedcircuit[row - 1] & flag))
				{
					transposedcnotlist[transposedtotalgates] = -(row);  //CNOT up gate
					transposedtotalgates++;
					transposedcircuit[row - 1] ^= transposedcircuit[row];
				}
				transposedcnotlist[transposedtotalgates] = (row - 1);   //CNOT down gate
				transposedtotalgates++;
				transposedcircuit[row] ^= transposedcircuit[row - 1];
			}
		}
		TransposeCircuit(N, transposedcircuit, circuit); 			
	}
	//Terminate both gate lists and combine
	cnotlist[totalgates] = INVALID;
	transposedcnotlist[transposedtotalgates] = INVALID;
	ReverseandTransposeCNOTList(transposedcnotlist, cnotlist + totalgates);

	return totalgates + transposedtotalgates;
}

// Linear Nearest Neighbor Alternating Elimination, solve for lower diagonal first
int LNNAE_L(int N, uint64_t * inputcircuit, int * cnotlist) //returns gate count
{
	int totalgates = 0, transposedtotalgates = 0, column, row, transposedcnotlist[CNOTLISTSIZE];
	uint64_t circuit[N], transposedcircuit[N], flag;
	
	for (int i = 0; i < N; i++) //use copy of circuit   
	circuit[i] = inputcircuit[i];
	for (flag = (uint64_t)1 << (N - 1), column = N - 1; column > 0; column--, flag >>= 1)    //first phase of Gaussian Elimination
	{
		for (row = 0; row < column; row++)
		{
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
			}
		} 
		//1. Transpose circuit
		//2. Use forward substitution and backwards elimination on column
		//3. Transpose back
		TransposeCircuit(N, circuit, transposedcircuit); 
		for (row = 0; row < column; row++)
		{
			if(transposedcircuit[row] & flag)
			{
				if(!(transposedcircuit[row + 1] & flag))
				{
					transposedcnotlist[transposedtotalgates] = row; //CNOT down gate
					transposedtotalgates++;
					transposedcircuit[row + 1] ^= transposedcircuit[row];
				}
				transposedcnotlist[transposedtotalgates] = -(row + 1); //CNOT up gate
				transposedtotalgates++;
				transposedcircuit[row] ^= transposedcircuit[row + 1];
			}
		} 
		TransposeCircuit(N, transposedcircuit, circuit); 
	}
	//Terminate both gate lists and combine
	cnotlist[totalgates] = INVALID;
	transposedcnotlist[transposedtotalgates] = INVALID;
	ReverseandTransposeCNOTList(transposedcnotlist, cnotlist + totalgates);
	return totalgates + transposedtotalgates;
}



// Linear Nearest Neighbor Alternating Elimination with Depth, solve for upper diagonal first
int LNNAED_U(int N, uint64_t * inputcircuit, int * cnotlist, int depth) //returns gate count
{
	int totalgates = 0, transposedtotalgates = 0, column, row, transposedcnotlist[CNOTLISTSIZE];
	uint64_t circuit[N], transposedcircuit[N], flag;
	if (depth == 0)
	return LNNAE_U(N, inputcircuit, cnotlist);

	for (int i = 0; i < N; i++) //use copy of circuit   
	circuit[i] = inputcircuit[i];
	for (flag = 1, column = 0; column < N - 1; column++, flag <<= 1)    //first phase of Gaussian Elimination
	{
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
					do
					{
						cnotlist[totalgates] = -(row); //CNOT up gate
						totalgates++;
						cnotlist[totalgates] = (row - 1); //CNOT down gate
						totalgates++;
						circuit[row - 1] ^= circuit[row];
						circuit[row] ^= circuit[row - 1];
					}
					while (--row > column);
				}
				else
				{//choose best heuristic and adjust row
					int minimumheuristic, mincnotdown = 0, cnotdown = 0, rown, cd, temp;
					//first set minimumheuristic to all CNOT up
					for (rown = row; rown - 1> rowabove; rown--)
					circuit[rown - 1] ^= circuit[rown];
					minimumheuristic = LNNAED_U(N, circuit, cnotlist + totalgates, depth - 1); 

					//compare against rest
					for (cnotdown = 1; cnotdown < row - rowabove; cnotdown++)
					{	// compute deltas, find cost, and ultimately restore
						circuit[rowabove + cnotdown] ^= circuit[rowabove + cnotdown + 1];
						circuit[rowabove + cnotdown] ^= circuit[rowabove + cnotdown - 1];
						temp = LNNAED_U(N, circuit, cnotlist + totalgates, depth - 1);
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
		//1. Transpose circuit
		//2. Use forward substitution and backwards elimination on column
		//3. Transpose back
		TransposeCircuit(N, circuit, transposedcircuit); 

		for (row = N-1; row > column; row--)
		if(transposedcircuit[row] & flag)
		{
			//now check if for another instance of this variable
			//on the row above, necessitating a CNOT down gate
			if(transposedcircuit[row - 1] & flag)
			{
				transposedcnotlist[transposedtotalgates] = (row - 1); //CNOT down gate
				transposedtotalgates++;
				transposedcircuit[row] ^= transposedcircuit[row - 1];						
			}
			else 
			{
				//Check for higher instance of variable in the same column
				//i.e. row[0] represents a wire that physically is higher
				//than row[1]
				int rowabove = row - 2, instancefound = FALSE;
				while (!instancefound && rowabove >= column) 
				if (transposedcircuit[rowabove] & flag) 
				instancefound = TRUE;							
				else
				rowabove--;
				if (!instancefound) 
				{
					do
					{
						transposedcnotlist[transposedtotalgates] = -(row); //CNOT up gate
						transposedtotalgates++;
						transposedcnotlist[transposedtotalgates] = (row - 1); //CNOT down gate
						transposedtotalgates++;
						transposedcircuit[row - 1] ^= transposedcircuit[row];
						transposedcircuit[row] ^= transposedcircuit[row - 1];
					}
					while (--row > column);
				}
				else
				{//choose best heuristic and adjust row
					int minimumheuristic, mincnotdown = 0, cnotdown = 0, rown, cd, temp;
					//first set minimumheuristic to all CNOT up
					for (rown = row; rown - 1> rowabove; rown--)
					transposedcircuit[rown - 1] ^= transposedcircuit[rown];
					//In order to keep all operations consistent recursive function
					//calls need to use the non-transposed circuit
					TransposeCircuit(N, transposedcircuit, circuit); 
					minimumheuristic = LNNAED_U(N, circuit, cnotlist + totalgates, depth - 1); 

					//compare against rest
					for (cnotdown = 1; cnotdown < row - rowabove; cnotdown++)
					{	// compute deltas, find cost, and ultimately restore
						transposedcircuit[rowabove + cnotdown] ^= transposedcircuit[rowabove + cnotdown + 1];
						transposedcircuit[rowabove + cnotdown] ^= transposedcircuit[rowabove + cnotdown - 1];
						TransposeCircuit(N, transposedcircuit, circuit); 
						temp = LNNAED_U(N, circuit, cnotlist + totalgates, depth - 1);
						if (temp < minimumheuristic){
							minimumheuristic = temp;
							mincnotdown = cnotdown;
						}
					}//restore transposedcircuit
					for (rown = row; rown - 1 > rowabove; rown--)
					transposedcircuit[rown - 1] ^= transposedcircuit[rown - 2];

					//choose best
					cnotdown = mincnotdown;
					for (cd = 0; cd < cnotdown; cd++)
					{
						transposedcnotlist[transposedtotalgates] = (rowabove + cd); //CNOT down gate
						transposedtotalgates++;
						transposedcircuit[rowabove + 1 + cd] ^= transposedcircuit[rowabove + cd];
					}
					for (rown = row; rown - 1> rowabove+cnotdown; rown--)
					{
						transposedcnotlist[transposedtotalgates] = -(rown); //CNOT up gate
						transposedtotalgates++;
						transposedcircuit[rown - 1] ^= transposedcircuit[rown];
					}

					row++; //adjustment so row calculation starts over
				}
			}
		} 
		TransposeCircuit(N, transposedcircuit, circuit); 					
	}
	//Terminate both gate lists and combine
	cnotlist[totalgates] = INVALID;
	transposedcnotlist[transposedtotalgates] = INVALID;
	ReverseandTransposeCNOTList(transposedcnotlist, cnotlist + totalgates);

	return totalgates + transposedtotalgates;
}

// Linear Nearest Neighbor Alternating Elimination with Depth, solve for lower diagonal first
int LNNAED_L(int N, uint64_t * inputcircuit, int * cnotlist, int depth) //returns gate count
{
	int totalgates = 0, transposedtotalgates = 0, column, row, transposedcnotlist[CNOTLISTSIZE];
	uint64_t circuit[N], transposedcircuit[N], flag;
	
	if (depth == 0)
	return LNNAE_L(N, inputcircuit, cnotlist);
	for (int i = 0; i < N; i++) //use copy of circuit   
	circuit[i] = inputcircuit[i];
	for (flag = (uint64_t)1 << (N - 1), column = N - 1; column > 0; column--, flag >>= 1)    //first phase of Gaussian Elimination
	{
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
					do
					{
						cnotlist[totalgates] = row; //CNOT down gate
						totalgates++;
						cnotlist[totalgates] = -(row+1); //CNOT up gate
						totalgates++;
						circuit[row+1] ^= circuit[row];
						circuit[row] ^= circuit[row+1];
					}
					while (++row < column);
				}
				else
				{//choose best heuristic and adjust row
					int minimumheuristic, mincnotup = 0, cnotup = 0, rown, cu, temp;
					//first set minimumheuristic to all CNOT down
					for (rown = row; rown + 1 < rowbelow; rown++)
					circuit[rown + 1] ^= circuit[rown];
					minimumheuristic = LNNAED_L(N, circuit, cnotlist + totalgates, depth - 1);
					
					//compare against rest
					for (cnotup = 1; cnotup < rowbelow - row; cnotup++)
					{
						// compute deltas, find cost, and eventually restore circuit
						circuit[rowbelow - cnotup] ^= circuit[rowbelow - cnotup + 1];
						circuit[rowbelow - cnotup] ^= circuit[rowbelow - cnotup - 1];
						temp = LNNAED_L(N, circuit, cnotlist + totalgates, depth - 1);
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
		//1. Transpose circuit
		//2. Use forward substitution and backwards elimination on column
		//3. Transpose back
		TransposeCircuit(N, circuit, transposedcircuit); 
		for (row = 0; row < column; row++)
		if(transposedcircuit[row] & flag)
		{
			//now check if for another instance of this variable
			//on the row above, necessitating a CNOT up gate
			if(transposedcircuit[row+1] & flag)
			{
				transposedcnotlist[transposedtotalgates] = -(row+1); //CNOT up gate
				transposedtotalgates++;
				transposedcircuit[row] ^= transposedcircuit[row+1];						
			}
			else 
			{
				//Check for lower instance of variable in the same column
				int rowbelow = row + 2, instancefound = FALSE;
				while (!instancefound && rowbelow <= column) 
				if (transposedcircuit[rowbelow] & flag) 
				instancefound = TRUE;							
				else
				rowbelow++;
				if (!instancefound)
				{
					do
					{
						transposedcnotlist[transposedtotalgates] = row; //CNOT down gate
						transposedtotalgates++;
						transposedcnotlist[transposedtotalgates] = -(row+1); //CNOT up gate
						transposedtotalgates++;
						transposedcircuit[row+1] ^= transposedcircuit[row];
						transposedcircuit[row] ^= transposedcircuit[row+1];
					}
					while (++row < column);
				}
				else
				{//choose best heuristic and adjust row
					int minimumheuristic, mincnotup = 0, cnotup = 0, rown, cu, temp;
					//first set minimumheuristic to all CNOT down
					for (rown = row; rown + 1 < rowbelow; rown++)
					transposedcircuit[rown + 1] ^= transposedcircuit[rown];
					TransposeCircuit(N, transposedcircuit, circuit); 
					minimumheuristic = LNNAED_L(N, circuit, cnotlist + totalgates, depth - 1); 
					
					//compare against rest
					for (cnotup = 1; cnotup < rowbelow - row; cnotup++)
					{
						// compute deltas, find cost, and eventually restore transposedcircuit
						transposedcircuit[rowbelow - cnotup] ^= transposedcircuit[rowbelow - cnotup + 1];
						transposedcircuit[rowbelow - cnotup] ^= transposedcircuit[rowbelow - cnotup - 1];
						TransposeCircuit(N, transposedcircuit, circuit); 
						temp = LNNAED_L(N, circuit, cnotlist + totalgates, depth - 1);
						if (temp < minimumheuristic)
						{
							minimumheuristic = temp;
							mincnotup = cnotup;
						}
					}

					//restore transposedcircuit
					for (rown = row; rown + 1 < rowbelow; rown++)
					transposedcircuit[rown + 1] ^= transposedcircuit[rown + 2];
					
					//choose best
					cnotup = mincnotup;
					for (cu = 0; cu < cnotup; cu++) 
					{
						transposedcnotlist[transposedtotalgates] = -(rowbelow - cu); //CNOT up gate
						transposedtotalgates++;
						transposedcircuit[rowbelow - 1 - cu] ^= transposedcircuit[rowbelow - cu];
					}
					for (rown = row; rown + 1 < rowbelow - cnotup; rown++)
					{
						transposedcnotlist[transposedtotalgates] = rown; //CNOT down gate
						transposedtotalgates++;
						transposedcircuit[rown+1] ^= transposedcircuit[rown];
					}

					row--;//adjustment so row calculation starts over
				}
			}
		} 
		TransposeCircuit(N, transposedcircuit, circuit); 
	}
	//Terminate both gate lists and combine
	cnotlist[totalgates] = INVALID;
	transposedcnotlist[transposedtotalgates] = INVALID;
	ReverseandTransposeCNOTList(transposedcnotlist, cnotlist + totalgates);
	return totalgates + transposedtotalgates;
}

