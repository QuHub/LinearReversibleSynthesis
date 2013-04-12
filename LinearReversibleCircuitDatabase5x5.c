//LinearReversibleCircuitDatabase5x5.c
//Notes: minimum maximum marker write_count smart

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

//work on 5x5 proof of concept
#define INCOMPLETE_A (64) //marker to indicate this circuit needs to be searched
#define INCOMPLETE_B (128) //marker to indicate this circuit needs to be searched
#define IDENTITY_CIRCUIT (0xf) //special value for terminal circuit in database searches
#define LOWEST_GATE (-4)
#define HIGHEST_GATE (3)
#define GATE_OFFSET (5) //gates are saved to the database with GATE_OFFSET added
#define NMAX (64)
#define INVALID (128)
#define ApplyCNOT(a, c) if (c < 0) a[-c-1] ^= a[-c]; else a[c+1] ^= a[c]
#define FALSE (0)
#define TRUE (1)
#define CNOTLISTSIZE (4*NMAX*NMAX)

static unsigned int LNNGED4optimaltotal = 0;
static unsigned int LNNAED4optimaltotal = 0;
static unsigned char *buffer;
static unsigned char *counts;
static unsigned control_mask[5] = {0x1f, 0x3e0, 0x7c00, 0xf8000, 0x1f00000};
static uint64_t identity64[64];
static unsigned int verification_counter = 0;
void Initialize64(void);

bool CircuitsAreEquivalent(int N, uint64_t * circuitA, uint64_t * circuitB);
//Returns a count of equivalent optimal circuits given problem "circuit"
//Returns 0 if circuit is not reversible
unsigned long long CountEquivalentOptimalCircuits5x5(unsigned circuit);

//Do not call this function directly as it is only used as a helper
//function by the public function
unsigned long long PrivateCountEquivalentOptimalCircuits5x5(unsigned circuit);

void _5x5_verify(int i);
void _5x5_comparisons(int i);

void TransposeCircuit(int N, uint64_t * source, uint64_t * destination);
void ReverseandTransposeCNOTList(int * source, int * destination); //Transposes a INVALID-terminated gate list
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
void CopyCircuit(int N, uint64_t * source, uint64_t * destination);

int main (void)
{
    unsigned long write_count, minimum_a, maximum_a, minimum_b, maximum_b;
    unsigned long identity, circuit, next_circuit, control;
    int gate; //follows CNOT gate encoding where negative values indicate
    //CNOT up and positive values indicate CNOT down
    int iteration = 1;
    //identity is a function of N, and in this case N=4
    //					5=(N+1) 10=2(N+1) 15=3(N+1)
    identity = (1<<0) + (1<<6) + (1<<12) + (1<<18) + (1<<24);
    buffer = malloc(1<<25);
    counts = malloc(1<<25);
    if (!buffer || !counts)
    {
        puts("memory allocation failure");
        return 0;
    }
    //mark the identity matrix so it gets searched
    buffer[identity] = IDENTITY_CIRCUIT | INCOMPLETE_A;
    counts[identity] = 0;
    //Set starting minimum and maximum to index of identity matrix
    minimum_a = identity;
    maximum_a = identity;

    do
    {
        write_count = 0;
        minimum_b = identity;//Reset minimum and maximum for recalculation
        maximum_b = identity;
        for (circuit = minimum_a; circuit <= maximum_a; circuit++)
        {
            if (buffer[circuit] & INCOMPLETE_A)//if circuit has not been searched
            {
                buffer[circuit] ^= INCOMPLETE_A;//clear incomplete flag
                //Search all nearby circuits except predecessor circuit
                for (gate = LOWEST_GATE; gate <= HIGHEST_GATE; gate++)
                    if (gate != (int)buffer[circuit] - GATE_OFFSET)
                    {
                        if (gate >= 0)
                        {
                            control = circuit & control_mask[gate];
                            next_circuit = circuit ^ control << 5;
                        }
                        else
                        {
                            control = circuit & control_mask[-gate];
                            next_circuit = circuit ^ control >> 5;
                        }
                        if (buffer[next_circuit] == 0)//if circuit is unknown
                        {
                            //Now that it is known that this circuit is reversible and needs to be marked as incomplete
                            buffer[next_circuit] = (gate + GATE_OFFSET) | INCOMPLETE_B;
                            counts[next_circuit] = iteration;
                            write_count++;
                            //adjust minimum and maximum if necessary
                            if (minimum_b > next_circuit)
                                minimum_b = next_circuit;
                            else if (maximum_b < next_circuit)
                                maximum_b = next_circuit;
                        }
                    }
            }
        }
        if (write_count == 0)//Break if search is complete
            break;
        printf ("iteration = %2d, write_count = %ld\n", iteration++, write_count);
        write_count=0;
        //Second iteration
        minimum_a = identity;//Reset minimum and maximum for recalculation
        maximum_a = identity;
        for (circuit = minimum_b; circuit <= maximum_b; circuit++)
        {
            if (buffer[circuit] & INCOMPLETE_B)//if circuit has not been searched
            {
                buffer[circuit] ^= INCOMPLETE_B;//clear incomplete flag
                //Search all nearby circuits except predecessor circuit
                for (gate = LOWEST_GATE; gate <= HIGHEST_GATE; gate++)
                    if (gate != (int)buffer[circuit] - GATE_OFFSET)
                    {
                        if (gate >= 0)
                        {
                            control = circuit & control_mask[gate];
                            next_circuit = circuit ^ control << 5;
                        }
                        else
                        {
                            control = circuit & control_mask[-gate];
                            next_circuit = circuit ^ control >> 5;
                        }
                        if (buffer[next_circuit] == 0)//if circuit is unknown
                        {
                            //Now that it is known that this circuit is reversible and needs to be marked as incomplete
                            buffer[next_circuit] = (gate + GATE_OFFSET) | INCOMPLETE_A;
                            write_count++;
                            counts[next_circuit] = iteration;
                            //adjust minimum and maximum if necessary
                            if (minimum_a > next_circuit)
                                minimum_a = next_circuit;
                            else if (maximum_a < next_circuit)
                                maximum_a = next_circuit;
                        }
                    }
            }
        }
        printf ("iteration = %2d, write_count = %ld\n", iteration++, write_count);
    }
    while (write_count > 0);

    for (int i = 0; i< 1<<25; i++)
        if (buffer[i])
            _5x5_verify(i);
    if (verification_counter == 9999360)//expected value from equation
        //(2^5-1)*(2^5-2)*(2^5-4)*(2^5-8)*(2^5-16)
    {
        puts("Database Verified.");

        //FILE * f = fopen("5x5LNN_LRC.dat","wb");
        //fwrite(buffer, 1<<25, 1, f);
        //fclose(f);
        int temp = 0;
        for (int i = 0; i< 1<<25; i++)
        {
            if (buffer[i])
            {
                _5x5_comparisons(i);
            }
        }
        printf("LNNGED depth = 4 optimal total: %lld\n", LNNGED4optimaltotal);
        printf("LNNAED depth = 4 optimal total: %lld\n", LNNAED4optimaltotal);
        printf("Optimal total 9999360\n");
    }
    else
        printf("%lld errors in database detected", 9999360 -
               verification_counter);
    return 0;
}

unsigned long long CountEquivalentOptimalCircuits5x5(unsigned circuit)
//Returns a count of equivalent optimal circuits given problem "circuit"
//Returns 0 if circuit is not reversible
{
    if (buffer[circuit] == 0)
        return 0;
    return PrivateCountEquivalentOptimalCircuits5x5(circuit);
}

unsigned long long PrivateCountEquivalentOptimalCircuits5x5(unsigned circuit)
{
    unsigned long long count = 0;
    unsigned control, next_circuit;
    int gate;
    if (buffer[circuit] == IDENTITY_CIRCUIT)
        return 1;
    //Recursively add up all counts of equivalent circuits
    for (gate = LOWEST_GATE; gate <= HIGHEST_GATE; gate++)
    {
        if (gate >= 0)
        {
            control = circuit & control_mask[gate];
            next_circuit = circuit ^ control << 5;
        }
        else
        {
            control = circuit & control_mask[-gate];
            next_circuit = circuit ^ control >> 5;
        }
        if (counts[circuit] - 1 == counts[next_circuit])
            count += PrivateCountEquivalentOptimalCircuits5x5(next_circuit);
    }
    return count;
}

void Initialize64(void) {
    uint64_t one = 1;
    for (int i = 0; i < NMAX; i++) {
        identity64[i] = one << i;
    }
}

int VerifyCNOTList(int N, uint64_t * inputcircuit, int * cnotlist) //true return means verified, false fails
{
    int i;
    uint64_t circuit[NMAX];
    Initialize64();
    for (i = 0; i < N; i++) //use copy of circuit
        circuit[i] = identity64[i];
    for (i = 0; cnotlist[i] != INVALID; i++)
        ApplyCNOT(circuit, cnotlist[i]);
    return CircuitsAreEquivalent(N, circuit, inputcircuit);
}

bool CircuitsAreEquivalent(int N, uint64_t * circuitA, uint64_t * circuitB)
{
    for (int i = 0; i < N; i++)
        if (circuitA[i] != circuitB[i])
            return false;
    return true;
}

void _5x5_verify(int i)
{
    int caution = 25;//maximum optimal CNOT list length
    while(buffer[i] != IDENTITY_CIRCUIT)
    {
        if (!buffer[i])//not reversible
        {
            printf("Encountered database error... exiting");
            exit(1);
        }
        if (caution-- == 0)//check for CNOT list getting longer than maximum for 5x5
            return;//this circuit will not be counted but subsequent tests can continue
        if (buffer[i] - GATE_OFFSET>= 0)
        {
            i ^= (i & control_mask[buffer[i] - GATE_OFFSET]) << 5;
        }
        else
        {
            i ^= (i & control_mask[-(buffer[i] - GATE_OFFSET)]) >> 5;
        }
    }
    verification_counter++;
}

void _5x5_comparisons(int i)
{
    int optimalcount = 0;//maximum optimal CNOT list length
    int LNNGEDcount = 0;
    int LNNAEDcount = 0;
    uint64_t circuit[5], transposedcircuit[5], inversetransposedcircuit[5];
    uint64_t inversecircuit[5] = {1, 2, 4, 8, 16};
    int cnotlist [100];
    int temp;

    circuit[0] = i & 0x1f;
    circuit[1] = (i >> 5) & 0x1f;
    circuit[2] = (i >> 10) & 0x1f;
    circuit[3] = (i >> 15) & 0x1f;
    circuit[4] = (i >> 20) & 0x1f;

    while(buffer[i] != IDENTITY_CIRCUIT)
    {
        optimalcount++;
        if (buffer[i] - GATE_OFFSET>= 0)
        {
            i ^= (i & control_mask[buffer[i] - GATE_OFFSET]) << 5;
        }
        else
        {
            i ^= (i & control_mask[-(buffer[i] - GATE_OFFSET)]) >> 5;
        }
    }

    //Prepare transposed matrix for future function calls
    TransposeCircuit(5, circuit, transposedcircuit);
    LNNGE_UTM(5, circuit, cnotlist);
    for(temp = 0; cnotlist[temp] != INVALID; temp++)
        ApplyCNOT(inversecircuit, cnotlist[temp]);
    TransposeCircuit(5, inversecircuit, inversetransposedcircuit);

    //LNNGED best of 8 approaches
    LNNGEDcount = LNNGED_UTM(5,circuit, cnotlist, 4);
    temp = LNNGED_LTM(5,circuit, cnotlist, 4);
    if (LNNGEDcount > temp)
        LNNGEDcount = temp;
    temp = LNNGED_UTM(5,transposedcircuit, cnotlist, 4);
    if (LNNGEDcount > temp)
        LNNGEDcount = temp;
    temp = LNNGED_LTM(5,transposedcircuit, cnotlist, 4);
    if (LNNGEDcount > temp)
        LNNGEDcount = temp;
    temp = LNNGED_UTM(5,inversecircuit, cnotlist, 4);
    if (LNNGEDcount > temp)
        LNNGEDcount = temp;
    temp = LNNGED_LTM(5,inversecircuit, cnotlist, 4);
    if (LNNGEDcount > temp)
        LNNGEDcount = temp;
    temp = LNNGED_UTM(5,inversetransposedcircuit, cnotlist, 4);
    if (LNNGEDcount > temp)
        LNNGEDcount = temp;
    temp = LNNGED_LTM(5,inversetransposedcircuit, cnotlist, 4);
    if (LNNGEDcount > temp)
        LNNGEDcount = temp;
    if (LNNGEDcount == optimalcount)
        LNNGED4optimaltotal++;

    //LNNAED best of 8 approaches
    LNNAEDcount = LNNAED_U(5,circuit, cnotlist, 4);
    //printf("%d,", LNNAEDcount);
    temp = LNNAED_L(5,circuit, cnotlist, 4);
    if (LNNAEDcount > temp)
        LNNAEDcount = temp;
    //printf("%d,", LNNAEDcount);
    temp = LNNAED_U(5,transposedcircuit, cnotlist, 4);
    if (LNNAEDcount > temp)
        LNNAEDcount = temp;
    temp = LNNAED_L(5,transposedcircuit, cnotlist, 4);
    if (LNNAEDcount > temp)
        LNNAEDcount = temp;
    //printf("%d,", LNNAEDcount);
    temp = LNNAED_U(5,inversecircuit, cnotlist, 4);
    if (LNNAEDcount > temp)
        LNNAEDcount = temp;
    //printf("%d,", LNNAEDcount);
    temp = LNNAED_L(5,inversecircuit, cnotlist, 4);
    if (LNNAEDcount > temp)
        LNNAEDcount = temp;
    //printf("%d,", LNNAEDcount);
    temp = LNNAED_U(5,inversetransposedcircuit, cnotlist, 4);
    if (LNNAEDcount > temp)
        LNNAEDcount = temp;
    temp = LNNAED_L(5,inversetransposedcircuit, cnotlist, 4);
    if (LNNAEDcount > temp)
        LNNAEDcount = temp;
    if (LNNAEDcount == optimalcount)
        LNNAED4optimaltotal++;
}

void TransposeCircuit(int N, uint64_t * source, uint64_t * destination) {
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

void ReverseandTransposeCNOTList(int * source, int * destination) {
    int lower = 0, higher = 0;
    while(source[higher] != INVALID)
        higher++;
    destination[higher] = INVALID;
    while(source[lower] != INVALID)
        destination[--higher] = -(source[lower++]+1);
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

int LNNGE_LTM(int N, uint64_t * inputcircuit, int * cnotlist) { //returns gate count
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
    for (; column < N-1; column++)
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
            for (; --row >= column;)
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

int LNNGED_UTM(int N, uint64_t * inputcircuit, int * cnotlist, int depth) { //returns gate count
    int totalgates = 0, column = 0, row;
    uint64_t circuit[N], flag;
    if (depth == 0)
        return LNNGE_UTM(N, inputcircuit, cnotlist);

    for (int i = 0; i < N; i++) //use copy of circuit
        circuit[i] = inputcircuit[i];
    for (; column < N-1; column++)
    {   //first phase of Gaussian Elimination
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
                    {   //choose best heuristic and adjust row
                        int minimumheuristic, mincnotdown = 0, cnotdown = 0, rown, cd, temp;
                        //first set minimumheuristic to all CNOT up
                        for (rown = row; rown - 1> rowabove; rown--)
                            circuit[rown - 1] ^= circuit[rown];
                        minimumheuristic = LNNGED_UTM(N, circuit, cnotlist+totalgates, depth - 1);

                        //compare against rest
                        for (cnotdown = 1; cnotdown < row - rowabove; cnotdown++)
                        {   // compute deltas, find cost, and ultimately restore
                            circuit[rowabove + cnotdown] ^= circuit[rowabove + cnotdown + 1];
                            circuit[rowabove + cnotdown] ^= circuit[rowabove + cnotdown - 1];
                            temp = LNNGED_UTM(N, circuit, cnotlist+totalgates, depth - 1);
                            if (temp < minimumheuristic) {
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

int LNNGED_LTM(int N, uint64_t * inputcircuit, int * cnotlist, int depth) { //returns gate count
    int totalgates = 0, column, row;
    uint64_t circuit[N], flag;
    if (depth == 0)
        return LNNGE_LTM(N, inputcircuit, cnotlist);

    for (int i = 0; i < N; i++) //use copy of circuit
        circuit[i] = inputcircuit[i];
    for (column = N-1; column > 0; column--)
    {   //first phase of Gaussian Elimination
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
                    {   //choose best heuristic and adjust row
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
    for (; column < N-1; column++)
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
            for (; --row >= column;)
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
    for (; i<N; ++i) {
        for(j=0; j<depth; j++)
            printf("    ");//indentation based on depth
        for(j=0; j<N; j++) {
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
                    {   //choose best heuristic and adjust row
                        int minimumheuristic, mincnotdown = 0, cnotdown = 0, rown, cd, temp;
                        //first set minimumheuristic to all CNOT up
                        for (rown = row; rown - 1> rowabove; rown--)
                            circuit[rown - 1] ^= circuit[rown];
                        minimumheuristic = LNNAED_U(N, circuit, cnotlist + totalgates, depth - 1);

                        //compare against rest
                        for (cnotdown = 1; cnotdown < row - rowabove; cnotdown++)
                        {   // compute deltas, find cost, and ultimately restore
                            circuit[rowabove + cnotdown] ^= circuit[rowabove + cnotdown + 1];
                            circuit[rowabove + cnotdown] ^= circuit[rowabove + cnotdown - 1];
                            temp = LNNAED_U(N, circuit, cnotlist + totalgates, depth - 1);
                            if (temp < minimumheuristic) {
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
                    {   //choose best heuristic and adjust row
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
                        {   // compute deltas, find cost, and ultimately restore
                            transposedcircuit[rowabove + cnotdown] ^= transposedcircuit[rowabove + cnotdown + 1];
                            transposedcircuit[rowabove + cnotdown] ^= transposedcircuit[rowabove + cnotdown - 1];
                            TransposeCircuit(N, transposedcircuit, circuit);
                            temp = LNNAED_U(N, circuit, cnotlist + totalgates, depth - 1);
                            if (temp < minimumheuristic) {
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
                    {   //choose best heuristic and adjust row
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
                    {   //choose best heuristic and adjust row
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

void CopyCircuit(int N, uint64_t * source, uint64_t * destination) {
    for (int i = 0; i < N; i++)
        destination[i] = source[i];
}

