//Controller behavior:
//When the controller starts up it enters the state IDLE and "count" is latched
//The external line "start" must be held low to maintain the IDLE state
//When a computation request is desired "start" is pulsed high for one cycle
//which causes the controller to enter the LOAD state and changes "count" to
//zero. "count" will be held low on the subsequent N	
//clock cycles as new rows are shifted into the 2D matrix. At this point the
//controller enters the CALCULATE state which takes 2(N-1)**2 cycles which
//this consists of (N-1) repetitions of (N-1) vertical computations followed
//(N-1) horizontal computations. During this state"count" is incremented
//for each elementary row operation performed in the LNNAE algorithm.


//`define DEBUG

import "DPI-C" function shortint unsigned random4by4();
import "DPI-C" function shortint unsigned run4by4(shortint unsigned array);
import "DPI-C" function void seedRNG();

localparam  N = 4;  // size of one dimension of the square array
localparam NUM_RANDOMIZED = 20000;  // numer of randomized tests to perform
localparam  HORIZONTAL = 1'b0;
localparam  VERTICAL = 1'b1;

module LNNAE2D_Systolic #(parameter N = 4) (input wire clock,  shiftdirection, loading, enable, input logic [0:N-1] datain, output  logic [7:0] totalCNOTsOut);

//Cells represent three-input, eight entry LUT, single FF FPGA cell
module TwoDShiftCell0(input wire left, top, shiftdirection, clock, enable, 
	output wire right, bottom);
logic out=0;
assign right = out;
assign bottom = out;
always_ff @(posedge clock iff enable)
  if (shiftdirection == HORIZONTAL) out <= left;
  else out <= top;
endmodule

module mux2to1(input wire d0, input wire d1, input wire s, output wire y);
logic out;
assign out = (s) ? d1:d0;
assign y=out;
endmodule

module combL(input wire low, input wire high, input logic control, output wire y);
assign y =  (control) ? low^high:low;
endmodule

module combH(input wire low, input wire high, input logic control, output wire y);
assign y =  (control) ? low^high:high;
endmodule


wire logic h[N][N+1], v[N+2][N];//Horizontal and vertical wires in matrix
wire logic hCombLout[N], hCombHout[N], vCombLout[N], vCombHout[N]; 

wire logic shifting; 
var logic [7:0] totalCNOTs;
var logic [1:0] CNOTs;
assign shifting = enable | loading;
assign totalCNOTsOut = totalCNOTs;

var logic columncombLcontrol, columncombHcontrol, rowcombLcontrol, rowcombHcontrol;

genvar i, j;
//matrix declarations
generate
	for (i = 0; i < N-1; i++) begin: rowvar
		for (j = 0; j < N-1; j++) begin: colvar
		    TwoDShiftCell0 a(h[i][j], v[i][j], shiftdirection, clock, shifting, h[i][j+1], v[i+1][j]);
		end
		//vertical wire glue 
		combL cvlow (  v[N-1][i], v[N][i], rowcombLcontrol, vCombLout[i]);
		TwoDShiftCell0 abottom(h[N-1][i], vCombLout[i], shiftdirection, clock, shifting, h[N-1][i+1], v[N][i]);
		combH cvhigh(  vCombLout[i], v[N][i], rowcombHcontrol, vCombHout[i]);
		mux2to1 m21(vCombHout[i], datain[i], loading, v[N+1][i]);
		assign v[0][i] = v[N+1][i];
		
		//horizontal wire glue 
		combL chlow (  h[i][N-1], h[i][N], columncombLcontrol, hCombLout[i]);
		TwoDShiftCell0 aend(hCombLout[i], v[i][N-1], shiftdirection, clock, shifting, h[i][N], v[i+1][N-1]);
		combH chhigh(  hCombLout[i], h[i][N], columncombHcontrol, hCombHout[i]);
		assign h[i][0] = hCombHout[i];		
		 
	end
endgenerate
//corner cell, wires + glue
TwoDShiftCell0 acorner(hCombLout[N-1], vCombLout[N-1], shiftdirection, clock, shifting, h[N-1][N], v[N][N-1]);
//vertical wire glue 
combL cvlow (  v[N-1][N-1], v[N][N-1], rowcombLcontrol, vCombLout[N-1]);
combH cvhigh(  vCombLout[N-1], v[N][N-1], rowcombHcontrol, vCombHout[N-1]);
mux2to1 m21end(vCombHout[N-1], datain[N-1], loading, v[N+1][N-1]);
assign v[0][N-1] = v[N+1][N-1];

//horizontal wire glue 
combL chlow (  h[N-1][N-1], h[N-1][N], columncombLcontrol, hCombLout[N-1]);
combH chhigh(  hCombLout[N-1], h[N-1][N], columncombHcontrol, hCombHout[N-1]);
assign h[N-1][0] = hCombHout[N-1];

assign rowcombLcontrol = ~h[N-2][1] & h[N-1][1] & ~loading;
assign rowcombHcontrol = h[N-1][1] & ~loading;
assign columncombLcontrol = ~h[N-1][N-1] & h[N-1][N];
assign columncombHcontrol = h[N-1][N];


always_comb
begin
if (shiftdirection == VERTICAL)
	case({h[N-2][1],h[N-1][1]}) // bits to check to count CNOT gates
	2'b00: CNOTs <= '0;
	2'b01: CNOTs <= 2'b10;
	2'b10: CNOTs <= '0;
	2'b11: CNOTs <= 2'b01;
	endcase
else 
	case({h[N-1][N-1],h[N-1][N]}) // bits to check to count CNOT gates
	2'b00: CNOTs <= '0;
	2'b01: CNOTs <= 2'b10;
	2'b10: CNOTs <= '0;
	2'b11: CNOTs <= 2'b01;
	endcase
end

always_ff@(posedge clock)
begin
if (loading)
	totalCNOTs <= '0;
else 
	begin 
	if (enable)
		begin
		totalCNOTs <= totalCNOTs + {6'b0, CNOTs};
		end
	else
		totalCNOTs <= totalCNOTs;
	end
end

`ifdef DEBUG
always @(posedge clock)
begin
    $display("Matrix Changed at time: %t, loading = %b, enable = %b, shifting = %b", $time, loading, enable, shifting );
    $display("Matrix Changed:  col1 col2  col3 col4, totalCNOTs = %d CNOTs = %d", totalCNOTs,CNOTs );
    $display("         Row 1:    %b   %b    %b   %b", h[0][1], h[0][2], h[0][3], h[0][4]);
    $display("         Row 2:    %b   %b    %b   %b", h[1][1], h[1][2], h[1][3], h[1][4]);
	$display("         Row 3:    %b   %b    %b   %b", h[2][1], h[2][2], h[2][3], h[2][4]);
    $display("         Row 4:    %b   %b    %b   %b", h[3][1], h[3][2], h[3][3], h[3][4]);
end
`endif

property vertAndEnable;
    ((shiftdirection == VERTICAL) &&  enable && ~loading) |=> (h[0][1] == 0);
endproperty    

property horizAndEnable;
    ((shiftdirection == HORIZONTAL) &&  enable && ~loading) |=> (h[N-1][1] == 0);
endproperty    

property load4CyclesOnly;
    $rose(loading) |=> loading ##1 loading ##1 loading ##1 ~loading;
endproperty

always @(posedge clock)
begin
    assert property (vertAndEnable);
    assert property (horizAndEnable);
    assert property (load4CyclesOnly);
    assert (! (loading && enable));
end

endmodule


// Behavioural Model is just a function that emulates the processing in the systolic array, we made this beforehand to verify the theory of operation

function automatic int behavioralLNNAE(input logic [0:N-1][0:N-1]a);
int row, column, i, j, k, count = 0;
logic [0:N-1][0:N-1]m;
logic [0:N-1] temp, L, H;
//emulate shifting in new data
for(i = N-1; i >= 0; i--) begin
	temp = a[i];
	for(j = N-1; j >= 1; j--) begin
		m[j] =m[j-1];
	end
	m[0] = temp;
end
`ifdef DEBUG
$display ("Processing new matrix:");
for(i = 0; i < N; i++) $display ("%b", m[i]);
`endif

// outer loop, iterate N-1 times
for(i = N; i > 1; i--) begin
	//row processing phase, iterate N-1 times
	for(row = N; row > 1; row--) begin
		if (m[N-1][0]) begin
			count++;
			if (!m[N-2][0]) 
				count++;
		end
		// first assign combinational logic for L and H
		if (!m[N-2][0] && m[N-1][0])
			L = m[N-2]^m[N-1];
		else
			L = m[N-2];
		if (m[N-1][0])
			H = L^m[N-1];
		else
			H = m[N-1];
		//perform new assignments
		m[N-1] = L;
		for(j = N - 2; j > 0; j--) begin
			m[j] =m[j-1];
		end
		m[0] = H;
	end
	
	for(column = N; column > 1; column--) begin
	//column processing phase, iterate N-1 times
		if (m[N-1][N-1]) begin
			count++;
			if (!m[N-1][N-2]) 
				count++;
		end
		if (!m[N-1][N-2] && m[N-1][N-1])
			L = {m[0][2]^m[0][3], m[1][2]^m[1][3], m[2][2]^m[2][3], m[3][2]^m[3][3]};
		else//I would prefer to make these statements based on N
			L = {m[0][2], m[1][2], m[2][2], m[3][2]};
		if (m[N-1][N-1])
			H = L^{m[0][3], m[1][3], m[2][3], m[3][3]};
		else
			H = {m[0][3], m[1][3], m[2][3], m[3][3]};
		//perform new assignments
		for(k = 0; k < N; k++) m[k][N-1] = L[k];
		for(j = N - 2; j > 0; j--) begin
			for(k = 0; k < N; k++) 
				m[k][j] =m[k][j-1];
		end
		for(k = 0; k < N; k++) m[k][0] = H[k];
	end	
end

return count;
endfunction


module Testing01();

var logic clock   = 0;
var logic enable  = 0;
var logic loading = 0;
var logic shiftdirection = VERTICAL;
var logic [0:N-1] datain;
var logic [7:0] totalCNOTs;

shortint unsigned randomArray, gateCount, behavioralGateCount, errorCount = 0;
int scoreBoard[int]; 
// in the case that any error is detected, we just pop it in an associative array to report at the end.
shortint unsigned errorBoard[shortint unsigned]; 
var logic [15:0] inArrayBits;
var logic [N-1:0][N-1:0] behavioralInput;

LNNAE2D_Systolic #(.N(N)) dut(.clock(clock), .shiftdirection(shiftdirection), .loading(loading), .enable(enable), .datain(datain), .totalCNOTsOut(totalCNOTs));

always #50 clock = ~clock;

task loadRow(input logic[N-1:0] value);
    shiftdirection=VERTICAL;
    enable = 0;
    loading = 1;
    datain = value; 
    @(posedge clock);
    loading = 0;
endtask

task resetMachine();
    loading = 0;
    enable = 0;
endtask

task runCalculation();
    enable = 1;
    shiftdirection=VERTICAL;
    @(posedge clock);
`ifdef DEBUG
    $display("end loading"); 
    $display("shiftdirection=VERTICAL"); 
`endif
    for(int i = N-1; i>0; i--) 
    begin
        repeat(N-2)@(posedge clock);
        shiftdirection=HORIZONTAL;
        @(posedge clock);
`ifdef DEBUG
        $display("shiftdirection=HORIZONTAL"); 
`endif
        repeat(N-2)@(posedge clock);
        shiftdirection=VERTICAL;
        if (i == 1) 
        begin
            enable = 0;
        end
        @(posedge clock);
        if (i != 1) 
        begin
`ifdef DEBUG
            $display("shiftdirection=VERTICAL"); 
`endif
        end
    end
endtask

task loadCalc1Array(shortint unsigned randomArray);
    inArrayBits = randomArray;
    gateCount = run4by4(randomArray);
    loadRow({<< {inArrayBits[3:0]}});
    loadRow({<< {inArrayBits[7:4]}});
    loadRow({<< {inArrayBits[11:8]}});
    loadRow({<< {inArrayBits[15:12]}});
    runCalculation();
    behavioralInput[0] = ({<< {inArrayBits[3:0]}});
    behavioralInput[1] = ({<< {inArrayBits[7:4]}});
    behavioralInput[2] = ({<< {inArrayBits[11:8]}});
    behavioralInput[3] = ({<< {inArrayBits[15:12]}});
    behavioralGateCount = behavioralLNNAE(behavioralInput);
`ifdef DEBUG
    $display("posting results of golden model:");
    $display(gateCount);
    $display("posting results of RTL model:");
    $display(totalCNOTs);
    $display("posting results of behavioral model:");
    $display(behavioralGateCount);
`endif

    if (gateCount != totalCNOTs)
    begin
        errorCount++;
        errorBoard[randomArray] = 1;
    end            
    else if (gateCount != behavioralGateCount)
    begin            
        errorCount++;
        errorBoard[randomArray] = 1;
    end            



    if (scoreBoard.exists(gateCount))
        scoreBoard[gateCount] = scoreBoard[gateCount] + 1;
    else
        scoreBoard[gateCount] = 1;
endtask

initial
begin
    shortint unsigned directArray;
    resetMachine();
// Do Directed testing of our 4 known arrays    
    directArray = (1 << 12) + (2 << 8) + (4 << 4) + (8);
    loadCalc1Array(directArray);
    directArray = (8 << 12) + (4 << 8) + (2 << 4) + (1);
    loadCalc1Array(directArray);
    directArray = (1 << 12) + (3 << 8) + (7 << 4) + (15);
    loadCalc1Array(directArray);
    directArray = (3 << 12) + (2 << 8) + (4 << 4) + (8);
    loadCalc1Array(directArray);

// Do Randomized testing    
    seedRNG();
    for (int i = 0; i < NUM_RANDOMIZED ; i++)
    begin
        randomArray = random4by4();
        loadCalc1Array(randomArray);
                        
    end        
    $display("Error Count:");
    $display(errorCount);
    if (errorCount > 0)
    begin
        $display("Error Details (list of erroring input arrays) :");
        foreach (errorBoard[errorArray])
            $display(errorBoard[errorArray]);
    end            
     
    $display("gateCount Coverage");
    for (int i = 0 ; i < 31 ; i++ )
    begin
        if (scoreBoard.exists(i))
            $display("gateCount = %d covered %d times", i, scoreBoard[i]);
    end

    $finish;
end
endmodule
