# Multi-Bit Flip-Flop (MBFF)


## Compile  
"make mbff": builds the MBFF binary     
"make rng": builds the random test generation binary  


## Usage
#### ./rng -B $b$ -N $n$ -M $m$ > test.in 
- Builds a set of $n$ flops (all constrained to the bounding box $[0 ... b] \times [0 ... b]$), and $m$ timing-critical paths.   
- Outputs the input file for the mbff binary to "test.in". 

Sample input file (test.in):
```txt
4 2 (#flops, #paths)
74.12 214.32 (x and y location of the 1st flop)
97.14 90.32
121.66 74.81
83.54 154.32
1 3 (path between flops 1 and 3)
2 4 
```



#### ./mbff $a$ $b$ test.in > test.out   
- Runs the MBFF clustering algorithm with $ALPHA = a$, $BETA = b$, and the input "test.in".     
- Outputs D, Z, W, and #trays used of each size after the algorithm is finished to "test.out".    

## Required Libraries    
- CPLEX (OR-Tools will be supported in The OpenROAD)  
- lemon (used for min-cost flow)   
- C++11 or higher (Makefile uses C++14)  
