------------------------------------------------------------------------
		CUFID_align network alignment algorithm
------------------------------------------------------------------------
CUFID_align finds the pairwise alignment (one-to-one mapping) of biological networks

------------------------------Release Note------------------------------
	10.31.2018
	1. A bug fixed in input networks read code
	2. Optimizing network read code  

--------------------------Matlab Implementaion--------------------------
 Input:
       Net_id_list  : List of the name for the PPI networks
       input_folder : File path for the input data
       id_flag      : If the nodes are in numeric format ("species id+number"),
                      set id_flag = 1
                      Otherwise, set id_flag = 0
       out_file     : Name of the output file
 Output:
       net_align    : Alignment results

 Example:
 Please type the following commands in the matlab command window

       Net_id_list = {'a', 'b'};
       input_folder = './test';
       id_flag = 1;
       out_file = './test/test_out.txt';
       net_align = CUFID_align(Net_id_list, input_folder, id_flag, out_file);

 Input file format
       - Network file: 
               It supports tab-delimited file format. 
               If two proteins have a interaction, 
               a1  a2
               a1  a3
               a4  a2
               Note that it also supports the edge weight:
               a1  a2  0.7
               a1  a3  0.5
               a4  a2  0.3
       - Similarity score file: 
               Elements in the first column indicate the nodes in the one network
               Elements in the second column indicate the nodes in another network
               Value in the third column indicates their node simialrity score. 

               As an example, a-b.sim may be as follows:
               a1  b1  100
               a1  b3  78
               a2  b2  80
               a2  b4  45
               a3  b4  60

 Output file format
       Each row in the output file is the aligned node pair.
               a1  b1
               a2  b2            
               a3  b4

 Reference:
 For more information on the algorithms, please see: 
       H. Jeong, X. Qian, and B.-J. Yoon (2016), Effective comparative 
       analysis of protein-protein interaction networks by measuring 
       the steady-state network flow using a Markov model, 
       BMC Bioinformatics, In Press.
 
 Contact: bjyoon@ece.tamu.edu
