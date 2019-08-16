% The program can produce the parity check matrix of DVB-S2 with the size
% 5400*2700, and 1/2 code rate.

% Copyright (C) Yang XIAO, Beijing Jiaotong University, March 28, 2009, E-Mail: yxiao@bjtu.edu.cn.
%
% In our following paper, we find that the shortten codes of DVB-S2 have
% some problems,this program can provide a good  DVB-S2 LDPC code without cascaded BCH codes.  

% [1] Yang Xiao, Kiseon Kim, "Alternative good LDPC codes for DVB-S2", 9th International Conference\
% on  Signal Processing, 2008 (ICSP 2008), Beijing, 26-29 Oct. 2008, page(s): 1959-1962, ISBN: 978-1-4244-2178-7,
% INSPEC Accession Number: 10411270, Digital Object Identifier: 10.1109/ICOSP.2008.4697527
% Current Version Published: 2008-12-08 
% Abstract
%  This paper reveals that the current LDPC codes in DVB-S2 are not good, since they have the codes
%  of minimum weight of 5, which leads to their BER can not satisfy the requirements of DVB-S2. 
%  This paper develops an algorithm to estimate the minimum weight of LDPC codes based on generator 
% matrices. Using it, we can find the minimum weight problem of the current LDPC codes in DVB-S2. 
% To enhance the minimum weight problem of the LDPC codes, we propose a new approach to construct 
% the encodable irregular QC codes, which can not produce new short girths and can have much larger
% minimum weight and minimum distance of the codes. Simulations verify the construction of alternative
% QC LDPC codes for DVB-S2 to be valid, since better BER performance than that of existing irregular 
% LDPC codes of DVB-S2 is obtained.  
% The paper can be downloaded from Web site of IEEE Explore.

load DVB_S2_5400 H -mat

% When you get H, you can use the well known LDPC simulator of 1/2 rate to get the good
% BER curve as we provided.