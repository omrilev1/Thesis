% buildLDPC.m 
% This matlab function calls GenerateLDPC_mex.C
% which is a mex interface for the PEG code written
% by Xiao-Yu Hu.
%
% This function builds an low-density parity-check 
% (LDPC) code  using PEG algorithm of Hu and Arnold.
% This code ONLY generates 'regular' LDPC codes
% with symbol node degree (number of 1's in each column) 
% equal to 'd'. The check node degree (number of 1's in each
% row) will be as uniform as possible. Refer to Paper 
% of Hu and Arnold "Regular and Irregular 
% Progressive Edge Growth Tanner Graphs"
% for more information.
% 
% Author: Hatef Monajemi (monajemi@stanford.edu)
% Date: 2013 Feb 23
% 
% 
% M: number of rows
% N: number of columns
% d: Symbol degree of regular PEG Tanner graph(number of 1's in
%    each column)
% LDPCmat = LDPC matrix

function LDPCmat = buildLDPCmatrix(M,N,d)

codeName = ['regTannerGraph-M',num2str(M),'-N',num2str(N),'-d',num2str(d),'.txt'];
GenerateLDPC_mex(M,N,codeName,d);
LDPCmat = ReadTannerGraph(codeName);

end

