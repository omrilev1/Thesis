% Program to test the girth4 for LDPC codes
% Copyright Yang XIAO, BJTU, July 22, 2007 
%
% Given a parity check marix H, the program can determine wether it has
% girth 4. Girth 4 of H will lead to a poor BER performance  
%
% Ref:
% [1] Y. Xiao, M -H Lee,Low complexity MIMO-LDPC CDMA systems over multipath channels, 
% IEICE Transactions on Communications, v E89-B, n 5, May, 2006, p 1713-1717
% [2] J. Fan, Y. Xiao, A method of counting the number of cycles in LDPC codes,
% 8th International Conference on Signal Processing, ICSP 2006,Volume: 3,  
% ISBN: 0-7803-9737-1, Digital Object Identifier: 10.1109/ICOSP.2006.345906
%
% The papers [1] and [2] can be downloaded from Web site of IEICE and IEEE Explore.
%

clear;

% If you have own H, you need not use the sub-program.
% Example. 

%Step 1:Input LDPC parity chech matrix
s=load('B1.txt');
M=252; N=504;
H(1:M,1:N)=0;
for j=1:N
   for i=1:3
      k=s(j,i);
       if k<=N;
       H(k,j)=1;
   end
end
end
  
%%%%%%% Girth4 Test %%%%%
% Here, we only concern whether a parity matrix H has girth4, instead of
% the number of girth4. If you need know the nember of girth 4, please
% refer Ref. [2].
% get the number of H
rows=size(H,1);
cols=size(H,2);

% We can test girth 4 by using the matrix O, see [1].

O=H*H';   
for i=1:rows
    O(i,i)=0;
end
for i=1:rows
girth(i)=max(O(i,:));
end
girth4=max(girth);

if girth4<2 
  fprintf('No girth 4')
else
   fprintf('The H matrix has girth 4')  % Provde the test result.
end    
% Display the matrice H and O
% If H matrix has no gith4, the O matrix in Fig.2 has no entry value to
% larger than 1.
figure(1)
mesh(H)
figure(2)
mesh(O)
