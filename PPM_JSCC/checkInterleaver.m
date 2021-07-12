% check interleaver

% option 1 
k = 3;
B = 3; 

table = zeros(k,B^k);

for i=1:k
    basic_vec = repmat(1:B^(i - 1),1,B);
    for j=1:1:B^(k - i)
        table(i,mod(1 + (j - 1)*B^i + (0:1:(B^i - 1)),B^k + 1)) = basic_vec + B^(i - 1) * (j - 1);
        
        disp(strcat('Indices = ',num2str(1 + (j - 1)*B^i + (0:1:(B^i - 1)))));
        
    end
end
