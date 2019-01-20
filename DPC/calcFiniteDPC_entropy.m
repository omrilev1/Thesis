function H = calcFiniteDPC_entropy(Delta,alpha,sigma)
% This function calculate H(Y|v) for the Finiter DPC Case with alpha MMSE :
% Y = alpha*(v-alpha*S+U) - U + alpha*Z + alpha*Z

nPoints = 1e2;
v = linspace(-Delta,Delta,nPoints);
H = 0;
for i = 1 : length(v)
    
    % generate sample sequence , for current v
    S = 2*Delta * (rand(1,1e6) - 0.5);
    U = 2*Delta * (rand(1,1e6) - 0.5);
    Z = sigma*randn(1,1e6);
    X = alpha*(mod(v(i) - alpha*S + U + Delta,2*Delta) - Delta) - U;
    Y = X + alpha*S + alpha*Z;
    
    currH = calc_entropy(Y);
    
    H = H + currH;
end
H = H/nPoints;
end

function [ emp_entropy ] = calc_entropy( samples )
%CALC_ENTROPY
[ w, e ] = histcounts( samples );
bin_size    = e(2)-e(1);
w = w ./ sum(w);
emp_entropy = 0.00;
for i = 1:length(w)
    if( w(i) > 0.0 )
        emp_entropy = emp_entropy + (w(i) * log2(w(i)/bin_size));
    end
end
emp_entropy = -1 * emp_entropy;
end