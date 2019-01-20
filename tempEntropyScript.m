
Delta = 1;
alpha = 0.1; sigma = 2;
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





function [ emp_entropy ] = calc_entropy( samples )
%CALC_ENTROPY 
    [ w, e ] = histcounts( samples );
    bin_size    = e(2)-e(1);
    w = w ./ sum(w);
    emp_entropy = 0.00;
    for i = 1:length(w)
        if( w(i) > 0.0 )
            emp_entropy = emp_entropy + (w(i) * log(w(i)/bin_size));
        end
    end
    emp_entropy = -1 * emp_entropy;
end