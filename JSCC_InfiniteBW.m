
% JSCC in high BW regime

P_N0 = -3:1:20;
P_N0_linear = 10.^(P_N0/10);

j = 1:1:1e3;
c = 

n_iter = 1e2;

SDR = zeros(length(P_N0),5);

for i=1:length(W)
    
    for j=1:length(P_N0)
        
        parfor iter = 1:n_iter
            
            beta = 1/sqrt(W(i));
            curr_source = 2*sqrt(12) * (rand - 0.5) % uniform source over [-1/sqrt(12),1/sqrt(12)]
            
            for k=1:n(i)
                
                Q = mod(beta*
        end
        
    end
    
end