 function [lambda] = AWGN_BPSK_LLR(Y,Sigma)
            
    lambda = (2.*Y).*(1/(Sigma^2));                
    
 end