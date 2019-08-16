 function [ll_p0,ll_p1] = AWGN_BPSK_LL(Y,Sigma)
 
    temp_2var = 2*Sigma^2;
 
    temp_add = - 0.5*log(pi()*temp_2var);
        
    ll_p0 = temp_add - ((Y-1).^2)./temp_2var;
    ll_p1 = temp_add - ((Y+1).^2)./temp_2var;      
        
 end