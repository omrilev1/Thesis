 function [ll_p0,ll_p1] = AWGN_BIPCM_LL(m,Constellation_Mapping_Array,Y,Sigma,future_matrix)
  
    Y = Y.';
 
    future_matrix = future_matrix.';
        
    N = length(Y);
                        
    ll_p0 = zeros(m,N);
    ll_p1 = zeros(m,N);
    
    temp_real = real(Y);
    temp_imag = imag(Y);
    
    temp_powers = (2.^(m-1:-1:0)).';
    
    temp_zeros = zeros(length(future_matrix(:,1)),1);
    
    exp_multiplier = (1/(2*Sigma^2));
    
    for k=1:1:m
            
        u = [future_matrix(:,1:k-1),temp_zeros,future_matrix(:,k:m-1)];
        t0 = u*temp_powers;
        t1 = t0+2^(m-k);

        X0 = Constellation_Mapping_Array(t0+1);
        X1 = Constellation_Mapping_Array(t1+1);

        p0 = sum(exp(-((bsxfun(@minus,temp_real,real(X0))).^2+(bsxfun(@minus,temp_imag,imag(X0))).^2)*exp_multiplier),2);
        p1 = sum(exp(-((bsxfun(@minus,temp_real,real(X1))).^2+(bsxfun(@minus,temp_imag,imag(X1))).^2)*exp_multiplier),2);

        ll_p0(k,:) = log(p0);
        ll_p1(k,:) = log(p1);
        
    end
        
    ll_p0 = ll_p0(:)';
    ll_p1 = ll_p1(:)';
    
    ll_p0(isinf(ll_p0)) = intmin();
    ll_p1(isinf(ll_p1)) = intmin();
    
 end