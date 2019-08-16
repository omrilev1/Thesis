 function [lambda] = AWGN_BIPCM_LLR(m,Constellation_Mapping_Array,Y,Sigma,future_matrix)
  
    Y = Y.';
 
    future_matrix = future_matrix.';
        
    N = length(Y);
                        
    lambda = zeros(m,N);
    
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

        lambda(k,:) = log(p0./p1);
        
    end
    
    lambda = lambda(:)';
    
 end