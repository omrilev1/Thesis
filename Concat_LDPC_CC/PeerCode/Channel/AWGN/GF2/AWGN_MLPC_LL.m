 function [ll_p0,ll_p1] = AWGN_MLPC_LL(m,Constellation_Mapping_Array,Y,Sigma,future_matrix,k,Estimated_X)
            
    Y = Y.';
    
    future_length = m-k;
    
    temp_future_size = size(future_matrix,1);
    
    future_matrix = future_matrix(temp_future_size-future_length+1:temp_future_size,1:2^future_length).';
        
    temp_future_size = size(future_matrix,1);
    
    temp_Estimated_X_size = size(Estimated_X,3);
    
    N = length(Y);
                            
    temp_real = real(Y);
    temp_imag = imag(Y);
    
    temp_powers = (2.^(m-1:-1:0)).';
                                        
    past = zeros(k-1,N);
    
    past(1:k-1,:) = Estimated_X(1,1:k-1,:);
    
    past = past.';
    
    u = [kron(past, ones(temp_future_size,1)), zeros(temp_Estimated_X_size*temp_future_size,1) ,repmat(future_matrix, [size(past, 1), 1])];
        
    t0 = u*temp_powers;
    t0 = reshape(t0,temp_future_size,temp_Estimated_X_size);

    t1 = t0+2^(m-k);

    X0 = Constellation_Mapping_Array(t0+1).';
    X1 = Constellation_Mapping_Array(t1+1).';
    
    exp_multiplier = (1/(2*Sigma^2));
    
    p0 = sum(exp(-((bsxfun(@minus,temp_real,real(X0))).^2+(bsxfun(@minus,temp_imag,imag(X0))).^2)*exp_multiplier),2);
    p1 = sum(exp(-((bsxfun(@minus,temp_real,real(X1))).^2+(bsxfun(@minus,temp_imag,imag(X1))).^2)*exp_multiplier),2);
    
    ll_p0 = log(p0).';
    ll_p1 = log(p1).';
    
    ll_p0(isinf(ll_p0)) = intmin();
    ll_p1(isinf(ll_p1)) = intmin();
        
 end