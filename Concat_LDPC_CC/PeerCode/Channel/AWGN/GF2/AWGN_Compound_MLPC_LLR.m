 function [lambda] = AWGN_Compound_MLPC_LLR(m,Constellation_Mapping_Array,Y,Sigma,future_matrix,k,Estimated_X,g0)
  
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
    
    temp_past = kron(past, ones(temp_future_size,1));
    temp_present_size = temp_Estimated_X_size*temp_future_size;
    temp_future = repmat(future_matrix, [size(past, 1), 1]);
    
    u0 = [temp_past, zeros(temp_present_size,1) ,temp_future];
    u1 = [temp_past, ones(temp_present_size,1) ,temp_future];
       
    x0 = mod(u0*g0,2);
    x1 = mod(u1*g0,2);
    
    t0 = x0*temp_powers;
    t0 = reshape(t0,temp_future_size,temp_Estimated_X_size);
    t1 = x1*temp_powers;
    t1 = reshape(t1,temp_future_size,temp_Estimated_X_size);

    X0 = Constellation_Mapping_Array(t0+1).';
    X1 = Constellation_Mapping_Array(t1+1).';
    
    exp_multiplier = (1/(2*Sigma^2));
    
    p0 = sum(exp(-((bsxfun(@minus,temp_real,real(X0))).^2+(bsxfun(@minus,temp_imag,imag(X0))).^2)*exp_multiplier),2);
    p1 = sum(exp(-((bsxfun(@minus,temp_real,real(X1))).^2+(bsxfun(@minus,temp_imag,imag(X1))).^2)*exp_multiplier),2);
    
    lambda = log(p0./p1).';
    
 end