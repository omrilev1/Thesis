 function [ll_p0,ll_p1,ll_p2,ll_p3] = AWGN_GF4_MLPC_LL_N_m(m,Constellation_Mapping_Array,Y,Sigma,future_matrix,k,Estimated_X)
            
    Y = Y.';
    
    future_length = m-(k+1);
    
    temp_future_size = size(future_matrix,1)-1;
    
    future_matrix = future_matrix(temp_future_size-future_length+1:temp_future_size,1:2:end).';
        
    temp_future_size = size(future_matrix,1);
    
    temp_Estimated_X_size = size(Estimated_X,3);
    
    N = length(Y);
                            
    temp_real = real(Y);
    temp_imag = imag(Y);
    
    temp_powers = (2.^(m-1:-1:0)).';
                                        
    past = zeros(k-1,N);
    
    past(1:k-1,:) = Estimated_X(1,1:k-1,:);
    
    past = past.';
    
    u = [kron(past, ones(temp_future_size,1)), zeros(temp_Estimated_X_size*temp_future_size,2) ,repmat(future_matrix, [size(past, 1), 1])];
        
    t0 = u*temp_powers;
    t0 = reshape(t0,temp_future_size,temp_Estimated_X_size);
    t_01 = 2^(m-(k+1));
    t_10 = 2^(m-k);
    t1 = t0+t_01;
    t2 = t0+t_10;
    t3 = t0+t_01+t_10;

    X0 = Constellation_Mapping_Array(t0+1).';
    X1 = Constellation_Mapping_Array(t1+1).';
    X2 = Constellation_Mapping_Array(t2+1).';
    X3 = Constellation_Mapping_Array(t3+1).';
    
    exp_multiplier = (1/(2*Sigma^2));
    
    p0 = sum(exp(-((temp_real-real(X0)).^2+(temp_imag-imag(X0)).^2)*exp_multiplier),1);
    p1 = sum(exp(-((temp_real-real(X1)).^2+(temp_imag-imag(X1)).^2)*exp_multiplier),1);
    p2 = sum(exp(-((temp_real-real(X2)).^2+(temp_imag-imag(X2)).^2)*exp_multiplier),1);
    p3 = sum(exp(-((temp_real-real(X3)).^2+(temp_imag-imag(X3)).^2)*exp_multiplier),1);
     
    ll_p0 = log(p0).';
    ll_p1 = log(p1).';
    ll_p2 = log(p2).';
    ll_p3 = log(p3).';
    
    ll_p0(isinf(ll_p0)) = intmin();
    ll_p1(isinf(ll_p1)) = intmin();
    ll_p2(isinf(ll_p2)) = intmin();
    ll_p3(isinf(ll_p3)) = intmin();
        
 end