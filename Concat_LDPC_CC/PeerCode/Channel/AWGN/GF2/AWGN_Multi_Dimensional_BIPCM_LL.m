 function [ll_p0,ll_p1] = AWGN_Multi_Dimensional_BIPCM_LL(m,Constellation_Mapping_Array,Y,Sigma,future_matrix)
  
    Y = Y.';
 
    future_matrix = future_matrix.';
        
    Constellation_2D_Dimension = size(Constellation_Mapping_Array,1);
    
    N = length(Y);
    
    U_N_div_m = N/Constellation_2D_Dimension;
                            
    ll_p0 = zeros(m,U_N_div_m);
    ll_p1 = zeros(m,U_N_div_m);
    
    temp_real = reshape(real(Y),Constellation_2D_Dimension,[]).';
    temp_imag = reshape(imag(Y),Constellation_2D_Dimension,[]).';
    
    temp_powers = (2.^(m-1:-1:0)).';
    
    temp_zeros = zeros(length(future_matrix(:,1)),1);
    
    exp_multiplier = (1/(2*Sigma^2));
    
    for k=1:1:m
            
        u = [future_matrix(:,1:k-1),temp_zeros,future_matrix(:,k:m-1)];
        t0 = u*temp_powers;
        t1 = t0+2^(m-k);
        
        p0 = 0;
        p1 = 0;
        for i=1:1:Constellation_2D_Dimension
            Constellation_Mapping_Array_temp = Constellation_Mapping_Array(i,:);
            X0 = Constellation_Mapping_Array_temp(t0+1);
            X1 = Constellation_Mapping_Array_temp(t1+1);
            p0 = p0+((bsxfun(@minus,temp_real(:,i),real(X0))).^2+(bsxfun(@minus,temp_imag(:,i),imag(X0))).^2);
            p1 = p1+((bsxfun(@minus,temp_real(:,i),real(X1))).^2+(bsxfun(@minus,temp_imag(:,i),imag(X1))).^2);
        end
        p0 = sum(exp(-p0*exp_multiplier),2);
        p1 = sum(exp(-p1*exp_multiplier),2);
        
        ll_p0(k,:) = log(p0);
        ll_p1(k,:) = log(p1);
        
    end
        
    ll_p0 = ll_p0(:)';
    ll_p1 = ll_p1(:)';
    
    ll_p0(isinf(ll_p0)) = intmin();
    ll_p1(isinf(ll_p1)) = intmin();
    
 end