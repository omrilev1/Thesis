 function [ll_p0,ll_p1,ll_p2,ll_p3] = AWGN_GF4_BIMLPCM_LL(m,BIPCM_m,Constellation_Mapping_Array,Y,Sigma,future_matrix,k,Estimated_X,g0)
    
    if(k>BIPCM_m)
 
        Y = Y.';
    
        future_length = m-(k+1);

        temp_future_size = size(future_matrix,1)-1;

        future_matrix = future_matrix(temp_future_size-future_length+1:temp_future_size,1:4:end).';

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

        u0 = [temp_past, zeros(temp_present_size,2), temp_future];
        u1 = [temp_past, zeros(temp_present_size,1), ones(temp_present_size,1), temp_future];
        u2 = [temp_past, ones(temp_present_size,1), zeros(temp_present_size,1), temp_future];
        u3 = [temp_past, ones(temp_present_size,2), temp_future];
        
        x0 = mod(u0*g0,2);
        x1 = mod(u1*g0,2);
        x2 = mod(u2*g0,2);
        x3 = mod(u3*g0,2);

        t0 = x0*temp_powers;
        t1 = x1*temp_powers;
        t2 = x2*temp_powers;
        t3 = x3*temp_powers;
        
        t0 = reshape(t0,temp_future_size,temp_Estimated_X_size);
        t1 = reshape(t1,temp_future_size,temp_Estimated_X_size);
        t2 = reshape(t2,temp_future_size,temp_Estimated_X_size);
        t3 = reshape(t3,temp_future_size,temp_Estimated_X_size);

        X0 = Constellation_Mapping_Array(t0+1).';
        X1 = Constellation_Mapping_Array(t1+1).';
        X2 = Constellation_Mapping_Array(t2+1).';
        X3 = Constellation_Mapping_Array(t3+1).';

        exp_multiplier = (1/(2*Sigma^2));
        
        ll_p0 = log(sum(exp(-((bsxfun(@minus,temp_real,real(X0))).^2+(bsxfun(@minus,temp_imag,imag(X0))).^2)*exp_multiplier),2)).';
        ll_p1 = log(sum(exp(-((bsxfun(@minus,temp_real,real(X1))).^2+(bsxfun(@minus,temp_imag,imag(X1))).^2)*exp_multiplier),2)).';
        ll_p2 = log(sum(exp(-((bsxfun(@minus,temp_real,real(X2))).^2+(bsxfun(@minus,temp_imag,imag(X2))).^2)*exp_multiplier),2)).';
        ll_p3 = log(sum(exp(-((bsxfun(@minus,temp_real,real(X3))).^2+(bsxfun(@minus,temp_imag,imag(X3))).^2)*exp_multiplier),2)).';
        
        ll_p0(isinf(ll_p0)) = intmin();
        ll_p1(isinf(ll_p1)) = intmin();
        ll_p2(isinf(ll_p2)) = intmin();
        ll_p3(isinf(ll_p3)) = intmin();
        
    else
        
        Y = Y.';
 
        future_matrix = future_matrix(1:end-2,1:4:end).';

        N = length(Y);

        U_m = m/2;
        
        ll_p0 = zeros(U_m,N);
        ll_p1 = zeros(U_m,N);
        ll_p2 = zeros(U_m,N);
        ll_p3 = zeros(U_m,N);

        temp_real = real(Y);
        temp_imag = imag(Y);

        temp_powers = (2.^(m-1:-1:0)).';

        temp_zeros = zeros(length(future_matrix(:,1)),2);

        exp_multiplier = (1/(2*Sigma^2));
        
        for k=1:2:m

            u = [future_matrix(:,1:k-1),temp_zeros,future_matrix(:,k:end)];
            t0 = u*temp_powers;
            t_01 = 2^(m-(k+1));
            t_10 = 2^(m-k);
            t1 = t0+t_01;
            t2 = t0+t_10;
            t3 = t0+t_10+t_01;

            X0 = Constellation_Mapping_Array(t0+1);
            X1 = Constellation_Mapping_Array(t1+1);
            X2 = Constellation_Mapping_Array(t2+1);
            X3 = Constellation_Mapping_Array(t3+1);

            p0 = sum(exp(-((bsxfun(@minus,temp_real,real(X0))).^2+(bsxfun(@minus,temp_imag,imag(X0))).^2)*exp_multiplier),2);
            p1 = sum(exp(-((bsxfun(@minus,temp_real,real(X1))).^2+(bsxfun(@minus,temp_imag,imag(X1))).^2)*exp_multiplier),2);
            p2 = sum(exp(-((bsxfun(@minus,temp_real,real(X2))).^2+(bsxfun(@minus,temp_imag,imag(X2))).^2)*exp_multiplier),2);
            p3 = sum(exp(-((bsxfun(@minus,temp_real,real(X3))).^2+(bsxfun(@minus,temp_imag,imag(X3))).^2)*exp_multiplier),2);

            index = (k+1)/2;
        
            ll_p0(index,:) = log(p0);
            ll_p1(index,:) = log(p1);
            ll_p2(index,:) = log(p2);
            ll_p3(index,:) = log(p3);
            
        end

        ll_p0 = ll_p0(:)';
        ll_p1 = ll_p1(:)';
        ll_p2 = ll_p2(:)';
        ll_p3 = ll_p3(:)';

        ll_p0(isinf(ll_p0)) = intmin();
        ll_p1(isinf(ll_p1)) = intmin();
        ll_p2(isinf(ll_p2)) = intmin();
        ll_p3(isinf(ll_p3)) = intmin();
        
    end
    
 end