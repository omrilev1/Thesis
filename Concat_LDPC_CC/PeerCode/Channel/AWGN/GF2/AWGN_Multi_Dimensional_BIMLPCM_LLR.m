 function [lambda] = AWGN_Multi_Dimensional_BIMLPCM_LLR(m,BIPCM_m,Constellation_Mapping_Array,Y,Sigma,future_matrix,k,Estimated_X,g0)
    
    if(k>BIPCM_m)
 
        Y = Y.';
    
        future_length = m-k;

        temp_future_size = size(future_matrix,1);

        future_matrix = future_matrix(temp_future_size-future_length+1:temp_future_size,1:2^future_length).';

        temp_future_size = size(future_matrix,1);

        temp_Estimated_X_size = size(Estimated_X,3);

        Constellation_2D_Dimension = size(Constellation_Mapping_Array,1);

        N = length(Y);

        U_N_div_m = N/Constellation_2D_Dimension;

        temp_real = reshape(real(Y),Constellation_2D_Dimension,[]).';
        temp_imag = reshape(imag(Y),Constellation_2D_Dimension,[]).';

        temp_powers = (2.^(m-1:-1:0)).';

        past = zeros(k-1,U_N_div_m);

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

        exp_multiplier = (1/(2*Sigma^2));
        p0 = 0;
        p1 = 0;
        for i=1:1:Constellation_2D_Dimension
            Constellation_Mapping_Array_temp = Constellation_Mapping_Array(i,:);
            X0 = Constellation_Mapping_Array_temp(t0+1).';
            X1 = Constellation_Mapping_Array_temp(t1+1).';
            p0 = p0+((bsxfun(@minus,temp_real(:,i),real(X0))).^2+(bsxfun(@minus,temp_imag(:,i),imag(X0))).^2);
            p1 = p1+((bsxfun(@minus,temp_real(:,i),real(X1))).^2+(bsxfun(@minus,temp_imag(:,i),imag(X1))).^2);
        end
        p0 = sum(exp(-p0*exp_multiplier),2);
        p1 = sum(exp(-p1*exp_multiplier),2);

        lambda = log(p0./p1).';
        
    else
        
        Y = Y.';
 
        future_matrix = future_matrix.';

        Constellation_2D_Dimension = size(Constellation_Mapping_Array,1);

        N = length(Y);

        U_N_div_m = N/Constellation_2D_Dimension;

        lambda = zeros(m,U_N_div_m);

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

            lambda(k,:) = log(p0./p1);

        end

        lambda = lambda(:)';
        
    end
    
 end