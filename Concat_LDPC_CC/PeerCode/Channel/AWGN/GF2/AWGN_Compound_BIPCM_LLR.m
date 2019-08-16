 function [lambda] = AWGN_Compound_BIPCM_LLR(m,Constellation_Mapping_Array,Y,Sigma,future_matrix,Is_Second_Stage,Estimated_X,g0)

    Y = Y.';
    
    half_m = m/2;
    
    if(Is_Second_Stage)
        future_length = half_m-1;
    else
        future_length = m-1;
    end
    
    temp_future_size = size(future_matrix,1);
    
    future_matrix = future_matrix(temp_future_size-future_length+1:temp_future_size,1:2^future_length).';
    
    temp_future_size = size(future_matrix,1);
    
    temp_Estimated_X_size = size(Estimated_X,2);
    
    N = length(Y);
    
    lambda = zeros(half_m,N);
    
    temp_real = real(Y);
    temp_imag = imag(Y);
    
    temp_powers = (2.^(m-1:-1:0)).';

    exp_multiplier = (1/(2*Sigma^2));
    
    if(Is_Second_Stage)
                
        for k=half_m+1:1:m
    
            past = Estimated_X.';
            temp_past = kron(past, ones(temp_future_size,1));
            temp_present_size = temp_Estimated_X_size*temp_future_size;
            temp_future = repmat(future_matrix, [size(past, 1), 1]);

            temp_u0 = [zeros(temp_present_size,1) ,temp_future];
            temp_u1 = [ones(temp_present_size,1) ,temp_future];

            u0 = zeros(size(temp_past,1),2*size(temp_past,2));
            u1 = zeros(size(temp_past,1),2*size(temp_past,2));

            u0(:,2:2:end) = [temp_u0(:,2:k-half_m),temp_u0(:,1),temp_u0(:,k-half_m+1:end)];
            u1(:,2:2:end) = [temp_u1(:,2:k-half_m),temp_u1(:,1),temp_u1(:,k-half_m+1:end)];

            u0(:,1:2:end) = temp_past;
            u1(:,1:2:end) = temp_past;

            temp_x0 = mod(u0*g0,2);
            temp_x1 = mod(u1*g0,2);
            
            x0 = [temp_x0(temp_x0(:,k)==0,:);temp_x1(temp_x1(:,k)==0,:)];
            x1 = [temp_x0(temp_x0(:,k)==1,:);temp_x1(temp_x1(:,k)==1,:)];

            t0 = x0*temp_powers;
            t0 = reshape(t0,temp_future_size,temp_Estimated_X_size);
            t1 = x1*temp_powers;
            t1 = reshape(t1,temp_future_size,temp_Estimated_X_size);

            X0 = Constellation_Mapping_Array(t0+1).';
            X1 = Constellation_Mapping_Array(t1+1).';
            
            p0 = sum(exp(-((bsxfun(@minus,temp_real,real(X0))).^2+(bsxfun(@minus,temp_imag,imag(X0))).^2)*exp_multiplier),2);
            p1 = sum(exp(-((bsxfun(@minus,temp_real,real(X1))).^2+(bsxfun(@minus,temp_imag,imag(X1))).^2)*exp_multiplier),2);
            
            lambda(k,:) = log(p0./p1);
            
        end

    else
        
        temp_zeros = zeros(length(future_matrix(:,1)),1);
        
        for k=1:1:m
            
            u = [future_matrix(:,1:k-1),temp_zeros,future_matrix(:,k:end)];
            t0 = u*temp_powers;
            t1 = t0+2^(m-k);

            X0 = Constellation_Mapping_Array(t0+1);
            X1 = Constellation_Mapping_Array(t1+1);

            p0 = sum(exp(-((bsxfun(@minus,temp_real,real(X0))).^2+(bsxfun(@minus,temp_imag,imag(X0))).^2)*exp_multiplier),2);
            p1 = sum(exp(-((bsxfun(@minus,temp_real,real(X1))).^2+(bsxfun(@minus,temp_imag,imag(X1))).^2)*exp_multiplier),2);

            lambda(k,:) = log(p0./p1);

        end
            
    end
    
    lambda = lambda(:)';
    
 end