function [Estimated_U,Estimated_X,Estimated_L] = gArikan_GF4_SC_Decoder(Lambda,U,N,Is_Frozen_bit_Index_Vec,Is_Genie_Aided)

    persistent gArikan_GF4_gf2

    gArikan_GF4_gf2 = [1     0     0     0;...
                        0     1     0     0;...
                        1     0     0     1;...
                        0     1     1     1];
            
    if(N==4)
           
        Estimated_U = zeros(1,4);
        Estimated_L = zeros(4,2);
        Estimated_X = zeros(1,4);
        
        r = [1:1:2];
        t_gf2 = [0,0;0,1;1,0;1,1];
        
        for l=1:1:2
            
            w = Create_All_Options(2*(2-l),1);

            v = mod(w*gArikan_GF4_gf2(2*l+1:4,:),2);

            lambda_temp = v;

            if(~size(lambda_temp,1))
                lambda_temp = zeros(1,4);
            end
            
            t_temp = mod(t_gf2*gArikan_GF4_gf2(2*l-1:2*l,:),2);
            
            temp_Add = mod(bsxfun(@plus,Estimated_X,t_temp),2);
            
            x = mod(repmat(lambda_temp,4,1) + kron(temp_Add,ones(size(lambda_temp,1),1)),2);
            x_gf4 = 2*x(:,1:2:end)+x(:,2:2:end);

            minus_arg_sum = reshape(-sum(Lambda(1+bsxfun(@plus,x_gf4,4*(r-1))),2),size(lambda_temp,1),[]);
            
            if(size(minus_arg_sum,1)>1)
                max_arg_sum = max(minus_arg_sum);
                Estimated_L(:,l) = (max_arg_sum+log(sum(exp(bsxfun(@minus,minus_arg_sum,max_arg_sum)))));
            else
                Estimated_L(:,l) = minus_arg_sum;
            end
            
            if(Is_Frozen_bit_Index_Vec(l))
                Estimated_U(2*l-1:2*l) = U(2*l-1:2*l);
            else
                [~,Estimated_U_gf4] = max(Estimated_L(:,l));
                Estimated_U(1,2*l-1:2*l) = [(Estimated_U_gf4 - 1 - mod(Estimated_U_gf4 - 1,2))/2,mod(Estimated_U_gf4 - 1,2)];
            end
            
            if(Is_Genie_Aided)
                Estimated_X = mod(U(1:2*l)*gArikan_GF4_gf2(1:2*l,:),2);
            else
                Estimated_X = mod(Estimated_U(1:2*l)*gArikan_GF4_gf2(1:2*l,:),2);
            end
            
        end
        
        Estimated_L = bsxfun(@minus,Estimated_L(1,:),Estimated_L);
        
    else
        
        N_gf4 = N/2;
        next_N = N/2;
        next_N_gf4 = N_gf4/2;
        
        Estimated_U = zeros(1,N);
        Estimated_L = zeros(4,N_gf4);
        Estimated_X = zeros(1,N);
        
        L = zeros(4,next_N_gf4);
        
        r = [1:1:2];
        
        for l=1:1:2

            w = Create_All_Options(2*(2-l),1);

            v = mod(w*gArikan_GF4_gf2(2*l+1:end,:),2);

            lambda_temp = v;

            if(~size(lambda_temp,1))
                lambda_temp = zeros(1,4);
            end

            for t_gf4=0:1:3

                t = [(t_gf4-mod(t_gf4,2))/2,mod(t_gf4,2)];

                t_temp = mod(t*gArikan_GF4_gf2(2*l-1:2*l,:),2);

                for j=1:1:next_N_gf4
                
                    x = mod(bsxfun(@plus,lambda_temp,Estimated_X(1,4*(j-1)+1:4*j)+t_temp),2);
                    x_gf4 = 2*x(:,1:2:end)+x(:,2:2:end);

                    minus_arg_sum = reshape(-sum(Lambda(1+bsxfun(@plus,x_gf4,4*((j-1)*2+r-1))),2),size(lambda_temp,1),[]);

                    max_arg_sum = max(minus_arg_sum);
                    L(1+t_gf4,j) = (max_arg_sum+log(sum(exp(minus_arg_sum-max_arg_sum))));

                end
                
            end

            L = bsxfun(@minus,L(1,:),L);
            
            [Estimated_U(1,(l-1)*next_N+1:l*next_N),Estimated_X_temp,Estimated_L(:,(l-1)*next_N_gf4+1:l*next_N_gf4)] = gArikan_GF4_SC_Decoder(L,U(1,(l-1)*next_N+1:l*next_N),next_N,Is_Frozen_bit_Index_Vec(1,(l-1)*next_N_gf4+1:l*next_N_gf4),Is_Genie_Aided);
            
            for j=1:1:next_N_gf4
            
                Estimated_X(1,4*(j-1)+1:4*j) = mod(Estimated_X(1,4*(j-1)+1:4*j) + Estimated_X_temp(2*j-1:2*j)*gArikan_GF4_gf2(2*l-1:2*l,:),2);
            
            end
            
        end
                    
    end 
    
end
