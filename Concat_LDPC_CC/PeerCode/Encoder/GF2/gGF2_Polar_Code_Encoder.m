function [X] = gGF2_Polar_Code_Encoder(U)
       
    g = [1,0;1,1];
    
    N = length(U);
    
    X = gArikan_Polar_Code_Encoder_Recurtion(U,g,N);
           
end

function [X] = gArikan_Polar_Code_Encoder_Recurtion(U,g,N)
   
    if(N==2)
        
        X = mod(U*g,2);
        
    else
    
        X = mod(reshape(((reshape(U,2,[]))'*g),1,[]),2);

        next_N = N/2;

        X(1:next_N) = gArikan_Polar_Code_Encoder_Recurtion(X(1:next_N),g,next_N);
        X(next_N+1:N) = gArikan_Polar_Code_Encoder_Recurtion(X(next_N+1:N),g,next_N);   
        
    end

end