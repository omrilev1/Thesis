function [X_after_compound, X_before_compound, X_after_encoding_gf] = gGF4_Compound_Polar_Code_Encoder(g,g0,U,Num_Of_Encoded_Symbols_In_Constellation)

    U_N = length(U);

    Bit_Level_U_N = U_N/Num_Of_Encoded_Symbols_In_Constellation;
    
    X_before_compound = zeros(Num_Of_Encoded_Symbols_In_Constellation,Bit_Level_U_N);

    for BitLevel = 1:Num_Of_Encoded_Symbols_In_Constellation
        X_before_compound(BitLevel,:) = mod(U((BitLevel-1)*Bit_Level_U_N+1:BitLevel*Bit_Level_U_N)*g,2);
    end
    
    X_before_compound_gf4 = reshape(2*X_before_compound(:,1:2:end) + X_before_compound(:,2:2:end),1,[]);
    
    X_before_compound_gf = [(X_before_compound_gf4-mod(X_before_compound_gf4,2))/2;mod(X_before_compound_gf4,2)];
    
    X_after_encoding_gf = reshape(X_before_compound_gf,2*Num_Of_Encoded_Symbols_In_Constellation,[]);
    
    X_before_compound_gf2 = reshape(X_before_compound_gf,size(g0,1),[]).';
    
    X_after_compound = mod(X_before_compound_gf2*g0,2);
   
    X_before_compound = reshape(X_before_compound,1,[]);
    
    X_after_compound = reshape(X_after_compound.',1,[]);
    
end

