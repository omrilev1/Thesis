function [ X_after_encoding, X_after_encoding_gf ] = gGF4_Separated_Polar_Code_Encoder(g,U,Num_Of_Encoded_Symbols_In_Constellation)

    U_N = length(U);

    Bit_Level_U_N = U_N/Num_Of_Encoded_Symbols_In_Constellation;
    
    X_after_encoding = zeros(Num_Of_Encoded_Symbols_In_Constellation,Bit_Level_U_N);

    for BitLevel = 1:Num_Of_Encoded_Symbols_In_Constellation
        X_after_encoding(BitLevel,:) = mod(U((BitLevel-1)*Bit_Level_U_N+1:BitLevel*Bit_Level_U_N)*g,2);
    end
        
    X_after_encoding_gf4 = reshape(2*X_after_encoding(:,1:2:end) + X_after_encoding(:,2:2:end),1,[]);
    
    X_after_encoding_gf = [(X_after_encoding_gf4-mod(X_after_encoding_gf4,2))/2;mod(X_after_encoding_gf4,2)];
        
    X_after_encoding_gf = reshape(X_after_encoding_gf,2*Num_Of_Encoded_Symbols_In_Constellation,[]);

    X_after_encoding = reshape(X_after_encoding_gf,1,[]);
    
end