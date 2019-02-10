function llr = calcLLR(y,m,noiseVar)
% This Function calculates llr of the recieved signal, according to the
% relevant constellations

switch m
    case 1
        llr = 4*real(y)/noiseVar;
    case 2
        temp_llr = [ 2*real((y(:)).')/noiseVar;2*imag((y(:)).')/noiseVar];
        llr = temp_llr(:);
end

end

