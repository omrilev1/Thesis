function [H,curr_LDPC_r,ldpc_dec] = ldpcGen(curr_ldpc_rate,Ensemble,curr_decMethod)

% build LDPC file
if strcmp(Ensemble,'Gallager')
    
    switch curr_ldpc_rate
        case 0.25
            ldpc_file = '.\Ensembles\Gallager\rate_0.25\MacKay_13928.3296';
        case 0.33
            ldpc_file = '.\Ensembles\Gallager\rate_0.33\MacKay_1920.1280.3.303';
        case 0.5
            ldpc_file = '.\Ensembles\Gallager\rate_0.5\MacKay_20000.r0.5';
        case 0.82
            ldpc_file = '.\Ensembles\Gallager\rate_0.82\MacKay_4095.738.4.102';
        case 0.87
            ldpc_file = '.\Ensembles\Gallager\rate_0.87\MacKay_16383.2131.4.104';
        case 0.93
            ldpc_file = '.\Ensembles\Gallager\rate_0.93\MacKay_32000.2240.3.105';
    end
elseif strcmp(Ensemble,'Irregular_BEC')
    switch curr_ldpc_rate
        case 0.33
            ldpc_file = '.\Ensembles\Irregular_optimized\BEC\rate_0.33\tenBrink__16383.v16.c8';
        case 0.5
            ldpc_file = '.\Ensembles\Irregular_optimized\BEC\rate_0.5\tenBrink__36000.v16.c8';
        case 0.66
            ldpc_file = '.\Ensembles\Irregular_optimized\BEC\rate_0.66\tenBrink_1500.v8.c10';
        case 0.82
            ldpc_file = '.\Ensembles\Irregular_optimized\BEC\rate_0.82\tenBrink__4095.v12.c20';
        case 0.87
            ldpc_file = '.\Ensembles\Irregular_optimized\BEC\rate_0.87\tenBrink_32000.v10.c30';
        case 0.93
            ldpc_file = '.\Ensembles\Irregular_optimized\BEC\rate_0.9\tenBrink_16384.v10.c30';
    end
elseif strcmp(Ensemble,'Irregular_AWGN')
    switch curr_ldpc_rate
        case 0.33
            ldpc_file = '.\Ensembles\Irregular_optimized\AWGN\rate_0.33\tenBrink__16383.v16';
        case 0.5
            ldpc_file = '.\Ensembles\Irregular_optimized\AWGN\rate_0.5\tenBrink__16383.v16.c8';
        case 0.66
            ldpc_file = '.\Ensembles\Irregular_optimized\AWGN\rate_0.66\tenBrink_32000.v15.c12';
        case 0.82
            ldpc_file = '.\Ensembles\Irregular_optimized\AWGN\rate_0.82\tenBrink__16384.v8';
        case 0.87
            ldpc_file = '.\Ensembles\Irregular_optimized\AWGN\rate_0.87\tenBrink_32000.v10.c30';
        case 0.93
            ldpc_file = '.\Ensembles\Irregular_optimized\AWGN\rate_0.9\tenBrink_16384.v10.c30';
    end
end

if strcmp(Ensemble,'DVB')
        H = dvbs2ldpc(curr_ldpc_rate);
else
    [H,~] = parse_alist(ldpc_file,'alist',0);
end
curr_LDPC_r = (size(H,2) - size(H,1))/size(H,2);
if curr_decMethod == 0
    ldpc_dec = comm.LDPCDecoder(H,'DecisionMethod','Hard decision','MaximumIterationCount',63,...
        'OutputValue','Whole codeword','NumIterationsOutputPort',1,'IterationTerminationCondition','Parity check satisfied');
else
    ldpc_dec = comm.LDPCDecoder(H,'DecisionMethod','Soft decision','MaximumIterationCount',63,...
        'OutputValue','Whole codeword','NumIterationsOutputPort',1,'IterationTerminationCondition','Parity check satisfied');
end

end

