function [trellis,conv_enc,conv_dec] = convCodeGen(v,r,rsc)
% convCode is a struct, where input contain the constraint length
% and rate
% the output include the trellis, and viterbi decoder memory length
if rsc
    
    switch r 
        
        case 1/4 
            v = 6;
            poly = [67,51,45,71];
            feedback = [67];
        case 1/3 
            v = 7;
            poly = [163,131,101];
            feedback = [163];
            
        case 1/2            
            v = 7;
            poly = [147,115];
            feedback = [147];
    end
    % generate convolutional code
    trellis = poly2trellis(v,poly,feedback);
else
    
    switch v
        
        case 5
            
            switch r
                
                case 1/4
                    poly = [37,33,27,25];
                case 1/3
                    poly = [37,33,25];
                case 2/5
                    poly = [21,27,0;37,33,25];
                case 1/2
                    poly = [31,27];
            end
            
        case 6
            
            switch r
                
                case 1/4
                    poly = [77,67,55,51];
                case 1/3
                    poly = [75,53,47];
                case 2/5
                    poly = [55,73,0;73,75,51];
                case 1/2
                    poly = [65,57];
            end
        case 7
            
            switch r
                
                case 1/4
                    poly = [171,155,127,117];
                case 1/3
                    poly = [155,127,117];
                case 2/5
                    poly = [157,165,0;157,127,111];
                case 1/2
                    poly = [155,117];
                    
            end
            
        case 8
            
            switch r
                
                case 1/4
                    poly = [235,275,313,357];
                case 1/3
                    poly = [225,331,367];
                case 1/2
                    poly = [247,371];
            end
            
        case 10
            
            switch r
                
                case 1/4
                    poly = [1117,1365,1633,1653];
                case 1/3
                    poly = [1117,1365,1633];
                case 1/2
                    poly = [1167,1545];
            end
        case 12
            
            switch r
                
                case 1/4
                    poly = [4767,5723,6265,7455];
                case 1/3
                    poly = [4767,5723,6265];
                case 1/2
                    poly = [4335,5723];
            end
        case 14
            
            switch r
                
                case 1/4
                    poly = [21113,23175,35577,35537];
                case 1/3
                    poly = [21645,35661,37133];
                case 1/2
                    poly = [21675,27123];
            end
    end
    % generate convolutional code
    trellis = poly2trellis(v,poly);
end
conv_enc = comm.ConvolutionalEncoder('TrellisStructure',trellis);
conv_dec = comm.APPDecoder(...
    'TrellisStructure',trellis,'CodedBitLLROutputPort',false,'Algorithm','True App');
% conv_dec = comm.ViterbiDecoder('TrellisStructure',trellis);

end

