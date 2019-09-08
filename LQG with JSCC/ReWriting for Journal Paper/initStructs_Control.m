function [Results,Debug,PredictErr,EstimErr] = initStructs_Control(Params)
% This function initialize the results and debug sturcts for the simulation

% Dictionary : 
% 1 = FullAccess
% 2 = ZeroAccess
% 3 = ZeroAccess with Kochman Encoder
% 4 = ZeroAccess with Tuncel Encoder
% 5 = partialAccess with modulo based feedback
% 6 = partialAccess with linear based feedback 

Results = zeros(Params.NumOfSchemes,length(Params.SNR),Params.T,Params.N_avg);
Debug = zeros(Params.NumOfSchemes,length(Params.SNR),Params.T,Params.N_avg);
PredictErr = zeros(Params.NumOfSchemes,length(Params.SNR),Params.T,Params.N_avg);
EstimErr = zeros(Params.NumOfSchemes,length(Params.SNR),Params.T,Params.N_avg);


end

