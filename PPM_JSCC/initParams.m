function [Params] = initParams()
% Parameters initialization tracking

%% model coefficients
Params.SNR = 5:0.25:15;  % 8              % Direct channel SNR
Params.snrLin = 10.^(Params.SNR/10);

Params.maxStages = 3; % number of PPM levels

%% Modulo parameters:
% zooming factor 
Params.etaLinear = [1.05*sqrt(12) 8 10];
Params.etaPPM = [1.2*sqrt(12) 9 10];

% interval 
Params.IntLin = [sqrt(12)/2.125 sqrt(12)/2.125 sqrt(12)/2.125];

%% Simulation parameters
Params.N_avg = 2^17;


end


