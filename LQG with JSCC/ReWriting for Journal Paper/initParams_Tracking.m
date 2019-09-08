function [Params] = initParams_Tracking()
% Parameters initialization tracking 

%% model coefficients
Params.lambda = 0.95;               % Gauss markov memory parameter
Params.W = 1;                       % input gaussian variance
Params.T = 2^8;                     % Horizon
Params.P = 1;                       % power constraint normalized to 1
Params.P_timesharing = Params.P;    % Power during uncoded transmission in time-sharing 

Params.SNR = 10;  % 14                % Direct channel SNR
Params.deltaSNR = -2;  % 8           % Delta between feedback path and direct path
Params.snrLin = 10.^(Params.SNR/10);
Params.snrFeedback_Lin = 10.^((Params.SNR + Params.deltaSNR)/10);

Params.statePower = zeros(1,Params.T);
for i=1:Params.T
    if i == 1
        Params.statePower(i) = Params.W;
    else
        Params.statePower(i) = (Params.lambda)^2 * Params.statePower(i-1) + Params.W;
    end
end

%% Transmission scheme parameters
Params.deltaKochman = 2.25; % sqrt(3*Params.P);
Params.gammaKochman = 1.05;

Params.sign_cheat = 0; % to transmit in the feedback path only the absolute value, with half of the power

% set timesharing parameters - we need to scale the powers to the overall
% avergae power will be P
Params.uncoded_timesharing_factor = 2^8;


Params.P_alias = 9*1e-3;        % Aliasing probability for Kochman- Encoder
Params.N_feedback = 4;  % 4        % Number of feedback iteration over 1 control samples

Params.deltaTuncel    = 3.2;% 3.15;            % Delta for Tuncel HDA Transmission scheme
partition             = (-6.5:1:6.5)*Params.deltaTuncel;    
Params.centersTuncel  = [partition(1) (partition(1:end-1) + partition(2:end))/2 partition(end)];
Params.codebookTuncel = generateGaussianCodeBook(partition);
Params.alpha          = 0.75;% 0.8;             %  Tuncel scheme parameters
Params.beta           = -1.18;% -1.15;           % Tuncel scheme parameters;



%% Simulation parameters
Params.N_avg = 200;
Params.NumOfSchemes = 8;
% Dictionary : 
% 1 = FullAccess
% 2 = ZeroAccess
% 3 = ZeroAccess with Kochman Encoder
% 4 = ZeroAccess with Tuncel Encoder
% 5 = partialAccess with modulo based feedback
% 6 = partialAccess with linear based feedback 
% 7 = partialAccess with modulo based reversed feedback
% 8 = partialAccess with linear based reversed feedback 

end

function [codebook] = generateGaussianCodeBook(partition)

dx = 0.001;
x = (-3 + partition(1)):dx:(partition(end) + 3);

pdf = (1/sqrt(2*pi))*exp(-0.5*x.^2);
codebook = zeros(1,length(partition) - 1);
for i=1:(length(partition) - 1)
    
    curIdx = find((x >= partition(i)) & (x <= partition(i+1)));
    curPDF = pdf(curIdx) / sum(pdf(curIdx));
    curX = x(curIdx);
    codebook(i) = sum(curX.*curPDF);
end

% generate first and last codewords

% first
curIdx = find(x <= partition(1));
curPDF = pdf(curIdx) / sum(pdf(curIdx));
curX = x(curIdx);
codebook = [sum(curX.*curPDF) codebook];

% last codeword
curIdx = find(x >= partition(end));
curPDF = pdf(curIdx) / sum(pdf(curIdx));
curX = x(curIdx);
codebook = [codebook sum(curX.*curPDF)];
end