function [Params] = initParams_Control()
% Parameters initialization tracking

%% model coefficients
Params.alpha = 2;                   % Gauss markov memory parameter
Params.W = 1;                       % input gaussian variance
Params.T = 2^8;                     % Horizon
Params.P = 1;                       % power constraint normalized to 1
Params.P_timesharing = Params.P;    % Power during uncoded transmission in time-sharing

Params.Q = 5;                       % LQG Cost parameters
Params.R = 1;                       % LQG Cost parameters

Params.ratio_up = 3;                % Control Cost Uncertainity - upper limit
Params.ratio_down = 1/3;            % Control Cost Uncertainity - upper limit
Params.scale = 'randomized';        % Scaling Method : 'randomized','worst_case'

% calculate LQG parameters
[Params.k_vec,Params.s_vec] = calcLQG(Params.Q,Params.R,Params.T,Params.alpha);
Params.k_vec_up = calcLQG(Params.Q,Params.Q*Params.ratio_up,Params.T,Params.alpha);
Params.k_vec_down = calcLQG(Params.Q,Params.Q*Params.ratio_down,Params.T,Params.alpha);
if strcmp(Params.scale,'randomized')
    Params.k_for_power = (Params.k_vec_up - Params.k_vec_down).*rand(1,Params.T-1) + Params.k_vec_down;
else
    Params.k_for_power = Params.k_vec_up;
end

Params.SNR = 6;  % 8              % Direct channel SNR
Params.deltaSNR = -1;  % 8          % Delta between feedback path and direct path
Params.snrLin = 10.^(Params.SNR/10);
Params.snrFeedback_Lin = 10.^((Params.SNR + Params.deltaSNR)/10);

Params.rho = 0.9; % correlation of receiver side informaiton

%% Transmission scheme parameters
Params.deltaKochman = 3.8; % sqrt(3*Params.P);
Params.gammaKochman = 1.025;

Params.sign_cheat = 0; % to transmit in the feedback path only the absolute value, with half of the power

% set timesharing parameters - we need to scale the powers to the overall
% avergae power will be P
Params.uncoded_timesharing_factor = 2;


Params.P_alias = 9*1e-3;        % Aliasing probability for Kochman- Encoder
Params.N_feedback = 4;  % 4        % Number of feedback iteration over 1 control samples

Params.deltaTuncel    = 3.2;% 3.15;            % Delta for Tuncel HDA Transmission scheme
partition             = (-6.5:1:6.5)*Params.deltaTuncel;
Params.centersTuncel  = [partition(1) (partition(1:end-1) + partition(2:end))/2 partition(end)];
Params.codebookTuncel = generateGaussianCodeBook(partition);
Params.alphaTuncel    = 0.75;% 0.75;             %  Tuncel scheme parameters
Params.betaTuncel     = -1.18;% -1.18;           % Tuncel scheme parameters;

%% Simulation parameters
Params.N_avg = 1024;
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

function [k_vec,s_vec] = calcLQG(Q,R,T,alpha)
% calculate k
k_vec = zeros(1,T-1);
s_vec = zeros(1,T-1);
for t=T:-1:2
    if t==T
        s_vec(t) = Q;
        k_vec(t-1) = alpha*s_vec(t)/(s_vec(t) + R);
    else
        s_vec(t) = ((alpha^2)*R*s_vec(t+1))/(s_vec(t+1) + R) + Q;
        k_vec(t-1) = alpha*s_vec(t)/(s_vec(t) + R);
    end
end
end
