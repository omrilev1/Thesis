function [Params] = initParams_LQGtracking()
% Parameters initialization tracking

%% model coefficients
Params.alpha = 2;                   % Gauss markov memory parameter
Params.W = 1;                       % input gaussian variance
Params.T = 2^10;                     % Horizon
Params.P = 1;                       % power constraint normalized to 1
Params.P_timesharing = Params.P;    % Power during uncoded transmission in time-sharing

Params.Q = 5;                       % LQG Cost parameters
Params.R = 1;                       % LQG Cost parameters

Params.ratio_up = 3;                % Control Cost Uncertainity - upper limit
Params.ratio_down = 1/3;              % Control Cost Uncertainity - upper limit
Params.scale = 'randomized';        % Scaling Method : 'randomized','worst_case'

% calculate LQG parameters
Qmat = Params.Q * [1 0;0 1]; Rmat = Params.R * [1 0;0 1]; 
Params.Amat = [Params.alpha 0;-1 1]; Params.Bmat = [1;0];
[Params.k_mat,Params.s_mat,Params.M_mat] = calcLQG_vector(Qmat,Params.R,Params.T,Params.Amat,Params.Bmat);
[Params.k_mat_up,~,~] = calcLQG_vector(Qmat,Params.R*Params.ratio_up,Params.T,Params.Amat,Params.Bmat);
[Params.k_mat_down,~,~] = calcLQG_vector(Qmat,Params.R*Params.ratio_down,Params.T,Params.Amat,Params.Bmat);
if strcmp(Params.scale,'randomized')
    Params.k_for_power = (Params.k_mat_up - Params.k_mat_down).*rand(1,2,Params.T-1) + Params.k_mat_down;
else
    Params.k_for_power = Params.k_mat_up;
end

Params.SNR = 9;  % 8              % Direct channel SNR
Params.deltaSNR = -1;  % 8          % Delta between feedback path and direct path
Params.snrLin = 10.^(Params.SNR/10);
Params.snrFeedback_Lin = 10.^((Params.SNR + Params.deltaSNR)/10);

Params.rho = 0.9; % correlation of receiver side informaiton

Params.maxRefInput = 0.4; % unknown step input is in the interval [-maxRefInput,maxRefInput]
%% Transmission scheme parameters

Params.sign_cheat = 0; % to transmit in the feedback path only the absolute value, with half of the power

Params.P_alias = 9*1e-3;        % Aliasing probability for Kochman- Encoder
Params.N_feedback = 4;  % 4        % Number of feedback iteration over 1 control samples

Params.deltaTuncel    = 3.2;% 3.15;            % Delta for Tuncel HDA Transmission scheme
partition             = (-6.5:1:6.5)*Params.deltaTuncel;
Params.centersTuncel  = [partition(1) (partition(1:end-1) + partition(2:end))/2 partition(end)];
Params.codebookTuncel = generateGaussianCodeBook(partition);
Params.alphaTuncel    = 0.75;% 0.75;             %  Tuncel scheme parameters
Params.betaTuncel     = -1.18;% -1.18;           % Tuncel scheme parameters;

%% Simulation parameters
Params.N_avg = 64;
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

function [k_mat,s_mat,M_mat] = calcLQG_vector(Q,R,T,A,B)
% calculate k
k_mat = zeros(1,2,T-1);
s_mat = zeros(2,2,T-1);
M_mat = zeros(2,2,T-1);

for t=T:-1:2
    if t==T
        s_mat(:,:,t) = Q;
        k_mat(:,:,t-1) = (B'*s_mat(:,:,t)*A) / (B'*s_mat(:,:,t)*B + R);
    else
        M_mat(:,:,t) = (A'*s_mat(:,:,t+1)*B) * (B'*s_mat(:,:,t+1)*A)/(R + B'*s_mat(:,:,t+1)*B);
        s_mat(:,:,t) = A'*s_mat(:,:,t+1)*A - (A'*s_mat(:,:,t+1)*B) * (B'*s_mat(:,:,t+1)*A)/(R + B'*s_mat(:,:,t+1)*B) + Q;
        k_mat(:,:,t-1) = (B'*s_mat(:,:,t)*A) / (B'*s_mat(:,:,t)*B + R);
    end
end
end
