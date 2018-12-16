function [InstCap] = PAMCapacity(SNR,probVec,PAMsize)
%% Capacity Calculation of PAM constellations with different a-priori probabilites
% The calculation is with the Gauss-Hermite quadrature rules , for
% approximating integrals of gaussian functionss :
% http://en.wikipedia.org/wiki/Gauss%E2%80%93Hermite_quadrature
% Input : 
%   o SNR     - snr for current capacity calculation
%   o probVec - a-priori probabilites of PAM symbols
%   o PAMsize - order of PAM constellation

% initiate approximation related values
x=[
-4.688738939
-3.869447905
-3.176999162
-2.546202158
-1.951787991
-1.380258539
-0.822951449
-0.273481046
0.273481046
0.822951449
1.380258539
1.951787991
2.546202158
3.176999162
3.869447905
4.688738939
];
w=[
2.65e-10
2.32e-07
2.71e-05
0.000932284
0.012880312
0.083810041
0.280647459
0.507929479
0.507929479
0.280647459
0.083810041
0.012880312
0.000932284
2.71e-05
2.32e-07
2.65e-10
];

%% SNR affect 
% We assume noise power is normalized. Thus the snr will take part in
% constellation values
SNRlin = 10^(SNR/10);
sigAmp = sqrt(SNRlin);

%% Generate PAM constellation
pamSymbols = zeros(1,PAMsize);
for i = 1 : PAMsize/2
    pamSymbols(i) = -1*(2*i - 1);
    pamSymbols(i + PAMsize/2) = 2*i - 1;
end
NormFactor = sqrt(PAMsize/sum(abs (pamSymbols).^2));
Constellation = NormFactor.*h.Constellation;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C1 = zeros(1,16);  % BPSK Gaussi Hermite
C2 = zeros(16,16); % Two dimensional signals
Capforx = zeros(1,QAMsize);
if (QAMsize==2) % BPSK
    SNRaffBPSK = sqrt(2).*SNRaff;  % see page 363 of Digital Communications 5th A/sigma
    for m=1:16
        C1(m)=1/2*(w(m)*(1/sqrt(pi))*log2(2/(1+exp(-2*(sqrt(2)*x(m)+SNRaffBPSK)*SNRaffBPSK)))+...
                       w(m)*(1/sqrt(pi))*log2(2/(1+exp( 2*(sqrt(2)*x(m)-SNRaffBPSK)*SNRaffBPSK))));
    end
    InstCap = sum(C1);

elseif (QAMsize==4)||(QAMsize==16)||(QAMsize==64)||(QAMsize==8)||(QAMsize==32)
    for xindex = 1:QAMsize
        for m1=1:16
            for m2=1:16
                sumoverxprime = 0;
                for xprimeindex=1:QAMsize
                    sumoverxprime = sumoverxprime + ...
                    exp(-abs(SNRaff.*(Constellation(xindex)-Constellation(xprimeindex))+x(m1)+sqrt(-1).*x(m2)).^2 ...
                    +x(m1).^2+x(m2).^2);
                end
                C2(m1,m2)=1/(pi)*w(m1)*w(m2)*log2(sumoverxprime);
            end
        end
        Capforx(xindex) = sum(sum(C2,1),2);
    end
    InstCap = log2(QAMsize) - 1/QAMsize *sum (Capforx);               
else
    error('the modulation size %d is not supported!',  QAMsize);
end

    
    
    