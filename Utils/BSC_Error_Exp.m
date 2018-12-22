function BSC_Error_Exp(p)
% This code calculates and draws the error exponents for BSC channel 
% the code was built using p=0.2

%% parameters
if (nargin<1) p=0.2; end
n=100;
Rate_Res = 1000;
q_Res = 1000;
q = 0:(1/q_Res):0.5;

% Capacity = 1-IT_Hber(p);
% Divergence = IT_Divergence(q,p);
% Rate = 0:(1/Rate_Res):Capacity;


%% Bouded Distance Decoding
% Ebd(R) = D(q*||p),
% q* is the q that satisfies D(q||p) = 1 - Hb(q) - R
% here the q* is observed as min( D(q||p) - (1 - Hb(q) - R) )

Ebd_q_magshim_arr = zeros(1,Rate_Res+1);
Ebd_Rate_arr = zeros(1,Rate_Res+1);

for ii = 0:Rate_Res
    Ebd_Rate = (ii/Rate_Res)*(1-IT_Hber(p)); % rate between 0 and capacity
    [M,I] = min(abs(IT_Divergence(q,p) - (1 - IT_Hber(q) - Ebd_Rate))) ;
    Ebd_q_magshim = (1/q_Res)*I;
    Ebd_q_magshim_arr(ii+1) = Ebd_q_magshim;
    Ebd_Rate_arr(ii+1) = Ebd_Rate;
end

Ebd_prob = IT_Divergence(Ebd_q_magshim_arr,p);
Ebd = 2.^((-1)*n*Ebd_prob);




%% Maximum Likelyhood Decoding
% Eml(R) = D(q*||p),
% here the q* is observed as min( D(q||p) + (1 - Hb(q) - R) )
% with the exception that (1 - Hb(q) - R) cannot be negative. 
% (see plus functions comment)

Eml_q_magshim_arr = zeros(1,Rate_Res+1);  % no need to build array here... just for debug
Eml_Rate_arr = zeros(1,Rate_Res+1);
Eml_prob =  zeros(1,Rate_Res+1);

for ii = 0:Rate_Res
    Eml_Rate = (ii/Rate_Res)*(1-IT_Hber(p));
    C_minus_R = (1 - IT_Hber(q) - Eml_Rate);
    C_minus_R(C_minus_R<0) = 0;  %% plus fundtion
    [M,I] = min(abs(IT_Divergence(q,p) + C_minus_R)) ;
    Eml_q_magshim = (1/q_Res)*I;   % no need to build array here... just for debug
    Eml_q_magshim_arr(ii+1) = Eml_q_magshim;   % no need to build array here... just for debug
    Eml_Rate_arr(ii+1) = Eml_Rate;
    Eml_prob(ii+1) =M;
end

Eml = 2.^((-1)*n*Eml_prob);

% R0 of Maximum Likelyhood Decoding
R0 = (-1)*log2(0.5*(sqrt(p)+sqrt(1-p))^2);
Eml_Rate_arr_R0 = Eml_Rate_arr(Eml_Rate_arr<=R0);


%% Plot

figure;  
plot(Eml_Rate_arr_R0,R0-Eml_Rate_arr_R0,'c.-','LineWidth',2)
hold all;
plot(Ebd_Rate_arr,Ebd_prob,'g.-','LineWidth',2)
plot(Eml_Rate_arr,Eml_prob,'r--','LineWidth',2)
grid on
xlabel('Rate')
ylabel('E(R)')
title('Error Exponent')
legend('R0-R','Bounded Distance','Maximum Likelihood');

end


%% Functions

function [D]=IT_Divergence(q,p)

element1 = q.*log2(q/p);
element1(isnan(element1)==1)=0;

element2 = (1-q).*log2((1-q)/(1-p));
element2(isnan(element2)==1)=0;

D = element1 + element2;

end


function [Hber]=IT_Hber(p)

Hber = (-1)* ( p.*log2(p) + (1-p).*log2(1-p) );
Hber(isnan(Hber)==1)=0;
end

