function [probs,cost] = findBestBeta(alpha,delta,codebook,rho,snr)

%% start with calculating input distribution related parameters : 

% Probabilites Pr{T = k*Delta}
pVec = [qfunc(3*delta/2), qfunc(delta/2) - qfunc(3*delta/2), ...
    1 - 2*qfunc(delta/2), qfunc(delta/2) - qfunc(3*delta/2), qfunc(3*delta/2)];

% quantization error power : E[S^2 | T = k*Delta]
sVec = zeros(size(pVec));
dx = 1e-3;
boundaries = [-3*delta/2,-1*delta/2,delta/2,3*delta/2];
for i=1:length(pVec)
    
    if i==1
        currX = -8:dx:boundaries(i);
    elseif i==5
        currX = boundaries(i-1):dx:8;
    else
        currX = boundaries(i-1):dx:boundaries(i);
    end
    
    currPDF = (1/(sqrt(2*pi)))*exp(-currX.^2 / 2);
    currPDF = currPDF/sum(currPDF);
    
    sVec(i) = sum((currX - codebook(i)).^2 .* currPDF);
    
end

% E[T^2]
intP = pVec*([codebook(1);codebook(2);codebook(3);codebook(4);codebook(5)].^2);

a = rho^2/(1-rho^2);
b = 1/(1-rho^2);
%% Build the optimization objective + constraint

f    = @(x)objectiveFun(x,pVec,sVec,delta,a,b,snr,codebook);
gfun = @(x) deal(constraintFun(x,pVec,sVec,alpha,intP),[]);

x0 =  [1.4 1.15 1.3 1.15 1.4];

options = optimoptions('fmincon','Algorithm','interior-point','Display','iter');
[probs,cost,exitflag,output] = fmincon(f,x0,[],[],[],[],[],[],gfun,options);

end

function y = objectiveFun(x,pVec,sVec,delta,a,b,snr,codebook)
y = pVec(1) * (sVec(1) + a + x(1)^2 * snr + (codebook(1))^2)/(b + x(1)^2 * snr)^2 + ...
    pVec(2) * (sVec(2) + a + x(2)^2 * snr + (codebook(2))^2)/(b + x(2)^2 * snr)^2 + ...
    pVec(3) * (sVec(3) + a + x(3)^2 * snr + (codebook(3))^2)/(b + x(3)^2 * snr)^2 + ...
    pVec(4) * (sVec(4) + a + x(4)^2 * snr + (codebook(4))^2)/(b + x(4)^2 * snr)^2 + ...
    pVec(5) * (sVec(5) + a + x(5)^2 * snr + (codebook(5))^2)/(b + x(5)^2 * snr)^2;
end

function y = constraintFun(x,pVec,sVec,alpha,intP)
y = pVec(1) * sVec(1)*(x(1))^2 + ...
    pVec(2) * sVec(2)*(x(2))^2 + ...
    pVec(3) * sVec(3)*(x(3))^2 + ...
    pVec(4) * sVec(4)*(x(4))^2 + ...
    pVec(5) * sVec(5)*(x(5))^2 + alpha^2*intP - 1;
end
