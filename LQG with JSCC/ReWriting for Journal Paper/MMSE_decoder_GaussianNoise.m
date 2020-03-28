function [X_hat] = MMSE_decoder_GaussianNoise(RxSig,y,rho,sigmaW,f)

dx = 0.001;
x = (-6:dx:6);

% Encoding function Improved
f_x = f;

expTerm1 = -1*(0.5*(x).^2);
expTerm2 = - (0.5/(1-rho^2))*(y - rho*x).^2;
expTerm3 = - 0.5*(1/sigmaW^2)*(RxSig - f_x).^2;
Integrand = exp(expTerm1 + expTerm2 + expTerm3);

num = sum(Integrand.*x,2);
den = sum(Integrand);
X_hat = num/den;

end