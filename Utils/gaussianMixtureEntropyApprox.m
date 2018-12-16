function [Entropy] = gaussianMixtureEntropyApprox(centers,probs,sigma,ord)
% This function approximate entropy of gaussian mixture distribution
% It gets the mixutre model probabilites , and approximate the integral.
% For now , we use the same variance for all gaussians
% The function uses the gauss-hermite approximation
% Input :
%  o centers : centers of gaussians , interpreted in rows
%  o probs   : probability of each gaussian , interpreted as rows
%  o sigma   : gaussians variance

% initiate approximation related values
if ord == 20
    x=[-5.38748089,-4.60368244,-3.94476404,-3.34785456,-2.78880605,-2.25497400,...
        -1.73853771,-1.23407621,-0.73747372,-0.2453407,0.2453407,0.73747372,...
        1.23407621,1.73853771,2.25497400,2.78880605,3.34785456,3.94476404,...
        4.60368244,5.38748089];
    w=[2.229e-13,4.399e-10,1.086e-07,7.802e-06,2.283e-04,0.00324377,0.02481052,...
        0.10901720,0.28667550,0.46224266,0.46224266,0.28667550,0.10901720,...
        0.02481052,0.00324377,2.283e-04,7.802e-06,1.086e-07,4.399e-10,2.229e-13];    
else
    x=[-4.688738939,-3.869447905,-3.176999162,-2.546202158,-1.951787991,-1.380258539,-0.822951449,-0.273481046,0.273481046 ...
        0.822951449,1.380258539,1.951787991,2.546202158,3.176999162,3.869447905,4.688738939];
    w=[2.65e-10,2.32e-07,2.71e-05,0.000932284,0.012880312,0.083810041,0.280647459,0.507929479,0.507929479 ...
        0.280647459,0.083810041,0.012880312,0.000932284,2.71e-05,2.32e-07,2.65e-10];
end
N = length(w);

% arrays initialization
Entropy = zeros(size(probs,1),1);

% approximation : iterate over gauss-hermite roots , and sum
% H = sum over (probs .* g)/sqrt(pi)

% iterate over gaussians probability
for i=1:size(probs,2)
    % iterate over gauss-hermite roots
    curr_g = zeros(size(probs,1),1);
    for j=1:N
        currInnerSum = zeros(size(probs,1),1);
        for k=1:size(probs,2)
            currInnerSum = currInnerSum + (probs(:,k)/(sqrt(2*pi)*sigma)).*exp(-0.5*((sqrt(2)*sigma*x(j) + centers(:,i) - centers(:,k)).^2)/(sigma^2));
        end
%         if currInnerSum ~= 0
            curr_g = curr_g + w(j)*log2(currInnerSum);
            if isnan(curr_g)
                disp('WTF');
            end
%         end
    end
    Entropy = Entropy + probs(:,i).*curr_g/sqrt(pi);
end



end

