function [Entropy] = complexMixtureEntropyApprox(sigma,Delta,alpha,xAxis)
% This function approximate entropy of mixture distribution
% where the "continous" part pdf is weighted sum of uniform and gaussian:
% N = (alpha-1)*X + alpha*Z
% and the overall pdf is K + N

% Input :
%  o centers : centers of gaussians , interpreted in rows
%  o probs   : probability of each gaussian , interpreted as rows
%  o sigma   : gaussians variance
%  o Delta   : uniform distribution is between [-Delta,Delta]


% arrays initialization
Entropy = zeros(size(probs,1),1);
dX = mean(diff(xAxis));
% iterate over gaussians probability , for each probability construct the
% pdf and calculate the entropy
for i=1:size(probs,1)
    currCenters = centers(i,:);
    currProbs = probs(i,:);
    
    %% generate pdf
    % generate composite noise
    alphaTimesNoise = (1/alpha)*(1/sqrt(2*pi*sigma^2)) * exp(-0.5*(xAxis.^2)/((alpha^2)*(sigma^2)));
    alphaTimesInput = zeros(size(xAxis));
    alphaTimesInput(xAxis >= -1*Delta*(1-alpha) & xAxis <= Delta*(1-alpha)) = 1/(2*Delta*(1-alpha));
    equivNoiseDPC = conv(alphaTimesNoise,alphaTimesInput) * dX;
    equivNoiseDPC = equivNoiseDPC(ceil(length(alphaTimesNoise)/2) : end - ceil(length(alphaTimesNoise)/2));
    
    % generate mixture
    intPart_pdf = zeros(size(xAxis));
    [~,idx1] = min(abs(xAxis - currCenters(1)));
    [~,idx2] = min(abs(xAxis - currCenters(2)));
    [~,idx3] = min(abs(xAxis - currCenters(3)));
    
    intPart_pdf(idx1) = currProbs(1);
    intPart_pdf(idx2) = currProbs(2);
    intPart_pdf(idx3) = currProbs(3);
    intPart_pdf = sparse(intPart_pdf);
    
    mixture_pdf = sparseconv(intPart_pdf,equivNoiseDPC);
    
    Entropy(i) = -1*sum(mixture_pdf(mixture_pdf > 0).*log2(mixture_pdf(mixture_pdf > 0))) * dX;
    
    if (mod(i,100) == 0)
        display(num2str(i));
    end
end



end

