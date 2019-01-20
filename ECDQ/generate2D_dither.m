function dither = generate2D_dither(Delta,nPoints,method,latticePoints)
% method 1 : take random uniform points over a big rectangular , and
% start eliminate only those that are inside the basic cell
% method 2 : using construction-A , where we take coarser lattice,
% quantize to the original lattice points and take the error

% generate uniform dither
x_cord = linspace(-2*sqrt(3) * Delta,2*sqrt(3) * Delta,sqrt(30*nPoints));
y_cord = linspace(-2*sqrt(3) * Delta,2*sqrt(3) * Delta,sqrt(30*nPoints));
[meshX,meshY] = meshgrid(x_cord,y_cord);
uniDither = meshX(:) + 1j*meshY(:);

switch method
    
    case 1
        %     uniDither = 4*sqrt(3) * Delta * ((rand(1,nPoints^2) - 0.5) + 1j*(rand(1,nPoints^2) - 0.5));
        
        % start eliminate points
        
        % take only points which lie between y = [-Delta,Delta]
        step1 = uniDither(abs(imag(uniDither)) < Delta);
        
        step2 = step1((imag(step1) < -2*sqrt(3) * real(step1) + 3*Delta) & ...
            (imag(step1) > -2*sqrt(3) * real(step1) - 3*Delta));
        
        step3 = step2((imag(step2) < 2*sqrt(3) * real(step2) + 3*Delta) & ...
            (imag(step2) > 2*sqrt(3) * real(step2) - 3*Delta));
        
        dither = ([real(step3(:)) imag(step3(:))]).';
        
    case 2
        uniDither = [real(uniDither.') ; imag(uniDither.')];
        error = sum(abs(uniDither - reshape(latticePoints,2,1,[])).^2,1);
        [~,minIdx] = min(error,[],3);
        
        % calculate the error relative to the lattice points
        dither = uniDither - latticePoints(:,minIdx);
        
        % take only nPoints
        tempIdx = randperm(size(dither,2),nPoints);
        dither = dither(:,tempIdx);
end
        
        
        
        
end