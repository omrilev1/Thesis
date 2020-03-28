clc; 
rho = 0.999;
sigma = 0.3;

X = linspace(-7,7,201);
indx0 = find(X==0);
dx = X(2)-X(1);
A = exp(-(rho^2)*(X'-X).^2/(2*sigma^2));
A = A - eye(size(A));
fx = 1/sqrt(2*pi)*exp(-X.^2/2);
B = dx^2*A.*fx';
B((B<10^-6)) = 0;
B = sparse(B);
C = B(indx0,:)+B(:,indx0)';
%f = rand(size(X))-0.5;
%f = f/sqrt(sum(f.^2.*fx*dx));
f = mod(X + 0.25,0.5) - 0.25;

%loss = @(f)
for j = 0:1000
for i = 0:100
    %difff = sum(B.*(-(f'-f)/sigma^2).*exp(-(f'-f).^2/(2*sigma^2)),2);
    %theSqrt = sqrt(sqrt(2*pi)/dx-sqrt(2*pi)*sum(fx.*f.^2));
    %difff = difff + C'.*exp(-(theSqrt-f').^2/(2*sigma^2)).*...
    %    ((-1/sigma^2)*(theSqrt-f')).*...
    %    ((-sqrt(2*pi)*fx'.*f')./theSqrt - 1);
    difff = sum(B.*(-(f'-f)/sigma^2).*exp(-(f'-f).^2/(2*sigma^2)),2)/dx;
    difff = difff + 1*(2*f'.*fx').*(sum(fx.*f.^2)*dx-1)/dx;
    f = f - 0.001*difff';
    %f(indx0) = sqrt(sqrt(2*pi)*(1/dx - sum(f.^2.*fx) + f(indx0)^2*fx(indx0)));
    
end
loss = lossFunction(B,f,sigma)
figure(2)
plot(X,f)
ylim([-3,3])
drawnow 
end
plot(X,f)

function loss = lossFunction(B,f,sigma)
loss = sum(sum(B.*exp(-(f'-f).^2/(2*sigma^2))));
end


