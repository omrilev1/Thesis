% Power after modulo

sigma = [1 1.5 2 3 4 5 6];

delta = linspace(0.1,10,100);

P = zeros(1,length(delta));

figure; hold all
for i=1:length(sigma)
    for j=1:length(delta)
        
        x = sqrt(sigma(i))*randn(1,1e6);
        
        x_mod = mod(delta(j) + x,2*delta(j)) - delta(j);
        
        P(i,j) = mean(x_mod.^2);
        
    end
    plot(delta,P(i,:),'LineWidth',1.5)
end
grid on; grid minor;
legend('\sigma^2 = 1','\sigma^2 = 1.5','\sigma^2 = 2','\sigma^2 = 3',...
    '\sigma^2 = 4','\sigma^2 = 5','\sigma^2 = 6')
xlabel('\Delta'); ylabel('P_{[x]_{\Delta}}')