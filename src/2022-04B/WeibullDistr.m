% clc,clear,close all;
% x = linspace(-2,2);
% y = zeros(10,length(x));
% n=10;
% for i=1:n
%     for j=1:length(x)
%         y(i,j) = probHermiteH(i,x(j));
%     end 
% end
% for i=1:n
%    figure;
%    plot(x,y(i,:))
% end
clc,clear,close all;
digits(64)

a = [0.702044, 0.5, 0.360655, 0.263116, 0.193938, 0.144293, 0.108282, ...
0.0819066, 0.0624142, 0.0478892];
a = a.*a;
pd=makedist('Weibull','a',sqrt(0.5/gamma(1.2)),'b',10);

vpa(a,64)

kappa = calcCumulantByMoment(a);

%theory result: by edgeworth expansion
snr = linspace(-20,5);
ew = zeros(1,length(snr));
n = 3;
miu = a(1);
r = 6; Rate = 1; k2 = 10^(-2);
sigma = sqrt(a(2) - a(1)*a(1));
for i=1:length(snr)
    xx = (sqrt((k2 +1/(10^(snr(i)/10)))*(2^Rate-1)) - n*miu)/(sqrt(n)*sigma);
    coef = 0;
    for k=3:r
        coef = coef + (n^(-(k-2)/2))*kappa(k)/(sigma^k)/(factorial(k)) * probHermiteH(k-1,xx);
    end
    ew(i) = normcdf(xx,0,1) - normpdf(xx,0,1) * coef;
end
plot(snr,ew);
axis([-20 6 -0.15 1.15])


%simulation
hold on;
snr = linspace(-20,5,20);
mcResult = zeros(1,length(snr));
simulationTime = 100000;
for i = 1:length(snr)
   xx = sqrt((k2 +1/(10^(snr(i)/10)))*(2^Rate-1));
   h = zeros(simulationTime,1);
   for j=1:n
        h = h + random(pd,simulationTime,1).*random(pd,simulationTime,1);
   end
   mcResult(i) = sum(h < xx)/simulationTime;
end
plot(snr, mcResult,'ro');


snr = linspace(-20,5);
ew = zeros(1,length(snr));
n = 4;
miu = a(1);
Rate = 1; k2 = 10^(-2);
sigma = sqrt(a(2) - a(1)*a(1));
for i=1:length(snr)
    xx = (sqrt((k2 +1/(10^(snr(i)/10)))*(2^Rate-1)) - n*miu)/(sqrt(n)*sigma);
    coef = 0;
    for k=3:r
        coef = coef + (n^(-(k-2)/2))*kappa(k)/(sigma^k)/(factorial(k)) * probHermiteH(k-1,xx);
    end
    ew(i) = normcdf(xx,0,1) - normpdf(xx,0,1) * coef;
end
plot(snr,ew);

%simulation
hold on;
snr = linspace(-20,5,20);
mcResult = zeros(1,length(snr));
simulationTime = 100000;
for i = 1:length(snr)
   xx = sqrt((k2 +1/(10^(snr(i)/10)))*(2^Rate-1));
   h = zeros(simulationTime,1);
   for j=1:n
        h = h + random(pd,simulationTime,1).*random(pd,simulationTime,1);
   end
   mcResult(i) = sum(h < xx)/simulationTime;
end
plot(snr, mcResult,'o');


snr = linspace(-20,5);
ew = zeros(1,length(snr));
n = 8;
miu = a(1);
Rate = 1; k2 = 10^(-2);
sigma = sqrt(a(2) - a(1)*a(1));
for i=1:length(snr)
    xx = (sqrt((k2 +1/(10^(snr(i)/10)))*(2^Rate-1)) - n*miu)/(sqrt(n)*sigma);
    coef = 0;
    for k=3:r
        coef = coef + (n^(-(k-2)/2))*kappa(k)/(sigma^k)/(factorial(k)) * probHermiteH(k-1,xx);
    end
    ew(i) = normcdf(xx,0,1) - normpdf(xx,0,1) * coef;
end
plot(snr,ew);

%simulation
hold on;
snr = linspace(-20,5,20);
mcResult = zeros(1,length(snr));
simulationTime = 100000;
for i = 1:length(snr)
   xx = sqrt((k2 +1/(10^(snr(i)/10)))*(2^Rate-1));
   h = zeros(simulationTime,1);
   for j=1:n
        h = h + random(pd,simulationTime,1).*random(pd,simulationTime,1);
   end
   mcResult(i) = sum(h < xx)/simulationTime;
end
plot(snr, mcResult,'o');

legend('theory:N=3','simulation:N=3','theory:N=4','simulation:N=4','theory:N=8','simulation:N=8');
xlabel('Average SNR[dB]')
ylabel('Average Outage Probability')