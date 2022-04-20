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

a = [0.785398163397448309615660845820, 1.00000000000000000000000000000, ...
1.76714586764425869663523690309, 4.00000000000000000000000000000, ...
11.0446616727766168539702306443, 36.0000000000000000000000000000, ...
135.297105491513556461135325393, 576.000000000000000000000000000, ...
2739.76638620314951833799033921, 14400.0000000000000000000000000, ...
82877.9331826452729297242077612, 518400.000000000000000000000000, ...
3.50159267696676278128084777791*10^6, ...
2.54016000000000000000000000000*10^7, ...
1.96964588079380406447047687507*10^8];
vpa(a,64)

kappa = calcCumulantByMoment(a);

%theory result: by edgeworth expansion
snr = linspace(-20,20);
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
axis([-20 20 -0.15 1.15])


%simulation
hold on;
snr = linspace(-20,20,15);
mcResult = zeros(1,length(snr));
simulationTime = 100000;
for i = 1:length(snr)
   xx = sqrt((k2 +1/(10^(snr(i)/10)))*(2^Rate-1));
   h = zeros(simulationTime,1);
   for j=1:n
        h = h + raylrnd(1/sqrt(2),100000,1).*raylrnd(1/sqrt(2),100000,1);
   end
   mcResult(i) = sum(h < xx)/simulationTime;
end
plot(snr, mcResult,'ro');


snr = linspace(-20,20);
ew = zeros(1,length(snr));
n = 4;
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

%simulation
hold on;
snr = linspace(-20,20,15);
mcResult = zeros(1,length(snr));
simulationTime = 100000;
for i = 1:length(snr)
   xx = sqrt((k2 +1/(10^(snr(i)/10)))*(2^Rate-1));
   h = zeros(simulationTime,1);
   for j=1:n
        h = h + raylrnd(1/sqrt(2),100000,1).*raylrnd(1/sqrt(2),100000,1);
   end
   mcResult(i) = sum(h < xx)/simulationTime;
end
plot(snr, mcResult,'o');


snr = linspace(-20,20);
ew = zeros(1,length(snr));
n = 8;
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

%simulation
hold on;
snr = linspace(-20,20,15);
mcResult = zeros(1,length(snr));
simulationTime = 100000;
for i = 1:length(snr)
   xx = sqrt((k2 +1/(10^(snr(i)/10)))*(2^Rate-1));
   h = zeros(simulationTime,1);
   for j=1:n
        h = h + raylrnd(1/sqrt(2),100000,1).*raylrnd(1/sqrt(2),100000,1);
   end
   mcResult(i) = sum(h < xx)/simulationTime;
end
plot(snr, mcResult,'o');

legend('theory:N=3','simulation:N=3','theory:N=4','simulation:N=4','theory:N=8','simulation:N=8');
xlabel('Average SNR[dB]')
ylabel('Average Outage Probability')