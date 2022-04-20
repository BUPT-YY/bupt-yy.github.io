%simulation and verify for thm 3.
clc,clear, close all;
digits(64);
figure; hold on;

pd=makedist('Weibull','a',sqrt(0.5/gamma(1.2)),'b',10);

n = 4; k2 = 10^(-4); simulationTime = 100000;
h = zeros(simulationTime,1);
for j=1:n
        h = h + random(pd,simulationTime,1).*random(pd,simulationTime,1);
end


alpha = 1; beta = 2;
snr = linspace(-20,5,15);
result = zeros(1,length(snr));
for i = 1:length(snr)
    for j = 1:simulationTime
        result(i) = result(i) + alpha * qfunc(sqrt(beta)*h(j)/sqrt(k2 +1/(10^(snr(i)/10)))); 
    end
end
result = result / simulationTime;
plot(snr,result,'o');


alpha = 2; beta = 1;
snr = linspace(-20,5,15);
result = zeros(1,length(snr));
for i = 1:length(snr)
    for j = 1:simulationTime
        result(i) = result(i) + alpha * qfunc(sqrt(beta)*h(j)/sqrt(k2 +1/(10^(snr(i)/10)))); 
    end
end
result = result / simulationTime;
plot(snr,result,'o');


alpha = 4-sqrt(2); beta = 3/7;
snr = linspace(-20,5,15);
result = zeros(1,length(snr));
for i = 1:length(snr)
    for j = 1:simulationTime
        result(i) = result(i) + alpha * qfunc(sqrt(beta)*h(j)/sqrt(k2 +1/(10^(snr(i)/10)))); 
    end
end
result = result / simulationTime; result = min(result,ones(size(result)));
plot(snr,result,'o');




%Theory results
a = [0.702044, 0.5, 0.360655, 0.263116, 0.193938, 0.144293, 0.108282, ...
0.0819066, 0.0624142, 0.0478892];
a = a.*a;
vpa(a,64);
kappa = calcCumulantByMoment(a);
snr = linspace(-20,5);
ew = zeros(1,length(snr));
n = 4;
miu = a(1); r = 6; Rate = 1; 
sigma = sqrt(a(2) - a(1)*a(1));
M = 40;
phim = zeros(1,M);
for m=1:M
   phim(m) =  cos((2*m-1)*pi/2/M);
end

ew = zeros(1,length(snr));
alpha = 1; beta = 2;
for i = 1:length(snr)
    snrReal = 10^(snr(i)/10);
    paxb = sqrt(snrReal/(k2*snrReal+1));
     for m=1:M
         xx1 = sqrt(log(2/(phim(m)+1))*2/beta)/paxb;
         xx = (xx1-n*miu)/(sqrt(n)*sigma);
         coef = 0;
         for k=3:r
             coef = coef + (n^(-(k-2)/2))*kappa(k)/(sigma^k)/(factorial(k)) * probHermiteH(k-1,xx);
         end
         yy = normcdf(xx,0,1) - normpdf(xx,0,1) * coef;
         ew(i) = ew(i) + yy * alpha*sqrt(pi)/4/M*sqrt(1-phim(m)^2)/sqrt(log(2/(phim(m)+1)));
     end
    
end
plot(snr,ew);

ew = zeros(1,length(snr));
alpha = 2; beta = 1;
for i = 1:length(snr)
    snrReal = 10^(snr(i)/10);
    paxb = sqrt(snrReal/(k2*snrReal+1));
     for m=1:M
         xx1 = sqrt(log(2/(phim(m)+1))*2/beta)/paxb;
         xx = (xx1-n*miu)/(sqrt(n)*sigma);
         coef = 0;
         for k=3:r
             coef = coef + (n^(-(k-2)/2))*kappa(k)/(sigma^k)/(factorial(k)) * probHermiteH(k-1,xx);
         end
         yy = normcdf(xx,0,1) - normpdf(xx,0,1) * coef;
         ew(i) = ew(i) + yy * alpha*sqrt(pi)/4/M*sqrt(1-phim(m)^2)/sqrt(log(2/(phim(m)+1)));
     end
    
end
plot(snr,ew);



ew = zeros(1,length(snr));
alpha = 4-sqrt(2); beta = 3/7;
for i = 1:length(snr)
    snrReal = 10^(snr(i)/10);
    paxb = sqrt(snrReal/(k2*snrReal+1));
     for m=1:M
         xx1 = sqrt(log(2/(phim(m)+1))*2/beta)/paxb;
         xx = (xx1-n*miu)/(sqrt(n)*sigma);
         coef = 0;
         for k=3:r
             coef = coef + (n^(-(k-2)/2))*kappa(k)/(sigma^k)/(factorial(k)) * probHermiteH(k-1,xx);
         end
         yy = normcdf(xx,0,1) - normpdf(xx,0,1) * coef;
         ew(i) = ew(i) + yy * alpha*sqrt(pi)/4/M*sqrt(1-phim(m)^2)/sqrt(log(2/(phim(m)+1)));
     end
    
end
ew = min(ew,1); plot(snr,ew); 

axis([-20 6 -0.15 1.15])







xlabel('Average SNR[dB]')
ylabel('Average Symbol Error Rate')
legend('BPSK-simulation', '4-QAM-simulation', '8-QAM-simulation','BPSK-theory','4QAM-theory','8QAM-theory');