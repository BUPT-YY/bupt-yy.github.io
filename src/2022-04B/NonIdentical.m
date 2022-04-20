clc,clear,close all;
figure; hold on;

%Rician Distr.s in different parameters.
pd{1} = makedist('Rician','s',sqrt(5/6),'sigma',1/sqrt(12));
pd{2} = makedist('Rician','s',sqrt(6/7),'sigma',1/sqrt(14));
pd{3} = makedist('Rician','s',sqrt(5*1.5/6),'sigma',sqrt(1.5)/sqrt(12));
pd{4} = makedist('Rician','s',sqrt(6*1.5/7),'sigma',sqrt(1.5)/sqrt(14));

snr = linspace(-20,5,20);
mcResult = zeros(1,length(snr));
simulationTime = 100000;
n=4; k2 = 0.01; Rate = 1;
for i = 1:length(snr)
   xx = sqrt((k2 +1/(10^(snr(i)/10)))*(2^Rate-1));
   h = zeros(simulationTime,1);
   for j=1:n
        h = h + random(pd{j},simulationTime,1).*random(pd{j},simulationTime,1);
   end
   mcResult(i) = sum(h < xx)/simulationTime;
end
plot(snr, mcResult,'o');



%Weibull Distr.s in different parameters.
pd{1} = makedist('Weibull','a',sqrt(0.5/gamma(1.2)),'b',10);
pd{2} = makedist('Weibull','a',sqrt(1/gamma(1.2)),'b',10);
pd{3} = makedist('Weibull','a',sqrt(0.5/gamma(1.2)),'b',12);
pd{4} = makedist('Weibull','a',sqrt(1/gamma(1.2)),'b',12);
snr = linspace(-20,5,20);
mcResult = zeros(1,length(snr));
simulationTime = 100000;
n=4; k2 = 0.01; Rate = 1;
for i = 1:length(snr)
   xx = sqrt((k2 +1/(10^(snr(i)/10)))*(2^Rate-1));
   h = zeros(simulationTime,1);
   for j=1:n
        h = h + random(pd{j},simulationTime,1).*random(pd{j},simulationTime,1);
   end
   mcResult(i) = sum(h < xx)/simulationTime;
end
plot(snr, mcResult,'o');


axis([-20 6 -0.15 1.15])




%Theoretical results
snr = linspace(-20,5);
ewRice = zeros(1,length(snr));
ewWeibull = zeros(1,length(snr));


%Moments for every Rician fading channels
aRice1 = [0.959930110752018935758544730598, 1.00000000000000000000000000000, ...
1.11177654821320601342971559345, 1.30555555555555555555555555556, ...
1.60766919039780783665966434452, 2.06481481481481481481481481481, ...
2.75437675388674387876672567442, 3.80324074074074074074074074074, ...
5.42076220057757990402010317735, 7.95653292181069958847736625514, ...
12.0025685773826563282439747515, 18.5762817215363511659807956104, ...
29.4523316201544209273410362094, 47.7723122427983539094650205761, ...
79.1793689528662987735713942951];
aRice2 = [0.965348853621578678468436675781, 1.00000000000000000000000000000, ...
1.09744183586807030791985098660, 1.26530612244897959183673469388, ...
1.52344680796155856179256299064, 1.90670553935860058309037900875, ...
2.47164327335909142379747392953, 3.30862140774677217825905872553, ...
4.56242484472779320714680711141, 6.46730528946272386505622657226, ...
9.40695523110186199731721529398, 14.0182406990284660303104998768, ...
21.3727144401187523907301675902, 33.2980208683699576099851495307, ...
52.9537992469912584758279819686];
aRice3 = [1.20622938159501243249301545976, 1.62500000000000000000000000000, ...
2.37738880211323171916436996616, 3.71875000000000000000000000000, ...
6.15616069054109032523352423414, 10.7070312500000000000000000000, ...
19.4570605400364835726455252789, 36.7832031250000000000000000000, ...
72.0875586055901300531654192039, 146.031494140625000000000000000, ...
305.030772974487382592160736924, 655.603332519531250000000000000, ...
1447.27322850374131073345102828, 3276.28597259521484375000000000, ...
7594.96424428136809428822668045];
aRice4 = [1.20768248136941820772760993571, 1.60714285714285714285714285714, ...
2.30418993114744874104463372372, 3.51275510204081632653061224490, ...
5.64431985837204106324952127070, 9.49772230320699708454810495627, ...
16.6545789419058004518302905351, 30.3149078508954602249062890462, ...
57.0954978837416872720988024201, 110.972691345124055453085024097, ...
222.085769673953677939859618479, 456.746281742991865634217035419, ...
963.713749535386683947772820254, 2083.04171304157577382116052228, ...
4606.34252885947196946098588435];
aRice1 = aRice1 .* aRice1;
aRice2 = aRice2 .* aRice2;
aRice3 = aRice3 .* aRice3;
aRice4 = aRice4 .* aRice4;
cumuRice = zeros(1,length(aRice1));
cumuRice = cumuRice + calcCumulantByMoment(aRice1);
cumuRice = cumuRice + calcCumulantByMoment(aRice2);
cumuRice = cumuRice + calcCumulantByMoment(aRice3);
cumuRice = cumuRice + calcCumulantByMoment(aRice4);
miuRice = cumuRice(1);
sigmaRice = sqrt(cumuRice(2));
for i=1:length(snr)
    xx = (sqrt((k2 +1/(10^(snr(i)/10)))*(2^Rate-1)) - miuRice)/(sigmaRice);
    coef = 0;
    coef = coef + cumuRice(3)/6*probHermiteH(2,xx)/(sigmaRice^3);
    coef = coef + cumuRice(4)/factorial(4)*probHermiteH(3,xx)/(sigmaRice^4);
    coef = coef + cumuRice(5)/factorial(5)*probHermiteH(4,xx)/(sigmaRice^5);
    coef = coef + (cumuRice(3)^2/72)*probHermiteH(5,xx)/(sigmaRice^6);%cumuRice(6)/factorial(6)
    coef = coef + (cumuRice(3)*cumuRice(5)/factorial(3)/factorial(5)+cumuRice(4)^2/1152)*probHermiteH(7,xx)/(sigmaRice^8);
    coef = coef + (cumuRice(3)^3)/(6^4)*probHermiteH(8,xx)/(sigmaRice^9);
    coef = coef + (cumuRice(3)^2)*(cumuRice(4)/1728)*probHermiteH(9,xx)/(sigmaRice^10);
    coef = coef + (cumuRice(3)^4)/31104*probHermiteH(11,xx)/(sigmaRice^12);
    ewRice(i) = normcdf(xx,0,1) - normpdf(xx,0,1) * coef;
end
plot(snr, ewRice);

%Moments for every Weibull Fading channels
aWeibull1 = [0.702044149057774030800908289979, 0.500000000000000000000000000000, ...
0.360654643286712918229296503396, 0.263116311561201047117190589102, ...
0.193938348041420992091624531730, 0.144292752835079773679807835278, ...
0.108282306375652196572320061603, 0.0819065566016783461030467590059, ...
0.0624141551844396101183375252769, 0.0478891912330610626279449383057];
aWeibull2 = [0.992840356982182751943700460599, 1.00000000000000000000000000000, ...
1.02008537573780025019615987133, 1.05246524624480418846876235641, ...
1.09708096825764456877801206893, 1.15434202268063818943846268222, ...
1.22507444993188796295398649030, 1.31050490562685353764874814409, ...
1.41227111593429632438080036437, 1.53245411945795400409423802578];
aWeibull3 = [0.707161730041398840218902541151, 0.505200891072220996652912384003, ...
0.364243939531454591711570573318, 0.264811289211632186430664005563, ...
0.193994191112687185662847614025, 0.143115753760794624581633087823, ...
0.106269325528056659760867479747, 0.0793880683106221290059797395640, ...
0.0596429081507161272271875125511, 0.0450472482805577177349263283005];
aWeibull4 = [1.00007770941576760146912534990, 1.01040178214444199330582476801, ...
1.03023743859517723916176660939, 1.05924515684652874572265602225, ...
1.09739686437264144171377283926, 1.14492603008635699665306470258, ...
1.20230017140815279145801049255, 1.27020909296995406409567583302, ...
1.34956495369784897382355371252, 1.44151194497784696751764250562];
aWeibull1 = aWeibull1 .* aWeibull1;
aWeibull2 = aWeibull2 .* aWeibull2;
aWeibull3 = aWeibull3 .* aWeibull3;
aWeibull4 = aWeibull4 .* aWeibull4;
cumuWeibull = zeros(1,length(aWeibull1));
cumuWeibull = cumuWeibull + calcCumulantByMoment(aWeibull1);
cumuWeibull = cumuWeibull + calcCumulantByMoment(aWeibull2);
cumuWeibull = cumuWeibull + calcCumulantByMoment(aWeibull3);
cumuWeibull = cumuWeibull + calcCumulantByMoment(aWeibull4);
miuWeibull = cumuWeibull(1);
sigmaWeibull = sqrt(cumuWeibull(2));
for i=1:length(snr)
    xx = (sqrt((k2 +1/(10^(snr(i)/10)))*(2^Rate-1)) - miuWeibull)/(sigmaWeibull);
    coef = 0;
    coef = coef + cumuWeibull(3)/6*probHermiteH(2,xx)/(sigmaWeibull^3);
    coef = coef + cumuWeibull(4)/factorial(4)*probHermiteH(3,xx)/(sigmaWeibull^4);
    coef = coef + cumuWeibull(5)/factorial(5)*probHermiteH(4,xx)/(sigmaWeibull^5);
    coef = coef + (cumuWeibull(3)^2/72)*probHermiteH(5,xx)/(sigmaWeibull^6);%cumuWeibull(6)/factorial(6)
    %coef = coef + (cumuWeibull(3)*cumuWeibull(5)/factorial(3)/factorial(5)+cumuWeibull(4)^2/1152)*probHermiteH(7,xx)/(sigmaWeibull^8);
    coef = coef + (cumuWeibull(3)^3)/(6^4)*probHermiteH(8,xx)/(sigmaWeibull^9);
    %coef = coef + (cumuWeibull(3)^2)*(cumuWeibull(4)/1728)*probHermiteH(9,xx)/(sigmaWeibull^10);
    %coef = coef + (cumuWeibull(3)^4)/31104*probHermiteH(11,xx)/(sigmaWeibull^12);
    ewWeibull(i) = normcdf(xx,0,1) - normpdf(xx,0,1) * coef;
end
plot(snr, ewWeibull);

legend('Rician N=4, simulation','Weibull N=4, simulation','Rician, Theory','Weibull, Theory');
xlabel('Average SNR[dB]')
ylabel('Average Outage Probability')