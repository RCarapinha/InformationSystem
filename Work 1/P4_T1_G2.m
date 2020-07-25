%% Begin
clear all;                                        % clears all variables
close all;                                        % closes all windows
clc

% TRANSMISSION SYSTEM VARIABLES
R = 2.5e9;                                        % transmission rate in bit/s
F = 0.05;                                         % raised cosine roll-off factor
A = 12;                                           % amplitude of isolated raised cosine pulses
NoisePower = 4;                                   % average noise power (W)

% SIMULATION VARIABLES
NSYM = 20000;                                     % number of simulated symbols
N = 16;                                           % number of samples per symbol

% DISPLAY VARIABLES
NT = 32;                                          % number of symbols to display in temporal diagrams

%% 3 -> Transmitter
%% 3.1 -> Signal Generation
%% 3.1.1 -> Delay Estimation

fa = N*R;                                         % sampling frequency
inf = randi([0 1],1,NSYM);                        % bit seq. with equal prob. for '0' and '1'
[s_tx0,t] = rcosflt(inf,R,fa,'fir/normal',F, 4);  % raised cosine pulses
s_tx = A*s_tx0;

figure
plot(t(1:N*NT)*R,s_tx(1:N*NT),'k');               % time diag. for the formatted signal
title('Sequences with raised cosine formatting')
hold on
x = 0:1:(NT-1);
stem(x,1.3*A*inf(1:NT),'r'); shg
hold off

%% 3.1.2 -> Delay Calculation

[corr,atraso] = correlacao(inf,s_tx,N,64,16);
figure
plot(0:length(corr)-1,corr);

%% 3.1.3 -> Delay Calculation

figure
plot(t(1:N*NT),s_tx(1+atraso:atraso+N*NT),'b');
title('Original and raised cosine sequences after delay adjustment ')
hold on;
x = 0:1/R:(NT-1)/R;
stem(x,1.3*A*inf(1:NT),'r');
hold off 

%% 3.2 -> Eye Diagram
%% 3.2.1 -> Delay Calculation

eyediagram(s_tx(1:256*N),2*N, 2);

%% 3.3 -> Power Spectral Density
%% 3.3.1

figure
windowL = floor(length(s_tx)/10);                                       % Segment length in Welch's method
windowT = hamming(windowL);                                             % Hamming window type (Welch's method)
[Pxx,f] = pwelch(s_tx, windowT, windowL/2, windowL, fa, 'onesided');
Pxx=Pxx/2;                                                              % To compensate for the fact that function pwelch performs
% unilateral PSD
plot(f,10*log10(Pxx/1e-3),'r');
hold on
FR = f(2)-f(1);
FR_theorical = fa/windowL;
xlabel('Hz'); ylabel('dBm/Hz'); 
xlim([0 2.5*R]); ylim([-130 -40]);

%% 3.3.2
T = 1/(R);
y = ones(length(f),1);

for i = 1:1:length(f)
    if abs(f(i)) < ((1-F)/(2*T))
        y(i) = A*T;
    elseif abs(f(i)) >= ((1-F)/(2*T)) && abs(f(i)) <= ((1+F)/(2*T))
        y(i) = A*T/2*(1-sin((pi*T/F)*(abs(f(i))-1/(2*T))));
    elseif abs(f(i)) > ((1+F)/(2*T))
        y(i) = 0;
    end
end
plot(f,10*log10((abs(y)).^2/(4*T)./(1*10^-3)),'b');
title('Experimental and Theoretical PSD F=0.05');
hold off

%% 3.3.3
    inf1 = zeros(1,NSYM);                        % bit seq. with equal prob. for '0' and '1'
    inf1(4)=1;
    %inf1(10)=1;
    %inf1(100)=1;
    [s_tx10,t1] = rcosflt(inf1,R,fa,'fir/normal',F, 4);  % raised cosine pulses
    s_tx1 = A*s_tx10;

    figure
    windowL1 = floor(length(s_tx1)/10);                                       % Segment length in Welch's method
    windowT1 = hamming(windowL1);                                               % Hamming window type (Welch's method)
    [Pxx1,f1] = pwelch(s_tx1, windowT, windowL1/2, windowL1, fa, 'onesided');
    Pxx1=Pxx1/2;                                                              % To compensate for the fact that function pwelch performs unilateral PSD
    plot(f1,10*log10(Pxx1/1e-3),'r');
    hold on
    FR = f1(2)-f1(1);
    xlabel('Hz'); ylabel('dBm/Hz'); 
    %xlim([0 2.5*R]); ylim([-130 -40]);

%% 3.3.4
%% 3.4 -> Impact of the raised cosine roll-off factor
%% 3.4.1

F = 1;
[s_tx_1,t] = rcosflt(inf,R,fa,'fir/normal',F, 4);  % raised cosine pulses
s_tx_2 = A*s_tx_1;
eyediagram(s_tx_2(1:256*N),2*N, 2);

%% 3.4.2
figure
hold on
F = 1;
t = linspace(-3*T,3*T,1000);
h1 = zeros(length(t));
for j = 1:length(t)
    if t(j) == -(3*T)/(2*F) || t(j) == -(2*T)/(2*F) || t(j) == -T/(2*F) || t(j) == T/(2*F) ||t(j) == (2*T)/(2*F) ||t(j) == (3*T)/(2*F)
        h1(j) = (pi/(4*T))*sinc(1/(2*F));
    else
        h1(j) = (1/T)*sinc(t(j)/T)*((cos((pi*F*t(j))/(T)))/(1-((2*F*t(j))/(T))^2));
    end
end
plot(t,h1,'r--');

F = 0.05;
h2 = zeros(length(t));
for j = 1:length(t)
    if t(j) == -(3*T)/(2*F) || t(j) == -(2*T)/(2*F) || t(j) == -T/(2*F) || t(j) == T/(2*F) ||t(j) == (2*T)/(2*F) ||t(j) == (3*T)/(2*F)
        h2(j) = (pi/(4*T))*sinc(1/(2*F));
    else
        h2(j) = (1/T)*sinc(t(j)/T)*((cos((pi*F*t(j))/(T)))/(1-((2*F*t(j))/(T))^2));
    end
end
plot(t,h2,'b');
grid on
legend('F = 1','F = 0.05')
xlabel('Time');
ylabel('h(t)');
hold off

%% 3.4.3
F = 1;
figure
hold on
windowL = floor(length(s_tx_2)/10);                                       % Segment length in Welch's method
windowT = hamming(windowL);                                               % Hamming window type (Welch's method)
[Pxx,f] = pwelch(s_tx_2, windowT, windowL/2, windowL, fa, 'onesided');
Pxx=Pxx/2;                                                              % To compensate for the fact that function pwelch performs
% unilateral PSD
plot(f,10*log10(Pxx/1e-3),'r');
FR1 = f(2)-f(1);
xlabel('Hz'); ylabel('dBm/Hz'); 
xlim([0 2.5*R]); ylim([-130 -40]);

T = 1/R;
y2 = ones(length(f),1);
for i = 1:1:length(f)
    if abs(f(i)) < ((1-F)/(2*T))
        y2(i) = A*T;
    elseif abs(f(i)) >= ((1-F)/(2*T)) && abs(f(i)) <= ((1+F)/(2*T))
        y2(i) = A*T/2*(1-sin((pi*T/F)*(abs(f(i))-1/(2*T))));
    elseif abs(f(i)) > ((1+F)/(2*T))
        y2(i) = 0;
    end
end
plot(f,10*log10((abs(y2)).^2/(4*T)./(1*10^-3)),'b');
title('Experimental and Theoretical PSD F=1');
hold off

%% 4 -> Noise
%% 4.1 -> Gausian Noise Generation

F = 0.05;
NoisePower_dB = 10*log10(NoisePower);
noise = wgn(length(s_tx),1,NoisePower_dB); % NoisePower_dB is the average power expressed in dBW (dBs relative the Watt). 

%% 4.2 -> Spectral Characterization of Noise
%% 4.2.1

figure
[Pxx,f1] = pwelch(noise, windowT, windowL/2, windowL, fa, 'onesided'); % windowL and windowT are defined above for the signal PSD
Pxx_dB = 10*log10(Pxx)+30-3;
plot(f1, Pxx_dB); 
FR1 = f1(2)-f1(1);
var1 = var(Pxx);
avg1 = mean(Pxx);

%% 4.2.2

avg1_theorical = 2*(NoisePower/fa);
ratio_theorical = fa*2*avg1_theorical;

%% 4.2.3

r1 = var1/avg1;

%% 4.2.4

NSYM = 10000;
F=0.05;
inf = randi([0 1],1,NSYM);                        % bit seq. with equal prob. for '0' and '1'
[s_tx0,t] = rcosflt(inf,R,fa,'fir/normal',F, 4);  % raised cosine pulses
s_tx = A*s_tx0;

NoisePower_dB = 10*log10(NoisePower);
noise = wgn(length(s_tx),1,NoisePower_dB); % NoisePower_dB is the average power expressed in dBW (dBs relative the Watt). 

figure
[Pxx1,f2] = pwelch(noise, windowT, windowL/2, windowL, fa, 'onesided'); % windowL and windowT are defined above for the signal PSD
Pxx1_dB = 10*log10(Pxx1)+30-3;
plot(f2, Pxx1_dB); 
FR2 = f2(2)-f2(1);
var2 = var(Pxx1);
avg2 = mean(Pxx1);

r2 = var2/avg2;