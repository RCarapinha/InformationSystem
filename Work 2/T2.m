%% Trabalho 2
close all
clear all
clc

% INPUT PARAMETERS FOR THE TRANSMISSION SYSTEM
R = 100e6;              % transmission rate in bit/s
F = 0.75;               % raise-cosine roll-off factor
NSYM = 40000;           % number of symbols to be simulated
NUM = 0.1;              % fraction of the symbols used to estimate the BER
N = 16;                 % number of samples per symbol
noisePSD = -66;         % noise power spectral density (bilateral) (dBm/Hz)
Nfilt = 2;              % receiver filter order
%Bfilt = 0.2*R;         % cut-off frequency of the receiver filter at -3 dB (Hz)
Bfilt = 0.7*R;         % cut-off frequency of the receiver filter at -3 dB (Hz)
%Bfilt = 2.0*R;         % cut-off frequency of the receiver filter at -3 dB (Hz)

% INPUT PARAMETERS FOR GRAPHICAL REPRESENTATION
NT = 32;                % nº of symbols to be visualized in the temporal diagrams
NEye = 256;             % nº of symbols to be visualized in the eye diagram
fs = N*R;                                       % Sampling Frequency
noisePower = noisePSD + 10*log10(fs) - 30;      % dBm to dBW
sigma2_n = 10^(noisePower/10);

%% 3.1.1
figure(1);
signal = randi([0 1],1,NSYM); % information signal
[signalTransmitted,t] = rcosflt(signal,R,N*R,'fir/normal',F);% formatted
noise = wgn(length(signalTransmitted),1,noisePower); % noisePower in dBW
[n, x] = hist(noise,1000);
dx = abs(x(2)-x(1));
bar(x,n/(sum(n)*dx),'r');
hold on
plot(x, 1/sqrt(2*pi*sigma2_n)*exp(-x.^2/(2*sigma2_n)),'k');% sigma2_n is the theoretical value of the noise variance (average power)
hold off;
title('Noise Normalized Histogram');
xlabel('Noise (W^{1/2})');
ylabel ('Normalized Relative Frequency')
view([90 -90]);

%% 3.1.3
M_Noise = mean(noise);
MP_Noise = mean(noise.^2);

%% 3.2
signalReceived = signalTransmitted + noise;
signalReceived_amost = zeros(1,NSYM);
[~,atraso1] = correlacao(signal,signalTransmitted,N,16,16);
for k=1:NSYM % Align and sampled the received signal plus noise
    signalReceived_amost(k) = signalReceived(N*(k-1)+atraso1);
end

%Zero conditioned histogram
signalReceived_amost_0 = signalReceived_amost(signal==0);% Samples corresponding to ZEROs
[n0, x0] = hist(signalReceived_amost_0,1000);
dx0 = abs(x0(2)-x0(1));
bar(x0,n0/(sum(n0)*dx0),'b');
hold on

%ONE conditioned histogram
signalReceived_amost_1 = signalReceived_amost(signal==1); % Samples corresponding to ONEs
[n1, x1] = hist(signalReceived_amost_1,1000);
dx1 = abs(x1(2)-x1(1));
bar(x1,n1/(sum(n1)*dx1),'r');
title('Histograms at input to receiver filter (conditioned to Zero and One)');
xlabel('tensão [W^{1/2}]');
ylabel ('Normalized relative frequency')
hold off
legend('Zeros','Ones')

%% 4.1.a
[b,a] = butter(Nfilt,Bfilt/(fs/2));
Nfft = 512;
[H,f,s] = freqz(b,a,Nfft,fs);
HH = 10*log10(H.*conj(H));
figure;
hold on;
plot(f, HH);
plot(Bfilt*ones(1,100),linspace(min(HH),-3),'g--');
plot(linspace(0, Bfilt), -3*ones(1,100),'r--');
xlim([0, 3*R])
hold off;
xlabel('Frequency');
ylabel('Power')
legend('Frequency Response of Filter','','-3dB')

%% 4.1.b
signalReceived = signalTransmitted + noise;
signalReceivedFiltered = filter(b,a,signalReceived);
noiseFiltered = filter(b,a,noise);
signalTransmittedFiltered = filter(b,a,signalTransmitted);

[~,atraso1] = correlacao(signal,signalTransmitted,N,16,16);
[corr,atraso2] = correlacao(signal,signalTransmittedFiltered,N,16,16);

figure;
subplot(2,1,1);
hold on;
plot(t(1:N*NT),signalTransmitted(1+atraso1:atraso1+N*NT),'g');
plot(t(1:N*NT),signalTransmittedFiltered(1+atraso2:atraso2+N*NT),'r');
plot(linspace(0,max(t(1:N*NT)),10),0.5*ones(1,10),'k--');
title ('Signal without noise (before and after the received filter)')
x = 0:1/R:(NT-1)/R;
stem(x,signal(1:NT),'b');
hold off;
axis([0 max(x) -0.5 1.5]);
legend('Before Filter','After Filter','0.5','Signal')

subplot(2,1,2);
hold on;
plot(t(1:N*NT),signalReceived(1+atraso1:atraso1+N*NT),'g');
plot(t(1:N*NT),signalReceivedFiltered(1+atraso2:atraso2+N*NT),'r');
plot(linspace(0,max(t(1:N*NT)),10),0.5*ones(1,10),'k--');
title ('Signal with noise (before and after the received filter)')
x = 0:1/R:(NT-1)/R;
stem(x,signal(1:NT),'b');
hold off;
axis([0 max(x) -0.5 1.5]);

%% 4.2
eyediagram(signalTransmittedFiltered(1+atraso2+N/2:NEye*N+atraso2+N/2),2*N, 2);
title('Filtered signal without noise');
grid on;

eyediagram(signalReceivedFiltered(1+atraso2+N/2:NEye*N+atraso2+N/2),2*N, 2);
title('Filtered signal with noise');
grid on;

%% 4.3
s_amost = zeros(1,NSYM);
for k=1:NSYM
    s_amost(k) = signalReceivedFiltered(N*(k-1)+1+atraso2);
end

%% 4.4
decisionThreshold = 0.5; % Tem de ser entre '0' e '1'
s_reg = s_amost > decisionThreshold;

%% 4.5.a
err = sum(s_reg~=signal);
realBER = err/NSYM;
fprintf('Nº de bits errados = %g, num total de %g\n',err,NSYM);
fprintf('BER (REAL) = %f\n', realBER);
err = sum(s_reg(1:NUM*NSYM)~=signal(1:NUM*NSYM));
Pe_mc = err/(NUM*NSYM);
fprintf('Nº de bits errados = %g, num total de %g\n',err,NUM*NSYM);
fprintf('BER (Monte-Carlo) = %f, Erro Associado = %3.0f%% \n',Pe_mc, abs((Pe_mc-realBER)/realBER)*100); % percentage difference

%% 4.5.c
Bn_filt = (pi/(2*Nfilt))/sin(pi/(2*Nfilt))*Bfilt; %
DEP_uni_ruido_wpHz = 2*10^((noisePSD-30)/10);
No = DEP_uni_ruido_wpHz*Bn_filt;
Pe_a = 0.5*erfc(0.5/sqrt(2*No));
fprintf('Potência de Ruído depois do Filtro (Analítico) = %f W \n', No);
fprintf('Potência de Ruído depois do Filtro (Numérico) = %f W\n', std(noiseFiltered)^2);
fprintf('BER (Analítico) = %f, Erro Associado = %3.0f%% \n',Pe_a, abs((Pe_a-realBER)/realBER)*100); % percentage difference

%% 4.5.d
N_0 = 0; % nº de simbolos '0'
M1_0 = 0; % momento de primeira ordem para o '0'
M2_0 = 0; % momento de segunda ordem para o '0'
N_1 = 0; % nº de simbolos '1'
M1_1 = 0; % momento de primeira ordem para o '1'
M2_1 = 0; % momento de segunda ordem para o '1'
for k = 1:NUM*NSYM
    if signal(k)==0
        N_0 = N_0+1;
        M1_0 = M1_0+s_amost(k);
        M2_0 = M2_0+s_amost(k)^2;
    else
        N_1 = N_1+1;
        M1_1 = M1_1+s_amost(k);
        M2_1 = M2_1+s_amost(k)^2;
    end
end
med0 = M1_0/N_0;
sig0 = sqrt(M2_0/N_0-med0^2);
med1 = M1_1/N_1;
sig1 = sqrt(M2_1/N_1-med1^2);
Qsg = (med1-med0)/(sig1+sig0);
Pe_g = 0.5*erfc(Qsg/sqrt(2));
fprintf('BER (Q-Factor) = %f, Erro Associado = %3.0f%% \n',Pe_g, abs((Pe_g - realBER)/realBER)*100);

%% 4.5.e
% Monte-Carlo
figure;
[n, x] = hist(s_amost,100);
dx = abs(x(2)-x(1));
bar(x,n/(sum(n)*dx)); 
hold on;
y0 = 1/2*1/sqrt(2*pi*No)*exp(-x.^2/(2*No));
y1 = 1/2*1/sqrt(2*pi*No)*exp(-(x-1).^2/(2*No));
plot(x, y0,'g');
plot(x, y1,'r');
plot(x, y0+y1, 'b');
hold off;
title('Histograma (Analitico)');
xlabel('tensão [W^{-1/2}]');
ylabel ('Frequência relativa')
% Q Factor (Gaussian approximation)
figure;
[n, x] = hist(s_amost,100);
dx = abs(x(2)-x(1));
bar(x,n/(sum(n)*dx)); hold on;
y0=1/2*1/sqrt(2*pi*sig0^2)*exp(-(x-med0).^2/(2*sig0^2));
y1=1/2*1/sqrt(2*pi*sig1^2)*exp(-(x-med1).^2/(2*sig1^2));
plot(x, y0,'g');
plot(x, y1,'r');
plot(x, y0+y1, 'b');
hold off;
title('Histograma (Q-Factor)');
xlabel('tensão [W^{-1/2}]');
ylabel ('Frequência relativa')