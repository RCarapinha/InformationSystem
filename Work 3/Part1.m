%% Begin
clear all
close all
clc

m=1000;       % number of symbols to simulate
n=100;      %number of samples per symbol
KI=8;       %number of levels of the baseband signal modulating the in-phase carrier
KQ=8;        %number of levels of the baseband signal modulating the quadrature carrier
d=2;        %difference between consecutive levels in the baseband modulating signals
Rc=5;       %ration between carrier frequency and symbol rate.
Tbit=6;
M=KI*KQ;
k=log2(M);

%% 2.1
xd=d/2*(mary(KI,m,n)); % signal modulating the in-phase carrier
xq=d/2*(mary(KQ,m,n)); % signal modulating the quadrature carrier
xd = reshape(xd,[m*n,1]);
xq = reshape(xq,[m*n,1]);

%% 2.2
%a)
figure
plot(xd,xq,'o')
title('64-QAM Constellation')
xlabel('In-phase')
ylabel('Quadrature')
hold on

s1=[xd(n) xq(n)];
s2=[xd(2*n) xq(2*n)];
s3=[xd(3*n) xq(3*n)];
s4=[xd(4*n) xq(4*n)];

text(s1(1),s1(2),'<- S1')
text(s2(1),s2(2),'<- S2')
text(s3(1),s3(2),'<- S3')
text(s4(1),s4(2),'<- S4')

line([0 0], [-8 8],'Color','black')
line([-8 8], [0 0],'Color','black')

%b)
figure
t=(0:1/n:k-1/n)'; % vetor tempo em bits
subplot(2,1,1)
stem(t,xd(1:length(t)))
title('In-phase modulating signal')
axis([0 4 min(xd)-1 max(xd)+1])

subplot(2,1,2)
stem(t(1:end),xq(1:length(t)))
title('Quadrature modulating signal')
axis([0 4 min(xq)-1 max(xq)+1])

%c)
t=(0:1/n:m-1/n)';
p1=sqrt(2/Tbit)*cos(2*pi*Rc*t);
p2=sqrt(2/Tbit)*sin(2*pi*Rc*t);

C=xd.*p1;
D=xq.*p2;
E=C+D;

figure
subplot(3,1,1)
plot(t,C(1:length(t)))
title('In-phase modulated carrier')
axis([0 4 min(C)-1 max(C)+1])

subplot(3,1,2)
plot(t,D(1:length(t)))
title('Quadrature modulated carrier')
axis([0 4 min(D)-1 max(D)+1])

subplot(3,1,3)
plot(t,E(1:length(t)))
hold on;
title('64-QAM output')
axis([0 4 min(E)-1 max(E)+1])

%e)
Amp=sqrt(s4(1).^2+s4(2).^2);
phase=atan2(s4(2),s4(1))*180/pi;

%f)
E1 = xd.*xd;
E2 = xq.*xq;

Em_pratico = (sum(E1(:))+sum(E2(:)))/(m*n);
Es = ((M-1)*d^2)/Tbit;

%g)
segmentL=m*n/10;
[Pxd,f] = pwelch(xd, hamming(segmentL), segmentL /2, segmentL, n/k,'onesided');
figure
Pxd=Pxd/2; 
fres = f(2)-f(1);
plot(f,10*log10(Pxd),'b');
title('Signal xd PSD');
ylabel('dBm/Hz');
grid on;
hold on
plot(f,10*log10(Es*(sin(pi.*f'*Tbit)./(pi.*f'*Tbit)).^2),'r');
legend('PSD modulated','PSD theoretical')

%h)
[PE,f] = pwelch(E, hamming(segmentL), segmentL /2, segmentL, n/k,'onesided');
figure
hold on
PE=PE/2; 
plot(f,10*log10(PE),'b');
title('Signal E PSD');
ylabel('dBm/Hz');
xlim([0.3 1.4])
grid on;
hold on
plot(f,10*log10(((Es/2)*(sin(pi.*(f-Rc)'*Tbit)./(pi.*(f-Rc)'*Tbit)).^2)+((Es/2)*(sin(pi.*(f+Rc)'*Tbit)./(pi.*(f+Rc)'*Tbit)).^2)),'r');
legend('PSD 64-QAM modulated','PSD 64-QAM theoretical')

%i)
KI4=2;
KQ4=2;
M4=KI4*KQ4;
k4=log2(M4);
d4=sqrt((Es*k4)/(M4-1));
xd4=d4/2*(mary(KI4,m,n)); % signal modulating the in-phase carrier
xq4=d4/2*(mary(KQ4,m,n)); % signal modulating the quadrature carrier
xd4 = reshape(xd4,[m*n,1]);
xq4 = reshape(xq4,[m*n,1]);
C4=xd4.*(sqrt(2/Tbit4)*cos(2*pi*Rc*t));
D4=xq4.*(sqrt(2/Tbit4)*sin(2*pi*Rc*t));
E4=C4+D4;
E1_4 = xd4.*xd4;
E2_4 = xq4.*xq4;
Em_pratico4 = (sum(E1_4(:))+sum(E2_4(:)))/(m*n);
Es4 = (M4-1)*d4^2/6;

[PE4,f] = pwelch(E4, hamming(segmentL), segmentL /2, segmentL, n/k,'onesided');
figure
hold on
PE4=PE4/2; 
plot(f,10*log10(PE4),'b');
title('Signal E PSD');
ylabel('dBm/Hz');
xlim([0.3 1.4])
grid on;
hold on
plot(f,10*log10(((Es4/2)*(sin(pi.*(f-Rc)'*Tbit)./(pi.*(f-Rc)'*Tbit)).^2)+((Es4/2)*(sin(pi.*(f+Rc)'*Tbit)./(pi.*(f+Rc)'*Tbit)).^2)),'r');
legend('PSD 4-QAM modulated','PSD 4-QAM theoretical')