%% Begin
clear all
close all
clc

m=1000;       % number of symbols to simulate
n=100;      %number of samples per symbol
KI=4;       %number of levels of the baseband signal modulating the in-phase carrier
KQ=4;        %number of levels of the baseband signal modulating the quadrature carrier
Rc=5;       %ration between carrier frequency and symbol rate.
M=KI*KQ;
k=log2(M);

%% 3.1
mapping_16_QAM=[3   3  0 0 0 0;
                3   1  0 0 0 1;
                3  -3  0 0 1 0;
                3  -1  0 0 1 1;
                1   3  0 1 0 0;
                1   1  0 1 0 1;
                1  -3  0 1 1 0;
                1  -1  0 1 1 1;
                -3   3  1 0 0 0;
                -3   1  1 0 0 1;
                -3  -3  1 0 1 0;
                -3  -1  1 0 1 1;
                -1   3  1 1 0 0;
                -1   1  1 1 0 1
                -1  -3  1 1 1 0
                -1  -1  1 1 1 1];
            
source=randsrc(k,m,[0,1])';
A = zeros(m,2);
for i=1:m
    for j=1:M
        if source(i,1:k)==mapping_16_QAM(j,3:k+2)
            A(i,1)=mapping_16_QAM(j,1); % signal modulating the in-phase carrier
            A(i,2)=mapping_16_QAM(j,2); % signal modulating the quadrature carrier
        end
    end
end

%% 3.2
noisePSD=17.8;
PES = 0.5;
aux = (1-sqrt(1-PES))/(1-(1/sqrt(M)));
Eb = ((erfcinv(aux)*sqrt(2*(M-1)*noisePSD))^2)/(3*log2(M));
Es = Eb*k;
d = sqrt((Es*6)/(M-1));

const=mapping_16_QAM*d/2;
A=A*d/2;
noise=wgn(m,2,noisePSD-3,'dBm'); % noise samples with average power noisePower
%%noisePSD-3 para unilateral, noisePSD para bilateral
B=A+noise;

dist=ones(m,M);
for i=1:m
    for j=1:M
        dist(i,j)=B(i,1)*const(j,1)+B(i,2)*const(j,2)-(const(j,1).^2+const(j,2).^2)/2;
    end
    [max_val,max_index]=max(dist(i,:));
    C(i,1)=const(max_index,1);
    C(i,2)=const(max_index,2);
end

%% 3.3
figure
plot(A(:,1),A(:,2),'ro')
plot(B(:,1),B(:,2),'bo')
title('16-QAM Constellation')
xlabel('In-Phase Component')
ylabel('Quadrature Component')
hold on
grid on

A1=[A(1,1) A(1,2)];
A2=[A(2,1) A(2,2)];
A3=[A(3,1) A(3,2)];
A4=[A(4,1) A(4,2)];

grid on
plot(A1(1),A1(2),'g*','MarkerSize',12)
plot(A2(1),A2(2),'g*','MarkerSize',12)
plot(A3(1),A3(2),'g*','MarkerSize',12)
plot(A4(1),A4(2),'g*','MarkerSize',12)

text(A1(1),A1(2),'<- A1')
text(A2(1),A2(2),'<- A2')
text(A3(1),A3(2),'<- A3')
text(A4(1),A4(2),'<- A4')

figure
hold on
grid on
plot(C(:,1),C(:,2),'bo')

%% 3.3
Eb = 4; %CORRIGIR
PES_Max = 0.5;
PES_Min = 10e-3;
INI = (((3*log2(M)*Eb)/(erfcinv((1-sqrt(1-PES_Min))/(1-(1/(sqrt(M)))))))^2)/(2*(M-1));
FIN = (((3*log2(M)*Eb)/(erfcinv((1-sqrt(1-PES_Max))/(1-(1/(sqrt(M)))))))^2)/(2*(M-1));