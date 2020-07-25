%input parameters should be:
% m=          % number of symbols to simulate
% n=100       %number of samples per symbol
% KI=        %number of levels of the baseband signal modulating the in-phase carrier
% KQ=         %number of levels of the baseband signal modulating the quadrature carrier
% d=          %difference between consecutive levels in the baseband modulating signals
% Rc =5       %ration between carrier frequency and symbol rate.

%% 2.1 QAM Transmitter

%xd=mary(KI,m,n);     %represents the signal modulating the in-phase carrier (with distance 2 between levels)

%xq=mary(KQ,m,n);         %represents the signal modulating the quadrature carrier (with distance 2 between levels)

%% 2.2 (a)
close all
clear
clc

m=100000;       % number of symbols to simulate
n=100;      %number of samples per symbol
KI=8;       %number of levels of the baseband signal modulating the in-phase carrier
KQ=8;        %number of levels of the baseband signal modulating the quadrature carrier
d=2;        %difference between consecutive levels in the baseband modulating signals
Rc=5;       %ration between carrier frequency and symbol rate.
Tbit =6;
fc = 1/(Tbit)*Rc;

%plot da constelação
xa= -(KI-1):2:KI-1; 
xb=repmat(xa,KI,1);
marca=4;

plot(xa,xb, 'ob')
grid on

xd=mary(KI,m,n);    %sinal A    %represents the signal modulating the in-phase carrier (with distance 2 between levels)
xq=mary(KQ,m,n);    %sinal B     %represents the signal modulating the quadrature carrier (with distance 2 between levels)

figure
hold on
for i=1:m
       
    if(i<=marca)
        plot(xd(1,i),xq(1,i), '*r')
        text(xd(1,i),xq(1,i), ['\leftarrow symbol ', num2str(i)])
    end
    
    plot(xd(1,i),xq(1,i), 'ob')
    
    
end
grid on
hold off

%% 2.2 (b)
figure
subplot(2,1,1)
hold on
t=(0:marca*n-1)*Tbit;
for i=1:marca
   if(i==1)
     plot(t(1:n)./n, xd(1:n,i))
   else
     plot(t(n*(i-1):n*i)./n, [xd(1,i-1); xd(1:n,i)])  
   end
end
hold off
title('A')
legend('symbol 1', 'symbol 2', 'symbol 3', 'symbol 4');
subplot(2,1,2)
hold on

for i=1:marca
   if(i==1)
     plot(t(1:n)./n, xq(1:n,i))
    else
     plot(t(n*(i-1):n*i)./n, [xq(1,i-1); xq(1:n,i)])  
   end
end
hold off
title('B')
legend('symbol 1', 'symbol 2', 'symbol 3', 'symbol 4');

%% 2.2 (d)
t=(0:marca*n-1)./n;
fb=1/Tbit;
wc=2*pi*Rc;

for i=1:marca
    inicio = n*(i-1)+1;
    fim = n*i;
    
    
    C(inicio:fim)=xd(1:n,i).*sqrt(2*fb).*cos(wc.*(t(inicio:fim))');
    D(inicio:fim)=xq(1:n,i)*sqrt(2*fb).*sin(wc.*(t(inicio:fim))');

end

figure
subplot(2,1,1)
hold on
for i=1:marca
   if(i==1)
     plot(t(1:n)*Tbit, C(1:n))
   else
     plot(t(n*(i-1):n*i)*Tbit, C(n*(i-1):i*n))  
   end

end
hold off
title('C')
legend('symbol 1', 'symbol 2', 'symbol 3', 'symbol 4');
subplot(2,1,2)
hold on
for i=1:marca
   if(i==1)
     plot(t(1:n)*Tbit, D(1:n))
   else
     plot(t(n*(i-1):n*i)*Tbit, D(n*(i-1):i*n))  
   end
end
hold off
title('D')
legend('symbol 1', 'symbol 2', 'symbol 3', 'symbol 4');

%% 2.2(e)

E = C + D; 
figure

subplot(2,1,1)
hold on
for i=1:marca
   if(i==1)
     plot(t(1:n)*Tbit, E(1:n))
   else
     plot(t(n*(i-1):n*i)*Tbit, E(n*(i-1):i*n))  
   end

end
hold off
title('E = C + D')
legend('symbol 1', 'symbol 2', 'symbol 3', 'symbol 4');

subplot(2,1,2)
hold on
for i=1:marca
    Cc(i) = sqrt(2*fb)*xd(1,i);
    Ss(i) =-sqrt(2*fb)*xq(1,i);
    R(i) = sqrt((Cc(i))^2+(Ss(i))^2);
    fi(i) = atan(Ss(i)/Cc(i));
    
     if(xd(1,i)<0)
         fi(i)=fi(i)+pi;
     end
    
  if(i==1)
     plot(t(1:n)*Tbit, R(i).*cos(wc.*(t(inicio:fim))'+fi(i)))
   else
     plot(t(n*(i-1)+1:n*i)*Tbit, R(i).*cos(wc.*(t(inicio:fim))'+fi(i)))
   end
end
title('E = R cos(\theta+\phi)')
hold off
legend('symbol 1', 'symbol 2', 'symbol 3', 'symbol 4');

%% 2.2 (f)
M=2^6;


E1 = xd.*xd;
E2 = xq.*xq;
Em_pratico = (sum(E1(:))+sum(E2(:)))/(m*n);

Em_teorico = (M-1)*d^2/6;


%% 2.2 (g)
yq = xq(:,i)';
for i=2:m
    yq =[yq xq(:,i)'];
end

figure
segmentL=n*m/10
[Pxx,f] = pwelch(yq, hamming(segmentL), segmentL /2, segmentL ,n, 'onesided');
xlim([0 0.5])
FR = f(2)-f(1);
teorico = Em_teorico*(sin(pi.*f'*Tbit)./(pi.*f'*Tbit)).^2;
plot(f/Tbit,10*log10(Pxx))
hold on
plot(f,10*log10(teorico))
hold off
xlim([0 0.5])
legend('Pratical PSD', 'Theoretical PSD')

%% 2.2 (h)
t=(0:m*n-1)./n;
for i=1:m
    inicio = n*(i-1)+1;
    fim = n*i;
    
    
    Yc(inicio:fim)=xd(1:n,i).*sqrt(2*fb).*cos(wc.*(t(inicio:fim))');
    Yd(inicio:fim)=xq(1:n,i)*sqrt(2*fb).*sin(wc.*(t(inicio:fim))');
    Ye(inicio:fim)= Yc(inicio:fim)+Yd(inicio:fim);
end



figure
segmentL=n*m/10
[Pxx_E,f_E] = pwelch(Ye, hamming(segmentL), segmentL /2, segmentL ,n, 'onesided');
xlim([0.3 1.4])

teorico = Em_teorico/2*((sin(pi.*(f_E'-fc)*Tbit))./(pi.*(f_E'-fc)*Tbit)).^2 + Em_teorico/2*((sin(pi.*(f_E'+fc)*Tbit))./(pi.*(f_E'+fc)*Tbit)).^2;
plot(f_E/Tbit,10*log10(Pxx_E))
hold on
plot(f_E,10*log10(teorico))
hold off
xlim([0.3 1.4])
legend('Pratical PSD', 'Theoretical PSD')


%% 2.2(i)
M_4Q=2^2;
d_4Q= sqrt(Em_teorico*6/(M_4Q-1));
KI_4Q=2;       %number of levels of the baseband signal modulating the in-phase carrier
KQ_4Q=2;        %number of levels of the baseband signal modulating the quadrature carrier
Tbit_4Q =2;
fc_4Q = 1/(Tbit_4Q)*Rc;
fb_4Q =1/Tbit_4Q;
xd_4Q=mary(KI_4Q,m,n).*d_4Q/2;    %sinal A    %represents the signal modulating the in-phase carrier (with distance 2 between levels)
xq_4Q=mary(KQ_4Q,m,n).*d_4Q/2;    %sinal B     %represents the signal modulating the quadrature carrier (with distance 2 between levels)

t_4Q=(0:m*n-1)./n;
for i=1:m
    inicio = n*(i-1)+1;
    fim = n*i;
    
    
    Yc_4Q(inicio:fim)=xd_4Q(1:n,i).*sqrt(2*fb_4Q).*cos(wc.*(t_4Q(inicio:fim))');
    Yd_4Q(inicio:fim)=xq_4Q(1:n,i)*sqrt(2*fb_4Q).*sin(wc.*(t_4Q(inicio:fim))');
    Ye_4Q(inicio:fim)= Yc_4Q(inicio:fim)+Yd_4Q(inicio:fim);
end

figure
segmentL=n*m/10
[Pxx_E_4Q,f_E_4Q] = pwelch(Ye_4Q, hamming(segmentL), segmentL /2, segmentL ,n, 'onesided');
xlim([0.3 1.4])

teorico_4Q = Em_teorico/2*((sin(pi.*(f_E_4Q'-fc_4Q)*Tbit_4Q))./(pi.*(f_E_4Q'-fc_4Q)*Tbit_4Q)).^2 + Em_teorico/2*((sin(pi.*(f_E_4Q'+fc_4Q)*Tbit_4Q))./(pi.*(f_E_4Q'+fc_4Q)*Tbit_4Q)).^2;
plot(f_E_4Q/Tbit_4Q,10*log10(Pxx_E))
hold on
plot(f_E_4Q,10*log10(teorico_4Q))
hold off
xlim([1 4])
legend('Pratical PSD', 'Theoretical PSD')
