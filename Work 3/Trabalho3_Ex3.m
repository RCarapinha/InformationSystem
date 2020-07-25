clear
close all

m=10000;       % number of symbols to simulate  => 10/Pe
n=100;      %number of samples per symbol
KI=4;       %number of levels of the baseband signal modulating the in-phase carrier
KQ=4;        %number of levels of the baseband signal modulating the quadrature carrier
n0=17.8;        %dBm/Hz
Tbit =2;
Rc=5;
M=2^KI;
fc = 1/(Tbit)*Rc;
marca =4;
%% 3.2(a)
Pes = 0.5;

d= erfcinv((1-sqrt(1-Pes))/(1-1/sqrt(M)))*2*sqrt(n0);



%% 3.2(b)
% A_x=mary(KI,m,n).*d/2;
% A_y=mary(KI,m,n).*d/2;
m=16;
n=4;
[A_x, A_y] = Mapping_16_QAM(m,d);
       

noise_x = wgn(1,m,n0/2,'dBm');  %specifies the units of power as 'dBW', 'dBm', or 'linear' in addition to the input arguments in any of the previous syntaxes.
noise_y = wgn(1,m,n0/2,'dBm');  %specifies the units of power as 'dBW', 'dBm', or 'linear' in addition to the input arguments in any of the previous syntaxes.

for i=1:m
    B_x(i) = A_x(i)+noise_x(i);
    B_y(i) = A_y(i)+noise_y(i);
end

for i=1:m
    [~,indice]=min((A_x(1:16)-B_x(1,i)).^2+(A_y(1:16)-B_y(1,i)).^2);
    C_x(i)=A_x(indice);
    C_y(i)=A_y(indice);
end
figure
hold on
for i=1:m
       
    if(i<=marca)
        plot(A_x(i),A_y(i), '.r')
        text(A_x(i),A_y(i), ['\leftarrow A', num2str(i)])
        
        plot(B_x(i),B_y(i), '*g')
        text(B_x(i),B_y(i), ['\leftarrow B', num2str(i)])
        
        plot(C_x(i),C_y(i), '+c')
        text(C_x(i),C_y(i), ['C', num2str(i), '\rightarrow' ], 'HorizontalAlignment', 'right')
        
    end
    
    plot(A_x(i),A_y(i), 'ob')
    
    
end
grid on
hold off

%% 3.3(a)
Es=4;
Pes_ini=10^-3;
Pes_fin=0.5;
INI=3*Es/(( erfcinv((1-sqrt(1-Pes_ini))/(1-1/sqrt(M))))^2*2*(M-1));
FIN=3*Es/(( erfcinv((1-sqrt(1-Pes_fin))/(1-1/sqrt(M))))^2*2*(M-1));

%% 3.3(b)
fprintf('\nTabela I\n');
var = FIN-INI;
Eb = Es/4;
for i=0:4
    n_16 = INI+(var*i)/4;
    Pes = 1-(1-(1-1/sqrt(M))*erfc(sqrt(3*log2(M)*Eb/(2*(M-1)*n_16))))^2;
    PeI = 1-(1-erfc(sqrt(6*Eb/(n_16*(M-1)))))^2;
    ratio = PeI/Pes;
    fprintf('PeI = %.4f, Pes = %.4f, ratio = %.4f para n0=%.4f dBm\n',PeI,Pes,ratio,10*log10(1000*n_16))
end



%% 3.4 (16-QAM) & 3.3 (c)(e)
fprintf('\n16-QAM\n');
M=2^4;
d= sqrt(6*Es/(M-1));

for i=0:4
    n_16 = INI+(var*i)/4;
    Pes = 1-(1-(1-1/sqrt(M))*erfc(sqrt(3*log2(M)*Eb/(2*(M-1)*n_16))))^2;
    mm(i+1)=round(10/Pes);
end
m=max(mm);
[A_x_16, A_y_16, A_bin_16] = Mapping_16_QAM(m,d);

erro_SER = zeros(5,5);
erro_BER = zeros(5,5);
for l=0:4
    
    n_16 = INI+(var*l)/4;
    Pes = 1-(1-(1-1/sqrt(M))*erfc(sqrt(3*log2(M)*Eb/(2*(M-1)*n_16))))^2;

    for j=1:5
        noise_x_16 = wgn(1,m,n_16/2,'linear');  %specifies the units of power as 'dBW', 'dBm', or 'linear' in addition to the input arguments in any of the previous syntaxes.
        noise_y_16 = wgn(1,m,n_16/2,'linear');  %specifies the units of power as 'dBW', 'dBm', or 'linear' in addition to the input arguments in any of the previous syntaxes.

        for i=1:m
           B_x_16(i) = A_x_16(i)+noise_x_16(i);
           B_y_16(i) = A_y_16(i)+noise_y_16(i);
        end
        
        for i=1:m
            [~,indice(i)]=min((A_x_16(1:16)-B_x_16(i)).^2+(A_y_16(1:16)-B_y_16(i)).^2);
             C_x_16(i)=A_x_16(indice(i));
             C_y_16(i)=A_y_16(indice(i));
             C_bin_16(i,:)=A_bin_16(indice(i),:);
             
             if (A_x_16(i)~=C_x_16(i) || A_y_16(i)~=C_y_16(i))
                erro_SER(l+1,j)=erro_SER(l+1,j)+1; 
             end
             erro_BER(l+1,j)=erro_BER(l+1,j)+sum(C_bin_16(i,:)~=A_bin_16(i,:));

        end
    
    end
    
    erro_medio_SER = (sum(erro_SER(l+1,:))/j)/m;
    variacao=abs(erro_medio_SER-Pes)/Pes*100;
    fprintf('3.3 - Numero simbolos = %.0f, Erro médio = %.4f, Erro teórico = %.4f, var= %.4f para n0=%.4f dBm\n',m, erro_medio_SER,Pes,variacao,10*log10(1000*n_16))


    erro_medio_BER = (sum(erro_BER(l+1,:))/j)/(m*4);
    BER_SER=erro_medio_BER/erro_medio_SER;
    fprintf('3.4 - Numero simbolos = %.0f, BER = %.4f, BER/SER = %.4f  para n0=%.4f dBm\n',m, erro_medio_BER,BER_SER,10*log10(1000*n_16))

    figure
    hold on
    for a=1:m
         plot(A_x_16(a),A_y_16(a), 'ob')
         plot(B_x_16(a),B_y_16(a), '*g')
         plot(C_x_16(a),C_y_16(a), '+c')
    end
    hold off
    title(['n0= ',num2str(10*log10(1000*n_16)), 'dBm/Hz'])
end
   fprintf('\n');
%    figure
%    hold on
%    for i=1:16
%        plot(A_x_16(i),A_y_16(i),'ob')
%        text(A_x_16(i),A_y_16(i), ['\leftarrow ', num2str(A_bin_16(i,:))])
%    end
%    hold off

%% 3.4 (64-QAM)  & 3.3 (f)
fprintf('64-QAM\n');
M=2^6;
d= sqrt(6*Es/(M-1));
Eb = Es/6;
for i=0:4
    n_64 = INI+(var*i)/4;
    Pes = 1-(1-(1-1/sqrt(M))*erfc(sqrt(3*log2(M)*Eb/(2*(M-1)*n_64))))^2;
    mm(i+1)=round(10/Pes);
end

m=max(mm);
[A_x_64, A_y_64, A_bin_64] = Mapping_64_QAM(m,d);

erro_SER = zeros(5,5);
erro_BER = zeros(5,5);
for l=0:4
    
    n_64 = INI+(var*l)/4;
    Pes = 1-(1-(1-1/sqrt(M))*erfc(sqrt(3*log2(M)*Eb/(2*(M-1)*n_64))))^2;
    for j=1:5
        noise_x_64 = wgn(1,m,n_64/2,'linear');  %specifies the units of power as 'dBW', 'dBm', or 'linear' in addition to the input arguments in any of the previous syntaxes.
        noise_y_64 = wgn(1,m,n_64/2,'linear');  %specifies the units of power as 'dBW', 'dBm', or 'linear' in addition to the input arguments in any of the previous syntaxes.

        for i=1:m
           B_x_64(i) = A_x_64(i)+noise_x_64(i);
           B_y_64(i) = A_y_64(i)+noise_y_64(i);
        end
        
        for i=1:m
            [~,indice(i)]=min((A_x_64(1:64)-B_x_64(i)).^2+(A_y_64(1:64)-B_y_64(i)).^2);
             C_x_64(i)=A_x_64(indice(i));
             C_y_64(i)=A_y_64(indice(i));
             C_bin_64(i,:)=A_bin_64(indice(i),:);
             
             if (A_x_64(i)~=C_x_64(i) || A_y_64(i)~=C_y_64(i))
                erro_SER(l+1,j)=erro_SER(l+1,j)+1; 
             end
             erro_BER(l+1,j)=erro_BER(l+1,j)+sum(C_bin_64(i,:)~=A_bin_64(i,:));

        end
    
    end
    
    erro_medio_SER = (sum(erro_SER(l+1,:))/j)/(m);
    variacao=abs(erro_medio_SER-Pes)/Pes*100;
    fprintf('3.3 - Numero simbolos = %.0f, Erro médio = %.4f, Erro teórico = %.4f, var= %.4f para n0=%.4f dBm\n',m, erro_medio_SER,Pes,variacao,10*log10(1000*n_64))


    erro_medio_BER = (sum(erro_BER(l+1,:))/j)/(m*6);
    BER_SER=erro_medio_BER/erro_medio_SER;
    fprintf('3.4 - Numero simbolos = %.0f, BER = %.4f, BER/SER = %.4f  para n0=%.4f dBm\n',m, erro_medio_BER,BER_SER,10*log10(1000*n_64))


end