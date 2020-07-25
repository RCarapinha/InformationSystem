function [z,atraso]=correlacao(inf,tx,Amost,Nmax,Mmax)
% function [z,atraso]=correlacao(inf,tx,Amost,Nmax,Mmax)
% where inf is the information sequence (an element / symbol),
% tx is the sampled sequence, Amost is the no. of samples per symbol,
% Nmax is the no. of information symbol to consider in the correlation
% and Mmax is the maximum number of symbols of evauated delay
% z is the correlation function between inf and tx in samples
% and atraso is the delay between inf and tx, also in samples
a = zeros(1,Amost*Mmax); % Preallocate matrix
for j = 1:Amost*Mmax % For each delay T=(j-1)*Ta
    for i = 1:Nmax % For each information sample
        k = j+(i-1)*Amost; % Position in the sampled signal
        a(j) = a(j)+inf(i)*tx(k);
    end
end
z=a./Nmax;

% Determination of the maximum position
max =0;
pos = 0;
for k = 1:length(a)
    if a(k) > max
        max=a(k);
        pos = k;
    end
end
atraso = pos-1;