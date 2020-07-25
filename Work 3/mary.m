function [y] = mary(levels,m,n)
%UNTITLED3 Summary of this function goes here
%   generates a sequence of m amplitudes drawn randomly

a=-(levels-1):2:levels-1;
x=randsrc(1,m,a);
y=repmat(x,n,1);

end

