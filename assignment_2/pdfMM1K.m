function [p,b] = pdfMM1K(rho,C)
%
% This function returns the theoretical pdf of the number of tokens in
% an MM1 system with rho = interServiceTime / interArrivalTime = arrivalRate / serviceRate;
% It calculates the pdf for buffer sizes up-to and including C, if tail = 0.
% If tail = 1, it returns for C the probability of finding a buffersize
% larger-than or equal to C
%
b = [0:C];
n = [1:C];
% p = (1-rho)*rho.^b;
p0 = (1-rho)/(1-rho.^((C+1)+1));
p = p0*rho.^n;
p = [p0 p];
end