function [dixy, diyx] = DI(X,Y,mem,K)
%DI      Directed Information.
%
%  DI takes as input multiple samples of vectors X and Y, memory, kNN size.
%
%  DI code can be run in parallel by using inside a parfor loop
%
%   Inputs,     X,     timeseries X, an n x m matrix, n = number of independent samples,
%                      m = length of the time series
%               Y,     timeseries Y, an n x m matrix, n = number of independent samples,
%                      m = length of the time series
%               mem,   Number of memory size samples of timeseries
%               K,     k neareast neighbor size, typically 3-4 provides
%               good estimates, can be increased when number of samples is
%               in the order of 10000. Ref to paper "Estimating Mutual
%               INformation" by Kroskov. et al, 2004.
%   Outputs:    dixy,     DI from X to Y in bits
%               diyx,     DI from Y to X in bits
%
%   Sudha Yellapantula, Rice University

nlen = size(X,2); % length of time series

%initialization
DIX_Y = 0; 
DIY_X = 0;

%count to accrue number of DI estimates for average DI estimate
count =1;

for i = mem+1 : nlen
    X1 = X(:,i-mem+1:i);
    X2 = Y(:,i);
    X3 = Y(:,i-mem:i-1);
    
    DIX_Y = DIX_Y+ ConditionalMI(X1,X2,X3,K);
    
    Y1 = Y(:,i-mem+1:i);
    Y2 = X(:,i);
    Y3 = X(:,i-mem:i-1);
    
    DIY_X = DIY_X + ConditionalMI(Y1,Y2,Y3,K);
    count = count+1;
end

%average DI across the time series window under consideration to obtain DI
%rate
dixy = DIX_Y/count;
diyx = DIY_X/count;
    
end