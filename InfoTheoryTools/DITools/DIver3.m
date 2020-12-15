function [dixy, diyx] = DIver3(Xorig,Yorig,mem,K, down_sample)
%DIver3     Directed Information.
%
%  DIver3 takes as input multiple samples of vectors X and Y, memory, 
%  kNN size and downsampling size, to be able to estimate higher memory
%  lengths
%
%  DI code can be run in parallel by using inside a parfor loop
%
%   Inputs,     X,     timeseries X, an n x m matrix, n = number of independent samples,
%                      m = length of the time series
%               Y,     timeseries Y, an n x m matrix, n = number of independent samples,
%                      m = length of the time series
%               mem,   Number of memory size samples of timeseries
%               K,     k neareast neighbor size, typically 3-4 provides
%                       good estimates when number samples are less than 1000, 
%                       can be increased when number of samples is
%                       in the order of 10000. (details on how to choose 
%                       a good k is in "Estimating Mutual Information" 
%                       by Kroskov. et al, 2004.)
%               down_sample,     the downsampled version of the time series 
%                                  is used, inorder to be able to use
%                                  higher memory orders.
%   Outputs:    dixy,     DI from X to Y in bits
%               diyx,     DI from Y to X in bits
%
%   Note: DIver3 differs from DIver2 in that multiple timeshifted 
%   downsampled versions of X and Y are used in estimations, to improve 
%   the average. Data is nt thrown away, like in the DIver2. 
%
%   Sudha Yellapantula, Rice University


nlen = size(Xorig,2);% length of time series

%initialization
DIX_Y = 0;
DIY_X = 0;

%count to accrue number of DI estimates  to report final average DI rate
count =1;
for shift_index = 1:down_sample % going through each data between between the downsamples
    X = [NaN(size(Xorig,1),shift_index-1) Xorig];
    Y = [NaN(size(Xorig,1),shift_index-1) Yorig]; 
    for i = mem+1 : floor(nlen/down_sample) %downsampling from different starting indices
        X1 = X(:,(i-mem+1:i)*down_sample);
        X2 = Y(:,i*down_sample);
        X3 = Y(:,(i-mem:i-1)*down_sample);
        
        DIX_Y = DIX_Y+ ConditionalMI(X1,X2,X3,K);
        
        
        Y1 = Y(:,(i-mem+1:i)*down_sample);
        Y2 = X(:,i*down_sample);
        Y3 = X(:,(i-mem:i-1)*down_sample);
        
        DIY_X = DIY_X + ConditionalMI(Y1,Y2,Y3,K);
        count = count+1;
    end
end

%average DI across the time series window to obtain DI rate
dixy = DIX_Y/count;
diyx = DIY_X/count;
    
end
