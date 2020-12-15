function [condMI] = ConditionalMI(x,y,z,K)
% ConditionalMI     Conditional Mutual Information.
%
%  ConditionalMI takes as input multiple samples of vectors x, y and z and
%  kNN parameter K
%
%
%   Inputs,     x,     timeseries x, an n x m1 matrix, n = number of independent samples,
%                      m1 = length of the time series
%               y,     timeseries y, an n x m2 matrix, n = number of independent samples,
%                      m2 = length of the time series
%               z,     timeseries z, an n x m3 matrix, n = number of independent samples,
%                      m3 = length of the time series
%               K,     k neareast neighbor size, typically 3-4 provides
%                       good estimates when number samples are less than 1000, 
%                       can be increased when number of samples is
%                       in the order of 10000. (details on how to choose 
%                       a good k is in "Estimating Mutual Information" 
%                       by Kroskov. et al, 2004.)
%
%   Outputs:    condMI,     Conditional Mutual Information between 
%                           x and y, conditioned on z. Joint entropies are
%                           estimated and CondMI is estimated as sums and
%                           differences of joint entropies.
%
%
%   Sudha Yellapantula, Rice University

% ensure data is in the correct format
    if (size(x,1) < size(x,2))
        x = x.';
    end
   
    if (size(y,1) < size(y,2))
        y = y.';
    end

    if (size(z,1) < size(z,2))
        z = z.';
    end

    data_HXZ = [x z];
    data_HYZ = [y z];
    data_HZ = [z];
    data_HXYZ = [x y z];
    
    dim_xz = size(data_HXZ,2);
    dim_yz = size(data_HYZ,2);
    dim_z = size(data_HZ,2);
    dim_xyz = size(data_HXYZ,2);
    
    NumPoints = size(y,1);
    
    dim_z = size(z,2);

    atria_xyz = nn_prepare(data_HXYZ,'maximum');
    [indexMat, distance_xyz] = nn_search(data_HXYZ,atria_xyz,1:NumPoints,K,0);

    H_XYZ = -psi(K) + psi(NumPoints) + (dim_xyz)*mean(log(2*distance_xyz(:,K)));

    %entropy of x z
    atria_xz = nn_prepare(data_HXZ,'maximum');
    ncnt_xz = range_search(data_HXZ,atria_xz,1:NumPoints,distance_xyz(1:NumPoints,K)-eps,0);
    H_XZ = -mean(psi(ncnt_xz+1)) + psi(NumPoints) + (dim_xz)*mean(log(2*distance_xyz(:,K)));

    %entropy of y z
    atria_yz = nn_prepare(data_HYZ,'maximum');
    ncnt_yz = range_search(data_HYZ,atria_yz,1:NumPoints,distance_xyz(1:NumPoints,K)-eps,0);
    H_YZ = -mean(psi(ncnt_yz+1)) + psi(NumPoints) + (dim_yz)*mean(log(2*distance_xyz(:,K)));

    %entropy of z
    atria_z = nn_prepare(data_HZ,'maximum');
    ncnt_z = range_search(data_HZ,atria_z,1:NumPoints,distance_xyz(1:NumPoints,K)-eps,0);
    H_Z = -mean(psi(ncnt_z+1)) + psi(NumPoints) + (dim_z)*mean(log(2*distance_xyz(:,K)));


    %MI KSG 1
    condMI = H_XZ +H_YZ -H_Z -H_XYZ;
    %MI = psi(K)  - mean(psi(ncnt_x+1) + psi(ncnt_y+1)) + psi(NumPoints); 
     

end