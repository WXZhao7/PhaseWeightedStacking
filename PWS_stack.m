function STACK=PWS_stack(DATA,dt)
 
% FUNCTION PWS_stack=PWS_stack(data)
% Phase weighted stacking following 
% Schimmel M., Stutzmann E., Gallart J.,2011. Using instantaneous phase coherence for signal 
% extraction from ambient noise data at a local to a global scale,Geophys. J. Int.,184(1),
% 494–506.10.1111/j.1365-246X.2010.04861.x
% DATA is the input data arranged such that columns are time samples 
% and rows are ecents. dt is the time sample rate.
% The stacked data is returned as a vector STACK. 

% load data
[a1 a2] = size(DATA);

% number of events to use
Nv = a1;

% time vector
tvec = [0:(a2-1)]*dt;
% power of stacking as in Schimmel et al. (2010; GJI)
pwr = 2;

% S-transform of linear stack
[stranl,fvec] = S_transform_FD_fullspec(mean(DATA(1:Nv,:)),dt);

% make phase weight
sumr = zeros(a2,a2);
for ii=1:Nv
    % S-transform of ii-th seismogram
    [stran,fvec] = S_transform_FD_fullspec(DATA(ii,:),dt);
    % Equation 6 in Schimmel et al. (2010; GJI)
    
    sss = abs(stran);
    %%%Deal with possible zero amplitude
    sss(sss==0) = 1e32;
    sss(sss==0) = min(min(sss));
    bbb = (stran./sss).*exp(i*2*pi*(fvec'*tvec));
    %%bbb(isnan(bbb)) = 0;
    sumr = sumr + bbb;
    %sumr = sumr + (stran./sss).*exp(i*2*pi*(fvec'*tvec));
    %find(isnan(sumr))
    ii
end
% Equation 6 in Schimmel et al. (2010; GJI)
sumr = abs(sumr/Nv).^pwr;

% apply phase weight to linear stack to get phase-weighted-stack
% Equation 7 in Schimmel et al. (2010; GJI)
stranpws = sumr.*stranl;

% transform back to the time domain
STACK = S_transform_inverse_fullspec(stranpws,fvec);