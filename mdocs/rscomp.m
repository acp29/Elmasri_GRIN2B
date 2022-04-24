function y = rscomp (y, sampInt, R_s_final, R_s, C_m, V_hold, V_rev)

% Function for offline series resistance compensation
%    
% This script is my adaptation of Erwin Neher's Igor function:
%  SeriesresistanceComp
%  This were freely available in in Proc02_Apr4Web.ipf from Erwin Neher's webpage:
%    http://www3.mpibpc.mpg.de/groups/neher/index.php?page=software
%    (last accessed: 01 July 2014)
%
% Instead of setting the fraction compensation, this function asks for the desired
% uncompensated series resistance. Whole-cell properties are calcuated using the wcp
% function.
%
% rsomp replaces current traces by their series-compensated version;
% the value at i is replaced by the average at i and i+1
% R_s is in ohms, C_m in Farads, fraction is the fraction to be compensated
% if R_s was 5 MOhm in the experiment and if it was 50% hardware compensated
% then R_s = 2.5e6 has to be entered and f=1 for a complete overall compensation
% The routine, similarly to that published by Traynelis J. Neurosc. Meth. 86:25,
% compensates the frequency response, assuming a single R_s*C_m - time constant
% (at constant V-hold) and a three component equivalent circuit for the pipette cell
% assembly with  R_s, C_m, R_m
% 
% Theory: V_h = R_s*I+U_m;
%  I_r is membrane resistive current,
%  I is total current
%  I_r = I-I_c = I-C_m*dU_m/dt = I+C_m*R_s*dI/dt (because R_s*I+U_m = const.)
%  G_m=I_r/ U_m = (I+C_m*R_s*dI/dt)/ (V_h-R_s*I)
% For complete correction (fraction = 1) : I_corr = V_h*(I+C_m*R_s*dI/dt)/ (V_h-R_s*I)
 
numpoints = size(y,1); 
numtraces = size(y,2);
if numtraces > 1
 error('This function only supports one trace/wave')
end

% Units:
%  y should be in A
%  SampInt should be in s
%  R_s should be in ohms
%  R_s_final should be in ohms
%  Cm should be in Farads
tau = R_s * C_m;
if nargin < 7
  V_rev = 0;
end
voltage = V_hold - V_rev;
fraction = max(1 - (R_s_final / R_s), 0);
tau_corr = R_s * C_m * (1 - fraction);

% First point: (we have to calculate this separately, because we need the value at i-1 below)
denominator = voltage - R_s * fraction * y(1);
if denominator ~= 0
  y(1) = y(1) * (voltage / denominator);
end

for i = 2:numpoints-1
  % this is the loop doing the correction for all other points
  % first calculate R_m for zero series resistance under the assumptions
  % that  U_m + U_Rs = const = voltage
  current = (y(i+1) + y(i)) / 2;  % The in between(mean) value
  derivative =  (y(i+1) - y(i)) / sampInt;
  denominator = current + tau * derivative;
  if denominator ~= 0
   R_m = (voltage - R_s * current) / denominator;  % calculate the true R_m
  else
    % Do nothing. Leave R_m as is
  end
  % Now calculate current for new series resitance
  denominator = (R_m + (1 - fraction) * R_s) * (1 + tau_corr / sampInt);
  if denominator ~= 0
   y(i) = tau_corr / (tau_corr + sampInt) * y(i-1) + voltage/denominator;
  else
   y(i) = y(i-1);  % old value
  end
end

