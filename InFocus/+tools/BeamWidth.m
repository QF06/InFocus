function [width, tlead, ttrail] = BeamWidth(x,y)
%%  FWHM calculate the full width of the two independent variables whose corresonding 
%  dependent variables equal to the half of the maximum dependent variable value.
%%  Form
% width = BeamWidth(x,y)
%%  Input
% x(1,:); Array of one raw
% y(1,:); Array of one raw
%% Description
% Full-Width at Half-Maximum (FWHM) of the waveform y(x) and its polarity.
% The FWHM result in 'width' will be in units of 'x'
%% Copyright
% Rev 1.0, Jan. 2020 (Qingfeng Li)

%%Initialize
if(nargin < 2)
  error('Input: error','Not enough inputs');
elseif(nargin > 2)
  error('Input: error','Too many inputs');
end

N     = length(y);
levE2 = (max(y) +min(y))/exp(1)^2;
if y(1) < levE2	% find index of center (max or min) of pulse
  [~,centerIndex] = max(y);
  Pol             = +1;
  %disp('Pulse Polarity = Positive')
  y               = y / max(y);
  levE2           = max(y)/exp(1)^2;
else
  [~,centerIndex] = min(y);
  Pol = -1; 
  %disp('Pulse Polarity = Negative')
  y               = y / min(y);
  levE2           = max(y)/exp(1)^2;
end
k   = 2;
while sign(y(k)-levE2) == sign(y(k-1)-levE2)
    k = k+1;
    if(k == N)
      error('This beam size is out of bundary.');
    end
end	%first crossing is between y(k-1) & y(k) 
slope  = (y(k)-y(k-1)) / (x(k)-x(k-1));
tlead  = x(k-1) + (levE2-y(k-1)) / slope;
k      = centerIndex+1;	%start search for next crossing at center
while ((sign(y(k)-levE2) == sign(y(k-1)-levE2)) && (k <= N-1))
	k = k+1;
end
if k ~= N
  Ptype  = 1;  
  %disp('Pulse is Impulse or Rectangular with 2 edges')
  slope  = (y(k)-y(k-1)) / (x(k)-x(k-1));
  ttrail = x(k-1) + (levE2-y(k-1)) / slope;
  width  = ttrail - tlead;
else
  Ptype  = 2; 
  %disp('Step-Like Pulse, no second edge')
  ttrail = NaN;
  width  = NaN;
end


