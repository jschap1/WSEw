% Exploring relationships between hydraulic parameters under different
% channel shapes
%
% Maybe use the channel shape proposed by Neal et al. 2015
%
% Change in depth, WP introduced by using a trapezoidal channel instead of rectangular,
% given the same cross sectional area and top width:
%

%%

% Wide, rectangular channel
yr = 10;
w = 200;
A = yr*w;
WPr = 2*yr+w;

% Equivalent trapezoidal channel
s = 1; % choose a side slope

yt1 = (w*s+sqrt(w*s*(w*s-4*yr)))/2; % positive sqrt
yt2 = (w*s-sqrt(w*s*(w*s-4*yr)))/2; % negative sqrt

% At = wt*yt-yt^2/s % check (the area is always equal, regardless of the choice of the sign of the square root.)
% However, the wetted perimeter may be negative, if the sign is chosen
% incorrectly. Require WPt>0

Ht1 = yt1*sqrt(1+1/s^2); % calculate hypotenuse
Ht2 = yt2*sqrt(1+1/s^2);

yt1<(w+2*Ht1)*s/2
yt2<(w+2*Ht2)*s/2

% yt = yt(yt>0 & yt<(w+2*Ht)*s/2);

WPt = w - 2*yt1/s + 2*Ht1;
WPt = w - 2*yt2/s + 2*Ht2;

At1 = w*yt1-yt1^2/s;
At2 = w*yt2-yt2^2/s;


