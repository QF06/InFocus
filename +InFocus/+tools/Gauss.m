function [Dat] = Gauss(Sc, Wid, ofst)
%Gauss Summary of this function goes here
%   Sc: coordimate info
%   Wid: waist radius for 1/e (amplitute); 1/e^2(intensity)
%   Wid(1) for x axis, Wid(2) for y axis
lw_x = Wid(1); lw_y = Wid(2);
dx = ofst(1); dy = ofst(2);
Dat = exp(-(Sc.x-dx).^2/lw_x^2-(Sc.y-dy).^2/lw_y^2);
end

