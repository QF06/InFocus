function [ Sc,Ph ] = Vortex(Sc,m,max)
%% vortex, this function add a vortex phase distribution on the input field
% num: topological charge of the vortex phase plate (m)
Tmp   = Sc.x+1i*Sc.y;
Azimuth = angle(Tmp);
Ph     = zeros(Sc.pix_w,Sc.pix_h);
full   = 2*pi/m;
for ki=1:Sc.pix_w
    for kj=1:Sc.pix_h
        Ph(ki,kj) = mod(Azimuth(ki,kj),full)/full*max;
    end
end

Sc.E_x = Sc.E_x.*exp(1i*Ph);
Sc.E_y = Sc.E_y.*exp(1i*Ph);
end