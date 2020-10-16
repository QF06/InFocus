function [Sc,Ph] = Bessel(Sc,Lsr,conic)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[X,Y] = meshgrid(1:Sc.pix_w,1:Sc.pix_h);
x = Sc.wid/Sc.pix_w .*(X-(Sc.pix_w-1)/2);
y = Sc.wid/Sc.pix_w .*(Y-(Sc.pix_h-1)/2);
r = sqrt(x.^2+y.^2);
r(r>12.5e-3)=0;
phi = atan2(y,x);
%Ph = mod(20*pi*r/5.2e-3 ,2*pi);
%Ph = mod(phi + 2*pi*r/5.2e-3 + 2*pi*r.*cos(phi)*4.5e3,2*pi);
%Ph = -Lsr.w0/Lsr.c * 0.1* sind(conic) .* r;
Ph = mod(Lsr.k0 * (1.51-1)*sind(conic) * r + Lsr.k0*sind(conic),2*pi);


% Ph1   = zeros(Sc.pix_w,Sc.pix_h);
% for ki=1:Sc.pix_w
%     for kj=1:Sc.pix_h
%         r = sqrt((ki-(Sc.pix_w-1)/2)^2+(kj-(Sc.pix_h-1)/2)^2)*Sc.wid/Sc.pix_w;
%         Ph1(ki,kj)=mod( Lsr.w0/Lsr.c * (sind(conic) * r ),2*pi);  
%     end
% end
% Ph2 = zeros(Sc.pix_w,Sc.pix_h);
% for ki=1:Sc.pix_w
%   r = sqrt((ki-(Sc.pix_w-1)/2)^2)*Sc.wid/Sc.pix_w;
%   Ph2(ki,:) = mod( Lsr.w0/Lsr.c * (sind(alpha) * r ),2*pi);  
% end
% Ph = mod(Ph1+Ph2,2*pi);
Sc.E_x=Sc.E_x.*exp(1i*Ph);
Sc.E_y=Sc.E_y.*exp(1i*Ph);
end
