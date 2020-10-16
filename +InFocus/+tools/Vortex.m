function [ Sc,Ph ] = Vortex(Sc,Lsr,num,max)
%VERTEX �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

angle=Angle(Sc)+2*pi/3;
Ph=zeros(Sc.pix_w,Sc.pix_h);
full=2*pi/num;
for ki=1:Sc.pix_w
    for kj=1:Sc.pix_h
        Ph(ki,kj)=mod(angle(ki,kj),full)/full*max;
    end
end
n = 12;
for j=1:n
 Ph(Ph>=(j-1)*pi/n & Ph<j*pi/n) = (j-1)*pi/n;
end

Sc.E_x=Sc.E_x.*exp(1i*Ph);
Sc.E_y=Sc.E_y.*exp(1i*Ph);
end