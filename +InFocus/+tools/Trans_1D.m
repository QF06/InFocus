function [Time, Freq] = Trans_1D(Data1, Data2, Lsr, Direction)
%Trans_1D Summary of this function goes here
%   Detailed explanation goes here
%
if Direction>0
    Time=Data1;
    Freq=zeros(1,Lsr.pix_f);
	for kw=1:(Lsr.pix_f-1)
        for kt=1:(Lsr.pix_t-1)
            t=Lsr.t(kt); w=Lsr.w(kw)-Lsr.w0;
            FF=Time(kt)*exp(-1i*w*t);
            Freq(kw)=Freq(kw)+FF*Lsr.dt;
        end
	end
elseif Direction<0
    Time=zeros(1,Lsr.pix_t);
    Freq=Data2;
	for kt=1:(Lsr.pix_t-1)
        for kw=1:(Lsr.pix_f-1)
            t=Lsr.t(kt); w=Lsr.w(kw);
            FF=Freq(kw)*exp(1i*w*t);
            Time(kt)=Time(kt)+FF*(2*pi*Lsr.c/Lsr.df);
        end
    end
Freq=Freq/max(Freq);
Time=Time/max(Time);
end
