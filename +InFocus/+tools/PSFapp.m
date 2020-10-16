function [  I_fout, E_fxout, E_fyout, E_fzout ] = PSFapp( c,polarization, ~, para1,~, para2, ~, para3, ~, para4, ~, para5, ~, para6, ~, para7)
%% Calculation PSF of input light field with complex amplitude c
%%  Form
% [intensityOut] = PSFaberration(c, polarization,...
%                               'laser',[lambda],...
%                               'obj',[NA,f,entrancePupil,Zernike],...
%                               'materials',[n1,n2,d],...
%                               'inXY', [pixelsIncidentPlane,rhoIncident],...
%                               'outZ', [pixelsOutputZ,rhoZ],...
%                               'outX', [pixelsOutputX,physicalStart,physicalEnd],...
%                               'outY', [pixelsOutputY,physicalStart,physicalEnd])
%%  Description
% Calculation PSF of input light field with complex amplitude g
%% Input
% c (m,m) m x m elements array
% polarization '.' string
% 'laser',[lambda],(1,1) scalar value
% 'obj',[NA,f,entrancePupil, Zernike] ,(1,4) two elements array
% 'materials',[n1,n2,d], (1,3) three elements array
% 'inXY', [pixelsIncidentPlane,rhoIncident],...
% 'outZ', [pixelsOutputZ,rhoZ],...
% 'outX', [pixelsOutputX,physicalStart,physicalEnd],...
% 'outY', [pixelsOutputY,physicalStart,physicalEnd]
%% Output
% I_fout, E_fxout, E_fyout, E_fzout (m,m,mz)


%% Initialize

lambda         = para1(1); % parameter 1
NA             = para2(1); f         = para2(2); entrancePupil = para2(3); 
zernike1       = para2(4); zernike4  = para2(5); zernike11     = para2(6); 
zernike22      = para2(7); zernike37 = para2(8); % parameter 2
n1             = para3(1); n2        = para3(2); d             = - para3(3)/n2; % parameter 3
rhoIncident    = para4(1); pixIn     = para4(2); % parameter 4
pixOutZ        = para5(1); rhoZ      = para5(2);
zOut1          = para5(3); zOut2     = para5(4); % parameter 5
pixOutX				 = para6(1); xOut1     = para6(2); xOut2   = para6(3); % parameter 6
pixOutY				 = para7(1); yOut1     = para7(2); yOut2   = para7(3); % parameter 7

%% Incident plane
rPupil       = entrancePupil;    % Entrance Pupil Radius
k0           = 2*pi/lambda;
k1           = k0*n1;
k2           = k0*n2;
xIn1         = -pixIn/2 * rhoIncident; 
xIn2         = pixIn/2 * rhoIncident; 
meshIn       = (xIn2-xIn1) * (0:pixIn-1) / pixIn+xIn1;
[XIN,YIN]    = meshgrid(meshIn);
Aperture  	 = double((sqrt(XIN.^2+YIN.^2)<=rPupil));
Eix          = zeros(size(c));
Eiy 				 = zeros(size(c));
Eiz 				 = zeros(size(c));
switch polarization
    case 'linear_x'
        Eix=c.*Aperture;
    case 'linear_y'
        Eiy=c.*Aperture;
    case 'circular'
        Eix=c.*Aperture;
        Eiy=c.*Aperture*1i;
    case 'radial'
        Eix=c.*Aperture.*cos(phi);
        Eiy=c.*Aperture.*sin(phi);
    case 'azimuthal'
        Eix=-c.*Aperture.*sin(phi);
        Eiy=c.*Aperture.*cos(phi);
end

%% parallel calculation
updateWaitbar = InFocus.tools.waitbarParfor(pixOutZ, "Calculation in progress...");
if (zernike1 == 0 && zernike4 == 0 && zernike11 == 0 && zernike22 == 0 && zernike37 ==0) % No lens aberration
  if (n1==n2 || d==0) % No lens aberration + No interface aberration
    phi       	= atan2(YIN,XIN);
    theta     	= asin(Aperture.*sqrt(XIN.^2+YIN.^2)*NA/(n1*rPupil));
    Eox         = 0.5.*(Eix.*((cos(theta)+1)+(cos(theta)-1).*cos(2.*phi))+...
                       Eiy.*(cos(theta)-1).*sin(2.*phi)+...
                       Eiz.*2.*sin(theta).*cos(phi))./sqrt(cos(theta));
    Eoy         = 0.5.*(Eix.*(cos(theta)-1).*sin(2.*phi)+...
                       Eiy.*((cos(theta)+1)-(cos(theta)-1).*cos(2.*phi))+...
                       Eiz.*2.*sin(theta).*sin(phi))./sqrt(cos(theta));
    Eoz         = (-Eix.*sin(theta).*cos(phi)-...
                    Eiy.*sin(theta).*sin(phi)+...
                    Eiz.*cos(theta))./sqrt(cos(theta));
    ax        = exp(1j*2*pi*xOut1*(xIn2-xIn1)*NA/(rPupil*lambda*pixIn));
    wx        = exp(-1j*2*pi*(xOut2-xOut1)*(xIn2-xIn1)*NA/(rPupil*lambda*pixIn*pixOutX));
    zOut      = (zOut1:rhoZ:zOut2); % sampled region
    factor_fx = exp(-1j*2*pi*xIn1*xOut1*NA/(lambda*rPupil))*exp(-1j*2*pi*(xOut2-xOut1)*xIn1*(0:pixOutX-1)'*NA/(lambda*rPupil*pixOutX));
    factor_fx = factor_fx(:,ones(1,pixOutY));
    ay        = exp(1j*2*pi*yOut1*(xIn2-xIn1)*NA/(rPupil*lambda*pixIn));
    wy        = exp(-1j*2*pi*(yOut2 - yOut1)*(xIn2-xIn1)*NA/(rPupil*lambda*pixIn*pixOutY));
    factor_fy = exp(-1j*2*pi*xIn1*yOut1*NA/(lambda*rPupil))*exp(-1j*2*pi*(yOut2-yOut1)*xIn1*(0:pixOutY-1)'*NA/(lambda*rPupil*pixOutY));
    factor_fy = factor_fy(:,ones(1,pixIn));
    E_fx      = zeros(pixOutY,pixOutX,pixOutZ);
    E_fy      = zeros(pixOutY,pixOutX,pixOutZ);
    E_fz      = zeros(pixOutY,pixOutX,pixOutZ);
    parfor zz=1:pixOutZ
        E_fx(:,:,zz)=(czt((czt(Eox.*exp(-1j*k1*cos(theta)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';
        E_fy(:,:,zz)=(czt((czt(Eoy.*exp(-1j*k1*cos(theta)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';
        E_fz(:,:,zz)=(czt((czt(Eoz.*exp(-1j*k1*cos(theta)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';
        updateWaitbar();
    end
  elseif (n1~=n2) % No lens aberration + Interface aberration
    phi       	= atan2(YIN,XIN);
    theta1    	= asin(Aperture.*(sqrt(XIN.^2+YIN.^2)*NA/(n1*rPupil)));
    theta2      = asin(n1.*sin(theta1)/n2);
    tp        	= 2*n1*cos(theta1) ./ (n2*cos(theta1)+n1*cos(theta2));
    ts       		= 2*n1*cos(theta1) ./ (n1*cos(theta1)+n2*cos(theta2));
    phiTheta2 	= exp(1j*k0*d.*(n2.*cos(theta2) - n1.*cos(theta1)));
    Eox 				= n2/n1*0.5*(Eix.*((tp.*cos(theta2)+ts)+(tp.*cos(theta2)-ts).*cos(2*phi))+...
                                Eiy.*(tp.*cos(theta2)-ts)./sin(2.*phi)+...
                                Eiz.*2.*tp.*sin(theta2).*cos(phi))./sqrt(cos(theta1)).*phiTheta2;
    Eoy 				= n2/n1*0.5*(Eix.*(tp.*cos(theta2)-ts).*sin(2*phi)+...
                                Eiy.*((tp.*cos(theta2)+ts)-(tp.*cos(theta2)-ts).*cos(2*phi))+...
                                Eiz.*2.*tp.*sin(theta2).*sin(phi))./sqrt(cos(theta1)).*phiTheta2;
    Eoz 				= n2/n1*(-Eix.*tp.*sin(theta2).*cos(phi)-...
                                Eiy.*tp.*sin(theta2).*sin(phi)+...
                                Eiz.*tp.*cos(theta2))./sqrt(cos(theta1)).*phiTheta2;
    ax        = exp(1j*2*pi*xOut1*(xIn2-xIn1)*NA/(rPupil*lambda*pixIn));
    wx        = exp(-1j*2*pi*(xOut2-xOut1)*(xIn2-xIn1)*NA/(rPupil*lambda*pixIn*pixOutX));
    zOut      = (zOut1:rhoZ:zOut2); % sampled region
    factor_fx = exp(-1j*2*pi*xIn1*xOut1*NA/(lambda*rPupil)) * exp(-1j*2*pi*(xOut2-xOut1)*xIn1*(0:pixOutX-1)'*NA/(lambda*rPupil*pixOutX));
    factor_fx = factor_fx(:,ones(1,pixOutY));
    ay        = exp(1j*2*pi*yOut1*(xIn2-xIn1)*NA/(rPupil*lambda*pixIn));
    wy        = exp(-1j*2*pi*(yOut2 - yOut1)*(xIn2-xIn1)*NA/(rPupil*lambda*pixIn*pixOutY));
    factor_fy = exp(-1j*2*pi*xIn1*yOut1*NA/(lambda*rPupil))*exp(-1j*2*pi*(yOut2-yOut1)*xIn1*(0:pixOutY-1)'*NA/(lambda*rPupil*pixOutY));
    factor_fy = factor_fy(:,ones(1,pixIn));
    E_fx      = zeros(pixOutY,pixOutX,pixOutZ);
    E_fy      = zeros(pixOutY,pixOutX,pixOutZ);
    E_fz      = zeros(pixOutY,pixOutX,pixOutZ);
    parfor zz=1:pixOutZ
        updateWaitbar();
        E_fx(:,:,zz)=(czt((czt(Eox.*exp(-1j*k2*cos(theta2)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';
        E_fy(:,:,zz)=(czt((czt(Eoy.*exp(-1j*k2*cos(theta2)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';
        E_fz(:,:,zz)=(czt((czt(Eoz.*exp(-1j*k2*cos(theta2)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';  
    end
  end

elseif (zernike1 ~= 0 || zernike4 ~= 0 || zernike11 ~= 0 || zernike22 ~= 0 || zernike37 ~=0) % % Lens aberration
  if (n1==n2 || d==0) % % Lens aberrationd + No interface aberration 
    phi       	= atan2(YIN,XIN);
    theta     	= asin(sqrt(XIN.^2+YIN.^2)*NA/(n1*rPupil));
    r           = sqrt(XIN.^2+YIN.^2).*Aperture/rPupil;
%     Eix         = Eix .* (exp(1j*2*pi*(W040 * r.^4))); % Spherical aberration
%     Eiy         = Eiy .* (exp(1j*2*pi*(W040 * r.^4)));
%     Eiz         = Eiz .* (exp(1j*2*pi*(W040 * r.^4)));
    Eix         = Eix .* (exp(1j*2*pi* (zernike1 * sqrt(1) * (1)+ ...
                                        zernike4 * sqrt(3) * (2 * r.^2 -1)+ ...
                                        zernike11 * sqrt(5) * (6 * r.^4 - 6 * r.^2 +1)+ ...
                                        zernike22 * sqrt(7) * (20 * r.^6 - 30 * r.^4 + 12 * r.^2 -1) + ...
                                        zernike37 * sqrt(9) * (70 * r.^8 - 140 * r.^6 + 90 * r.^4 - 20 * r.^2 +1)))); % Spherical aberration
    Eiy         = Eiy .* (exp(1j*2*pi* (zernike1 * sqrt(1) * (1)+ ...
                                    zernike4 * sqrt(3) * (2 * r.^2 -1)+ ...
                                    zernike11 * sqrt(5) * (6 * r.^4 - 6 * r.^2 +1)+ ...
                                    zernike22 * sqrt(7) * (20 * r.^6 - 30 * r.^4 + 12 * r.^2 -1) + ...
                                    zernike37 * sqrt(9) * (70 * r.^8 - 140 * r.^6 + 90 * r.^4 - 20 * r.^2 +1)))); % Spherical aberration
    Eiz         = Eiz .* (exp(1j*2*pi* (zernike1 * sqrt(1) * (1)+ ...
                                    zernike4 * sqrt(3) * (2 * r.^2 -1)+ ...
                                    zernike11 * sqrt(5) * (6 * r.^4 - 6 * r.^2 +1)+ ...
                                    zernike22 * sqrt(7) * (20 * r.^6 - 30 * r.^4 + 12 * r.^2 -1) + ...
                                    zernike37 * sqrt(9) * (70 * r.^8 - 140 * r.^6 + 90 * r.^4 - 20 * r.^2 +1)))); % Spherical aberration

    Eox         = 0.5.*(Eix.*((cos(theta)+1)+(cos(theta)-1).*cos(2.*phi))+...
                       Eiy.*(cos(theta)-1).*sin(2.*phi)+...
                       Eiz.*2.*sin(theta).*cos(phi))./sqrt(cos(theta));
    Eoy         = 0.5.*(Eix.*(cos(theta)-1).*sin(2.*phi)+...
                       Eiy.*((cos(theta)+1)-(cos(theta)-1).*cos(2.*phi))+...
                       Eiz.*2.*sin(theta).*sin(phi))./sqrt(cos(theta));
    Eoz         = (-Eix.*sin(theta).*cos(phi)-...
                    Eiy.*sin(theta).*sin(phi)+...
                    Eiz.*cos(theta))./sqrt(cos(theta));
    ax        = exp(1j*2*pi*xOut1*(xIn2-xIn1)*NA/(rPupil*lambda*pixIn));
    wx        = exp(-1j*2*pi*(xOut2-xOut1)*(xIn2-xIn1)*NA/(rPupil*lambda*pixIn*pixOutX));
    zOut      = (zOut1:rhoZ:zOut2); % sampled region
    factor_fx = exp(-1j*2*pi*xIn1*xOut1*NA/(lambda*rPupil))*exp(-1j*2*pi*(xOut2-xOut1)*xIn1*(0:pixOutX-1)'*NA/(lambda*rPupil*pixOutX));
    factor_fx = factor_fx(:,ones(1,pixOutY));
    ay        = exp(1j*2*pi*yOut1*(xIn2-xIn1)*NA/(rPupil*lambda*pixIn));
    wy        = exp(-1j*2*pi*(yOut2 - yOut1)*(xIn2-xIn1)*NA/(rPupil*lambda*pixIn*pixOutY));
    factor_fy = exp(-1j*2*pi*xIn1*yOut1*NA/(lambda*rPupil))*exp(-1j*2*pi*(yOut2-yOut1)*xIn1*(0:pixOutY-1)'*NA/(lambda*rPupil*pixOutY));
    factor_fy = factor_fy(:,ones(1,pixIn));
    E_fx      = zeros(pixOutY,pixOutX,pixOutZ);
    E_fy      = zeros(pixOutY,pixOutX,pixOutZ);
    E_fz      = zeros(pixOutY,pixOutX,pixOutZ);
    parfor zz=1:pixOutZ
        updateWaitbar();
        E_fx(:,:,zz)=(czt((czt(Eox.*exp(-1j*k1*cos(theta)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';
        E_fy(:,:,zz)=(czt((czt(Eoy.*exp(-1j*k1*cos(theta)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';
        E_fz(:,:,zz)=(czt((czt(Eoz.*exp(-1j*k1*cos(theta)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';
    end
  elseif (n1~=n2) % Interface aberration  considered
    phi       	= atan2(YIN,XIN);
    r           = sqrt(XIN.^2+YIN.^2).*Aperture/rPupil;
%     Eix         = Eix .* (exp(1j*2*pi*(W040 * r.^4))); % Spherical aberration
%     Eiy         = Eiy .* (exp(1j*2*pi*(W040 * r.^4)));
%     Eiz         = Eiz .* (exp(1j*2*pi*(W040 * r.^4)));
    Eix         = Eix .* (exp(1j*2*pi* (zernike1 * sqrt(1) * (1)+ ...
                                        zernike4 * sqrt(3) * (2 * r.^2 -1)+ ...
                                        zernike11 * sqrt(5) * (6 * r.^4 - 6 * r.^2 +1)+ ...
                                        zernike22 * sqrt(7) * (20 * r.^6 - 30 * r.^4 + 12 * r.^2 -1) + ...
                                        zernike37 * sqrt(9) * (70 * r.^8 - 140 * r.^6 + 90 * r.^4 - 20 * r.^2 +1)))); % Spherical aberration
    Eiy         = Eiy .* (exp(1j*2*pi* (zernike1 * sqrt(1) * (1)+ ...
                                    zernike4 * sqrt(3) * (2 * r.^2 -1)+ ...
                                    zernike11 * sqrt(5) * (6 * r.^4 - 6 * r.^2 +1)+ ...
                                    zernike22 * sqrt(7) * (20 * r.^6 - 30 * r.^4 + 12 * r.^2 -1) + ...
                                    zernike37 * sqrt(9) * (70 * r.^8 - 140 * r.^6 + 90 * r.^4 - 20 * r.^2 +1)))); % Spherical aberration
    Eiz         = Eiz .* (exp(1j*2*pi* (zernike1 * sqrt(1) * (1)+ ...
                                    zernike4 * sqrt(3) * (2 * r.^2 -1)+ ...
                                    zernike11 * sqrt(5) * (6 * r.^4 - 6 * r.^2 +1)+ ...
                                    zernike22 * sqrt(7) * (20 * r.^6 - 30 * r.^4 + 12 * r.^2 -1) + ...
                                    zernike37 * sqrt(9) * (70 * r.^8 - 140 * r.^6 + 90 * r.^4 - 20 * r.^2 +1)))); % Spherical aberration

    theta1    	= asin(Aperture.*(sqrt(XIN.^2+YIN.^2)*NA/(n1*rPupil)));
    theta2      = asin(n1.*sin(theta1)/n2);
    tp        	= 2*n1*cos(theta1) ./ (n2*cos(theta1)+n1*cos(theta2));
    ts       		= 2*n1*cos(theta1) ./ (n1*cos(theta1)+n2*cos(theta2));
    phiTheta2 	= exp(1j*k0*d.*(n2.*cos(theta2) - n1.*cos(theta1)));
    Eox 				= n2/n1*0.5*(Eix.*((tp.*cos(theta2)+ts)+(tp.*cos(theta2)-ts).*cos(2*phi))+...
                                Eiy.*(tp.*cos(theta2)-ts)./sin(2.*phi)+...
                                Eiz.*2.*tp.*sin(theta2).*cos(phi))./sqrt(cos(theta1)).*phiTheta2;
    Eoy 				= n2/n1*0.5*(Eix.*(tp.*cos(theta2)-ts).*sin(2*phi)+...
                                Eiy.*((tp.*cos(theta2)+ts)-(tp.*cos(theta2)-ts).*cos(2*phi))+...
                                Eiz.*2.*tp.*sin(theta2).*sin(phi))./sqrt(cos(theta1)).*phiTheta2;
    Eoz 				= n2/n1*(-Eix.*tp.*sin(theta2).*cos(phi)-...
                                Eiy.*tp.*sin(theta2).*sin(phi)+...
                                Eiz.*tp.*cos(theta2))./sqrt(cos(theta1)).*phiTheta2;
    ax        = exp(1j*2*pi*xOut1*(xIn2-xIn1)*NA/(rPupil*lambda*pixIn));
    wx        = exp(-1j*2*pi*(xOut2-xOut1)*(xIn2-xIn1)*NA/(rPupil*lambda*pixIn*pixOutX));
    zOut      = (zOut1:rhoZ:zOut2); % sampled region
    factor_fx = exp(-1j*2*pi*xIn1*xOut1*NA/(lambda*rPupil))*exp(-1j*2*pi*(xOut2-xOut1)*xIn1*(0:pixOutX-1)'*NA/(lambda*rPupil*pixOutX));
    factor_fx = factor_fx(:,ones(1,pixOutY));
    ay        = exp(1j*2*pi*yOut1*(xIn2-xIn1)*NA/(rPupil*lambda*pixIn));
    wy        = exp(-1j*2*pi*(yOut2 - yOut1)*(xIn2-xIn1)*NA/(rPupil*lambda*pixIn*pixOutY));
    factor_fy = exp(-1j*2*pi*xIn1*yOut1*NA/(lambda*rPupil))*exp(-1j*2*pi*(yOut2-yOut1)*xIn1*(0:pixOutY-1)'*NA/(lambda*rPupil*pixOutY));
    factor_fy = factor_fy(:,ones(1,pixIn));
    E_fx      = zeros(pixOutY,pixOutX,pixOutZ);
    E_fy      = zeros(pixOutY,pixOutX,pixOutZ);
    E_fz      = zeros(pixOutY,pixOutX,pixOutZ);
    parfor zz=1:pixOutZ
        updateWaitbar();
        E_fx(:,:,zz)=(czt((czt(Eox.*exp(-1j*k2*cos(theta2)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';
        E_fy(:,:,zz)=(czt((czt(Eoy.*exp(-1j*k2*cos(theta2)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';
        E_fz(:,:,zz)=(czt((czt(Eoz.*exp(-1j*k2*cos(theta2)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';
    end
  end
end
%% Output    
I_f = abs(E_fx).^2+abs(E_fy).^2+abs(E_fz).^2;
if (nargout==1)
    I_fout = I_f;
elseif (nargout==4)
    E_fxout = E_fx;
    E_fyout = E_fy;
    E_fzout = E_fz;
    I_fout  = I_f;
end
end