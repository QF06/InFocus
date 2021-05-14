%% This demonstration calculate the intensity and field distributions near the focus.
% An x-polarized gaussian beam with 5.2 mm radius is focused by a Mitutoyo 50x NIR
% objective lens.
% The plotted images show the intensity as well as the field distribution
% of each components (Ix, Iy, Iz, Ex, Ey, Ez).

clc; close all; clear;     
%% Initiate the laser and the beam profile
wavelength       = 1555 * 1e-9;
laser            = InFocus.tools.Laser(1555e-9); % Input: wavelength
Polarization     = 'linear_x';
BeamProfile      = 'Gauss'; % 'Gauss', 'Vortex', etc
Wx               = 5.2 * 1e-3; % beam radius [m]
Wy               = 5.2 * 1e-3;
physicalWH       = 20 * 1e-3; % sampling area [m] 
pixelWH          = 801; % sampling pixels 
sectionIn        = InFocus.tools.Section('pixelWH',[pixelWH, pixelWH],...
       'physicalWH',[physicalWH, physicalWH],...
       'vector',[0,0,1],'offset',[0,0,0]);
intensityIn.all  = zeros(sectionIn.pix_w,sectionIn.pix_h);
intensityIn.freq = zeros(sectionIn.pix_w,sectionIn.pix_h);
sectionIn.E_x    = InFocus.tools.Gauss(sectionIn,[Wx, Wy],[0,0]);
intensityIn.freq(:,:,1) = abs((sectionIn.E_x).^2);
intensityIn.all         = intensityIn.all+intensityIn.freq(:,:,1);
%% Plot intensity profile of the incidenct beam
XTickIn          = 1 : (sectionIn.pix_w-1)/5 : sectionIn.pix_w;
XTickLabelIn     = round((- sectionIn.wid/2 : sectionIn.wid/5 : sectionIn.wid/2) * 1e3);
YTickIn          = 1 : (sectionIn.pix_h-1)/5 : sectionIn.pix_h;
YTickLabelIn     = round((- sectionIn.hei/2 : sectionIn.hei/5 : sectionIn.hei/2) * 1e3);
window=[800,600];
figure('position',[100,100,window(1),window(2)]);
ax1 = subplot(5,3,1); 
PlotIncident(ax1,intensityIn.all/max(max(max(intensityIn.all))),'Incident I_{In}',...
              XTickIn, XTickLabelIn, YTickIn, YTickLabelIn, 'hot','linear','x[mm]','y[mm]',[0,1]);
axis equal tight;
colorbar;

%% Initiate focus condition and the set the sampling
c             = sectionIn.E_x; % complex amplitude of input light field 
lambda        = laser.lambda0;
f             = 4.0 * 1e-3; % focal length [mm]
NA            = 0.65; % numerical aperature
entrancePupil = 2.6 * 1e-3; % 
Z1            = 0;
Z4            = 0;
Z11           = 0;
Z22           = 0;
Z37           = 0;
Zernike       = [Z1, Z4, Z11, Z22, Z37]; % Zernike standard coefficient
n1            = 1;
n2            = 3.475;
d             = 5 * 1e-3;
kt            = 2*pi/lambda*n2;
%-------sampling on the incident plane-------------%
rhoIncident   = sectionIn.wid/(sectionIn.pix_w-1); % pixel size on the incident plane
pixIn         = sectionIn.pix_w; % sampling number SLM plane
%-------sampling on the focal plane----------------%
pixOutX       = 201; % sampling at xy near focus
pixOutY       = pixOutX;
xOut1         = -20 * 1e-6; % AOI at xy near focus min
xOut2         = 20 * 1e-6; % AOI at xy near focus max
yOut1         = xOut1;
yOut2         = xOut2;
%-------sampling along z near the focus------------%
pixOutZ       = 251; % sampling along z
zOut1         = 3000 * 1e-6; % AOI at along z 
zOut2         = 5500 * 1e-6; % AOI at along z 
rhoZ          = (zOut2-zOut1)/(pixOutZ-1); %sampling interval along the z-direction in the focal region unit[m/pixel]
rhoXY         = (xOut2-xOut1)/(pixOutX-1);    
rho           = [rhoXY, rhoZ];

%% Calcuate the field near the focus
tic
[  I_fout, E_fxout, E_fyout, E_fzout ] = InFocus.tools.PSFapp( c, Polarization,...
                        'laser',lambda,...
                        'obj',[NA,f, entrancePupil, Zernike],...
                        'material',[n1,n2,d],...
                        'inXY',[rhoIncident,pixIn],...
                        'outZ',[pixOutZ,rhoZ,zOut1,zOut2],...
                        'outX',[pixOutX, xOut1,xOut2],...
                        'outY',[pixOutY,yOut1,yOut2]);
toc
%% Plot E_x E_y E_z and I at the focal plane
xmid          = fix(pixOutX/2);
ymid          = fix(pixOutY/2);
intensityZ    = squeeze(I_fout(xmid,ymid,:));  
[maxI, maxP]  = max(intensityZ);
XTickOutxyf  = 1 : (pixOutX-1)/5 : pixOutX;
XTickLabelxyf = (xOut1 : (xOut2-xOut1)/5 : xOut2) * 1e6;
YTickOutxyf   = 1 : (pixOutY-1)/5 : pixOutY;
YTickLabelxyf = (yOut1 : (yOut2-yOut1)/5 : yOut2) * 1e6;
subplot(5,3,4);
PlotIncident(gca, I_fout(:,:,maxP)/maxI, 'focus I_f',...
              XTickOutxyf, XTickLabelxyf, YTickOutxyf, YTickLabelxyf, 'hot','linear','x[\mum]','y[\mum]',[0,1]);
axis equal tight;
colorbar;
subplot(5,3,7);
maxEx          = max(max(max(abs(E_fxout))));
PlotIncident(gca, abs(E_fxout(:,:,maxP)/maxEx), 'focus Ex_f',...
              XTickOutxyf, XTickLabelxyf, YTickOutxyf, YTickLabelxyf, 'hot','linear','x[\mum]','y[\mum]',[0,1]);
axis equal tight;
colorbar;
subplot(5,3,10);
maxEy          = max(max(max(abs(E_fyout))));
PlotIncident(gca, abs(E_fyout(:,:,maxP)/maxEy), 'focus Ey_f',...
              XTickOutxyf, XTickLabelxyf, YTickOutxyf, YTickLabelxyf, 'hot','linear','x[\mum]','y[\mum]',[0,1]);
axis equal tight;
colorbar;
subplot(5,3,13);
maxEz          = max(max(max(abs(E_fzout))));
PlotIncident(gca, abs(E_fzout(:,:,maxP)/maxEz), 'focus Ez_f',...
              XTickOutxyf, XTickLabelxyf, YTickOutxyf, YTickLabelxyf, 'hot','linear','x[\mum]','y[\mum]',[0,1]);
axis equal tight;
colorbar;

%% Plot E_x E_y E_z and I along x-z and y-z plane 
XTickOutfz      = 1 : (pixOutZ-1)/5 : pixOutZ;
XTickLabelOutfz = (zOut1:(zOut2-zOut1)/5:zOut2) * 1e6;
YTickOutfx      = 1 : (pixOutX-1)/5 : pixOutX;
YTickLabelOutfx = (xOut1 : (xOut2-xOut1)/5 : xOut2) * 1e6;
YTickOutfy      = 1 : (pixOutY-1)/5 : pixOutY;
YTickLabelOutfy = (yOut1 : (yOut2-yOut1)/5 : yOut2) * 1e6;
%------------I_fout----------------------------%
I_foutXZC       = abs(squeeze(I_fout(:, ymid, :)));
I_foutYZC       = abs(squeeze(I_fout(xmid, :, :)));
subplot(5,3,5);
PlotIncident(gca, I_foutXZC/maxI, 'focus I_f',...
              XTickOutfz, XTickLabelOutfz, YTickOutfx, YTickLabelOutfx, 'hot','linear','z[\mum]','y[\mum]',[0,1]);
subplot(5,3,6);
PlotIncident(gca, I_foutYZC/maxI, 'focus I_f',...
              XTickOutfz, XTickLabelOutfz, YTickOutfy, YTickLabelOutfy, 'hot','linear','z[\mum]','x[\mum]',[0,1]);   
%------------E_fxout----------------------------%            
E_fxoutXZC     = abs(squeeze(E_fxout(:, ymid, :)));
E_fxoutYZC     = abs(squeeze(E_fxout(xmid, :, :)));
subplot(5,3,8);
PlotIncident(gca, E_fxoutXZC/maxEx, 'focus E_fx',...
              XTickOutfz, XTickLabelOutfz, YTickOutfx, YTickLabelOutfx, 'hot','linear','z[\mum]','y[\mum]',[0,1]);
subplot(5,3,9);
PlotIncident(gca, E_fxoutYZC/maxEx, 'focus E_fx',...
              XTickOutfz, XTickLabelOutfz, YTickOutfy, YTickLabelOutfy, 'hot','linear','z[\mum]','x[\mum]',[0,1]); 
%------------E_fyout----------------------------%               
E_fyoutXZC     = abs(squeeze(E_fyout(:, ymid, :)));
E_fyoutYZC     = abs(squeeze(E_fyout(xmid, :, :)));
subplot(5,3,11);
PlotIncident(gca, E_fyoutXZC/maxEy, 'focus E_fy',...
              XTickOutfz, XTickLabelOutfz, YTickOutfx, YTickLabelOutfx, 'hot','linear','z[\mum]','y[\mum]',[0,1]);
subplot(5,3,12);
PlotIncident(gca, E_fyoutYZC/maxEy, 'focus E_fy',...
              XTickOutfz, XTickLabelOutfz, YTickOutfy, YTickLabelOutfy, 'hot','linear','z[\mum]','x[\mum]',[0,1]); 
%------------E_fzout----------------------------%              
E_fzoutXZC     = abs(squeeze(E_fzout(:, ymid, :)));
E_fzoutYZC     = abs(squeeze(E_fzout(xmid, :, :)));
subplot(5,3,14);
PlotIncident(gca, E_fzoutXZC/maxEz, 'focus E_fz',...
              XTickOutfz, XTickLabelOutfz, YTickOutfx, YTickLabelOutfx, 'hot','linear','z[\mum]','y[\mum]',[0,1]);
subplot(5,3,15);
PlotIncident(gca, E_fzoutYZC/maxEz, 'focus E_fz',...
              XTickOutfz, XTickLabelOutfz, YTickOutfy, YTickLabelOutfy, 'hot','linear','z[\mum]','x[\mum]',[0,1]); 

           
function [] = PlotIncident(ax,intensity,name,XTick,XTickLabel,YTick,YTickLabel,map,scale,XLabel,YLabel,Caxis)

imagesc(intensity,'Parent',ax);
set(ax,'YDir','normal');
set(ax,'XTick',XTick);
set(ax,'XTickLabel',XTickLabel);
set(ax,'YTick',YTick);
set(ax,'YTickLabel',YTickLabel);
set(ax,'Colorscale',scale);
title(ax, name);
colormap(ax,map);
xlabel(XLabel);
ylabel(YLabel);
caxis(ax,Caxis);
end
