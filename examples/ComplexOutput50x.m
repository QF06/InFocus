clc; close all; clear;     
%% Initiate the laser and the beam profile
list             = SquareWave(101,300e-15, 70e-15); 
wavelength       = 1555 * 1e-9;
laser            = Laser('lambda',[101,300e-9, wavelength ,35e-9],...
                         'time',[101,9800e-15,0,900e-15],...
                         'F_trans',0,list);
Polarization     = 'linear_x';
BeamProfile      = 'Gauss'; % 'Gauss', 'Vortex', etc
Wx               = 5.2 * 1e-3; % beam wasit [m]
Wy               = 5.2 * 1e-3;
physicalWH       = 20 * 1e-3; % sampling area [m] 
pixelWH          = 801; % sampling pixels 
sectionIn        = Section('pixelWH',[pixelWH, pixelWH],...
       'physicalWH',[physicalWH, physicalWH],...
       'vector',[0,0,1],'offset',[0,0,0]);
intensityIn.all  = zeros(sectionIn.pix_w,sectionIn.pix_h);
intensityIn.freq = zeros(sectionIn.pix_w,sectionIn.pix_h, laser.pixelWavelength);
sectionIn.E_x    = Gauss(sectionIn,[Wx, Wy],[0,0]);
intensityIn.freq(:,:,1) = abs((sectionIn.E_x).^2);
intensityIn.all         = intensityIn.all+intensityIn.freq(:,:,1);
%% Plot intensity profile of the incidenct beam
XTickIn          = 1 : (sectionIn.pix_w-1)/5 : sectionIn.pix_w;
XTickLabelIn     = round((- sectionIn.wid/2 : sectionIn.wid/5 : sectionIn.wid/2) * 1e3);
YTickIn          = 1 : (sectionIn.pix_h-1)/5 : sectionIn.pix_h;
YTickLabelIn     = round((- sectionIn.hei/2 : sectionIn.hei/5 : sectionIn.hei/2) * 1e3);
window=[1000,800];
figure('position',[100,50,window(1),window(2)]);
ax1 = subplot(5,3,1); 
PlotIncidentF(ax1,intensityIn.all/max(max(max(intensityIn.all))),'Incident I_{In}',...
              XTickIn, XTickLabelIn, YTickIn, YTickLabelIn, 'fire','linear','x[mm]','y[mm]',[0,1]);
axis equal tight;
colorbar;

%% Initiate focus condition and the set the sampling
c             = sectionIn.E_x; % complex amplitude of input light field 
lambda        = laser.lambda0;
f             = 4.0 * 1e-3; % focal length [mm]
NA            = 0.65; % numerical aperature
entrancePupil = 2.6 * 1e-3; % 
W040          = 0; %seidel coefficient
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
pixOutZ       = 2501; % sampling along z
zOut1         = 3000 * 1e-6; % AOI at along z 
zOut2         = 5500 * 1e-6; % AOI at along z 
rhoZ          = (zOut2-zOut1)/(pixOutZ-1); %sampling interval along the z-direction in the focal region unit[m/pixel]
rhoXY         = (xOut2-xOut1)/(pixOutX-1);    
rho           = [rhoXY, rhoZ];

%% Calcuate the field near the focus
[  I_fout, E_fxout, E_fyout, E_fzout ] = PSFapp( c, Polarization,...
                        'laser',lambda,...
                        'obj',[NA,f, entrancePupil, W040],...
                        'material',[n1,n2,d],...
                        'inXY',[rhoIncident,pixIn],...
                        'outZ',[pixOutZ,rhoZ,zOut1,zOut2],...
                        'outX',[pixOutX, xOut1,xOut2],...
                        'outY',[pixOutY,yOut1,yOut2]);
% E_fx = E_fxout(:,:,1);     
% E_fy = E_fyout(:,:,1);  
% E_fz = E_fzout(:,:,1);  
% save('E_fx.mat','E_fx');
% save('E_fy.mat','E_fy');
% save('E_fz.mat','E_fz');
% fid = fopen('pixelRatio.txt','wt');
% fprintf(fid,'rhoXY = %0.3f\n', rho(1)*1e6);
% fprintf(fid,'rhoZ = %0.3f\n', rho(2)*1e6);
% fprintf(fid,'unit: um/pixel');

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
PlotIncidentF(gca, I_fout(:,:,maxP)/maxI, 'focus I_f',...
              XTickOutxyf, XTickLabelxyf, YTickOutxyf, YTickLabelxyf, 'fire','linear','x[\mum]','y[\mum]',[0,1]);
axis equal tight;
colorbar;
subplot(5,3,7);
maxEx          = max(max(max(abs(E_fxout))));
PlotIncidentF(gca, abs(E_fxout(:,:,maxP)/maxEx), 'focus Ex_f',...
              XTickOutxyf, XTickLabelxyf, YTickOutxyf, YTickLabelxyf, 'fire','linear','x[\mum]','y[\mum]',[0,1]);
axis equal tight;
colorbar;
subplot(5,3,10);
maxEy          = max(max(max(abs(E_fyout))));
PlotIncidentF(gca, abs(E_fyout(:,:,maxP)/maxEy), 'focus Ey_f',...
              XTickOutxyf, XTickLabelxyf, YTickOutxyf, YTickLabelxyf, 'fire','linear','x[\mum]','y[\mum]',[0,1]);
axis equal tight;
colorbar;
subplot(5,3,13);
maxEz          = max(max(max(abs(E_fzout))));
PlotIncidentF(gca, abs(E_fzout(:,:,maxP)/maxEz), 'focus Ez_f',...
              XTickOutxyf, XTickLabelxyf, YTickOutxyf, YTickLabelxyf, 'fire','linear','x[\mum]','y[\mum]',[0,1]);
axis equal tight;
colorbar;

%% Plot the E_x, E_y, 
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
PlotIncidentF(gca, I_foutXZC/maxI, 'focus I_f',...
              XTickOutfz, XTickLabelOutfz, YTickOutfx, YTickLabelOutfx, 'fire','linear','z[\mum]','y[\mum]',[0,1]);
subplot(5,3,6);
PlotIncidentF(gca, I_foutYZC/maxI, 'focus I_f',...
              XTickOutfz, XTickLabelOutfz, YTickOutfy, YTickLabelOutfy, 'fire','linear','z[\mum]','x[\mum]',[0,1]);   
%------------E_fxout----------------------------%            
E_fxoutXZC     = abs(squeeze(E_fxout(:, ymid, :)));
E_fxoutYZC     = abs(squeeze(E_fxout(xmid, :, :)));
subplot(5,3,8);
PlotIncidentF(gca, E_fxoutXZC/maxEx, 'focus E_fx',...
              XTickOutfz, XTickLabelOutfz, YTickOutfx, YTickLabelOutfx, 'fire','linear','z[\mum]','y[\mum]',[0,1]);
subplot(5,3,9);
PlotIncidentF(gca, E_fxoutYZC/maxEx, 'focus E_fx',...
              XTickOutfz, XTickLabelOutfz, YTickOutfy, YTickLabelOutfy, 'fire','linear','z[\mum]','x[\mum]',[0,1]); 
%------------E_fyout----------------------------%               
E_fyoutXZC     = abs(squeeze(E_fyout(:, ymid, :)));
E_fyoutYZC     = abs(squeeze(E_fyout(xmid, :, :)));
subplot(5,3,11);
PlotIncidentF(gca, E_fyoutXZC/maxEy, 'focus E_fy',...
              XTickOutfz, XTickLabelOutfz, YTickOutfx, YTickLabelOutfx, 'fire','linear','z[\mum]','y[\mum]',[0,1]);
subplot(5,3,12);
PlotIncidentF(gca, E_fyoutYZC/maxEy, 'focus E_fy',...
              XTickOutfz, XTickLabelOutfz, YTickOutfy, YTickLabelOutfy, 'fire','linear','z[\mum]','x[\mum]',[0,1]); 
%------------E_fzout----------------------------%              
E_fzoutXZC     = abs(squeeze(E_fzout(:, ymid, :)));
E_fzoutYZC     = abs(squeeze(E_fzout(xmid, :, :)));
subplot(5,3,14);
PlotIncidentF(gca, E_fzoutXZC/maxEz, 'focus E_fz',...
              XTickOutfz, XTickLabelOutfz, YTickOutfx, YTickLabelOutfx, 'fire','linear','z[\mum]','y[\mum]',[0,1]);
subplot(5,3,15);
PlotIncidentF(gca, E_fzoutYZC/maxEz, 'focus E_fz',...
              XTickOutfz, XTickLabelOutfz, YTickOutfy, YTickLabelOutfy, 'fire','linear','z[\mum]','x[\mum]',[0,1]); 
%% Plot polarization compare 
window=[1000,300];
figure('position',[100,50,window(1),window(2)]);
subplot(2,4,1);
maxEx          = max(max(max(abs(E_fxout))));
PlotIncidentF(gca, abs(E_fxout(:,:,maxP)/maxEx), 'focus Ex_f',...
              XTickOutxyf, XTickLabelxyf, YTickOutxyf, YTickLabelxyf, 'fire','linear','x[\mum]','y[\mum]',[0,1]);
axis equal tight;
subplot(2,4,2);
maxEy          = max(max(max(abs(E_fyout))));
PlotIncidentF(gca, 100*abs(E_fyout(:,:,maxP)/maxEx), 'focus Ey_f',...
              XTickOutxyf, XTickLabelxyf, YTickOutxyf, YTickLabelxyf, 'fire','linear','x[\mum]','y[\mum]',[0,1]);
axis equal tight;
subplot(2,4,3);
maxEz          = max(max(max(abs(E_fzout))));
PlotIncidentF(gca, 20*abs(E_fzout(:,:,maxP)/maxEx), 'focus Ez_f',...
              XTickOutxyf, XTickLabelxyf, YTickOutxyf, YTickLabelxyf, 'fire','linear','x[\mum]','y[\mum]',[0,1]);
axis equal tight;
subplot(2,4,4);
PlotIncidentF(gca, abs(E_fxout(:,:,maxP)/maxEx), 'focus Ez_f',...
              XTickOutxyf, XTickLabelxyf, YTickOutxyf, YTickLabelxyf, 'fire','linear','x[\mum]','y[\mum]',[0,1]);
axis equal tight;
colorbar;
%------Plot the tail------------%
intensityEz    = E_fzoutXZC(xmid,:);  
[maxIz, maxPz]  = max(intensityEz);
subplot(2,4,5);
PlotIncidentF(gca, abs(E_fxout(:,:,maxPz)/maxEx), 'focus Ex_f',...
              XTickOutxyf, XTickLabelxyf, YTickOutxyf, YTickLabelxyf, 'fire','linear','x[\mum]','y[\mum]',[0,1]);
axis equal tight;
subplot(2,4,6);
PlotIncidentF(gca, 100*abs(E_fyout(:,:,maxPz)/maxEx), 'focus Ey_f',...
              XTickOutxyf, XTickLabelxyf, YTickOutxyf, YTickLabelxyf, 'fire','linear','x[\mum]','y[\mum]',[0,1]);
axis equal tight;
subplot(2,4,7);
PlotIncidentF(gca, 20*abs(E_fzout(:,:,maxPz)/maxEx), 'focus Ez_f',...
              XTickOutxyf, XTickLabelxyf, YTickOutxyf, YTickLabelxyf, 'fire','linear','x[\mum]','y[\mum]',[0,1]);
axis equal tight;
%colorbar;
            
%% functions
function p = fire(m)
%FIRE   Blue-Purple Hot colormap
% 
% FIRE(M) returns an M-by-3 matrix containing a "fire" colormap.
% FIRE, by itself, is the same length as the current figure's
% colormap. If no figure exists, MATLAB creates one.
%
% To add this colormap as a default map, use 'addpath' with the 
% directory containing 'fire.m'.
%
% To reset the colormap of the current figure use 'colormap(fire)'.
%
% see also:  HSV, GRAY, HOT, COOL, BONE, COPPER, FLAG, PINK, COLORMAP,
% RGBPLOT.
%
% To create any custom colormap, see the directions on line 23 of this
% m-file.
if nargin < 1
    m = size(get(gcf,'colormap'),1); 
end
%You can replace this M x 3 matrix with any matrix whose values range
%between 0 and 1 to create a new colormap file.  Use copy / paste to create
%a matrix like the one below, you do not have to add these values
%manually.  To create a new colormap, change 'cmap_mat' to the desired
%matrix, rename the function *and* the m-file from 'fire' to your desired
%colormap name.
cmap_mat=[
         0         0         0
         0         0    0.0275
         0         0    0.0588
         0         0    0.0863
         0         0    0.1176
         0         0    0.1490
         0         0    0.1765
         0         0    0.2078
         0         0    0.2392
         0         0    0.2549
         0         0    0.2706
         0         0    0.2902
         0         0    0.3059
         0         0    0.3216
         0         0    0.3412
         0         0    0.3569
    0.0039         0    0.3765
    0.0157         0    0.3922
    0.0275         0    0.4078
    0.0392         0    0.4235
    0.0510         0    0.4431
    0.0627         0    0.4588
    0.0745         0    0.4745
    0.0863         0    0.4902
    0.0980         0    0.5098
    0.1098         0    0.5255
    0.1216         0    0.5412
    0.1333         0    0.5608
    0.1451         0    0.5765
    0.1569         0    0.5922
    0.1686         0    0.6118
    0.1804         0    0.6275
    0.1922         0    0.6471
    0.2039         0    0.6588
    0.2157         0    0.6706
    0.2275         0    0.6863
    0.2392         0    0.6980
    0.2510         0    0.7098
    0.2627         0    0.7255
    0.2745         0    0.7373
    0.2863         0    0.7529
    0.2980         0    0.7647
    0.3098         0    0.7804
    0.3216         0    0.7922
    0.3333         0    0.8078
    0.3451         0    0.8196
    0.3569         0    0.8353
    0.3686         0    0.8471
    0.3843         0    0.8627
    0.3961         0    0.8627
    0.4078         0    0.8667
    0.4196         0    0.8706
    0.4314         0    0.8745
    0.4431         0    0.8784
    0.4549         0    0.8824
    0.4667         0    0.8863
    0.4784         0    0.8902
    0.4902         0    0.8784
    0.5020         0    0.8706
    0.5137         0    0.8627
    0.5255         0    0.8549
    0.5373         0    0.8471
    0.5490         0    0.8392
    0.5608         0    0.8314
    0.5725         0    0.8235
    0.5804         0    0.8078
    0.5882         0    0.7922
    0.5961         0    0.7804
    0.6039         0    0.7647
    0.6118         0    0.7490
    0.6196         0    0.7373
    0.6275         0    0.7216
    0.6353         0    0.7098
    0.6392         0    0.6941
    0.6431         0    0.6784
    0.6510         0    0.6627
    0.6549         0    0.6510
    0.6588         0    0.6353
    0.6667         0    0.6196
    0.6706         0    0.6039
    0.6784         0    0.5922
    0.6824         0    0.5765
    0.6863         0    0.5608
    0.6941         0    0.5490
    0.6980         0    0.5333
    0.7020         0    0.5176
    0.7098         0    0.5059
    0.7137         0    0.4902
    0.7216         0    0.4784
    0.7255         0    0.4627
    0.7294         0    0.4471
    0.7373         0    0.4353
    0.7412         0    0.4196
    0.7451         0    0.4039
    0.7529         0    0.3922
    0.7569         0    0.3765
    0.7647         0    0.3647
    0.7686    0.0039    0.3490
    0.7765    0.0118    0.3333
    0.7804    0.0196    0.3216
    0.7882    0.0275    0.3059
    0.7922    0.0314    0.2902
    0.8000    0.0392    0.2784
    0.8039    0.0471    0.2627
    0.8118    0.0549    0.2510
    0.8157    0.0627    0.2353
    0.8196    0.0745    0.2196
    0.8235    0.0824    0.2078
    0.8314    0.0941    0.1922
    0.8353    0.1059    0.1765
    0.8392    0.1137    0.1647
    0.8431    0.1255    0.1490
    0.8510    0.1373    0.1373
    0.8549    0.1451    0.1216
    0.8627    0.1569    0.1059
    0.8667    0.1686    0.0902
    0.8745    0.1804    0.0784
    0.8784    0.1882    0.0627
    0.8863    0.2000    0.0471
    0.8902    0.2118    0.0314
    0.8980    0.2235    0.0196
    0.9020    0.2314    0.0157
    0.9059    0.2431    0.0118
    0.9137    0.2549    0.0118
    0.9176    0.2667    0.0078
    0.9216    0.2745    0.0039
    0.9294    0.2863    0.0039
    0.9333    0.2980         0
    0.9412    0.3098         0
    0.9451    0.3176         0
    0.9529    0.3294         0
    0.9569    0.3412         0
    0.9647    0.3529         0
    0.9686    0.3608         0
    0.9765    0.3725         0
    0.9804    0.3843         0
    0.9882    0.3961         0
    0.9882    0.4039         0
    0.9882    0.4118         0
    0.9922    0.4196         0
    0.9922    0.4275         0
    0.9922    0.4353         0
    0.9961    0.4431         0
    0.9961    0.4510         0
    1.0000    0.4588         0
    1.0000    0.4667         0
    1.0000    0.4745         0
    1.0000    0.4824         0
    1.0000    0.4902         0
    1.0000    0.4980         0
    1.0000    0.5059         0
    1.0000    0.5137         0
    1.0000    0.5216         0
    1.0000    0.5255         0
    1.0000    0.5333         0
    1.0000    0.5412         0
    1.0000    0.5490         0
    1.0000    0.5529         0
    1.0000    0.5608         0
    1.0000    0.5686         0
    1.0000    0.5765         0
    1.0000    0.5804         0
    1.0000    0.5882         0
    1.0000    0.5961         0
    1.0000    0.6039         0
    1.0000    0.6078         0
    1.0000    0.6157         0
    1.0000    0.6235         0
    1.0000    0.6314         0
    1.0000    0.6353         0
    1.0000    0.6431         0
    1.0000    0.6510         0
    1.0000    0.6588         0
    1.0000    0.6627         0
    1.0000    0.6706         0
    1.0000    0.6784         0
    1.0000    0.6863         0
    1.0000    0.6902         0
    1.0000    0.6980         0
    1.0000    0.7059         0
    1.0000    0.7137         0
    1.0000    0.7216         0
    1.0000    0.7294         0
    1.0000    0.7373         0
    1.0000    0.7451         0
    1.0000    0.7490         0
    1.0000    0.7569         0
    1.0000    0.7647         0
    1.0000    0.7725         0
    1.0000    0.7804         0
    1.0000    0.7882         0
    1.0000    0.7961         0
    1.0000    0.8039         0
    1.0000    0.8078         0
    1.0000    0.8157         0
    1.0000    0.8235         0
    1.0000    0.8314         0
    1.0000    0.8353         0
    1.0000    0.8431         0
    1.0000    0.8510         0
    1.0000    0.8588         0
    1.0000    0.8627         0
    1.0000    0.8706         0
    1.0000    0.8784         0
    1.0000    0.8863         0
    1.0000    0.8941         0
    1.0000    0.9020         0
    1.0000    0.9098         0
    1.0000    0.9176         0
    1.0000    0.9216    0.0157
    1.0000    0.9294    0.0314
    1.0000    0.9373    0.0510
    1.0000    0.9451    0.0667
    1.0000    0.9490    0.0824
    1.0000    0.9569    0.1020
    1.0000    0.9647    0.1176
    1.0000    0.9725    0.1373
    1.0000    0.9725    0.1647
    1.0000    0.9765    0.1961
    1.0000    0.9804    0.2275
    1.0000    0.9843    0.2588
    1.0000    0.9882    0.2902
    1.0000    0.9922    0.3216
    1.0000    0.9961    0.3529
    1.0000    1.0000    0.3843
    1.0000    1.0000    0.4118
    1.0000    1.0000    0.4431
    1.0000    1.0000    0.4745
    1.0000    1.0000    0.5059
    1.0000    1.0000    0.5333
    1.0000    1.0000    0.5647
    1.0000    1.0000    0.5961
    1.0000    1.0000    0.6275
    1.0000    1.0000    0.6549
    1.0000    1.0000    0.6863
    1.0000    1.0000    0.7176
    1.0000    1.0000    0.7490
    1.0000    1.0000    0.7804
    1.0000    1.0000    0.8118
    1.0000    1.0000    0.8431
    1.0000    1.0000    0.8745
    1.0000    1.0000    0.8902
    1.0000    1.0000    0.9059
    1.0000    1.0000    0.9216
    1.0000    1.0000    0.9373
    1.0000    1.0000    0.9529
    1.0000    1.0000    0.9686
    1.0000    1.0000    0.9843
    1.0000    1.0000    1.0000
    1.0000    1.0000    1.0000
    1.0000    1.0000    1.0000
    1.0000    1.0000    1.0000
    1.0000    1.0000    1.0000
    1.0000    1.0000    1.0000
    1.0000    1.0000    1.0000
    1.0000    1.0000    1.0000
    ];
%interpolate values
xin=linspace(0,1,m)';
xorg=linspace(0,1,size(cmap_mat,1));
p(:,1)=interp1(xorg,cmap_mat(:,1),xin,'linear');
p(:,2)=interp1(xorg,cmap_mat(:,2),xin,'linear');
p(:,3)=interp1(xorg,cmap_mat(:,3),xin,'linear');
end

function [] = PlotIncidentF(ax,intensity,name,XTick,XTickLabel,YTick,YTickLabel,map,scale,XLabel,YLabel,Caxis)

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

function [Dat] = Gauss(Sc, Wid, ofst)
%Gauss Summary of this function goes here
%   Sc: coordimate info
%   Wid: waist radius for 1/e (amplitute); 1/e^2(intensity)
%   Wid(1) for x axis, Wid(2) for y axis
lw_x = Wid(1); lw_y = Wid(2);
dx = ofst(1); dy = ofst(2);
Dat = exp(-(Sc.x-dx).^2/lw_x^2-(Sc.y-dy).^2/lw_y^2);
end

function [  I_fout, E_fxout, E_fyout, E_fzout ] = PSFapp( c,polarization, ~, para1,~, para2, ~, para3, ~, para4, ~, para5, ~, para6, ~, para7)
%% Calculation PSF of input light field with complex amplitude c
%%  Form
% [intensityOut] = PSFaberration(c, polarization,...
%                               'laser',[lambda],...
%                               'obj',[NA,f,entrancePupil,W040],...
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
% 'obj',[NA,f,entrancePupil, W040] ,(1,4) two elements array
% 'materials',[n1,n2,d], (1,3) three elements array
% 'inXY', [pixelsIncidentPlane,rhoIncident],...
% 'outZ', [pixelsOutputZ,rhoZ],...
% 'outX', [pixelsOutputX,physicalStart,physicalEnd],...
% 'outY', [pixelsOutputY,physicalStart,physicalEnd]
%% Output
% I_fout, E_fxout, E_fyout, E_fzout (m,m,mz)


%% Initialize
lambda         = para1(1); 
NA             = para2(1); f      = para2(2); entrancePupil = para2(3);
W040           = para2(4);
n1             = para3(1); n2     = para3(2); d       = - para3(3)/n2;
rhoIncident    = para4(1); pixIn  = para4(2);
pixOutZ        = para5(1); rhoZ   = para5(2); zOut1   = para5(3); zOut2 = para5(4);
pixOutX				 = para6(1); xOut1  = para6(2); xOut2   = para6(3);
pixOutY				 = para7(1); yOut1  = para7(2); yOut2   = para7(3);

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
ppm = ParforProgressbar(pixOutZ);
if (W040 == 0) % No lens aberration
  if (n1==n2 || d==0) % No lens aberration + No interface aberration
    phi       	= atan2(YIN,XIN);
    theta     	= asin(sqrt(XIN.^2+YIN.^2)*NA/(n1*rPupil));
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
        ppm.increment();
        E_fx(:,:,zz)=(czt((czt(Eox.*exp(-1j*k1*cos(theta)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';
        E_fy(:,:,zz)=(czt((czt(Eoy.*exp(-1j*k1*cos(theta)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';
        E_fz(:,:,zz)=(czt((czt(Eoz.*exp(-1j*k1*cos(theta)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';
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
        ppm.increment();
        E_fx(:,:,zz)=(czt((czt(Eox.*exp(-1j*k2*cos(theta2)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';
        E_fy(:,:,zz)=(czt((czt(Eoy.*exp(-1j*k2*cos(theta2)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';
        E_fz(:,:,zz)=(czt((czt(Eoz.*exp(-1j*k2*cos(theta2)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';
    end
  end

elseif (W040 ~=0) % % Lens aberration
  if (n1==n2 || d==0) % % Lens aberrationd + No interface aberration 
    phi       	= atan2(YIN,XIN);
    theta     	= asin(sqrt(XIN.^2+YIN.^2)*NA/(n1*rPupil));
    r           = sqrt(XIN.^2+YIN.^2).*Aperture/rPupil;
    Eix         = Eix .* (exp(1j*2*pi*(W040 * r.^4))); % Spherical aberration
    Eiy         = Eiy .* (exp(1j*2*pi*(W040 * r.^4)));
    Eiz         = Eiz .* (exp(1j*2*pi*(W040 * r.^4)));
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
        ppm.increment();
        E_fx(:,:,zz)=(czt((czt(Eox.*exp(-1j*k1*cos(theta)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';
        E_fy(:,:,zz)=(czt((czt(Eoy.*exp(-1j*k1*cos(theta)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';
        E_fz(:,:,zz)=(czt((czt(Eoz.*exp(-1j*k1*cos(theta)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';
    end
  elseif (n1~=n2) % Interface aberration  considered
    phi       	= atan2(YIN,XIN);
    r           = sqrt(XIN.^2+YIN.^2).*Aperture/rPupil;
    Eix         = Eix .* (exp(1j*2*pi*(W040 * r.^4))); % Spherical aberration
    Eiy         = Eiy .* (exp(1j*2*pi*(W040 * r.^4)));
    Eiz         = Eiz .* (exp(1j*2*pi*(W040 * r.^4)));
%     Eix         = Eix .* (exp(1j*2*pi* zernike11 * sqrt(5) * (6 * r.^4 - 6 * r.^2 +1))+ ...
%                           exp(1j*2*pi* zernike22 * sqrt(7) * (20 * r.^6 - 30 * r.^4 + 12 * r.^2 -1)) + ...
%                           exp(1j*2*pi* zernike37 * sqrt(9) * (70 * r.^8 - 140 * r.^6 + 90 * r.^4 - 20 * r.^2 +1))); % Spherical aberration
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
        ppm.increment();
        E_fx(:,:,zz)=(czt((czt(Eox.*exp(-1j*k2*cos(theta2)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';
        E_fy(:,:,zz)=(czt((czt(Eoy.*exp(-1j*k2*cos(theta2)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';
        E_fz(:,:,zz)=(czt((czt(Eoz.*exp(-1j*k2*cos(theta2)*zOut(zz)),pixOutY,wy,ay).*factor_fy).',pixOutX,wx,ax).*factor_fx).';
    end
  end
end
delete(ppm); 
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
