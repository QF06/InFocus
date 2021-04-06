function [ laser ] = Laser(lambda0)
%%	Laser initiate the data of the laser with all necessary informations
%%  Form
% laser = Laser(lambda0);
%% Description
%  All the laser parameters are saved in here built laser structure.
%% Input
%   lambda0 is the central wavelength of the laser.   
%% Output
% laser (.) Structure data, the initial data contain the basic info of the
% laser and it is easy to extend the data structure.

%% Initialize  
laser.lambda0          = lambda0;    
laser.c                = 299792458;% c
laser.k0               = 2*pi/laser.lambda0;  % wavenumber for the centrol wavelength
laser.w0               = 2*pi*laser.c/laser.lambda0;  % frequence of the centrol wavelength
end

