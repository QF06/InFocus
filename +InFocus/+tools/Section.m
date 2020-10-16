function [section] = Section(~, para1, ~, para2, ~, para3, ~, para4)
%%	Section initiate the coordinate of the simulation section
%% Form
% section = Section('pixelWH',[pixelW,pixelH],...
%               'physicalWH',[physicalW,physicalH],...
%               'vector', [0,0,1],...
%               'offset', [0,0,0]);
%% Description
%  Initiate the coordinate of the simulation section, which includes
%  defining the simulation plane, plane dimensions (pixel & physical)and
%  offset
%% Input
% para1 (1,2) pixel dimension
% para2 (1,2) physical dimension
% para3 (1,3) normal vector of the simulated plane 
% para4 (1,3) offset
%% Output
% Sc (.) Structure data
% Sc.x Sc.y Sc.z are the coordinate of the simulation section.
% Sc.Ex Sc.Ey are the initialized electromagnetic field distribution.

%% Initialize
section.pix_w  = para1(1);
section.pix_h  = para1(2);
section.wid    = para2(1);
section.hei    = para2(2);
section.vector = para3;
section.offset = para4;

%% Calculation
section.dx     = section.wid / (section.pix_w-1); % sampling rate [m/pixel]
section.dy     = section.hei / (section.pix_h-1);
section.ds     = section.dx * section.dy;
section.x      = zeros(section.pix_w, section.pix_h); % pixel space
section.y      = zeros(section.pix_w, section.pix_h);
section.z      = zeros(section.pix_w, section.pix_h);
if (section.vector == [0,0,1])
    [section.x,section.y] = meshgrid(-section.wid/2:section.dx:section.wid/2,...
                            -section.hei/2:section.dy:section.hei/2); %coordinate in real space [m]
elseif (section.vector==[0,1,0])
    [section.x,section.z] = meshgrid(-section.wid/2:section.dx:section.wid/2,...
                            -section.hei/2:section.dy:section.hei/2);
    section.y             = zeros(section.pix_w, section.pix_h);
elseif (section.vector == [1,0,0])
    [section.y,section.z] = meshgrid(-section.wid/2:section.dx:section.wid/2,...
                            -section.hei/2:section.dy:section.hei/2);
    section.x             = zeros(section.pix_w, section.pix_h);
end
section.x                 = section.x + section.offset(1);
section.y                 = section.y + section.offset(2);
section.z                 = section.z + section.offset(3);
section.E_x               = zeros(section.pix_w, section.pix_h);
section.E_y               = zeros(section.pix_w, section.pix_h);
end

