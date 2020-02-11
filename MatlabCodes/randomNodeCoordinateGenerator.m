function [X,Y,angle]=randomNodeCoordinateGenerator(disk_center,disk_radius,noOfPoints,angleLim1_in,angleLim2_in)
% Author: KB
% Purpose: For generating random points within a certain angle range around
% a circle. E.g -30(-pi/6) to 30(pi/6) degrees, -45(-pi/2) to 45(pi/2) degrees etc.

%%---input description-----
%disk_center: for the center of the disk
%disk_radius: for required radius of the disk
%noOfPoints: for the number of points to generate within the disk
%angleLim1_in: the starting angle of the sector of interest within the
%disk(in radians)
%angleLim2_in: the stopping angle of the sector of interest within the
%disk(in radians)
%%--------------------------

%%---output description-----
%X: for x-coordinates
%Y: for y-coordinates
%angle: angles in radians
%%--------------------------

% Reference paper:Network Beamforming Using Relays With Perfect Channel
% Information.



n=noOfPoints; % number of relays points that you want (i.e Relay positions)
center = [1 ,2]; % center coordinates of the circle [x0,y0]
radius = 10; % radius of the circle

center=disk_center;
radius=disk_radius;

%Between angleLim1 and angleLim2 (using the understanding of Affine
%functions)

%angleLim1=pi/3;%first angular limit. E.g pi/3
angleLim1=angleLim1_in;%first angular limit. E.g pi/3
angleLim1_ratio=rand(n,1);%first angular limit's ratio

%angleLim2=-pi/3;%second angular limit. Eg -pi/3
angleLim2=angleLim2_in;%second angular limit. Eg -pi/3
angleLim2_ratio=1-angleLim1_ratio;%sedond angular limit's ratio

%Affine fuction definition
angle = (angleLim1*angleLim1_ratio)+(angleLim2*angleLim2_ratio);%angle in degrees

r = radius*sqrt(rand(n,1));
X = r.*cos(angle)- center(1);
Y = r.*sin(angle)- center(2);

angle=atan(Y./X);
%[X,Y,angle];%first column will be the x values, second column the y values and 3 column the angles
end