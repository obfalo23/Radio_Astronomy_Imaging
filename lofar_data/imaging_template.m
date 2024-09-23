% Template script for radio astronomy image formation
% Millad Sardarabadi, Sept 2015
%
clear all
close all
clc
%% Loading Data
% The data from LOFAR core will be loaded in this cell.
%
% Rh: A p x p sample covariance matrix, where p = number of antennas
% poslocal: A p x 3 matrix containing the earth-bound coordinates of the
%           receiving elements.
% freq: Central frequency for the measured channel.
% t_obs_matlab: Observation time.
%               (datestr(t_obs_matlab): 11-Jul-2012 22:22:20) 
% lon: Longitude of the LOFAR core
% lat: Latitude of the LOFAR core

load lofar_DSP_data_1.mat

figure
scatter(poslocal(:,1),poslocal(:,2))

%% Initialization
c = 299792458; % speed of light
p = size(Rh,1);% number of receiving elements
lambda = c / freq; % wavelength at the observation frequency

% define a distance matrix of all baselines (=vectors between pairs of antennas)
dist = sqrt(  (meshgrid(poslocal(:,1)) - meshgrid(poslocal(:,1)).').^2 ...
            + (meshgrid(poslocal(:,2)) - meshgrid(poslocal(:,2)).').^2);

D = max(dist(:)); % Maximum Baseline

%% Coordinate system
% In astronomy it is common to use the (l,m,n) coordinates which 
% are the components of the more familiar direction vector in spherical 
% coordinates:
%    s = [cos(theta) * cos(phi)
%        cos(theta) * sin(phi)
%        sin(theta)]


% define a grid of image coordinates (l,m) at the right resolution
FracAngle = 4;% fraction of theoretical limit for angle resolution
                 % 0.25 means ~4x4 pixels in the main lobe
dl = FracAngle * lambda / D; % width/height of a pixel
l = 0:dl:1; % First image coordinate := cos(theta)cos(phi)
l = [-fliplr(l(2:end)),l]';
m = l;  % Second image coordinate := cos(theta)sin(phi)
        % m = l means an aspect ratio 1:1

%% Your Imaging Algorithm should come hereafter...
% Calculate dirty beam pattern

% Define normalized location vector
z = (pi*freq/c) * poslocal;


% Calculate the dirty_beam per direction vector in for loops (slow)
dirty_beam = ones(size(l,1), size(m,1))
sum = 0
for pos_l=1:size(l,1)
    for pos_m=1:size(m,1)
        n = sqrt(1 - l(pos_l)^2 - m(pos_m).^2);
        direction_vector = [l(pos_l),m(pos_m),n];
        for x = 1:p
            for y = 1:p
                baseline = (z(x,:)-z(y,:));
                component = baseline*transpose(direction_vector);
                sum = sum +  exp(1j*component);
            end
        end
        dirty_beam(pos_l,pos_m) = sum;
        sum = 0;
    end
    disp(pos_l)
end

% Plot dirty beam
figure;
imagesc(abs(dirty_beam));
axis equal;        % Make axes equal for proper aspect ratio
colormap('jet');   % Use the 'jet' colormap for colors
colorbar;          % Display a color bar to the right
% Set the color axis limits (optional, for consistent color scaling)
caxis([0, 5000]); % 5000 looks cool

