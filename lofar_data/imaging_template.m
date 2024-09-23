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
FracAngle = 0.25;% fraction of theoretical limit for angle resolution
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

% Calculate the baseline vector for x, y and z
baseline_x = (meshgrid(poslocal(:,1)) - meshgrid(poslocal(:,1)).');
baseline_y = (meshgrid(poslocal(:,2)) - meshgrid(poslocal(:,2)).');
baseline_z = (meshgrid(poslocal(:,3)) - meshgrid(poslocal(:,3)).');
baseline_vector = ones(p,p, 3);
baseline_vector(:,:,1) = baseline_x;
baseline_vector(:,:,2) = baseline_y;
baseline_vector(:,:,3) = baseline_z;

% define the direction matrix
% Create 2D grids of l and m
[L, M] = meshgrid(l, m);

% Pre-allocate a 3D matrix for the direction vectors (513x513x3)
direction_matrix = zeros(size(l,1), size(l,1), 3);

% Calculate n for each element where l^2 + m^2 <= 1
valid_points = L.^2 + M.^2 <= 1;  % Logical mask for valid points within the unit circle

% Compute n only for valid points, otherwise leave as zero
N = zeros(size(L));  % Initialize n to zero
N(valid_points) = sqrt(1 - L(valid_points).^2 - M(valid_points).^2);  % Calculate n

% Store the l, m, and n components in the direction_matrix
direction_matrix(:, :, 1) = L;  % l-coordinates
direction_matrix(:, :, 2) = M;  % m-coordinates
direction_matrix(:, :, 3) = N;  % n-coordinates

% Display the resulting direction_matrix
disp('Direction vector matrix created with size:');
disp(size(direction_matrix));

% Matrix exponential method
% Sum all the Multiplies with on element of the baseline vector and one
% element of the direction vector, do this for all direction vectors.
% Demensions of dirty beam
dirty_beam = zeros(size(l,1), size(m,1));
dirtyImage = zeros(size(l,1), size(m,1));
% Reshape baseline_vector to be (p*p, 3) for easy matrix multiplication
baseline_vector_reshaped = reshape(baseline_vector, p*p, 3);
RhReshaped = reshape(Rh, p*p, 1);
% Reshape direction_matrix to be (513*513, 3) for easy matrix multiplication
%direction_matrix_reshaped = reshape(direction_matrix, 513*513, 3);
%direction_matrix_reshaped_transpose = direction_matrix_reshaped.';
%display(size(direction_matrix_reshaped_transpose))
% Compute the matrix exponential in a vectorized manner
% Compute the dot product between each baseline and each direction vector
% We use the transpose of direction_matrix_reshaped to multiply with baseline_vector_reshaped
% Resulting matrix size will be (p*p, 513*513)

%components = sum(baseline_vector_reshaped * direction_matrix_reshaped_transpose, 2);

% Now exponentiate the result (element-wise) and sum over all the baselines
% First, exponentiate with 1j
%exp_components = exp(1j * components);

% Sum over all baseline components (first dimension) and reshape back to 513x513
%dirty_beam = reshape(sum(exp_components, 1), [513, 513]);

% Calculate the dirty beam using a matrix exponential (permute?)
% dirty_beam = exp(1j*baseline_vector*permute(direction_matrix,[3 1 2])),[1 2])

% % Calculate the dirty_beam per direction vector in for loops (slow)
% dirty_beam = ones(size(l,1), size(m,1))

for pos_l=1:size(l,1)
    for pos_m=1:size(m,1)

        n = sqrt(1 - l(pos_l)^2 - m(pos_m).^2);

        direction_vector = [l(pos_l),m(pos_m),n];

        component = baseline_vector_reshaped*transpose(direction_vector);

        dirty_beam(pos_l,pos_m) = sum(exp(1j*component),1);
        dirtyImage(pos_l,pos_m) = sum(RhReshaped.*exp(1j*component),1);

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


% Plot dirty Image
figure;
imagesc(abs(dirtyImage));
axis equal;        % Make axes equal for proper aspect ratio
colormap('jet');   % Use the 'jet' colormap for colors
colorbar;  
caxis([100, 300]);
% save any matrixes in xlsx files
% filename = 'dirty_beam_513.xlsx';
% writematrix(dirty_beam,filename);