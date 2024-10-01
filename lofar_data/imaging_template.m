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
res_l = size(l,1);
res_m = size(m,1);

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

% Pre-allocate a 3D matrix for the direction vectors (max resolution: 513x513x3)
direction_matrix = zeros(res_l, res_m, 3);

% Calculate n for each element where l^2 + m^2 <= 1
valid_points = L.^2 + M.^2 <= 1;  % Logical mask for valid points within the unit circle

% Compute n only for valid points, otherwise leave as zero
N = zeros(res_l);  % Initialize n to zero
N(valid_points) = sqrt(1 - L(valid_points).^2 - M(valid_points).^2);  % Calculate n

% Store the l, m, and n components in the direction_matrix
direction_matrix(:, :, 1) = L;  % l-coordinates
direction_matrix(:, :, 2) = M;  % m-coordinates
direction_matrix(:, :, 3) = N;  % n-coordinates

% Display the resulting direction_matrix
disp('Direction vector matrix created with size:');
disp(size(direction_matrix));

% Demensions of dirty beam
dirty_beam = zeros(res_l, res_m);
dirtyImage = zeros(res_l, res_m);

% Reshape baseline_vector to be (p*p, 3) for easy matrix multiplication
baseline_vector_reshaped = reshape(baseline_vector, p*p, 3);
RhReshaped = reshape(Rh, p*p, 1);

% Reshape direction_matrix to be (res_l,res_m, 3) for easy matrix multiplication
direction_matrix_reshaped = reshape(direction_matrix, res_l*res_m, 3);
direction_matrix_reshaped_transpose = direction_matrix_reshaped.';

divide_num = int16(res_l);
[pos_l_list, pos_m_list] = ind2sub([res_l, res_m], 1:(res_l * res_m));
for pos_l_m=1:(res_l*res_m)
    pos_l = pos_l_list(pos_l_m); % idivide(pos_l_m,divide_num)+1;
    pos_m = pos_m_list(pos_l_m); % pos_m = mod(pos_l_m,divide_num)+1;
    if l(pos_l)^2 + m(pos_m)^2 > 1
        dirty_beam(pos_l,pos_m) = 0;
        dirtyImage(pos_l,pos_m) = 0;
    else
        component = exp(1j .* baseline_vector_reshaped * direction_matrix_reshaped_transpose(:,pos_l_m));
        dirty_beam(pos_l,pos_m) = sum(component);
        dirtyImage(pos_l,pos_m) = sum(RhReshaped.*component);
    end
    % if pos_m == res_l
    %     disp(pos_l_m)
    %     disp(pos_l)
    %     disp(pos_m)
    % end
end

% Plot dirty beam
figure;
imagesc(abs(dirty_beam));
axis equal;        % Make axes equal for proper aspect ratio
colormap('jet');   % Use the 'jet' colormap for colors
colorbar;          % Display a color bar to the right
% Set the color axis limits (optional, for consistent color scaling)
caxis([100, 5000]); % 5000 looks cool

%% Clean Algorithm
dirtyImageCleanA = dirtyImage;
q = 0;
gamma = 0.01;
varianceTreshold = 2000;
sumOfVariances = sum(sum(abs(corrcoef(dirtyImageCleanA))))
P_q = uint8.empty;
VarianceQ = double.empty;
while varianceTreshold <= sumOfVariances
    q = q+1;
    [maxPeak , I] = max((dirtyImageCleanA),[] ,'all');
    [Pq, Pp] = ind2sub(size(dirtyImageCleanA),I);
    P_q  = [P_q ; Pq , Pp];
    VarianceQ = [VarianceQ ; maxPeak/dirty_beam(floor(pos_l/2),floor(pos_m/2))];
    %dirtyImageCleanA = dirtyImageCleanA - gamma*(VarianceQ(end))*minus(dirty_beam,dirty_beam(P_q(end)));
    dirtyImageCleanA = dirtyImageCleanA - (gamma*(VarianceQ(end))*shift_center(dirty_beam,P_q(end,:)));
    sumOfVariances = sum(sum((corrcoef(dirtyImageCleanA))));
end    

function new_matrix = shift_center(matrix, new_center)
    % Get the dimensions of the input matrix
    [m, n] = size(matrix);
    
    % Calculate the original center
    original_center = int16([ceil(m / 2), ceil(n / 2)]);
    
    % Calculate the shifts
    row_shift = (int16(new_center(1)) - original_center(1));
    col_shift = (int16(new_center(2)) - original_center(2));
    
    % Create the new matrix initialized with zeros
    new_matrix = zeros(size(matrix));
    
    % Create row and column indices for the original matrix
    [rows, cols] = ndgrid(1:m, 1:n);
    
    % Calculate new positions
    new_rows = int16(rows) + row_shift;
    new_cols = int16(cols) + col_shift;
    
    % Find valid positions within the matrix bounds
    valid_mask = (new_rows >= 1 & new_rows <= m) & (new_cols >= 1 & new_cols <= n);
    
    % Assign values to the new matrix at valid positions
    new_matrix(new_rows(valid_mask), new_cols(valid_mask)) = matrix(rows(valid_mask), cols(valid_mask));
end

function [beam] = Bsynth(image_size, beam_width)
    % Generate and display synthetic beam (Gaussian bell curve)
    % Arguments:
    % image_size : Size of the output image (e.g., 100 for a 100x100 image)
    % beam_width : Standard deviation (spread) of the Gaussian bell curve
    
    % Create a grid of (l, m) coordinates
    l = linspace(-1, 1, image_size);  % l coordinates from -1 to 1
    m = linspace(-1, 1, image_size);  % m coordinates from -1 to 1
    [L, M] = meshgrid(l, m);          % 2D grid of (l, m) values
    
    % Calculate the distance from the center (0, 0)
    d = sqrt(L.^2 + M.^2);
    
    % Gaussian beam shape
    beam = exp(-d.^2 / (2 * beam_width^2));
end

beam_synth = size(size(l,1),2);
disp(beam_synth);
beam_synth = Bsynth(size(l,1), 2/size(l,1)*8);

source_num = q;
sum_beam_synth = zeros(res_l,res_m)
for q=1:source_num
    %sum_beam_synth = sum_beam_synth + gamma*(VarianceQ(q))*minus(beam_synth, beam_synth(P_q(q)));
    sum_beam_synth = sum_beam_synth + gamma*(VarianceQ(q))*shift_center(beam_synth,P_q(q,:));
end    
CleanAlgImage = dirtyImageCleanA + sum_beam_synth;
%CleanAlgImage = dirtyImageCleanA + sum(gamma*VarianceQ*minus(beam_synth, beam_synth(P_q)));

close all

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

% Plot Clean Algorithm Image
figure;
imagesc(abs(CleanAlgImage));
axis equal;        % Make axes equal for proper aspect ratio
colormap('jet');   % Use the 'jet' colormap for colors
colorbar;  
caxis([100, 300]);

% Plot Clean Algorithm Image
figure;
imagesc(abs(dirtyImageCleanA));
axis equal;        % Make axes equal for proper aspect ratio
colormap('jet');   % Use the 'jet' colormap for colors
colorbar;  
%caxis([100, 300]);