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

% Dimensions of dirty beam
dirty_beam = zeros(res_l, res_m);
dirtyImage = zeros(res_l, res_m);

% Reshape baseline_vector to be (p*p, 3) for easy matrix multiplication
baseline_vector_reshaped = reshape(baseline_vector, p*p, 3);
RhReshaped = reshape(Rh, p*p, 1);

% Reshape direction_matrix to be (res_l,res_m, 3) for easy matrix multiplication
direction_matrix_reshaped = reshape(direction_matrix, res_l*res_m, 3);
direction_matrix_reshaped_transpose = direction_matrix_reshaped.';

% Calculate beamformer/array response for a certain direction.
function [dirty_beam, dirtyImage] = beam_former(centre, res_l, res_m, baseline_vector_reshaped, RhReshaped, direction_matrix_reshaped_transpose)

    %check if centre is viable
    if centre(1)^2 + centre(2)^2 > 1;
        disp("Centre is not a viable point")
        return;
    end

    % Dimensions of dirty beam
    dirty_beam = zeros(res_l, res_m);
    dirtyImage = zeros(res_l, res_m);

    % Define a vector for the target location
    N_centre = sqrt(1 - centre(1)^2 - centre(2)^2);
    L_centre = centre(1);
    M_centre = centre(2);
    P_q = transpose([L_centre,M_centre,N_centre]);

    divide_num = int16(res_l);
    [pos_l_list, pos_m_list] = ind2sub([res_l, res_m], 1:(res_l * res_m));
    for pos_l_m=1:(res_l*res_m)
        pos_l = pos_l_list(pos_l_m); % idivide(pos_l_m,divide_num)+1;
        pos_m = pos_m_list(pos_l_m); % pos_m = mod(pos_l_m,divide_num)+1;
        if direction_matrix_reshaped_transpose(3,pos_l_m) == 0 %l(pos_l)^2 + m(pos_m)^2 > 1
            dirty_beam(pos_l,pos_m) = 0;
            dirtyImage(pos_l,pos_m) = 0;
        else
            component = exp(1j .* baseline_vector_reshaped * (direction_matrix_reshaped_transpose(:,pos_l_m) - P_q));
            dirty_beam(pos_l,pos_m) = sum(component);
            dirtyImage(pos_l,pos_m) = sum(RhReshaped.*component);
        end
        if pos_l == res_l
            disp(pos_m)
        end
    end
end

[dirty_beam, dirtyImage] = beam_former([0,0], res_l, res_m, baseline_vector_reshaped, RhReshaped, direction_matrix_reshaped_transpose);

% Plot dirty beam
figure;
imagesc(abs(dirty_beam));
axis equal;        % Make axes equal for proper aspect ratio
colormap('jet');   % Use the 'jet' colormap for colors
colorbar;          % Display a color bar to the right
title('dirty_beam');
% Set the color axis limits (optional, for consistent color scaling)
%caxis([100, 5000]); % 5000 looks cool

% Plot dirty Image
figure;
imagesc(abs(dirtyImage));
axis equal;        % Make axes equal for proper aspect ratio
colormap('jet');   % Use the 'jet' colormap for colors
colorbar;
title('dirtyImage');
%caxis([100, 300]);
% save any matrixes in xlsx files
% filename = 'dirty_beam_513.xlsx';
% writematrix(dirty_beam,filename);
%% Normalize dirty_beam and dirtyImage
dirty_beam_normalized = (dirty_beam - min(dirty_beam(:)))/(max(dirty_beam(:)) - min(dirty_beam(:)));
dirtyImage_normalized = (dirtyImage - min(dirtyImage(:)))/(max(dirtyImage(:)) - min(dirtyImage(:)));
%% Clean Algorithm
close all
dirtyImageCleanA = dirtyImage_normalized;
dirtyBeamCleanA = dirty_beam_normalized;
q = 0;
gamma = 0.5;
sumOfVariances = sum(sum(abs(corrcoef(abs(dirtyImageCleanA)))))
varianceTreshold = sumOfVariances - 0.5*sumOfVariances;
P_q = int16.empty;
VarianceQ = double.empty;
while varianceTreshold <= sumOfVariances && q < 4
    q = q + 1;
    [maxPeak , Index] = max(abs(dirtyImageCleanA(:))); % Find the maximum of the dirtyImage
    if maxPeak < 1e-6  % If the maximum value is very small, break
            break;
    end

    [row, col] = ind2sub(size(dirtyImageCleanA),Index);
    P_q  = [P_q ; row , col];
    l = direction_matrix(P_q(end,1),P_q(end,2),1);
    m = direction_matrix(P_q(end,1),P_q(end,2),2);
    VarianceQ = [VarianceQ ; maxPeak/dirtyBeamCleanA(ceil(res_l/2),ceil(res_m/2))];
    
    [shifted_dirty_beam, shifted_dirty_image] = beam_former([l,m], res_l, res_m, baseline_vector_reshaped, RhReshaped, direction_matrix_reshaped_transpose);
    shifted_dirty_beam_norm = (shifted_dirty_beam - min(shifted_dirty_beam(:)))/(max(shifted_dirty_beam(:)) - min(shifted_dirty_beam(:)));
    dirtyImageCleanA = dirtyImageCleanA - (gamma*(VarianceQ(end))*shifted_dirty_beam_norm);
    sumOfVariances = sum(sum((corrcoef(abs(dirtyImageCleanA)))));

    % Plot dirty beam
    figure;
    imagesc(abs(shifted_dirty_beam_norm));
    axis equal;        % Make axes equal for proper aspect ratio
    colormap('jet');   % Use the 'jet' colormap for colors
    colorbar;          % Display a color bar to the right
    title('shifted_dirty_beam_norm');
    % Set the color axis limits (optional, for consistent color scaling)
    %caxis([100, 5000]); % 5000 looks cool
    
    % Plot dirty Image
    figure;
    imagesc(abs(dirtyImageCleanA));
    axis equal;        % Make axes equal for proper aspect ratio
    colormap('jet');   % Use the 'jet' colormap for colors
    colorbar;
    title('dirtyImageCleanA');
    %caxis([100, 300]);
end


function [beam] = Bsynth(image_size, beam_width, centre)
    % Generate and display synthetic beam (Gaussian bell curve)
    % Arguments:
    % image_size : Size of the output image (e.g., 100 for a 100x100 image)
    % beam_width : Standard deviation (spread) of the Gaussian bell curve
    % centre     : Coordinates of the centre of the bell curve
    
    % Create a grid of (l, m) coordinates
    l = linspace(-1+centre(1), 1+centre(1), image_size);  % l coordinates from -1 to 1
    m = linspace(-1+centre(2), 1+centre(2), image_size);  % m coordinates from -1 to 1
    [L, M] = meshgrid(l, m);          % 2D grid of (l, m) values
    
    % Calculate the distance from the center (0, 0)
    d = sqrt(L.^2 + M.^2);
    
    % Gaussian beam shape
    beam = exp(-d.^2 / (2 * beam_width^2));
end 

beam_synth = size(res_l,2);
bell_width = res_l/10;
gain = 300;
beam_synth = gain * Bsynth(res_l, bell_width, [0,0]);

% Plot synth beam
figure;
imagesc(beam_synth);
axis equal;        % Make axes equal for proper aspect ratio
colormap('jet');   % Use the 'jet' colormap for colors
colorbar;
title('example Bell shape in centre');
%caxis([100, 300]);

source_num = q;
sum_beam_synth = zeros(res_l,res_m);
for q=1:source_num
    l = direction_matrix(P_q(q,1),P_q(q,2),1);
    m = direction_matrix(P_q(q,1),P_q(q,2),2);
    sum_beam_synth = sum_beam_synth + gain*gamma*VarianceQ(q)*Bsynth(res_l, bell_width, [l,m]);
end
CleanAlgImage = dirtyImageCleanA + sum_beam_synth;

% Plot Clean Algorithm Image
figure;
imagesc(abs(dirtyImageCleanA));
axis equal;        % Make axes equal for proper aspect ratio
colormap('jet');   % Use the 'jet' colormap for colors
colorbar;  
title('dirtyImageCleanA');
%caxis([100, 300]);

% Plot Clean Algorithm Image
figure;
imagesc(abs(CleanAlgImage));
axis equal;        % Make axes equal for proper aspect ratio
colormap('jet');   % Use the 'jet' colormap for colors
colorbar;  
title('CleanAlgImage');
%caxis([100, 300]);