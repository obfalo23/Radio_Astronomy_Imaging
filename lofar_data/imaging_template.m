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
FracAngle = 1 % 0.0625;% fraction of theoretical limit for angle resolution
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
z = poslocal; %(pi*freq/c) * poslocal;

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

% Reshape baseline_vector to be (p*p, 3) for easy matrix multiplication
baseline_vector_reshaped = reshape(baseline_vector, p*p, 3);
RhReshaped = reshape(Rh, p*p, 1);

% Reshape direction_matrix to be (res_l,res_m, 3) for easy matrix multiplication
direction_matrix_reshaped = reshape(direction_matrix, res_l*res_m, 3);
direction_matrix_reshaped_transpose = direction_matrix_reshaped.';

% Calculate beam former function for a certain direction.
function [dirty_beam_internal, beam_former_image, MVDR_image_internal, AAR_image_internal] = beam_former(centre, res_l, res_m, normalized_location_vector, Rh, direction_matrix_reshaped, make_MVDR)

    %check if centre is viable
    if centre(1)^2 + centre(2)^2 > 1;
        disp("Centre is not a viable point")
        return;
    end

    % predefine arrays
    dirty_beam_internal = zeros(res_l,res_m);
    beam_former_image = zeros(res_l,res_m);
    MVDR_image_internal = zeros(res_l,res_m);
    AAR_image_internal = zeros(res_l,res_m);

    % Pre invert the covariance matrix
    Rh_inv = inv(Rh);

    direction_matrix_reshaped_transpose = direction_matrix_reshaped.';

    % Define a vector for the target location
    N_centre = sqrt(1 - centre(1)^2 - centre(2)^2);
    L_centre = centre(1);
    M_centre = centre(2);
    P_q = [L_centre,M_centre,N_centre];
    P_q_t = P_q.';

    [pos_l_list, pos_m_list] = ind2sub([res_l, res_m], 1:(res_l * res_m));
    for pos_l_m=1:(res_l*res_m)
        pos_l = pos_l_list(pos_l_m);
        pos_m = pos_m_list(pos_l_m);
        if direction_matrix_reshaped_transpose(3,pos_l_m) == 0 %l(pos_l)^2 + m(pos_m)^2 > 1
            a_p(pos_l,pos_m) = 0;
            a_p_h(pos_l,pos_m) = 0;
        else
            component = exp(1j * normalized_location_vector * (direction_matrix_reshaped_transpose(:,pos_l_m) - P_q_t));
            component_H = ctranspose(component);
            
            beam_former_image(pos_l,pos_m) = component_H*Rh*component;
            if make_MVDR == 1
                MVDR_component = (component_H*Rh_inv*component);
                MVDR_image_internal(pos_l,pos_m) = 1/MVDR_component;
                AAR_image_internal(pos_l,pos_m) = MVDR_component/((component_H*(Rh_inv^2)*component)^2);
            end
            dirty_beam_internal(pos_l,pos_m) = sum(component);
        end
        if pos_l == res_l
            disp(pos_m)
        end
    end
end
make_MVDR = 0;
[dirty_beam, dirty_image, MVDR_image, AAR_image] = beam_former([0,0], res_l, res_m, z, Rh, direction_matrix_reshaped, make_MVDR);

%% 
close all;

abs_dirty_beam = abs(dirty_beam);
abs_dirty_image = abs(dirty_image);
abs_MVDR_image = abs(MVDR_image);
abs_AAR_image = abs(AAR_image);

% Plot dirty beam
fig1 = figure;
imagesc(abs_dirty_beam);
axis equal;        % Make axes equal for proper aspect ratio
colormap('jet');   % Use the 'jet' colormap for colors
colorbar;
title('dirty beam');
max_color = max(abs_dirty_beam,[],"all");
abs_dirty_beam(abs_dirty_beam == 0) = inf;
min_color = min(abs_dirty_beam,[],"all");
caxis([min_color, max_color]);

% Plot dirty Image
fig2 = figure;
imagesc(abs_dirty_image);
axis equal;        % Make axes equal for proper aspect ratio
colormap('jet');   % Use the 'jet' colormap for colors
colorbar;
title('beam former image (dirty)');
max_color = max(abs_dirty_image,[],"all");
abs_dirty_image(abs_dirty_image == 0) = inf;
min_color = min(abs_dirty_image,[],"all");
caxis([min_color, max_color]);

if make_MVDR == 1
    % Plot dirty Image using MVDR
    fig3 = figure;
    imagesc(abs_MVDR_image);
    axis equal;        % Make axes equal for proper aspect ratio
    colormap('jet');   % Use the 'jet' colormap for colors
    colorbar;
    title('MVDR image');
    max_color = max(abs_MVDR_image,[],"all");
    abs_MVDR_image(abs_MVDR_image == 0) = inf;
    min_color = min(abs_MVDR_image,[],"all");
    caxis([min_color, max_color]);
    
    % Plot dirty Image using MVDR
    fig4 = figure;
    imagesc(abs_AAR_image);
    axis equal;        % Make axes equal for proper aspect ratio
    colormap('jet');   % Use the 'jet' colormap for colors
    colorbar;
    title('AAR image');
    max_color = max(abs_AAR_image,[],"all");
    abs_AAR_image(abs_AAR_image == 0) = inf;
    min_color = min(abs_AAR_image,[],"all");
    caxis([min_color, max_color]);
end

% Align all figures, you need the image processing toolbox for this!
iptwindowalign(fig1,"right",fig2,"left");
if make_MVDR == 1
    iptwindowalign(fig1,"bottom",fig3,"top");
    iptwindowalign(fig3,"right",fig4,"left");
end

% Celebrate
disp("yay done with this")

%% Normalize dirty_beam and dirty_image
close all;
function normalized_image = normalize_image(image)
    abs_image = abs(image);
    max_abs_image = max(abs_image,[],"all");
    abs_image(abs_image == 0) = inf;
    min_abs_image = min(abs_image,[],"all");
    abs_image(abs_image == inf) = min_abs_image;
    % Take absolute of original image again, so that the outside POV is 0
    % again.
    normalized_image = (abs_image - min_abs_image)/(max_abs_image - min_abs_image);
end

dirty_beam_normalized = normalize_image(dirty_beam);
dirty_image_normalized = normalize_image(dirty_image);

if make_MVDR == 1
    MVDR_image_normalized = normalize_image(MVDR_image);
    AAR_image_normalized = normalize_image(AAR_image);
end

% Plot dirty beam
fig1 = figure;
imagesc(dirty_beam_normalized);
axis equal;        % Make axes equal for proper aspect ratio
colormap('jet');   % Use the 'jet' colormap for colors
colorbar;
title('dirty beam normalized');

% Plot dirty Image
fig2 = figure;
imagesc(dirty_image_normalized);
axis equal;        % Make axes equal for proper aspect ratio
colormap('jet');   % Use the 'jet' colormap for colors
colorbar;
title('beam former image (dirty)');

if make_MVDR == 1
    % Plot dirty Image using MVDR
    fig3 = figure;
    imagesc(MVDR_image_normalized);
    axis equal;        % Make axes equal for proper aspect ratio
    colormap('jet');   % Use the 'jet' colormap for colors
    colorbar;
    title('MVDR image');
    
    % Plot dirty Image using MVDR
    fig4 = figure;
    imagesc(AAR_image_normalized);
    axis equal;        % Make axes equal for proper aspect ratio
    colormap('jet');   % Use the 'jet' colormap for colors
    colorbar;
    title('AAR image');
end

% Align all figures, you need the image processing toolbox for this!
iptwindowalign(fig1,"right",fig2,"left");
if make_MVDR == 1
    iptwindowalign(fig1,"bottom",fig3,"top");
    iptwindowalign(fig3,"right",fig4,"left");
end

disp("Yaay also done with this")

%% Clean Algorithm
close all
% If you want to use the MVDR image replace by MVDR_image_normalized
% If you want to use the Dirty image replace by dirty_image_normalized
if make_MVDR == 1
    dirty_image_CLEAN = MVDR_image_normalized;
else
    dirty_image_CLEAN = dirty_image_normalized;
end

dirty_beam_CLEAN = dirty_beam_normalized;
abs_DI_CLEAN = abs(dirty_image_CLEAN);

% Initializing Clean algorithm
q = 0;
gamma = 0.5;
% sumOfVariances = sum(sum(abs(corrcoef(abs_DI_CLEAN))))
sumOfVariances = norm(abs_DI_CLEAN)
varianceTreshold = sumOfVariances - 0.9*sumOfVariances
P_q = int16.empty;
VarianceQ = double.empty;
while varianceTreshold <= sumOfVariances && q < 5
    q = q + 1
    [maxPeak , Index] = max(dirty_image_CLEAN,[],"all"); % Find the maximum of the dirty_image
    disp(Index)
    if maxPeak < 0.25  % If the maximum value is very small, break
        break;
    end

    [row, col] = ind2sub(size(dirty_image_CLEAN),Index);
    P_q  = [P_q ; row , col];
    l = direction_matrix(P_q(end,1),P_q(end,2),1);
    m = direction_matrix(P_q(end,1),P_q(end,2),2);
    VarianceQ = [VarianceQ ; maxPeak/dirty_beam_CLEAN(ceil(res_l/2),ceil(res_m/2))];
    
    [shifted_dirty_beam, shifted_dirty_image] = beam_former([l,m], res_l, res_m, z, Rh, direction_matrix_reshaped, make_MVDR);
    shifted_dirty_beam_norm = normalize_image(shifted_dirty_beam);
    dirty_image_CLEAN = dirty_image_CLEAN - (gamma*(VarianceQ(end))*shifted_dirty_beam_norm);
    abs_DI_CLEAN = abs(dirty_image_CLEAN);
    sumOfVariances = norm(abs_DI_CLEAN)

    % Plot dirty beam
    fig1 = figure;
    abs_shifted_DB_norm = abs(shifted_dirty_beam_norm);
    imagesc(abs_shifted_DB_norm);
    axis equal;        % Make axes equal for proper aspect ratio
    colormap('jet');   % Use the 'jet' colormap for colors
    colorbar;          % Display a color bar to the right
    title('shifted dirty beam norm');
    max_color = max(abs_shifted_DB_norm,[],"all");
    abs_shifted_DB_norm(abs_shifted_DB_norm == 0) = inf;
    min_color = min(abs_shifted_DB_norm,[],"all");
    caxis([min_color, max_color]);
    
    % Plot dirty Image
    fig2 = figure;
    imagesc(abs_DI_CLEAN);
    axis equal;        % Make axes equal for proper aspect ratio
    colormap('jet');   % Use the 'jet' colormap for colors
    colorbar;
    title('Absolute dirty image clean');
    max_color = max(abs_DI_CLEAN,[],"all");
    abs_DI_CLEAN(abs_DI_CLEAN == 0) = inf;
    min_color = min(abs_DI_CLEAN,[],"all");
    caxis([min_color, max_color]);
    
    iptwindowalign(fig1,"right",fig2,"left");
end
source_num = q

%% Synthetic beam forming
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

close all;

% Plot dirty Image
fig2 = figure;
imagesc(abs_DI_CLEAN);
axis equal;        % Make axes equal for proper aspect ratio
colormap('jet');   % Use the 'jet' colormap for colors
colorbar;
title('Absolute dirty image clean');
max_color = max(abs_DI_CLEAN,[],"all");
abs_DI_CLEAN(abs_DI_CLEAN == 0) = inf;
min_color = min(abs_DI_CLEAN,[],"all");
caxis([min_color, max_color]);

bell_width = 0.04;
gain = 0.5;
beam_synth = gain * Bsynth(res_l, bell_width, [0,0]);

% Plot synth beam
figure;
imagesc(beam_synth);
axis equal;        % Make axes equal for proper aspect ratio
colormap('jet');   % Use the 'jet' colormap for colors
colorbar;
title('example Bell shape in centre');
% caxis([0, 1]);
disp(direction_matrix(64,64,1))
disp(direction_matrix(64,64,2))
sum_beam_synth = zeros(res_l,res_m);
for q=1:source_num - 1
    disp(P_q(q,1))
    disp(P_q(q,2))
    l = direction_matrix(res_l-P_q(q,2),res_l-P_q(q,1),2)
    m = direction_matrix(res_l-P_q(q,2),res_l-P_q(q,1),1)
    steered_synth_beam = gain*gamma*VarianceQ(q)*Bsynth(res_l, bell_width, [l,m]);
    sum_beam_synth = sum_beam_synth + steered_synth_beam;
end
CleanAlgImage = dirty_image_CLEAN + sum_beam_synth;

%% Plot results
close all;
% Plot all 4 or six images
abs_DB = abs(dirty_beam);
abs_DI = abs(dirty_image);
abs_DI_CLEAN = abs(dirty_image_CLEAN);
abs_CLEAN = abs(CleanAlgImage);

if make_MVDR == 1
    abs_MVDR = abs(MVDR_image_normalized);
    abs_AAR = abs(AAR_image_normalized);
end

fig1 = figure;
imagesc(abs_DB);
axis equal;        % Make axes equal for proper aspect ratio
colormap('jet');   % Use the 'jet' colormap for colors
colorbar;  
title('dirty beam');
max_color = max(abs_DB,[],"all");
abs_DB(abs_DB == 0) = inf;
min_color = min(abs_DB,[],"all");
caxis([min_color, max_color]);

fig2 = figure;
imagesc(abs_DI);
axis equal;        % Make axes equal for proper aspect ratio
colormap('jet');   % Use the 'jet' colormap for colors
colorbar;  
title('dirty Image');
max_color = max(abs_DI,[],"all");
abs_DI(abs_DI == 0) = inf;
min_color = min(abs_DI,[],"all");
caxis([min_color, max_color]);

fig3 = figure;
imagesc(abs_DI_CLEAN);
axis equal;        % Make axes equal for proper aspect ratio
colormap('jet');   % Use the 'jet' colormap for colors
colorbar;  
title('Clean Image without beamsynth');
max_color = max(abs_DI_CLEAN,[],"all");
abs_DI_CLEAN(abs_DI_CLEAN == 0) = inf;
min_color = min(abs_DI_CLEAN,[],"all");
caxis([min_color, max_color]);

fig4 = figure;
imagesc(abs_CLEAN);
axis equal;        % Make axes equal for proper aspect ratio
colormap('jet');   % Use the 'jet' colormap for colors
colorbar;  
title('Clean Image');
max_color = max(abs_CLEAN,[],"all");
abs_CLEAN(abs_CLEAN == 0) = inf;
min_color = min(abs_CLEAN,[],"all");
caxis([min_color, max_color]);

if make_MVDR == 1
    fig5 = figure;
    imagesc(abs_MVDR);
    axis equal;        % Make axes equal for proper aspect ratio
    colormap('jet');   % Use the 'jet' colormap for colors
    colorbar;  
    title('MVDR Image');
    max_color = max(abs_MVDR,[],"all");
    abs_MVDR(abs_MVDR == 0) = inf;
    min_color = min(abs_MVDR,[],"all");
    caxis([min_color, max_color]);
    
    fig6 = figure;
    imagesc(abs_AAR); % NEEDS TO BE CHANGED TO ARR
    axis equal;        % Make axes equal for proper aspect ratio
    colormap('jet');   % Use the 'jet' colormap for colors
    colorbar;  
    title('AAR Image');
    max_color = max(abs_AAR,[],"all");
    abs_AAR(abs_AAR == 0) = inf;
    min_color = min(abs_AAR,[],"all");
    caxis([min_color, max_color]);
end

% Align all figures, you need the image processing toolbox for this!
iptwindowalign(fig1,"right",fig2,"left");
iptwindowalign(fig1,"bottom",fig3,"top");
iptwindowalign(fig3,"right",fig4,"left");
if make_MVDR == 1
    iptwindowalign(fig2,"right",fig5,"left");
    iptwindowalign(fig4,"right",fig6,"left");
end

