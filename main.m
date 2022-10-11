% This code tests the color image segmentation using the CIELAB transform
% and the k-means clustering
clc, clear, close all

% Reading the image

[fn, fpath] = uigetfile('*.jpg; *.png');
I = imread(fullfile(fpath, fn));

[m, n, ~] = size(I);

% Removing the blue dimension marker and upper white portion
blue_mark_indx = I(:, :, 1) == 0 & I(:, :, 2) == 0 & I(:, :, 3) == 255;

white_portion_indx = I(:, :, 1) > 240 & I(:, :, 2) > 240 & I(:, :, 3) > 240;
white_portion_indx = bwareafilt(~bwfill(~white_portion_indx, 'holes'), ...
    [1e4, inf]);
I = I .* uint8(repmat(~or(blue_mark_indx, white_portion_indx), [1, 1, 3]));

% Performing the preprocssing
Iycr = rgb2ycbcr(I);
Iy = Iycr(:, :, 1);
% Iypre = medfilt2(imadjust(adapthisteq(Iy)), [5, 5]);
Iypre = medfilt2(histeq(imadjust(Iy)), [5, 5]);


Irec1 = Iycr; Irec1(:, :, 1) = Iypre;
Ipre = ycbcr2rgb(Irec1);

Ipre = Ipre .* uint8(repmat(~or(blue_mark_indx, white_portion_indx), [1, 1, 3]));
figure, imshow(Ipre), impixelinfo
%%

% Define number of desired color clusters 
nColors = 4;
 
% Perform the color image clustering k-means
[output_images, bw_output_images] = img_k_means_fcn(im2double(Ipre), nColors);
red_points = rescale(Iycr(:, :, 3)) > 0.8;

red_counter = zeros(1, nColors);
%%
for i = 1:nColors
    
%     figure, imshow(output_images{i}), impixelinfo
    red_mask = red_points .* bw_output_images{i};
    red_counter(i) = sum(red_mask(:));
end
%%
[~, col_max_idx] = max(red_counter);
collagen_layer = output_images{col_max_idx};
collagen_mean_profile =mean(bw_output_images{col_max_idx}(:, :, 1), 2);
collagen_mean_profile([1:round(0.25*m), round(0.75*m):end]) = 0;
collagen_mean_smooth = movmean(collagen_mean_profile, 50);
collagen_thresh = collagen_mean_smooth > (0.75 * max(collagen_mean_smooth));

collagen_points = find(collagen_thresh);
collagen_endpoints = [collagen_points(1), collagen_points(end)];
collagen_mask = zeros(m, n);
collagen_mask(collagen_endpoints(1):collagen_endpoints(2), :) = 1;
% plot(collagen_mean_profile)
% % disp(collagen_endpoints)
collagen_layer_masked = imdilate(collagen_mask, strel('disk', 50)).* collagen_layer;
figure, imshow(im2uint8(collagen_layer_masked)) , impixelinfo
title('Collagen Layer')

musc_mask = zeros(m, n);
musc_mask(collagen_endpoints(end):end, :) = 1;

musc_counter = zeros(1, nColors);
%%
for i = 1:nColors
    this_musc_mask = musc_mask .* bw_output_images{i};
    musc_counter(i) = sum(this_musc_mask(:));
end
[~, musc_max_idx] = max(musc_counter);
musc_layer = output_images{musc_max_idx};
musc_layer_masked = musc_mask .* musc_layer;
figure, imshow(musc_layer_masked)
title('Muscle Layer')
%    mean_profile = mean(bw_output_images{i}(:, :, 1), 2);
%    figure, plot(mean_profile)
% end

mucin_mask = zeros(m, n);
mucin_mask(1:collagen_endpoints(1), :) = 1;

inds = 1:nColors;
inds([musc_max_idx, col_max_idx]) = [];
mucin_counter = zeros(1, length(inds));

%%
for i = 1:length(inds)
%     figure, imshow(rgb2gray( mucin_mask .* output_images{inds(i)}) > 0.5)
    this_mucin_mask = mucin_mask .* rgb2gray(output_images{inds(i)}) > 0.5;
    mucin_counter(i) = sum(this_mucin_mask(:));
%     title(mucin_counter(i))
end
%%
[~, mucin_max_idx] = max(mucin_counter); 
mucin_layer = output_images{inds(mucin_max_idx)};
mucin_layer_masked = mucin_mask .* mucin_layer;

mucin_mean_profile =mean(mucin_layer_masked(:, :, 1), 2);
% mucin_mean_profile([1:round(0.25*m), round(0.75*m):end]) = 0;
mucin_mean_smooth = movmean(mucin_mean_profile, 50);
mucin_thresh = mucin_mean_smooth > (0.5 * max(mucin_mean_smooth));

mucin_points = find(mucin_thresh);
mucin_endpoints = [mucin_points(1), mucin_points(end)];
mucin_mask2 = zeros(m, n);
mucin_mask2(mucin_endpoints(1):mucin_endpoints(2), :) = 1;
% plot(collagen_mean_profile)
% disp(mucin_endpoints)
mucin_layer_masked2 = imdilate(mucin_mask2, strel('disk', 30)).* mucin_layer_masked;
figure, imshow(im2uint8(mucin_layer_masked2)) , impixelinfo

figure, imshow(mucin_layer_masked2)

%%
color_points = [1, collagen_endpoints(1), m, collagen_endpoints(1);
                1, collagen_endpoints(end), m, collagen_endpoints(end);
                1, mucin_endpoints(end), m, mucin_endpoints(end)];
            
I_disp = insertShape(Ipre, 'line', color_points, 'Opacity', 1, ...
    'LineWidth', 15, 'color', 'yellow');
% I_% disp(color_points, :, 1) = 0;

figure, imshow(I_disp), title('Layer Bounaries')