clearvars

% get parameters
load('D:\wlwhite\image_processing\20240201_array_density_range_5min_fresh_arrays/multiseg_params.mat','p');

% preprocess images
disp('Collecting and background-subtracting images...')
collate_images_from_pos_files(p)      

%% run segmentation
disp('Segmenting cells and arrays...')
multisegment(p)

%% extract relevant information about cells
disp('Extracting data...')
% extract_seg_overlap_data(p)
extract_multiseg_data_with_random(p)
