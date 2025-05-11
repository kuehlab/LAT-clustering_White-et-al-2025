clearvars

%p.path = '/Volumes/tinytim/wlwhite/image_processing/20211124_full_array_density_range/';
p.path = 'D:\wlwhite\image_processing\20240201_array_density_range_5min_fresh_arrays/';

p.basename = '20240201';
p.pos_files = {'raw_data'};
p.out_prefix = 'multiseg_';

p.flatten = [2,3,4,5]; %channels to apply background subtraction to
p.flatfield_correction = 'D:\wlwhite\image_processing\20220801_flatfield_correction\correction_vals.mat';
p.cellsize = 100; %radius of average cell in pixels

p.im_X = 1080;
p.im_Y = 1080;

%cell segmentation
p.seg_CTV.segfun = 'cellseg_v2';
p.seg_CTV.mindist = 5;   % minimum value for the watershed minima
p.seg_CTV.min_val = 250;   % pixel intensity quantile  of background
p.seg_CTV.max_val = 7000;   %pixel initensity quantile of max foreground value
p.seg_CTV.minsize = 1500;   % minimum size of a cell (area in pixels)
p.seg_CTV.maxsize = 60000; % maximum size of a cell (area in pixels)
p.seg_CTV.maxecc = 3;    % the maximum eccentricity of a cell
p.seg_CTV.cellsize = 100; %radius of average cell in pixels
p.seg_CTV.unsharp_size = 7;  % 'size' of the unsharp filter
p.seg_CTV.unsharp_alpha = 0.8;  % degree of unsharp filtering between 0.2 and 0.9
p.seg_CTV.laplace_cutoff = -0.00007; %maximum acceptible value of laplacian
p.seg_CTV.gauss_size = 21; %size of gaussian blur
p.seg_CTV.gauss_std = 4; %stdev of gaussian blur
p.seg_CTV.close1 = 8; %size of diamond for first close operation (after laplacian threshholding)
p.seg_CTV.min_hole_brightness = 1400; %minimum intensity of area in order to be filled in during the fill holes step
p.seg_CTV.open1 = 12; %size of disk to use for open operation to remove debris
p.seg_CTV.max_hull_area = 20; %maximum area of indent that will not be filled in by convhull
p.seg_CTV.min_hull_ratio = 0.2; %minimum exposed perimeter to total perimeter ratio of indent that will not be filled by convhull
p.seg_CTV.dilate = 2; %size of disk to dilate cell borders by at end of segmentation
p.seg_CTV.min_brightness = 2000; % the mininmum required mean pixel intensity in a cell
p.seg_CTV.close3D = 1; %radius to use for close operation on full z-stack after all 2D segmentation is done
p.seg_CTV.min_n_slices = 6; %niminum number of slices for a cell to be counted
p.seg_CTV.min_grad = 1900; %minimum gradient magnitude (to drop out of focus cells)

%array segmentation
p.seg_GFP.segfun = 'arrayseg';
p.seg_GFP.array_tophat_radius = 20; %radius for second tophat applied to array chanel only
p.seg_GFP.fringe_std = 15; %stdev of gaussian blur simulated fringe
p.seg_GFP.fringe_size = 201; %size of gaussian kernel for simulated fringe
p.seg_GFP.fringe_frac = 1.2; %fraction of true signal in fringe
p.seg_GFP.array_gauss_size = 31; %size of blurring element
p.seg_GFP.array_gauss_std = 4; %stdev of blurring element
p.seg_GFP.array_min_quantile = 0.3; %lower cutoff (after rescaling) for pixel intensity quantile
p.seg_GFP.array_rescale_val = 800; %raw pixel intensity to achieve half-saturation in hyperbolic rescaling
p.seg_GFP.array_laplace_cutoff = -0.0004; %cutoff for laplace to be considered truly negative (was 0.003)
p.seg_GFP.array_unsharp_alpha = 0.3;  % degree of unsharp filtering between 0.2 and 0.9
p.seg_GFP.array_open = 1; %size of diamond for open operation to remove noise/background
p.seg_GFP.array_close = 3; %size of diamond for close operation to fill in missing space due to noisy laplace cutoff
p.seg_GFP.array_dilate = 1; %size to dilate arrays by after all other segmentation steps
p.seg_GFP.array_mindist = 2; %minimum brightness-adjusted distance for watershed
p.seg_GFP.array_max_hull_area = inf; %indents smaller than this will not be filled in
p.seg_GFP.array_min_hull_ratio = 0.25; %indents with higher fractional exposed perimeter than this will not be filled in
p.seg_GFP.array_3D_open = true; %whether or not to require arrays to be present in two connsecutive slices to be considered real

%pLAT segmentation
p.seg_pLAT.segfun = 'arrayseg';
p.seg_pLAT.array_tophat_radius = 20; %radius for second tophat applied to array chanel only
p.seg_pLAT.fringe_std = 15; %stdev of gaussian blur simulated fringe
p.seg_pLAT.fringe_size = 201; %size of gaussian kernel for simulated fringe
p.seg_pLAT.fringe_frac = 1.0; %fraction of true signal in fringe
p.seg_pLAT.array_gauss_size = 31; %size of blurring element
p.seg_pLAT.array_gauss_std = 4.5; %stdev of blurring element
p.seg_pLAT.array_min_quantile = 0.3; %lower cutoff (after rescaling) for pixel intensity quantile
p.seg_pLAT.array_rescale_val = 400; %raw pixel intensity to achieve half-saturation in hyperbolic rescaling
p.seg_pLAT.array_laplace_cutoff = -0.0004; %cutoff for laplace to be considered truly negative (was 0.003)
p.seg_pLAT.array_unsharp_alpha = 0.3;  % degree of unsharp filtering between 0.2 and 0.9
p.seg_pLAT.array_open = 1; %size of diamond for open operation to remove noise/background
p.seg_pLAT.array_close = 3; %size of diamond for close operation to fill in missing space due to noisy laplace cutoff
p.seg_pLAT.array_dilate = 1; %size to dilate arrays by after all other segmentation steps
p.seg_pLAT.array_mindist = 2; %minimum brightness-adjusted distance for watershed
p.seg_pLAT.array_max_hull_area = inf; %indents smaller than this will not be filled in
p.seg_pLAT.array_min_hull_ratio = 0.25; %indents with higher fractional exposed perimeter than this will not be filled in
p.seg_pLAT.array_3D_open = true; %whether or not to require arrays to be present in two connsecutive slices to be considered real

%pCD3z segmentation
p.seg_pCD3z.segfun = 'arrayseg';
p.seg_pCD3z.array_tophat_radius = 20; %radius for second tophat applied to array chanel only
p.seg_pCD3z.fringe_std = 15; %stdev of gaussian blur simulated fringe
p.seg_pCD3z.fringe_size = 201; %size of gaussian kernel for simulated fringe
p.seg_pCD3z.fringe_frac = 1.2; %fraction of true signal in fringe
p.seg_pCD3z.array_gauss_size = 31; %size of blurring element
p.seg_pCD3z.array_gauss_std = 4; %stdev of blurring element
p.seg_pCD3z.array_min_quantile = 0.3; %lower cutoff (after rescaling) for pixel intensity quantile
p.seg_pCD3z.array_rescale_val = 800; %raw pixel intensity to achieve half-saturation in hyperbolic rescaling
p.seg_pCD3z.array_laplace_cutoff = -0.0004; %cutoff for laplace to be considered truly negative (was 0.003)
p.seg_pCD3z.array_unsharp_alpha = 0.3;  % degree of unsharp filtering between 0.2 and 0.9
p.seg_pCD3z.array_open = 1; %size of diamond for open operation to remove noise/background
p.seg_pCD3z.array_close = 3; %size of diamond for close operation to fill in missing space due to noisy laplace cutoff
p.seg_pCD3z.array_dilate = 1; %size to dilate arrays by after all other segmentation steps
p.seg_pCD3z.array_mindist = 2; %minimum brightness-adjusted distance for watershed
p.seg_pCD3z.array_max_hull_area = inf; %indents smaller than this will not be filled in
p.seg_pCD3z.array_min_hull_ratio = 0.25; %indents with higher fractional exposed perimeter than this will not be filled in
p.seg_pCD3z.array_3D_open = true; %whether or not to require arrays to be present in two connsecutive slices to be considered real

p.background_dilate_size = 5; %number of pixels around array/cell to take for background
p.background_sub_quantile = 0.5; %quantile of cellular overlap with ring around array to take as background
p.discount_from_bkgnd = [3,4,5]; %channel numbers whose segmentations should be removed from the background donut before calculating background value

p.min_num_placements = 100; %minimum nunmber of allowed random placements of an array to inculde a random placement of it in the background distribution
p.max_cell_bkgrnd_dist = 200; %max distance from a cell centroid for an array to use cell background as it's own
p.frac_keep_not_in_cell = 1; %fraction of the arrays not touching a cell to keep when calculating summary stats

% p.array_clump_max = 10.^2.7; %maximum median intensity in array to include it's volume in the cell-based stats

p.channels = {'DIC';'CTV';'GFP';'pLAT';'pCD3z'};
p.channel_names = {'Camera DIC'; 'LDI 405'; 'LDI 470';'LDI 555'; 'LDI 640'};

%data extraction parameters
p.stain_channels = [3,4,5]; %which chennels are activation stains (these will get real and randomized data)
p.summary_funcs = {
    'med', @(x) median(x);
    'max', @(x) max(x);
    'tot', @(x) sum(x)
};

%not really good variable names but need to use them to get segmentation to run
p.time_names = {'B','C','D'}; %these are the row names for each time point
% p.time_plates = {'1-';'1-';'2-';'2-'}; %these are the plate names corresponding to each row/time point
p.density_names = {'2','3','4','5','6','7','8','9'}; %these are the column names for each array density

p.grad_cutoff = 0;
p.grad_border_width = 2;
p.grad_quantile = 0.8;

p.slices = 1:6; %which z-slices to use for quantification (must be continuous and start from 1)

p.omit_brighter_than = [inf, inf, inf, inf, inf]; %omit images with pixels brighter than this in the corresponding channel (gets rid of brirght Ab-fluorophore blobs)
p.blob_bkgnd_cut = [inf,inf,inf,300,800]; %segmented regions with background values higher than this are trigger the area around them to be removed from analysis
p.blob_tot_int_cut = [inf,inf,inf,1000,10000]; %segmented regions with median+backgrond values higher than this are trigger the area around them to be removed from analysis
p.blob_clear_radius = 100; %how many pixels around bright blobs to clear
p.blob_max_number = 20; %ignore images with more than this many blobs

p.save_wells_separately = true;

p.omit = {
    }; %which images to omit from quantification

save([p.path p.out_prefix 'params.mat'],'p')