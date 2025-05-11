load('D:\wlwhite\image_processing\20211112_full_array_density_range\params.mat','p')

pos = 'E7-2-3';
slice = 1;

load([p.path '\corrected_data\' pos '.mat'])
seg = zeros(1080,1080,31);

p.segment_arrays = true;
p.array_channel = 3; %fluoescent channel for GFP/array
p.array_tophat_radius = 20; %radius for second tophat aplied to array chanel only
p.fringe_std = 50; %stdev of gaussian blur simulated fringe
p.fringe_size = 201; %size of gaussian kernel for simulated fringe
p.fringe_frac = 1.5; %fraction of true signal in fringe
p.array_gauss_size = 27; %size of blurring element
p.array_gauss_std = 5; %stdev of blurring element
p.array_min_quantile = 0.7; %lower cutoff (after rescaling) for pixel intensity quantile
p.array_rescale_val = 400; %raw pixel intensity to achieve half-saturation in hyperbolic rescaling
p.array_laplace_cutoff = -0.0003; %cutoff for laplace to be considered truly negative (was 0.003)
p.array_unsharp_alpha = 0.3;  % degree of unsharp filtering between 0.2 and 0.9
p.array_open = 1; %size of diamond for open operation to remove noise/background
p.array_close = 3; %size of diamond for close operation to fill in missing space due to noisy laplace cutoff
p.array_dilate = 0; %size to dilate arrays by after all other segmentation steps
p.array_mindist = 2; %minimum brightness-adjusted distance for watershed
p.array_max_hull_area = inf; %indents smaller than this will not be filled in
p.array_min_hull_ratio = 0.25; %indents with higher fractional exposed perimeter than this will not be filled in
p.array_3D_open = true; %whether or not to require arrays to be present in two connsecutive slices to be considered real

for i = p.slices
    if i == slice
        seg(:,:,i) = arrayseg(current_im(p.array_channel).im(:,:,i),p,true);
    else
        seg(:,:,i) = arrayseg(current_im(p.array_channel).im(:,:,i),p);
    end
end

%clean up 3D array segmentation to get rid of image noise
%[x,y,z] = ndgrid(-p.array_close3D:p.array_close3D);
%se = strel(sqrt(x.^2 + y.^2 + z.^2) <= p.array_close3D);

if p.array_3D_open
    base_grid = zeros([3,3,3]);
    base_grid(:,2,2) = 1;
    base_grid(2,:,2) = 1;

    se_grid1 = base_grid;
    se_grid1(2,2,1) = 1;
    se1 = strel(se_grid1);
    se_grid2 = base_grid;
    se_grid2(2,2,3) = 1;
    se2 = strel(se_grid2);
    seg = imopen(seg,se1) | imopen(seg,se2);
end


figure()
imshow(adjust(current_im(p.array_channel).im(:,:,slice),p.array_min_quantile,p.array_rescale_val))
hold on
perims = bwboundaries(seg(:,:,slice));
for i = 1:length(perims)
    perim = perims{i};
    plot(perim(:,2), perim(:,1), 'Color', 'r', 'LineWidth', 1)
end


function out = adjust(im,min_quantile,rescale)
    out = im./(im + rescale);
    
    low = quantile(out(:),min_quantile);
    high = 1;
    out = imadjust(out./high,[low/high; 1],[0 ; 1]);
end