function seg = arrayseg(im,im_CTV,p,varargin)

%% load and unpack parameters
tophat_r = p.array_tophat_radius;
min_quantile = p.array_min_quantile;
gauss_size = p.array_gauss_size;
gauss_std = p.array_gauss_std;
fringe_std = p.fringe_std;
fringe_size = p.fringe_size;
fringe_frac = p.fringe_frac;
laplace_cutoff = p.array_laplace_cutoff;
unsharp_alpha = p.array_unsharp_alpha;
close = p.array_close;
open = p.array_open;
dilate = p.array_dilate;
mindist = p.array_mindist;
rescale_val = p.array_rescale_val;
max_hull_area = p.array_max_hull_area;
min_hull_ratio = p.array_min_hull_ratio;


%% define image analysis filters
h2 = fspecial('laplacian', unsharp_alpha);  % what is this shape factor?  worry about this later.
se_close = strel('diamond',close);
se_open = strel('disk', open);   % removing debris
se_dilate = strel('disk',dilate);

%% start processing routine
im = double(im);  % convert to double
original = im; %save unprocessed image for regionprops calcs

%make list of processing functions for use in for loop
funcs = {@(im) im - im_CTV, 'CTV Subtraction';
         @(im) imtophat(im,strel('disk',tophat_r)), 'Tophat';
         @(im) max(median(im(:)), im - (imgaussfilt(im,fringe_std,'FilterSize',fringe_size).*fringe_frac)), 'Remove Fringe';
         @(im) rescale(im,rescale_val), 'Rescale';
         @(im) imgaussfilt(im,gauss_std,'FilterSize',gauss_size,'FilterDomain','frequency'), 'Gaussian Filter';
         @(im) contrast_adjust(im,min_quantile), 'Contrast Adjust';
         @(im) imfilter(im,h2), 'Laplacian';
         @(im) (im < laplace_cutoff), 'Laplace < cutoff';
         @(im) imclose(im, se_close), 'Close';
         @(im) imfill(im,'holes'), 'Fill Holes';
         @(im) imopen(im,se_open), 'Open';
         @(im) semiconvhull(im,max_hull_area,min_hull_ratio), 'Hull';
         @(im) imdilate(im,se_dilate), 'Dilate';
         @(im) watershed_plus(im,mindist), 'Watershed';
         };
     
N_steps = size(funcs,1);

show_process = false;     
if ~isempty(varargin)
    show_process = varargin{1};
end

if show_process
    figure()
    sub_cols = ceil(sqrt(N_steps+1));
    sub_rows = ceil((N_steps+1)/sub_cols);
    subplot(sub_rows,sub_cols,1)
    imshow(im,[min(im(:)),max(im(:))])
    title('Original Image')
end

for i = 1:size(funcs,1)
    fun = funcs{i,1};
    im = fun(im);
    if show_process
        subplot(sub_rows,sub_cols,i+1)
        imshow(im,[double(min(im(:))),double(max(im(:)))])
        title(funcs{i,2})
    end
end
seg = im; %this is now the fully segmented image
end

function im_out = rescale(im,rescale_val)
    im_out = im./(im+rescale_val);
end

% contrast adjustment on the gaussian filtered image
function im_out = contrast_adjust(im,min_quantile)
    low_in = quantile(im(:),min_quantile);
    high_in = 1;%quantile(im(:),max_quantile);
    im_out = imadjust(im./high_in,[low_in/high_in; 1],[0 ; 1]);
end

function seg = watershed_plus(seg,mindist)
%watershed algorithm
D = bwdist(~seg);
D = -D;
D = imhmin(D,mindist);
D(~seg) = Inf;
L = watershed(D);
seg = seg.*(L>0);
end

function im_out = semiconvhull(im,max_area,min_pr)
%get regions with small differences from hull
hull = bwconvhull(im,'objects');
diff = hull - im; %image of regions that got painted over by hull
diff_label = bwlabel(diff);
diff_props = regionprops(diff_label,'Area');
area = [diff_props.Area];
select = find(area < max_area); %get the small regions
small_regions = ismember(diff_label,select); %will subtract off of hull image later

%find regions that will not affect perimeter much if removed from hull
    %defined as regions of diff that have a low internal:total perimeter ratio
big_regions = ismember(diff_label,find(area >= max_area));
big_perim = bwperim(big_regions);
hull_perim = bwperim(hull);
big_label = bwlabel(big_perim);
%get info about the section of the perimeter of the concavities that is on the perimeter of the hull
    %area is actually perimeter because these regions are perimeters
    %mean intensity comes from hull perimeter image, so only counts pixels in the hull perimeter
    %want total pixels in union of hull perimeter and diff perimeter --> multiply mean intensity by area
surface_props = regionprops(big_label,hull_perim,'Area','MeanIntensity');
surface_perim = [surface_props.Area].*[surface_props.MeanIntensity];
ratio = surface_perim./[surface_props.Area];
select = find(ratio > min_pr); %get all diff regions that have high surface perimeter
high_exposure = ismember(big_label,select);
high_exposure = imfill(high_exposure,'holes');

im_out = hull - small_regions - high_exposure; %undo hull for regions that are small or have high exposed perimeter
end