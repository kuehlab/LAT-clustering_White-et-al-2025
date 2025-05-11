function seg = cellseg_v2(im,im_GFP,p,varargin)

%% load and unpack parameters
min_val = p.min_val; %for contrast adjustment
max_val = p.max_val; %for contrast adjustment
mindist = p.mindist;   % minimum value for the watershed minima
minsize = p.minsize;   % minimum size of a cell
maxsize = p.maxsize; % maximum size of a cell
min_brightness = p.min_brightness; %minimum mean brightnenss of a cell
if isfield(p,'max_autofluor_div') && isfield(p,'max_autofluor_mult')
    max_autofluor_div = p.max_autofluor_div; %the GFP autofluorescence cutiff for a dead cell to be thrown out
    max_autofluor_mult = p.max_autofluor_mult; %the GFP autofluorescence cutiff for a dead cell to be thrown out
else
    max_autofluor_div = inf;
    max_autofluor_mult = inf;
end
if isfield(p,'min_grad')
    min_grad = p.min_grad; %minimum gradient magnitude (to drop out of focus cells)
else
    min_grad = 0;
end

maxecc = p.maxecc;    % the maximum eccentricity of a cell
unsharp_size = p.unsharp_size;  % 'size' of the unsharp filter
unsharp_alpha = p.unsharp_alpha;  % degree of unsharp filtering between 0.2 and 0.9
laplace_cutoff = p.laplace_cutoff; %maximum acceptible value of laplacian
gauss_size = p.gauss_size; %size of gaussian blur
gauss_std = p.gauss_std; %stdev of gaussian blur
close1 = p.close1; %size of diamond for first close operation (after laplacian threshholding)
min_hole_brightness = p.min_hole_brightness; %minimum intensity of area in order to be filled in during the fill holes step
open1 = p.open1; %size of disk to use for open operation to remove debris
max_hull_area = p.max_hull_area;
min_hull_ratio = p.min_hull_ratio;
dilate = p.dilate;

%% define image analysis filters
h1 = fspecial('gaussian',gauss_size, gauss_std); % 2D guassian filter
h2 = fspecial('laplacian', unsharp_alpha);  % what is this shape factor?  worry about this later.
se1 = strel('diamond',close1);
se2 = strel('disk', open1);   % removing debris
se3 = strel('disk',dilate);

%% start processing routine
im = double(im);  % convert to double
original = im; %save unprocessed image for regionprops calcs

%make list of processing functions for use in for loop
funcs = {@(im) imfilter(im,h1,'replicate'), 'Gaussian Filter';
         @(im) imtophat(im,strel('disk',2*p.cellsize)), 'Tophat';
         @(im) contrast_adjust(im,min_val,max_val), 'Contrast Adjust'; %see function defined below
         @(im) imfilter(im,h2), 'Laplacian';
         @(im) (im < laplace_cutoff), 'Laplace < cutoff';
         %@(im) padarray(im(2:end-1,2:end-1),[1 1],0,'both'), 'Clear Edge Pixels';
         @(im) imclose(im, se1), 'Close';
         @(im) fill_real_holes(im,original,min_hole_brightness), 'Fill Holes';
         @(im) imopen(im,se2), 'Open';
         @(im) watershed_plus(im,mindist), 'Watershed Separation'; %see function defined below
         @(im) semiconvhull(im,max_hull_area,min_hull_ratio), 'Hull'; %see function defined below
         @(im) watershed_plus(im,mindist), 'Watershed Separation'; %see function defined below
         @(im) imdilate(im,se3), 'Dilate';
         @(im) watershed_plus(im,mindist), 'Watershed Separation'; %see function defined below
         @(im) imclearborder(im,8), 'Clear Border';
         @(im) bwlabel(im), 'Label Objects';
         @(im) filter_obj(im,original,im_GFP,minsize,maxsize,maxecc,min_brightness,max_autofluor_div,max_autofluor_mult,min_grad), 'Objects Filtered'; %see function defined below
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

% contrast adjustment on the gaussian filtered image
function im_out = contrast_adjust(im,min_val,max_val)
    high_in = max_val; %quantile(im(:),max_quantile);
    low_in = min_val; %quantile(im(:),min_quantile);
    im_out = imadjust(im./high_in,[low_in/high_in; 1],[0 ; 1]);
end

function seg = filter_obj(im,original,im_GFP,minsize,maxsize,maxecc,min_brightness,max_autofluor_div,max_autofluor_mult,min_grad)
props = regionprops(im, original, 'Area','Perimeter','MajorAxisLength','MinorAxisLength','MeanIntensity','PixelIdxList'); % properties of objects

[im_grad,~] = imgradient(original);
grad_props = regionprops(im, im_grad,'PixelValues');
grad_vals = cellfun(@(x) quantile(x,0.9),{grad_props.PixelValues});

pixels = {props.PixelIdxList};
autofluor_div_im = im_GFP./original;
autofluor_mult_im = im_GFP.*original;
autofluor_div = cellfun(@(x) quantile(autofluor_div_im(x),0.9), pixels);
autofluor_mult = cellfun(@(x) quantile(autofluor_mult_im(x),0.9), pixels);

area = [props.Area]; % area of each object
select = intersect(find(area > minsize), find(area < maxsize));   % select objects of a certain size

% ratio = [props.Perimeter].^(2)./[props.Area]; % perimeter squared/area ratio of each object
% [w] = wobble(i9,im); % wobble of each object
select = intersect(select, find([props.MajorAxisLength]./[props.MinorAxisLength] < maxecc));  % select objects that are not too eccentric
select = intersect(select, find([props.MeanIntensity] > min_brightness));
%select = intersect(select, find(ratio < maxratio)); % select objects that are below a certain perimeter squared/area ratio 
% select = intersect(select, find(w < maxwobble)); % select objects without too much wobble

%filter out dead cells based on GFP autofluorescence
select = intersect(select, find(autofluor_div < max_autofluor_div));
select = intersect(select, find(autofluor_mult < max_autofluor_mult));

%filter out out of focus cells
select = intersect(select, find(grad_vals > min_grad));

seg = ismember(im, select); % select objects within parameters from labeled objects060117 to output
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

function filled = fill_real_holes(im,original,min_brightness)
    filled = imfill(im,'holes');
    holes = filled & ~im;
    
    hole_props = regionprops(holes,original,'MeanIntensity','PixelIdxList');
    dim_holes = [hole_props.MeanIntensity] < min_brightness;
    
    if any(dim_holes)
        pxls_to_unfill = vertcat(hole_props(dim_holes).PixelIdxList);
        filled(pxls_to_unfill) = 0;
    end
end