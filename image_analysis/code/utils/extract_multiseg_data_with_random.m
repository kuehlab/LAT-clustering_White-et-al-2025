function extract_multiseg_data_with_random(p)
    if ~isfield(p, 'in_path')
        in_path = p.path;
        out_path = p.path;
    else
        in_path = p.out_path; %use out path here because it assumes only the raw data is in the in_path location
        out_path = p.out_path;
    end
    
    if ~isfield(p, 'out_prefix')
        p.out_prefix = '';
    end
    
    if ~isfield(p, 'frac_keep_not_in_cell')
        p.frac_keep_not_in_cell = 1;
    end 

    rows = p.time_names;
    cols = p.density_names;
    
    if isfield(p, 'time_plates')
        plates = p.time_plates;
    else
        plates = repmat({''},size(rows));
    end
    
    if ~isfield(p, 'omit_brighter_than')
        p.omit_brighter_than = ones([1,len(p.channels)]) * inf;
    end
    
    if ~isfield(p, 'discount_from_bkgnd')
        p.discount_from_bkgnd = [];
    end
    
    mkdir([out_path p.out_prefix 'analyzed_data_array_background']);
        
    se_size = p.background_dilate_size;
    [x,y,z] = ndgrid(-se_size:se_size);
    se1 = strel(sqrt(x.^2 + y.^2 + z.^2) <= se_size);
    
    rowcol = cell(length(rows).*length(cols),2);
    for r = 1:length(rows)
        for c = 1:length(cols)
            idx = c + (r-1)*length(cols);
            rowcol{idx,1} = [plates{r} rows{r}];
            rowcol{idx,2} = cols{c};
        end
    end

    parfor idx = 1:length(rowcol) %parallel
        row = rowcol{idx,1};
        col = rowcol{idx,2};
            
        well_data = struct();
        for ch = p.stain_channels
            well_data.(p.channels{ch}) = [];
        end

        file_info = dir([in_path  p.out_prefix 'corrected_data/' row col '*.mat']);
        for f = 1:length(file_info)
            file = file_info(f).name;

            if any(strcmp(p.omit,file))
                disp(['Skipping ' file])
                continue
            end

            disp({file, datetime('now')})
            loaded_im = load([in_path p.out_prefix 'corrected_data/' file],'current_im');
            current_im = loaded_im.current_im;
            loaded_seg = load([in_path p.out_prefix 'segmentations/' file],'seg');
            seg = loaded_seg.seg;

            %check for bright antibody blobs and ignore images with them
            skip = false;
            for ch = 1:length(p.channels)
                if quantile(current_im(ch).im, 0.999, 'all') > p.omit_brighter_than(ch)
                    skip = true;
                    break
                end
            end
            if skip
                disp({file , 'Skipped due to bright blob'})
                continue %skip this image
            end

            if isfield(seg,'CTV')
                cells = logical(seg.CTV(:,:,p.slices));
                cell_channel = find(strcmp('CTV',p.channels));
            else
                cells = false(p.im_X,p.im_Y,length(p.slices));
            end

            if isfield(p,'grad_cutoff')
                %if requested, leave out the slices for each cell that are out of focus based on the gradient of
                %the cell image
                corrected_cells = cells;

                props = table2struct(regionprops3(cells,'VoxelIdxList'));

                grad_stack = zeros([p.im_X,p.im_Y,max(p.slices)]);
                for z = p.slices
                    [gradient,~] = imgradient(current_im(cell_channel).im(:,:,z));
                    grad_stack(:,:,z) = gradient;
                end

                for i = 1:length(props)
                    bw_3D = zeros(size(cells));
                    bw_3D(props(i).VoxelIdxList) = 1;
                    for z = p.slices
                        bw = bw_3D(:,:,z);
                        border = logical(imdilate(bw,strel('disk',p.grad_border_width)) - ...
                            imerode(bw,strel('disk',p.grad_border_width)));
                        grad = grad_stack(:,:,z);
                        if quantile(grad(border),p.grad_quantile) < p.grad_cutoff
                            %if the cell is below the cutoff in this slice, delete this slice of the cell from the
                            %segmentation prior to analysis
                            corrected_cells_slice = corrected_cells(:,:,z);
                            corrected_cells_slice(logical(bw)) = 0;
                            corrected_cells(:,:,z) = corrected_cells_slice;
                        end
                    end
                end
                cells = corrected_cells;
            end

            %% extract data for each indicated channel's segmentation
            all_im_data = struct();
            for ch = p.stain_channels

                seg_ch_name = p.channels{ch};
                ch_seg_bw = seg.(seg_ch_name);

                regions = table2struct(regionprops3(ch_seg_bw,'VoxelIdxList','VoxelList','Centroid','BoundingBox'));

                im_data = struct('fname',{},'centroid',{},'volume',{},'n_slices',{},'frac_in_cell',{},...
                                 'edge_dist',{},'n_bkgnd',{},'rand_n_bkgnd',{},'quant_pixels',{});
                for chq = p.stain_channels
                    quant_ch_name = p.channels{chq};
                    for d = 1:length(p.summary_funcs)
                        func_name = p.summary_funcs{d,1};
                        [im_data(:).([quant_ch_name '_' func_name])] = deal(0);
                        [im_data(:).(['rand_' quant_ch_name '_' func_name])] = deal(0);
                    end
                    [im_data(:).([quant_ch_name '_bkgnd'])] = deal(0);
                    [im_data(:).(['rand_' quant_ch_name '_bkgnd'])] = deal(0);
                end

                for i = 1:length(regions) %parallel

                    %section out the region of the image around the array to improve speed
                    bb = round(regions(i).BoundingBox);
                    max_i = size(current_im(1).im,2);
                    max_j = size(current_im(1).im,1);
                    max_k = length(p.slices);
                    start_i = max(1,bb(2) - 2*se_size);
                    start_j = max(1,bb(1) - 2*se_size);
                    start_k = max(1,bb(3) - 2*se_size);
                    end_i = min(max_i, bb(2) + bb(5) + 2*se_size);
                    end_j = min(max_j, bb(1) + bb(4) + 2*se_size);
                    end_k = min(max_k ,bb(3) + bb(6) + 2*se_size);

                    seg_bw = zeros(size(current_im(1).im));
                    seg_bw(regions(i).VoxelIdxList) = 1;
                    seg_bw = seg_bw(start_i:end_i,start_j:end_j,start_k:end_k);
                    big_seg = imdilate(seg_bw,se1);

                    cells_bb = cells(start_i:end_i,start_j:end_j,start_k:end_k); 

                    if any(cells(regions(i).VoxelIdxList),'all') %for arrays touching a cell
                        in_cell = true;
                        %background is area around speficic array that overlaps a cell, minus any area covered by another array
                        background_bw = big_seg & ~ch_seg_bw(start_i:end_i,start_j:end_j,start_k:end_k) & cells_bb;
                        %area of interest is the array overlap with a cell
                        quant_pixels = cells_bb & seg_bw;
                    else
                        in_cell = false;
                        %background is area around speficic array that overlaps a cell, minus any area covered by another array
                        background_bw = big_seg & ~ch_seg_bw(start_i:end_i,start_j:end_j,start_k:end_k) & ~cells_bb;
                        %area of interest is the array not overlapping a cell
                        quant_pixels = ~cells_bb & seg_bw;
                    end

                    %remove pixels in other channels' segmentations if requested
                    for ch2 = p.discount_from_bkgnd
                        ch2_name = p.channels{ch2};
                        ch2_seg = seg.(ch2_name);
                        background_bw = background_bw & ~ch2_seg(start_i:end_i,start_j:end_j,start_k:end_k);
                    end
                    n_background = sum(background_bw,'all');

                    for chq = p.stain_channels

                        quant_ch_name = p.channels{chq};
                        im = current_im(chq).im(start_i:end_i,start_j:end_j,start_k:end_k);

                        background = quantile(im(background_bw),p.background_sub_quantile);
                        im = im - background;
                        im_data(i).([quant_ch_name '_bkgnd']) = background;

                        for d = 1:length(p.summary_funcs)
                            func_name = p.summary_funcs{d,1};
                            func = p.summary_funcs{d,2};
                            im_data(i).([quant_ch_name '_' func_name]) = func(im(quant_pixels));
                        end
                    end

                    quant_props = regionprops(quant_pixels,'BoundingBox');
                    bbs = cell2mat({quant_props.BoundingBox}');
                    if length(bbs(1,:)) == 4 %there is only 2D data
                        n_slices = 1;
                    else
                        n_slices = max(bbs(:,3)+bbs(:,6)) - min(bbs(:,3));
                    end

                    %calculate the distance to the edge of the cell
                    slice = int16(regions(i).Centroid(3));
                    cell_bdrys = bwboundaries(cells(:,:,slice));
                    bdry_coords = vertcat(cell_bdrys{:});
                    if ~isempty(bdry_coords)
                        dists = sqrt(sum((bdry_coords - regions(i).Centroid(1:2)).^2,2));
                    else
                        dists = -1; %default value if no cells
                    end

                    im_data(i).edge_dist = min(dists,[],1);
                    im_data(i).n_bkgnd = n_background;
                    im_data(i).n_slices = n_slices;
                    im_data(i).fname = file;
                    im_data(i).centroid = regions(i).Centroid;
                    im_data(i).volume = sum(quant_pixels,'all');
                    im_data(i).frac_in_cell = sum(seg_bw & cells_bb,'all')/sum(seg_bw,'all');
                    im_data(i).quant_pixels = quant_pixels;
                end

%                     im_data = im_data(cellfun(@(x) ~isempty(x),{im_data.in_cell}));
                all_im_data.(seg_ch_name) = im_data; %need a place to score all data for this image for all quantified channels so can remove things from it later
            end

            %% remove the requested fraction of segmentations not in cells
            for ch = p.stain_channels
                seg_ch_name = p.channels{ch};
                ch_data = all_im_data.(seg_ch_name);

                keep_idx = ([ch_data.frac_in_cell] > 0) | (rand(1,length(ch_data)) < p.frac_keep_not_in_cell);
                all_im_data.(seg_ch_name) = ch_data(keep_idx);
            end

            %% remove segmentations close to bright blobs to avoid signal contamination
            %find bright blob locations
            bad_locs = [];
            for blob_ch = p.stain_channels
                blob_ch_name = p.channels{blob_ch};
                blob_data = all_im_data.(blob_ch_name);
                background = [blob_data.([blob_ch_name '_bkgnd'])];
                med = [blob_data.([blob_ch_name '_med'])];
                bad_idx = (background > p.blob_bkgnd_cut(blob_ch)) & ((background+med) > p.blob_tot_int_cut(blob_ch));

                centroids = {blob_data(bad_idx).centroid};
                bad_locs = [bad_locs; cat(1,centroids{:})];
            end

            if length(bad_locs) > p.blob_max_number
                disp(['Skipping ' file '. Too many bright blobs.'])
                continue
            end

            disp(['Removing objects near ' num2str(length(bad_locs)) ' bright blobs. ' file])

            %remove regions near bright blobs
            if ~isempty(bad_locs)
                for ch = p.stain_channels
                    seg_ch_name = p.channels{ch};
                    ch_data = all_im_data.(seg_ch_name);

                    %determine closest blob to each segmented region
                    locs = {ch_data.centroid};
                    locs = cat(1,locs{:});
                    dists = sqrt((locs(:,1) - bad_locs(:,1)').^2 + (locs(:,2) - bad_locs(:,2)').^2);
                    dists = min(dists,[],2);

                    %remove data too close to blobs
                    good_idx = dists > p.blob_clear_radius;
                    all_im_data.(seg_ch_name) = ch_data(good_idx);
                end
            end

            %% calculate data for randomized locations for remaining points
            for ch = p.stain_channels
                seg_ch_name = p.channels{ch};
                ch_data = all_im_data.(seg_ch_name);

                for i = 1:length(ch_data)
                    [rand_quant_pixels, rand_background_bw,...
                    rand_start_i, rand_end_i,...
                    rand_start_j, rand_end_j,...
                    rand_start_k, rand_end_k,...
                    skip] = get_bbox_for_random_placement(ch_data(i).quant_pixels, seg, (ch_data(i).frac_in_cell > 0), bad_locs, p, file, out_path);

                    for chq = p.stain_channels
                        
                        quant_ch_name = p.channels{chq};
                        
                        if ~skip %don't calculate random props if the array cannot be randomly placed
                            rand_im = current_im(chq).im(rand_start_i:rand_end_i,rand_start_j:rand_end_j,rand_start_k:rand_end_k);
                            rand_bkgnd = quantile(rand_im(rand_background_bw),p.background_sub_quantile);
                            rand_im = rand_im - rand_bkgnd;
                            ch_data(i).(['rand_' quant_ch_name '_bkgnd']) = rand_bkgnd;
                        else
                            ch_data(i).(['rand_' quant_ch_name '_bkgnd']) = NaN;
                        end

                        for d = 1:length(p.summary_funcs)
                            func_name = p.summary_funcs{d,1};
                            func = p.summary_funcs{d,2};

                            if ~skip %calculate random props if the seg can be placed randomly
                                ch_data(i).(['rand_' quant_ch_name '_' func_name]) = func(rand_im(rand_quant_pixels));
                            else %otherwise fill in NaN
                                ch_data(i).(['rand_' quant_ch_name '_' func_name]) = NaN;
                            end
                        end
                    end
                    if ~skip
                        ch_data(i).rand_n_bkgnd = sum(rand_background_bw,'all');
                    else
                        ch_data(i).rand_n_bkgnd = NaN;
                    end
                end
                %now that the data for this channel fo this image will not be modified anymore, can add it to the pool of data for the whole well
                ch_data = rmfield(ch_data,'quant_pixels');
                well_data.(seg_ch_name) = cat(2,well_data.(seg_ch_name),ch_data);
            end
        end

        parsave(p, out_path, row, col, well_data)

    end
end

function [rand_quant_pixels, rand_background_bw,...
    rand_start_i, rand_end_i,...
    rand_start_j, rand_end_j,...
    rand_start_k, rand_end_k,...
    skip] = get_bbox_for_random_placement(quant_pixels, seg, in_cell, blob_locs, p, fname, out_path)
    
    if isfield(seg,'CTV')
        cells = seg.CTV;
    else
        f = fields(seg);
        f = f{1};
        cells = zeros(size(seg.(f)));
    end
    
    if ~isfield(p,'min_num_placements')
        p.min_num_placements = 2;
    end
        
    %determine size/shape of array
    quant_props = table2struct(regionprops3(quant_pixels,'VoxelList','BoundingBox'));
    
    voxels = {quant_props.VoxelList};
    voxels = cat(1,voxels{:});
    all_bbox = {quant_props.BoundingBox};
    all_bbox = cat(1,all_bbox{:});
    bounding_coords = [min(all_bbox(:,1:3),[],1),...
        max(all_bbox(:,1:3)+all_bbox(:,4:6),[],1)-min(all_bbox(:,1:3),[],1)];
    
    %set the allowed region for random location sampling (in or out of cell)
    radius = ceil(sqrt(sum(bounding_coords(4:5).^2)));
    if in_cell
        allowed = max(imerode(cells,strel('disk',radius)),[],3);
    else
        allowed = max(~imdilate(cells,strel('disk',radius)),[],3);
    end
    %have to use whole width instead of half in case array was on the edge
    %(and is therefore on the edge of the bbox)
    bbox_width_x = ceil(bounding_coords(4)./2)+1;
    bbox_width_y = ceil(bounding_coords(5)./2)+1;
    
    %also remove pixels from allowed region if they are near bright blobs
    if ~isempty(blob_locs)
        X = ones(p.im_Y,p.im_X).*[1:p.im_X];
        Y = ones(p.im_Y,p.im_X).*[1:p.im_Y]';
        blob_X = reshape(blob_locs(:,1),[1,1,size(blob_locs,1)]);
        blob_Y = reshape(blob_locs(:,2),[1,1,size(blob_locs,1)]);
        near_blob = min(sqrt((blob_X - X).^2 + (blob_Y - Y).^2),[],3) < p.blob_clear_radius;
        allowed = allowed & ~near_blob;
    end
    
    allowed(:,1:bbox_width_x) = 0;
    allowed(1:bbox_width_y,:) = 0;
    allowed(:,end-bbox_width_x:end) = 0;
    allowed(end-bbox_width_y:end,:) = 0;
    
    
    if sum(allowed,'all') < p.min_num_placements
        rand_quant_pixels = 0;
        rand_background_bw = 0;
        rand_start_i = 0;
        rand_end_i = 0;
        rand_start_j = 0;
        rand_end_j = 0;
        rand_start_k = 0;
        rand_end_k = 0;
        skip = true;
        return
    else
        skip = false;
    end

    %get random x and y coordinates for the random array placement, but leave z alone
    [max_i,max_j,max_k] = size(cells);
    bb_center = bounding_coords(1:2) + bounding_coords(4:5)/2;
    [rand_r,rand_c] = ind2sub([max_i,max_j],randsample(find(allowed),1));
    rand_Y = int16(voxels(:,2)-bb_center(2)+rand_r);
    rand_X = int16(voxels(:,1)-bb_center(1)+rand_c);
    
    try
        rand_idxs = sub2ind([max_i,max_j,max(p.slices)],rand_Y,rand_X,voxels(:,3));
    catch
        disp({fname, max_i ,max_j ,rand_r, rand_c, min(rand_X), max(rand_X), min(rand_Y), max(rand_Y),...
            min(voxels(:,3)), max(voxels(:,3)),...
            bounding_coords(1), bounding_coords(2) ,bounding_coords(3), bounding_coords(4), bounding_coords(5), bounding_coords(6)})
        save([out_path 'crash_report_allowed.mat'],'allowed')
        return
    end

    %get shifted bbox based on random coordinates
    se_size = p.background_dilate_size;
    rand_start_i = floor(max(1,bounding_coords(2) - bb_center(2) + rand_r - 2*se_size));
    rand_start_j = floor(max(1,bounding_coords(1) - bb_center(1) + rand_c - 2*se_size));
    rand_start_k = floor(max(1,bounding_coords(3) - 2*se_size));
    rand_end_i = ceil(min(max_i, bounding_coords(2) + bounding_coords(5) - bb_center(2) + rand_r + 2*se_size));
    rand_end_j = ceil(min(max_j, bounding_coords(1) + bounding_coords(4) - bb_center(1) + rand_c + 2*se_size));
    rand_end_k = ceil(min(max_k ,bounding_coords(3) + bounding_coords(6) + 2*se_size));

    rand_cells_bb = cells(rand_start_i:rand_end_i,rand_start_j:rand_end_j,rand_start_k:rand_end_k);
    
    %get array seg and ring around it for background, in bbox region
    rand_array_bw = zeros(size(cells));
    rand_array_bw(rand_idxs) = 1;
    rand_array_bw = rand_array_bw(rand_start_i:rand_end_i,rand_start_j:rand_end_j,rand_start_k:rand_end_k);
    
    [x,y,z] = ndgrid(-se_size:se_size);
    se1 = strel(sqrt(x.^2 + y.^2 + z.^2) <= se_size);
    rand_big_array = imdilate(rand_array_bw,se1);
    
    if in_cell
        rand_background_bw = rand_big_array & ~rand_array_bw & rand_cells_bb;
        rand_quant_pixels = rand_array_bw & rand_cells_bb;
    else
        rand_background_bw = rand_big_array & ~rand_array_bw & ~rand_cells_bb;
        rand_quant_pixels = rand_array_bw & ~rand_cells_bb;
    end
    
    %remove pixels in other channels' segmentations if requested
    for ch2 = p.discount_from_bkgnd
        ch2_name = p.channels{ch2};
        ch2_seg = seg.(ch2_name);
        rand_background_bw = rand_background_bw & ~ch2_seg(rand_start_i:rand_end_i,rand_start_j:rand_end_j,rand_start_k:rand_end_k);
    end
end

function parsave(p, out_path, row, col, var)
    save([out_path p.out_prefix 'analyzed_data_array_background/' row col '.mat'],'var','-v7.3')
end