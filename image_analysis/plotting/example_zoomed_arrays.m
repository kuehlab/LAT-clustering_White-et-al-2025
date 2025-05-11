%% Example array images
mkdir('figures/example_arrays_zoomed')

clearvars

load('array_stats_array_bkgnd.mat');
load('multiseg_params.mat')
load('../20230624_array_density_range_5-10min_and_std_cndl/GFP_counting_coefs.mat')
%%
pMHC_ratio = [1/3, 1/6, 1/12, 1/24, 1/48, 1/96, 1/192, 0];
pMHC_ratio_for_plot = [1/3, 1/6, 1/12, 1/24, 1/48, 1/96, 1/192, 1/384]; %last one is really 0
time = [5, 10];
time_names = {'5min', '10min',};
ratio_names = {'1:3', '1:6', '1:12', '1:24', '1:48', '1:96', '1:192', '0'};
rows = 'BCDE';

alpha = 0.1;
lines_cmap = colormap('lines');
parula_colormap = parula(11);
GFP_cut = 10^2.5;
min_n_slice = 3;
min_frac_in_cell = 0.99;
%%
colors = [
        1,1,1;
        1,0,1;
        0,1,0;
        1,0,0;
        0,1,1];
contrast = [
    0,0;
    100,8000;
    250,2000;
    100,1000;
    250,2000;
    ];

n_quant = 4;
GFP_quant_range = [0.33, 0.66, 1];
window = 12;

CD3_quant_low = [0.1,0.65];
CD3_quant_high = [0.35,0.9];

for r = 1:2
    for c = 1:8
        
        well_data = data{r,c};
        
        centroids = {well_data.centroid};
        centroids = cat(1,centroids{:});
        is_neigh = sqrt((centroids(:,1) - centroids(:,1)').^2 + (centroids(:,2) - centroids(:,2)').^2) < window*sqrt(2)*2;
        well_num = cellfun(@(x) str2double([num2str(double(x(1))), num2str(x(2)), num2str(x(4)), num2str(x(6))]), {well_data.fname});
        well_match = well_num == well_num';
        is_neigh = is_neigh & well_match;
        brightest_neigh = max(is_neigh.*repmat([well_data.GFP_med],length(well_data),1),[],1);
        n_neigh = sum(is_neigh,1);
        for i = 1:length(well_data)
            [well_data(i).n_neigh] = n_neigh(i);
            [well_data(i).brightest_neigh] = brightest_neigh(i);
        end
        
        in_cell = [well_data.frac_in_cell] > min_frac_in_cell;
        good_GFP = [well_data.GFP_med] < GFP_cut;
        good_slices = [well_data.n_slices] >= min_n_slice;
        N_GFP = ([well_data.GFP_tot] - coef_2(2))./coef_2(1);
        good_size = N_GFP > 10;
        good_bkgnd = [well_data.GFP_bkgnd] < 10^2.7;
        good_idx = in_cell & good_GFP & good_slices & good_size & good_bkgnd;
        
        well_data = well_data(good_idx);
        GFP_tot = [well_data.GFP_tot];
        neigh = [well_data.n_neigh];
        bright = [well_data.brightest_neigh];
        
        for d = GFP_quant_range
            
            figure('Position',[0,0,80*n_quant + 20,250])
            ha = tight_subplot(3,n_quant,0.02,0.001,0.001);
            
            low = quantile(GFP_tot,d-0.1);
            high = quantile(GFP_tot,d);
            
            for q = 1:length(CD3_quant_high)
                in_subset = (GFP_tot > low) & (GFP_tot < high) & (neigh <= 1) &(bright <= GFP_cut);

                subset = well_data(in_subset);                        
                CD3 = [subset.pCD3z_tot];
                subset = subset((CD3 > quantile(CD3,CD3_quant_low(q))) & (CD3 < quantile(CD3,CD3_quant_high(q))));
                subset = struct2table(subset);
                subset = sortrows(subset,'pLAT_med');
                subset = table2struct(subset);

                for direction = [-1,1]
                    if direction == -1
                        order = 1:length(subset);
                    else
                        order = length(subset):-1:1;
                    end
                    found = false;
                    for s = order
                        
                        fname = subset(s).fname;
                        load(['multiseg_corrected_data/' fname])
                        load(['multiseg_segmentations/' fname])

                        centroid = subset(s).centroid;
                        x_min = floor(max(1,centroid(1)-window));
                        y_min = floor(max(1,centroid(2)-window));
                        x_max = ceil(min(p.im_X,centroid(1)+window));
                        y_max = ceil(min(p.im_Y,centroid(2)+window));

                        array_bdry = bwboundaries(max(seg.GFP,[],3));
                        cell_bdry = bwboundaries(max(seg.CTV,[],3));

                        center = [mean([x_min,x_max],'all') mean([y_min,y_max],'all')];
                        is_target = cellfun(@(x) inpolygon(center(2),center(1),x(:,1),x(:,2)), array_bdry);

                        if sum(is_target) == 0
                            continue
                        end

                        target = array_bdry{is_target};
                        ellipse = fit_ellipse(target(:,1),target(:,2));
                        ellipticity = ellipse.long_axis/ellipse.short_axis;

                        if (quantile(current_im(4).im(y_min:y_max,x_min:x_max,:),0.95,'all') < contrast(4,2)) &&...
                                (quantile(current_im(5).im(y_min:y_max,x_min:x_max,:),0.95,'all') < contrast(5,2)) &&...
                                (ellipticity < 1.2)
                            
                            for ch = 3:5
                                
                                if ch == 3
                                   ch_num = 0;
                                elseif ch == 4
                                    ch_num = 2;
                                else
                                    ch_num = 1;
                                end
                                
                                axes(ha(n_quant*ch_num+2*(q-1)+(direction+1)/2+1))
                                im = max(current_im(ch).im,[],3);
                                im = imadjust(im./contrast(ch,2),[contrast(ch,1)./contrast(ch,2);1],[0;1]);
                                im = cat(3,im*colors(ch,1),im*colors(ch,2),im*colors(ch,3));
                                imshow(im)
                                hold on
        %                         for b = 1:length(cell_bdry)
        %                             bdry = cell_bdry{b};
        %                             if any(inpolygon(bdry(:,2),bdry(:,1),[x_min,x_min,x_max,x_max],[y_min,y_max,y_max,y_min]),'all')
        %                                 plot(gca,bdry(:,2), bdry(:,1), 'Color','white','LineWidth',2,'LineStyle','--');
        %                             end
        %                         end
                                for b = 1:length(array_bdry)
                                    bdry = array_bdry{b};
                                    if any(inpolygon(bdry(:,2),bdry(:,1),[x_min,x_min,x_max,x_max],[y_min,y_max,y_max,y_min]),'all')
                                        plot(gca,bdry(:,2), bdry(:,1), 'Color','white','LineWidth',1,'LineStyle','-');
                                    end
                                end
                                xlim([x_min,x_max])
                                ylim([y_min,y_max]) 
                            end
                            break
                        end
                    end
                end
            end
            
            export_fig(['figures/example_arrays_zoomed/' strrep(ratio_names{c},':','-') '_' num2str(time(r)) '_' num2str(d) '.png'])
            axis off
            close all
        end
    end
end