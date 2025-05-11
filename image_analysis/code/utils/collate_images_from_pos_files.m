function collate_images_from_pos_files(p)
    disp('Preprocessing...')
    
    if ~isfield(p, 'in_path')
        in_path = [p.path 'raw_data/'];
        out_path = p.path;
    else
        in_path = [p.in_path 'raw_data/'];
        out_path = p.out_path;
    end
    
    if ~isfield(p, 'out_prefix')
        p.out_prefix = '';
    end
    
    if ~isfield(p, 'start_pos')
        p.start_pos = 1;
    end 
    
    mkdir([out_path p.out_prefix 'corrected_data'])
    
    if ~isfield(p, 'plate_prefix')
        p.plate_prefix = cell(size(p.pos_files));
    end
    
    if ~isfield(p, 'flatfield_correction')
        for ch = 1:length(p.channel_names)
            flatfield_correction(1).(p.channel_names{ch}) = ones([p.im_Y,p.im_X,length(p.slices)]);
        end
    else
        load(p.flatfield_correction,'flatfield_correction')
    end
    
    for f = 1:length(p.pos_files)
        fname = p.pos_files{f};
        lines = readcell([in_path fname '.pos'],'FileType','text');
        
        if iscell(p.basename)
            basename = p.basename{f};
        else
            basename = p.basename;
        end
        
        parfor l = 1:length(lines) %parallel
            well_pos = [p.plate_prefix{f} lines{l,1}];
            
            current_im = struct();
            for ch = 1:length(p.channels)
                current_im(ch).im = zeros(1080,1080,max(p.slices));
                
                for z = p.slices
                    z_im = double(imread([in_path fname '/' basename '_w' num2str(ch) p.channel_names{ch} '_s' num2str(l + p.start_pos - 1) '.TIF'],z));
                    current_im(ch).im(:,:,z) = z_im;
                end
            end
            current_im = background_subtract(current_im,flatfield_correction,p);
            save_im(current_im,[out_path p.out_prefix 'corrected_data/' well_pos '.mat'])
            disp(well_pos)
        end
    end
end

function save_im(current_im,fname)
    save(fname,'current_im','-v7.3')
end