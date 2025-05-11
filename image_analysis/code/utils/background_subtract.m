function corrected = background_subtract(stacks,flatfield_correction,p)

    se = strel('disk',p.cellsize*2);

    for ch = 1:length(stacks)
        if ismember(ch,p.flatten)
            corrected(ch).im = zeros(size(stacks(ch).im));
            for z = 1:size(stacks(ch).im,3)
                corrected(ch).im(:,:,z) = imtophat(stacks(ch).im(:,:,z),se)./flatfield_correction(1).(strrep(p.channel_names{ch},' ','_'));
            end
        else
            corrected(ch).im = stacks(ch).im;
        end
    end
end