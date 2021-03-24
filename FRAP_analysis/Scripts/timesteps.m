function timestep = timestep(fnamepre)
image = Tiff(fnamepre, 'r');
y = strsplit(getTag(image, 'ImageDescription'));
    for i = y
        if contains(i, 'finterval')
            l = char(i);
        end
    end
timestep = strsplit(l, '=');
timestep = round(str2double(timestep(2)), 2);
end





