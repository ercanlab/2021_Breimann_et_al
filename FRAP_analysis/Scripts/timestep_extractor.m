t = uigetfile('.tif', 'Choose');
image = Tiff(t, 'r');
des = strsplit(getTag(image, 'ImageDescription'));
for i = des
    if contains(i, 'finterval')
        l = char(i);
        timesteps = strsplit(l, '=');
        timesteps = str2double(timesteps(2));
    end
end



