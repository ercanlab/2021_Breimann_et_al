
%% 0 choose data and results directory

% If just want to quick run through a data set run the script: run_script_through_dataset.m
% Only uncomment the following lines if you use this script without the
% run-script
% 
%disp('choose data_tif directory')
%data_dir = uigetdir(pwd, 'choose data_tif directory')
% 
%disp('Choose matlap_results directory')
%results_dir = uigetdir(pwd,'Choose matlap_results directory') %choose output directory

%disp('load pre stack')
%[fnamepreBase, fnamepreFolder] = uigetfile(fullfile(data_dir,'*.tiff'),'load pre stack'); %adds folder path to the filename
%fnamepre = fullfile(fnamepreFolder, fnamepreBase);
%filename_new = char(...
%                 extractBefore(fnamepreBase,'_FRAP_P'))

% disp('load post stack')
% [fnamepostBase, fnamepostFolder] = uigetfile(fullfile(data_dir,'*.tiff'),'load post stack'); %adds folder path to the filename
% fnamepost = fullfile(fnamepostFolder, fnamepostBase);

%% 1st step, read in tiff files 

%clearvars -except data_dir results_dir files
%   

disp('process pre stack')
[stack, img_read] = tiffread2(fnamepre); 
Fpre = cat(3, stack.data);
timePre = timesteps(fnamepre); %correct timesteps between frames [s]
image = Tiff(fnamepre, 'r');
y = strsplit(getTag(image, 'ImageDescription'));
    for i = y
        if contains(i, 'finterval')
            l = char(i);
        end
    end
timestep = strsplit(l, '=');
timestep = round(str2double(timestep(2)), 2);


disp('process post stack')
[stack, img_read] = tiffread2(fnamepost); %relative file path
Fpost = cat(3, stack.data);
timePos = timesteps(fnamepost); %correct timesteps between frames [s]


%% 2nd step, filters images, highlights the nucleus by edge preserving filtering and ROI drawing

clear timestepsPre timestepsPost timecombined
disp('Second step, filters images and check which nucleus got bleached (for ROI drawing)')
% get number of time points
dim3_pre = size(Fpre, 3);
dim3_post = size(Fpost, 3);

% get timesteps between each frame in Pre or Post-bleaching session [s]

timestepsPre = timestepsArray(dim3_pre, timePre);
timestepsPost = timestepsArray(dim3_post, timePos);

tmp_timesteps = timestepsPost+timestepsPre(length(timestepsPre));
timecombined = [timestepsPre tmp_timesteps];
timecombined = timecombined';


% filter .... just a bit
Fpre_filt = imguidedfilter(Fpre); 
Fpost_filt = imguidedfilter(Fpost);

% compare pre and post filter
%matVis(cat(2,Fpost, Fpost_filt))


%% 2.1. Cropping -> choose ROI by drawing
disp('cropping - choose ROI by drawing')
tmp_im = mean(Fpost_filt,3);
imagesc(tmp_im); %opens the image-stack with color scale overlay
test_roi = roipoly; %activates ROI drawing on opened image and saves ROI to variable

test_roi = uint16(test_roi);
close
imagesc(test_roi) %colormap destinguish ROI mask and background

clear Fpre_filt_crop Fpost_filt_crop

for i = 1:dim3_pre
    tmp = Fpre_filt(:,:,i);
    Fpre_filt_crop(:,:,i) = uint16(tmp).*test_roi; %combines the filtered image with the ROI mask and new grey scale
end
for i = 1:dim3_post
    tmp = Fpost_filt(:,:,i);
    Fpost_filt_crop(:,:,i) = uint16(tmp).*test_roi;
end

Fpre_filt = Fpre_filt_crop;
Fpost_filt = Fpost_filt_crop;

%take a look at cropped part
%cropped images are always used
%matVis(Fpost_filt_crop)

%take the first slice of the pre-bleach stack and save as csv for further
%histogram plotting of the grey values

first_preBleach_img = Fpre_filt_crop(:,:,1);
filepath_1stBleach = fullfile(results_dir,strcat(filename_new, '_greyValues_1st_preBleach.csv'));
csvwrite(filepath_1stBleach ,first_preBleach_img);

close 

%% Step 3, pre image drift correction 

disp('Step 3, pre image drift correction')

clear Fim_reg Greg Fpre_filt_reg
dim3_pre = size(Fpre,3);
reg_temp = mean(Fpre_filt(:,:,dim3_pre-2:dim3_pre),3);

for i = 1:dim3_pre
    tmp = double(Fpre_filt(:,:,i));
    [output, Greg] = dftregistration(fft2(reg_temp),fft2(tmp)); %fft2 does a fourier transformation
    Fpre_filt_reg(:,:,i) = abs(ifft2(Greg));
end

%compare pre and pre-drift corrected
%matVis(cat(2,Fpre_filt, Fpre_filt_reg));

Fpre_filt = Fpre_filt_reg; 
%% Step 4, post image drift correction
    
disp('Step 4, post image drift correction')

clear Fim_reg Greg Fpost_filt_reg

reg_temp = mean(Fpost_filt(:,:,1:3),3);
for i = 1:dim3_post
    tmp = double(Fpost_filt(:,:,i));
    [output, Greg] = dftregistration(fft2(reg_temp),fft2(tmp));
    Fpost_filt_reg(:,:,i) = abs(ifft2(Greg));
end

%compare post and post-driftcorrected
%matVis(cat(2,Fpost_filt, Fpost_filt_reg));

Fpost_filt = Fpost_filt_reg; 

%% Step 5, Save graph as mask.tiff in a new folder
disp('Step 5, Save graph as mask.tiff in a new folder')

% arbitrary threshold for bleaching ROI 0 to 1, change till happy with size
bl_thresh = 0.6; 

Fpre_mean = mean(Fpre_filt,3);
Fpre_mn = Fpre_mean./max(Fpre_mean(:)); %mean_norm - normalization of the background
Fpost_mean = mean(Fpost_filt(:,:,1:5),3); %the mean of the post-bleach images from the first 5 images

Fbl_filt = Fpre_mean - Fpost_mean; %difference = bleachpoint

%pre, post images as well as the determined difference at the bleach point
%are displayed with a colormap 
cmax =  max(Fpre_mean(:))*1.1;
subplot(2,3,1);
imagesc(Fpre_mean); colormap('hot'); caxis([0 cmax]); axis image; axis off;
title('F_{pre}');
subplot(2,3,2);
imagesc(Fpost_mean); colormap('hot'); caxis([0 cmax]); axis image; axis off;
title('F_{post}');
subplot(2,3,3);
imagesc(Fbl_filt); colormap('hot'); caxis([0 cmax]); axis image; axis off;
title('F_{diff}');

%displaying the ROI mask
bw = graythresh(Fpre_mn); %"Otsu's method"
Fnuc_mask = zeros(size(Fpre_mn));
Fnuc_mask(Fpre_mn>bw) = 1;
Fnuc_mask = logical(Fnuc_mask);
se = strel('disk',2);
Fnuc_mask = imerode(Fnuc_mask,se); % erode an dilate are responsible for drawing the right shape around the detected bleachpoint
% changed this from 'square' 10 to 'disk' and smaller number
se = strel('square',10);
Fnuc_mask = imdilate(Fnuc_mask,se);
subplot(2,3,4)
imagesc(Fnuc_mask); axis image; axis off;
title('mask_{nuc}');

ROI = zeros(size(Fbl_filt));
ROI(Fbl_filt > max(Fbl_filt(:))*bl_thresh)=1; 
se = strel('disk',3);
tmp_BW1 = imerode(ROI,se);
se = strel('disk',5);
tmp_BW2 = imdilate(tmp_BW1,se);
subplot(2,3,5)
imagesc(tmp_BW2); axis image; axis off;
title('mask_b_l');

%filename1 = strcat(strtok(fullfile(results_dir,fnamepreBase),'.'),'_mask')  %save in extra results dir
%print('-f1',(filename1), '-dtiff')

%% 
% request if mask-file is good enough, if not adjust the threshold and do
% it again

correct_bleach_threshold = 0;
while ~correct_bleach_threshold 
    mask_request = menu('Does the bleach point mask looks ok?','Yes','No');
if  mask_request == 1
    correct_bleach_threshold = 1;
else
    prompt = 'Enter a refined threshold (default: 0.6, up->more stringent): ';
    new_bl_thresh = input(prompt);
    
    bl_thresh = new_bl_thresh; 
    Fpre_mean = mean(Fpre_filt,3);
    Fpre_mn = Fpre_mean./max(Fpre_mean(:));
    Fpost_mean = mean(Fpost_filt(:,:,1:5),3);

    Fbl_filt = Fpre_mean - Fpost_mean;
    cmax =  max(Fpre_mean(:))*1.1;
    subplot(2,3,1);
    imagesc(Fpre_mean); colormap('hot'); caxis([0 cmax]); axis image; axis off;
    title('F_{pre}');
    subplot(2,3,2);
    imagesc(Fpost_mean); colormap('hot'); caxis([0 cmax]); axis image; axis off;
    title('F_{post}');
    subplot(2,3,3);
    imagesc(Fbl_filt); colormap('hot'); caxis([0 cmax]); axis image; axis off;
    title('F_{diff}');

    bw = graythresh(Fpre_mn);
    Fnuc_mask = zeros(size(Fpre_mn));
    Fnuc_mask(Fpre_mn>bw) = 1;
    Fnuc_mask = logical(Fnuc_mask);
    se = strel('disk',3);
    Fnuc_mask = imerode(Fnuc_mask,se);
    se = strel('disk',5);
    Fnuc_mask = imdilate(Fnuc_mask,se);
    subplot(2,3,4)
    imagesc(Fnuc_mask); axis image; axis off;
    title('mask_{nuc}');

    ROI = zeros(size(Fbl_filt));
    ROI(Fbl_filt > max(Fbl_filt(:))*bl_thresh)=1;
    se = strel('disk',2);
    tmp_BW1 = imerode(ROI,se);
    se = strel('square',10);
    tmp_BW2 = imdilate(tmp_BW1,se);
    subplot(2,3,5)
    imagesc(tmp_BW2); axis image; axis off;
    title('mask_{bl}');

end
end

filename1 = strcat(strtok(fullfile(results_dir,filename_new),'.'),'_mask')  %save in extra results dir 
print('-f1',(filename1), '-dtiff')

close 

%% Step 6, exponential fitting,  Save graph as rec1.tiff in a new folder
disp('Step 6, exponential fitting,  Save graph as rec1.tiff in a new folder')
Fcomp = (cat(3,Fpre_filt,Fpost_filt));
for i = 1:size(Fcomp,3)
    tmp = Fcomp(:,:,i);
    Fbl_time(i) = mean(tmp(logical(tmp_BW2)));
    Fnuc_time(i) = mean(tmp(logical(Fnuc_mask)));
end

% fit prebleach of whole nucleus

time = dim3_pre+1:dim3_post;
tmp = Fnuc_time(dim3_pre+1:dim3_post);


[xData, yData] = prepareCurveData( time, tmp );
ft = fittype( 'exp1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
[fitresult, gof] = fit( xData, yData, ft, opts );

Facq = fitresult(1:length(Fnuc_time));


subplot(2,2,1);
plot(timecombined, Fnuc_time); hold on; 
plot(timecombined, Facq);
box off; xlabel('time [s]'); ylabel('Intensity');

subplot(2,2,3);
Fnuc_norm = Fnuc_time'./Facq; 
Fnuc_norm = Fnuc_norm./mean(Fnuc_norm(dim3_pre-5:dim3_pre)); 
plot(timecombined, Fnuc_norm);
ylim([0 1.1]);
box off; ylim([0 1.5]); 
xlabel('time [s]'); ylabel('Norm. Int.');

totblfrac = mean(Fnuc_norm(length(Fnuc_norm)-dim3_pre:length(Fnuc_norm)));

subplot(2,2,2);
plot(timecombined, Fbl_time); 
box off; 
xlabel('time [s]'); ylabel('Intensity');

subplot(2,2,4);
% tmp = Fbl_time./mean(Fbl_time(1:2));
tmp = Fbl_time./mean(Fbl_time(dim3_pre-5:dim3_pre));
Fbl_norm = tmp'./(Facq./Facq(1));
plot(timecombined, Fbl_norm); 
box off; ylim([0 1.5]); 
xlabel('time [s]'); ylabel('Norm. Int.');

filename2 = strcat(strtok(fullfile(results_dir,filename_new),'.'),'_rec1')  
print('-f1',(filename2), '-dtiff')

close


%% Step 7, recovery fitting, Save rec2
disp('Step 7 recovery fitting,  Save graph as rec2.tiff in a new folder')
tmp = Fbl_norm(dim3_pre+1:length(Fbl_norm));
time = 1:length(tmp);
time = time.*timePos;

ft = fittype( '1-a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.72828324812449 0.146626886594027 0.83559941025035];

[fitresult, gof] = fit( time', tmp, ft, opts);

% calculate recovery half time
fit_var = fitresult(time); 
fit_half_amp = (fit_var(length(fit_var))-fit_var(1))/2+fit_var(1);
fit_half_amp_2 = ones(size(fit_var));
fit_half_amp_2 = fit_half_amp_2.*fit_half_amp;
fit_50_diff = abs(fit_var-fit_half_amp);
t_half = timestepsPost(find(fit_50_diff==min(fit_50_diff)));


plot(timecombined, Fnuc_norm,'.r'); hold on;
plot(timecombined, Fbl_norm,'.k'); 
plot(tmp_timesteps, fitresult(time));
ylim([0 1.1]);
line([0 timecombined(length(timecombined))],[totblfrac totblfrac],'Color',[0.6 0.6 0.6],'LineStyle','--');
line([0 timecombined(length(timecombined))],[Fbl_norm(dim3_pre+1) Fbl_norm(dim3_pre+1)],'Color',[0.6 0.6 0.6],'LineStyle','--');
ylabel('norm. Int. (a.u.)');
xlabel('time (s)');
box off;

legend('whole nucleus','bleached ROI','recovery fit');
legend('boxoff');


f_mo = (fitresult(time(length(time)))- Fbl_norm(dim3_pre+1)) /(totblfrac-Fbl_norm(dim3_pre+1))
f_im = (totblfrac-fitresult(time(length(time)))) /(totblfrac-Fbl_norm(dim3_pre+1))

text(timecombined(length(timecombined)),fitresult(timecombined(length(timecombined))),['f_m_o = ',num2str(f_mo,3),'']);
text(timecombined(length(timecombined)),(fitresult(timecombined(length(timecombined)))+totblfrac)/2,['f_i_m = ',num2str(f_im,3),'']);
text(dim3_pre+5,(fitresult(timecombined(length(timecombined)))+Fbl_norm(dim3_pre+1))/2.5,['\tau = ',num2str(1/fitresult.b,3),' s']);
text(dim3_pre+5,(fitresult(timecombined(length(timecombined)))+Fbl_norm(dim3_pre+1))/2.5-0.1,['t_0_._5 = ',num2str(t_half),' s']);

tau = 1/fitresult.b;



filename3 = strcat(strtok(fullfile(results_dir,filename_new),'.'),'_rec2')  
print('-f1',(filename3), '-dtiff')
close

%%%%%%%%



%% Step 8 get t-half value from raw data (not fit like tau) 

% find the position (index) of the minimum intensity value
min_ind = find(Fbl_norm == min(Fbl_norm)); 

% define F_rec as the recovery 
Frec = Fbl_norm(min_ind:length(Fbl_norm));

% lets 'smoothen' the recovery trace a bit to make sure we only cross the
% half recovery point once.
Frec_raw = Frec; %keep raw trace for comparison
Frec = smooth(Frec,20); % moving average of 20 points


Frec_max = mean(Frec(length(Frec)-20:length(Frec))); % max recovery is ave of last 20 points
Frec_min = (Frec(1)); % min value is the first point
Frec_half = Frec_min + 0.5*(Frec_max - Frec_min); % half recovery level


% lets look at what we have
plot(Frec_raw); hold on; plot(Frec); 
line([0 length(Frec)],[Frec_max Frec_max],'Color',[0.8 0.8 0.8],'LineStyle','--');
line([0 length(Frec)],[Frec_min Frec_min],'Color',[0.8 0.8 0.8],'LineStyle','--');
line([0 length(Frec)],[Frec_half Frec_half],'Color',[0.8 0.8 0.8],'LineStyle','--');

% rather than just a scalar (single value) lets make this a vector or a
% function of time 
Frec_half_v = ones(size(Frec)).*Frec_half;
Frec_ms = abs(Frec - Frec_half_v); 
Frec_half_idx = find(Frec_ms == min(Frec_ms)); %find the time (frame #) when Frec_ms is min
Frec_half_time = Frec_half_idx .* timestep % absolute time is frame # x timestep

line([Frec_half_idx Frec_half_idx],[Frec_min Frec_half],'Color',[0.8 0.8 0.8],'LineStyle','--');
text([Frec_half_idx+15 Frec_half_idx+15],[Frec_half Frec_half],['t_0_._5 =',num2str(Frec_half_idx),' frames']);

xlabel('time (# frames)'); 
ylabel('F_r_e_c_o_v_e_r_y');
legend('F_r_a_w', 'F_a_v_e _2_0','Orientation','horizontal');
legend('boxoff');
box off;


filename4 = strcat(strtok(fullfile(results_dir,filename_new),'.'),'_thalf')  
print('-f1',(filename4), '-dtiff')

close


%% Step 9 export workspace for reanalysis

%20200213_vic - quick and dirty way - creates reduntant data, workspace not
%completely filtered for needed and not needed variables 

%define filename for saved roi_workspace
roi_workspace = join(strcat(strtok(fullfile(results_dir,filename_new),'.'), '_roi.mat'));
%save('test', '-regexp', '^(?!(image|Facq|Fbl_filt|Fpl_time|Fcomp|file|files|fit_50_diff|fit_half_amp_2|fit_var|fitresult|)$).')
save(roi_workspace, '-regexp', '^(?!(files|filefolder|filename|filename_new|results_dir|data_dir|test_roi|image|Fpost|Fpre|Fpost_filt_crop|Fpost_filt_reg|Fpre_filt_crop|Fpre_filt_reg|first_preBleach_img)$).')


%% Step 10 export values to txt and csv formats


%export gof values to a csv file 

gof_filename = join(strcat(strtok(fullfile(results_dir,filename_new),'.'), '_gof.csv'));
struct2csv(gof,gof_filename)

%export values for Python analysis
fileID = fopen(strcat(strtok(fullfile(results_dir,filename_new),'.'),'_pyan.txt'),'wt');
fprintf(fileID, '%f \n', tau);
fprintf(fileID, '%f \n',f_mo);
fprintf(fileID, '%f \n',f_im);
fprintf(fileID, '%f \n',Fbl_norm);
fprintf(fileID, '%f \n',timecombined)
fclose(fileID);




%quick and dirty version only numbers under each other

fileID_t_half = fopen(strcat(strtok(fullfile(results_dir,filename_new),'.'),'_t_half_values.txt'),'wt');
fprintf(fileID_t_half, '%f \n', t_half);
fclose(fileID_t_half);




fileID_t_half_2 = fopen(strcat(strtok(fullfile(results_dir,filename_new),'.'),'_t_half_values_2.txt'),'wt');
fprintf(fileID_t_half_2, '%f \n', Frec_half_time);
fclose(fileID_t_half_2);



%% Step 10 gather t_half values
% only needed to be uncommented if you analyse multiple files without the
% runscript


%I have a python function I can call here
% create python function in some python editor and add it to the Matlab path: (in this case gather_t_half.py)

if count(py.sys.path,'') == 0 % for some reason this has to be done when the matlab src directory is the pwd
   insert(py.sys.path,int32(0),'');
end
%py.gather_t_half.gather_t_half(results_dir);





