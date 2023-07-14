%Droplet size Code v1
%Written by Sumon Sahu - 08/22/22
%Last Updated by Sumon Sahu - 05/12/23

%This code measures the droplet diameter values from an image stack in .tiff. It
%saves the variable 'diameter_array' into 'Stack_diameters.mat' which has all diameter values.
%The diameter values from each frame can be found in variable 'diameter'. 
%It also saves number of droplets per image into 'Stack_num_drops.mat'.

%These settings are for YSB microscope 100x mag
%Filemanes should be 'Stack%d.tif' format, %d = integer.

% close all windows, clear all variables from workspcae
clear all;
clc;
close all;

%This code uses otsu's threshold to make binary images

% check scale manually for each day data
scale = 0.166; % pixel size for 512*512 image only 

%Path to your .tif stack
path = './folder/';
file = 'Stack11.tif';
file_path = [path file];

%extract number from filename
num = regexp(file, '\d+', 'match');

%create filenames to save your data
save_name1 = [path 'Stack' num{1} '_diameters.mat']; % saves diameter of the droplets
save_name2 = [path 'Stack' num{1} '_num_drops.mat'];% saves number of the droplets per image


%get info about the dataset stack
info = imfinfo(file_path);                                              % info about file
num_images = numel(info);

%set frame numbers for running the dataset stack
start_frame =1;
end_frame = num_images;

%create empty arrays to store the results
diameter_array = [];
diameter = {};

% loop through all images in  the image stack
for loop = start_frame:end_frame
    
        % read the image
        im = imread(file_path,loop);
        image_size = size(im);
                
        % Clear droplets at boundary by imclearborder 
        im2 = imclearborder(im);
               
        % apply median filter to reduce random noise
        imFilt = medfilt2(im2, [4 4]);
        
        % scales the intensity range between 0 and 1
        imNorm = mat2gray(imFilt); 
        
        %Find intensity threshold using Otsu's global thresholding method
         thresh = graythresh(imNorm); 
         imThresh = imNorm > thresh;
        
        
        %alternatively: use manual threshold
%        imThresh = imNorm > 0.3;
       
       
        % Ignore smaller areas with bwareafilter
        % Adjust lower limit according to the scale. 
        % For 512*512 image, lower limit is 2.5
        imThresh2 = bwareafilt(imThresh, [2.5 50000]); 

% Uncomment this portion to see which droplets are selected and the section of droplets that is selected after thresholding        
%======================check the impact of threshold value==================================
      % location of pixels   
%         iii=find(imThresh2);
%         [y,x] = ind2sub(size(imThresh2),iii);
% 
%         % show raw image  
%         figure
%         imshow(imNorm,[])
%         impixelinfo 
%         hold on        
%         p=scatter(x,y,'ro','filled')        
%         hold off
%         alpha(p,.5)
%         pause(0.5)
%======================check the impact of threshold value==================================


        % use regionprops function to get all of the information about the binarized droplets          
        s = regionprops('table',imThresh2, 'Area', 'PixelIdxList', 'MajorAxisLength', 'MinorAxisLength', 'Centroid','Orientation','Eccentricity');
        
        % remove elonngated objects using eccentricity parameter.
        % eccentricity ~ 0 is circle and eccentricity ~ 1 is straight line
        toDelete = s.Eccentricity>0.75; % use 0.75 for diameter analysis. Although, choosing no eccentricity threshold does not change the statistical distribution significantly.
        s(toDelete,:) = [];

        % estimate diameter values from 'MajorAxisLength' and 'MinorAxisLength'
        diameter{loop} = scale*(([s.MajorAxisLength]+[s.MinorAxisLength])/2); 
        
        % store the number of droplets from this image to this variable
        num_drops(loop) = numel(diameter{loop});
        
% Uncomment this section to see the ellipse fit to measure the diameter         
%======================check the ellipse fit==================================        
%         figure
%         imshow(imNorm,'InitialMagnification','fit')
% 
%         t = linspace(0,2*pi,50);
% 
%         hold on
%         for k = 1:size(s,1)
%             a = s.MajorAxisLength(k)/2;
%             b = s.MinorAxisLength(k)/2;
%             Xc = s.Centroid(k,1);
%             Yc = s.Centroid(k,2);
%             phi = deg2rad(-s.Orientation(k));
%             x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
%             y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
%             plot(x,y,'r','Linewidth',2)
%             plot(Xc,Yc,'r*')
%         end
%         hold off        
%======================check the ellipse fit==================================
        
     
    % add diameter values from this current image to the diamater_array which will diameter values from all images finally 
     diameter_array = vertcat(diameter_array, diameter{loop});
             
    % show in the command window which image is running
    fprintf('Image number %d is running!\n',loop)  
            
end
%take transpose of the array to report
num_dropt = num_drops';

% save variables to the files
save(save_name1,'diameter_array')
save(save_name2,'num_dropt')
