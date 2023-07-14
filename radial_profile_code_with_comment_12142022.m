%Radial profile plot code - Code v1
%Written by Sumon Sahu - 10/04/22
%Last Updated by Sumon Sahu - 05/12/23

%Input - Any image stack in .tif format. 
%This code picks up round droplets and calculates the radial intensity
%distribution for droplet diameter between 2.5 um and 4 um. Based on the
%diameter, this code also calculates normalized r array for individual droplets.
%This code uses otsu's threshold to make binary images

%These settings are for YSB microscope 100x mag
%Filemanes should be 'Stack%d_%d.tif' format, %d = integer.

% close all windows, clear all variables from workspace
clear all;
clc;
close all;

% check scale manually in FIJI for the data
scale = 0.166 % for 512*512 YSB image


%Path to your .tif stack
path = './folder/';
file = 'Stack2_561.tif';
file_path=[path file];

num = regexp(file, '\d+', 'match');

%create filenames to save your data
save_name1= [path 'profile_array_Stack' num{1} '_' num{2} '.mat']; % saves radially averaged intensity profiles of the droplets
save_name2= [path 'norm_r_array_Stack' num{1} '_' num{2} '.mat']; % saves normalized r array


%get info about the dataset stack
info = imfinfo(file_path);                                              
num_images = numel(info);

%set frame numbers for running the dataset stack
start_frame =1;
end_frame =num_images;

%create empty arrays to store the results
profile_array=[];
norm_r_array=[];

%define the mask size for radial averaging
distance_window = 30; %for YSB images

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
%   imThresh = imNorm > 0.3;


    
    % Ignore smaller areas with bwareafilter
    % Adjust lower limit according to the scale. 
    % For 512*512 YSB image, lower limit is 2.5

    imThresh2 = bwareafilt(imThresh, [2.5 50000]);

    %label droplets in the binary image
    L = bwlabel(imThresh2);
    
    % use regionprops function to get all of the information about the binarized droplets          
    s = regionprops('table',imThresh2, 'Area', 'PixelIdxList', 'MajorAxisLength', 'MinorAxisLength', 'Centroid','Orientation','Eccentricity');
    
    diameter = [];
    % estimate diameter values from 'MajorAxisLength' and 'MinorAxisLength'
    diameter = scale*(([s.MajorAxisLength]+[s.MinorAxisLength])/2); %scaled with scale bar 0.166 um/px
    
    % use this to add diameter values to the s table
    s.Diameter = diameter;
    
    %sort the table s according to the area values in descending order
    [s1,index] = sortrows(s,'Area','descend'); 
    
%=======================Uncomment if you want to see the droplet number assignment =====================    
%     figure
%     imshow(imThresh2)
%     hold on
%     for k = 1:size(s1,1)
%         c(1) = s1.Centroid(k,1);
%         c(2) = s1.Centroid(k,2);
%         text(c(1), c(2), sprintf('%d', k), ...
%             'HorizontalAlignment', 'center', ...
%             'VerticalAlignment', 'middle', 'Color', 'r');
%     end
%     hold off
%=======================================================================================================    

    r=(1:distance_window)*scale; %define radial distance array
    
    % loop through the droplets in a single image
    for i=1:length(diameter) % length(diameter) is the number of droplets in the image
        
        %reintroduce variables to store numbers to arrays. This clears out
        %the previous array for previous droplet
        profile1=[];
        norm_r=[];
        
        % filter between the droplet size with diameter 2.5 um and 4 um.
        % Choose round-ish droplets with eccentricity <=0.5
        if (s1.Diameter(i) >= 2.5) & (s1.Diameter(i) <= 4) & (s1.Eccentricity(i) <= 0.5)
            
            % get the centroids and round them up
            centers1(1) = round(s1.Centroid(i,1));
            centers1(2) = round(s1.Centroid(i,2));
            
            %get the radial average intensity profile from the 'radialAverage' function
            profile1 = radialAverage(im, centers1(1), centers1(2), distance_window);
            
            
            % normalize the r with radius of each droplet
            norm_r = 2*r/s1.Diameter(i);
            
            % append the values from this loop to the main variable array
            profile_array = horzcat(profile_array,profile1');
            norm_r_array = horzcat(norm_r_array,norm_r');
        else
            continue;
        end    
    end
    
    fprintf('Image number %d is running!\n',loop)      

end

%save variables to the files
save(save_name1,'profile_array')
save(save_name2,'norm_r_array')
%%

%Make plot for intensity distribution, intensity vs r plot
figure
for loop=1:size(profile_array,2) 
     plot(r,profile_array(:,loop),'o-')
     hold on
end
hold off
xlim([0 4])
ylim([0 60000])
xlabel('r (\mum)')
ylabel('Intensity')
set(gca,'FontSize', 18)

%Make plot for intensity distribution, intensity vs r_norm plot
figure
for loop=1:size(profile_array,2) 
    plot(norm_r_array(:,loop),profile_array(:,loop),'o-')
     hold on
end
hold off
xlim([0 4])
ylim([0 60000])
xlabel('r_{norm}')
ylabel('Intensity')
set(gca,'FontSize', 18)

%% function space

function profile = radialAverage(G, cx, cy, rmax)
    % computes the radial average of the image G around the cx,cy point in
    % XY co-ordinates    
    % rmax is the mask distance from origin in (r,theta)
    profile = [];
    
    %get the size of the image
    [a,b] = size(G);
    
    % form a XY grid at the center of the selected droplet
    [COLUMN, ROW] = meshgrid( (1:b)-cx, (1:a)-cy);
    
    % convert ROW-COLUMN to (R,theta) co-ordinate
    R = sqrt(ROW.^2 + COLUMN.^2);
    
    % for each i, make rings with 1 px width to calculate intensity
    for i = 1:rmax % radius of the circle
        % create the ring mask
        mask = (i-1<R & R<i+1);
        
        %get the values from the mask
        values = [];
        values = G(mask); % without smooth

        % get the mean intensity from the radial ring mask
        profile(i) = mean(values(:));
    end
end

