%Radial profile binning and plotting code - Code v1
%Written by Sumon Sahu - 10/04/22
%Last Updated by Sumon Sahu - 05/12/23

%Input - excel files with intensity values as y and normalized r as x. 
%This code makes equal width bins of normalized r values from 0 to 2.5. 
%The intensity values are put into different bins accordingly.
%Finally, this code calculates mean and standard deviation of intensities in each bin.


% close all windows, clear all variables from workspace
clear all;
clc;
close all;


% To select cells from excel, choose from one cell above the topmost left
% cell to right most bottom cell!
% Example excel data import
unmod_561 = xlsread("**.xlsx", 'Sheet1','*:*'); % does not include header
unmod_561_r = xlsread("**.xlsx", 'Sheet1','*:*'); % does not include header

%for plotting central intensity of the droplets use the command below
%droplet_center_intensity  = unmod_561(1,:)'

%% 
%total number of droplets in the variable 'unmod_561'
num_drops = sprintf('droplet number = %d',size(unmod_561,2))

% Define variables to store data
Yn_array = [];
Yn_array_mean = [];
Yn_array_std = [];

% Call function bin_data to calculate the binned y values
[Xedges, dx, Yn_array] = bin_data(unmod_561_r,unmod_561,30);

% Get bin mid-value from bin edges
X_mid = Xedges+(dx/2);
X_mid2 = X_mid(1:end-1); % representative x-values of bins (mid value)

% Calculate mean and standard deviation of y value excluding NaN values
Yn_array_mean = mean(Yn_array,2,'omitnan');
Yn_array_std = std(Yn_array,[],2,'omitnan');

% Get the output r and intensity (excluding origin)
r_norm = X_mid2(2:end)
Yn_array_mean_final = Yn_array_mean(2:end)
Yn_array_std_final = Yn_array_std(2:end)


figure
%errorbar(X_mid2, Yn_array_mean, Yn_array_std, 'k', 'LineWidth', 2)
errorbar(r_norm, Yn_array_mean_final, Yn_array_std_final, 'k', 'LineWidth', 2)
xlim([0 2.5])
ylim([0 30000])
xlabel('r_{norm}')
ylabel('Intensity')
set(gca,'FontSize', 18)
legend(num_drops)

%% function space
 
function [Xedges, dx, Yn_array] = bin_data(r_normvalues,y_values,norm_r_range)
    % store the binned y-data in this array
    Yn_array = [];
    
    % loop over the droplets of interest
    for loop=1:size(y_values,2)
        
        % Define NaN array to store data
        Yn=NaN(1,norm_r_range)
        
        % Get X-Y values for each droplets
        X = r_normvalues(:,loop);   
        Y = y_values(:,loop); 

        % Range of normalized r
        Xmin = 0; Xmax = 2.5;
        
        % no. of bins we want
        N = norm_r_range;  
        
        % bin width
        dx = (Xmax - Xmin)/ (N-1);
        
        % bin edges
        Xedges = Xmin - dx/2 : dx : Xmax + dx/2;
        Xedges = Xedges';

        % Get the map of what values would be in which bin
        [bins,E] = discretize(X,Xedges); 

        % Start putting values in bins.
        for i=1:size(X,1)
            if isnan(bins(i))==0
                Yn(bins(i)) = Y(i);
            else
                continue;
            end
        end
        
        %finally store each droplet values into main variable array
        Yn_array = horzcat(Yn_array,Yn');

    end
end 



