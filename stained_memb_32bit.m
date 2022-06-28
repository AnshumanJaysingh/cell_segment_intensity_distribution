close all
clear all
clc

currentfolder=string(pwd)
srcFile = dir(currentfolder+'/*.tif');
image_folder = currentfolder+'/'

golgi_mean=[];
memb_mean=[];
file_name_arr=[];

%% SET THE HISTOGRAM PLOTTING INFORMATION BEFORE RUNNING THE CODE
hist_min = 1 %lower end of histogram
hist_max = 100000 %higher end of histogram
binw = 10 %Bin width for histogram

%%
for n=1:length(srcFile)        
        FileName = strcat(image_folder,srcFile(n).name);
        file_name_arr=[file_name_arr, string(srcFile(n).name)];
        path=strcat(image_folder, srcFile(n).name);

        % read image and split channels
        I1=imread(FileName,1); % channel 1 = HRAS
        I1=uint32(I1);
        figure
        imshow(I1,[])
        I1=medfilt2(I1);
        imshow(I1,[])
        close

        I2=imread(FileName,2); % channel 2 = Golgi
        I2=uint32(I2);
        imshow(I2,[])
        I2=medfilt2(I2);
        imshow(I2,[])
        close
     
        I3=imread(FileName,2); % channel 3 = Membrane
        I3=uint32(I3);
        imshow(I3,[])
        I3=medfilt2(I3);
        imshow(I3,[])
        close
        
%% GOLGI SEGMENTATION
        seg_golgi = imbinarize(I2,'adaptive'); %global adaptive thresholding using Otsu (GOLGI CHANNEL)
        
        figure;
        imshowpair(I2, seg_golgi, 'montage')
        close
        
        golgi = immultiply(uint32(seg_golgi),I1);


%% MEMBRANE SEGMENTATION

        seg_memb = imbinarize(I3,'adaptive'); %global adaptive thresholding using Otsu (MEMB CHANNEL)
        
        figure;
        imshowpair(I3, seg_memb, 'montage')
        close
        
        memb = immultiply(uint32(seg_memb),I1);

%%  MEMBRANE HISTOGRAM 

        figure 
        h_memb=histogram(memb,'BinWidth',binw, 'BinLimits',[hist_min, hist_max]); %plot histogram of membrane
        title(string(srcFile(n).name) + ' histogram of membrane')
        data_memb = h_memb.Values;
        counts_memb = data_memb;
        centerBinGrayLevels_memb = (h_memb.BinEdges(1:end-1) + h_memb.BinEdges(2:end)) / 2;
        meanBinnedGrayLevel_memb = sum(centerBinGrayLevels_memb .* counts_memb) / sum(counts_memb);
        memb_mean=[memb_mean,meanBinnedGrayLevel_memb];
        xline(meanBinnedGrayLevel_memb, 'Color', 'g', 'LineWidth', 2);
        saveas(gcf, string(path)+'_memb_mean.png')
        close

%% GOLGI HISTOGRAM 
        figure
        h_golgi=histogram(golgi,'BinWidth',binw, 'BinLimits',[hist_min, hist_max]); %plot histogram of golgi
        title(string(srcFile(n).name) +' histogram of golgi')
        data_golgi = h_golgi.Values;
        counts_golgi = data_golgi;
        centerBinGrayLevels_golgi = (h_golgi.BinEdges(1:end-1) + h_golgi.BinEdges(2:end)) / 2;
        meanBinnedGrayLevel_golgi = sum(centerBinGrayLevels_golgi .* counts_golgi) / sum(counts_golgi);
        golgi_mean=[golgi_mean,meanBinnedGrayLevel_golgi];
        xline(meanBinnedGrayLevel_golgi, 'Color', 'r', 'LineWidth', 2);
        saveas(gcf, string(path) + '_golgi_mean.png')
        close

%% Uncomment to save histogram data
%         data_golgi_filename=string(path_golgi)+'.csv';
%         data_memb_filename=string(path_memb)+'.csv';

%         csvwrite(string(path)+'_golgi_data.csv', transpose(data_golgi));
%         csvwrite(string(path)+'_memb_data.csv', transpose(data_memb));

end

%% Save R_mean data (for Method 1)

csvwrite('memb_mean_data_32bit.csv', transpose(memb_mean));
csvwrite('golgi_mean_data_32bit.csv', transpose(golgi_mean));
writematrix(transpose(file_name_arr), 'file_names.csv');

% GOLGI/MEMBRANE
mean_ratio=golgi_mean./memb_mean;
csvwrite('mean_ratio_32bit.csv', transpose(mean_ratio));
