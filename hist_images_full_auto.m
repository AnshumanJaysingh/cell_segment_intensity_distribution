close all
clear all
clc

% mkdir golgi_images
% mkdir membrane_images

currentfolder=string(pwd)

srcFile = dir(currentfolder+'/*.tif');
image_folder = currentfolder+'/'

% specify folder for saving
% memb_save_folder = ''; 
% golgi_save_folder = '';

golgi_mean=[];
memb_mean=[];
file_name_arr=[];

%%
for n=1:length(srcFile)        
        FileName = strcat(image_folder,srcFile(n).name);
        file_name_arr=[file_name_arr, string(srcFile(n).name)];

        % read image and split channels
        I1=imread(FileName,1);
        figure
        imshow(I1,[])
        I1=medfilt2(I1);
        imshow(I1,[])
        close

        I2=imread(FileName,2);
        imshow(I2,[])
        I2=medfilt2(I2);
        imshow(I2,[])
        close
        
%         maxI1=max(I1, [], 'all')
%         maxI2=max(I2, [], 'all')
%         max_absolute=max(maxI1,maxI2)


%% GOLGI SEGMENTATION
        [level,metric] = multithresh(I2); %multi level thresholding of red (golgi) channel
        seg_I2 = imquantize(I2,level);
        RGB_golgi = label2rgb(seg_I2); 	 

        figure;
        imshow(seg_I2,[])
        imshow(RGB_golgi)
%         saveas(gcf, string(path)+'_golgi_segment.png')
        close
        imwrite(RGB_golgi, FileName+'_golgi_segment.png','png'); 
        seg=seg_I2-1;
        golgi = immultiply(uint16(seg),I1);
%         max_golgi=max(golgi, [], 'all')

        % max_golgi=max(var_golgi, [], 'all')
        % golgi=var_golgi./max_golgi;
        
%         figure
%         imshow(golgi,[])
        % file=golgi;


        % imwrite(file,baseFileName+'golgi','tiff'); %write segmented golgi (red) channel image to a tiff file


%% MEMBRANE SEGMENTATION

        Golgi_rem_Image=I1-(golgi);

        % Get a binary image.
        mask = Golgi_rem_Image > 100;
        % Take only the largest blob.  Use 4-connectivity.
        mask = bwareafilt(mask, 2, 4);

        % To find the average length and width we need to get the Euclidean distance transform and the skeleton of it.
        edtImage = bwdist(~mask);
        edtImage = edtImage./max(edtImage,[],'all');
%         mean_edt = mean(edtImage(:))
        threshold = 0.3; % change threshold to select more or less of the membrane
        edtImage = edtImage > threshold;
        
        figure
        imshow(edtImage, [])
        close       
      
        memb = immultiply(uint16(edtImage),Golgi_rem_Image); 

        figure
        imshow(memb,[]);
        RGB_memb = label2rgb(memb);
        imwrite(RGB_memb, FileName+'_memb_segment.png','png'); 
        close
        
        path=strcat(image_folder, srcFile(n).name);
%         path_memb=strcat(memb_save_folder, srcFile(n).name)
%         path_golgi=strcat(golgi_save_folder, srcFile(n).name)


%%  MEMBRANE HISTOGRAM 
        figure 
        h_memb=histogram(memb,'BinWidth',5, 'BinLimits',[1,4090]); %plot histogram of membrane
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
        h_golgi=histogram(golgi,'BinWidth',5, 'BinLimits',[1,4090]); %plot histogram of golgi
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

%% Save mean data (Method 1)

csvwrite('memb_mean_data_4090.csv', transpose(memb_mean));
csvwrite('golgi_mean_data_4090.csv', transpose(golgi_mean));
writematrix(transpose(file_name_arr), 'file_names.csv');

% GOLGI/MEMBRANE
mean_ratio=golgi_mean./memb_mean;
csvwrite('mean_ratio_4090.csv', transpose(mean_ratio));

% %% Merge all csv files 
% 
% csv1 = csvread('file_names.csv')
% csv2 = csvread('golgi_mean_data_4090.csv')
% csv3 = csvread('memb_mean_data_4090.csv')
% csv3 = csvread('mean_ratio_4090.csv')
% 
% % Concatenate vertically
% allCsv = [csv1; csv2; csv3; csv3]; 
% csvwrite('Tabulated_data.csv', allCsv);
