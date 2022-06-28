close all
clear all


srcFile = dir('D:\swarang\golgi_membrane\Analysis-Cells-16022021\Doxorubicin-20-nM-120-hours-18-hours-Cycloheximide-+ve\*.tif');
image_folder ='D:\swarang\golgi_membrane\Analysis-Cells-16022021\Doxorubicin-20-nM-120-hours-18-hours-Cycloheximide-+ve\'

%specify folder for saving
% memb_save_folder = 'D:\swarang\golgi_membrane\yo\Analysis-Cells-16012021\trial_run\memb_data'; 
% golgi_save_folder = 'D:\swarang\golgi_membrane\yo\Analysis-Cells-16012021\trial_run\golgi_data';


for n=1:length(srcFile)
        
%         baseFileName="Cell4" 
%         FileName=baseFileName+'.tif'
        % change this to the file name of image file

        FileName = strcat(image_folder,srcFile(n).name);

        % read image and split channels
        I1=imread(FileName,1);
        I2=imread(FileName,2);


%         maxI1=max(I1, [], 'all')
%         maxI2=max(I2, [], 'all')
%         max_absolute=max(maxI1,maxI2)




        % --------------------------------------------------------------------------------------------------
        %%%%%%%%%%%%%%%%%%     GOLGI SECTION     %%%%%%%%%%%%%%%%%%
        % --------------------------------------------------------------------------------------------------


        [level,metric] = multithresh(I2); %multi level thresholding of red (golgi) channel
        seg_I2 = imquantize(I2,level);


        RGB = label2rgb(seg_I2); 	 

%         figure;
%         % imshow(seg_I2,[])
%         imshow(RGB)

        seg=seg_I2-1;

        % var_golgi = (immultiply(uint16(seg_I2),I1))./2; %multiply sugmented red channel with original green channel
        golgi = immultiply(uint16(seg),I1);
        max_golgi=max(golgi, [], 'all')

        % max_golgi=max(var_golgi, [], 'all')
        % golgi=var_golgi./max_golgi;
        % figure
        % imshow(golgi,[])
        % file=golgi;


        % imwrite(file,baseFileName+'golgi','tiff'); %write segmented golgi (red) channel image to a tiff file


        % --------------------------------------------------------------------------------------------------
        %%%%%%%%%%%%%%%%     MEMBRANE SECTION     %%%%%%%%%%%%%%%%%%
        % --------------------------------------------------------------------------------------------------



        Golgi_rem_Image=I1-(golgi);
%         figure;
%         imshow(Golgi_rem_Image,[])

        % Get the green channel.
        % imhist(grayImage);


        % Get a binary image.
        mask = Golgi_rem_Image > 100;
        % Take only the largest blob.  Use 4-connectivity.
        mask = bwareafilt(mask, 1, 4);

        % To find the average length and width we need to get the Euclidean distance transform and the skeleton of it.
        edtImage = bwdist(~mask);
        edtImage=edtImage./max(edtImage,[],'all');


        % var_memb = (immultiply(uint16(edtImage),I1))./(2^4); 
        memb = immultiply(uint16(edtImage),Golgi_rem_Image); 
        max_memb=max(memb, [], 'all')

        % max_memb=max(var_memb, [], 'all')
        % memb=var_memb./max_memb;

        % figure
        % imshow(memb,[]);
        
        path=strcat(image_folder, srcFile(n).name)
%         path_memb=strcat(memb_save_folder, srcFile(n).name)
%         path_golgi=strcat(golgi_save_folder, srcFile(n).name)

        mean_memb_arr=[]

        % %%%%%%%%%% MEMBRANE HISTOGRAM %%%%%%%%%%
        figure 
        h_memb=histogram(memb,'BinWidth',5, 'BinLimits',[1,4090]) %plot histogram of membrane
        title(string(srcFile(n).name) + ' histogram of membrane')
        data_memb = h_memb.Values;
        counts_memb = data_memb;
        centerBinGrayLevels_memb = (h_memb.BinEdges(1:end-1) + h_memb.BinEdges(2:end)) / 2;
        meanBinnedGrayLevel_memb = sum(centerBinGrayLevels_memb .* counts_memb) / sum(counts_memb)
        mean_memb_arr = [mean_memb_arr, meanBinnedGrayLevel_memb];
        xline(meanBinnedGrayLevel_memb, 'Color', 'g', 'LineWidth', 2);
        saveas(gcf, string(path)+'_memb_mean.png')
        close


        % %%%%%%%%%% GOLGI HISTOGRAM %%%%%%%%%%
        figure
        h_golgi=histogram(golgi,'BinWidth',5, 'BinLimits',[1,4090]) %plot histogram of golgi
        title(string(srcFile(n).name) +' histogram of golgi')
        data_golgi = h_golgi.Values;
        counts_golgi = data_golgi;
        centerBinGrayLevels_golgi = (h_golgi.BinEdges(1:end-1) + h_golgi.BinEdges(2:end)) / 2;
        meanBinnedGrayLevel_golgi = sum(centerBinGrayLevels_golgi .* counts_golgi) / sum(counts_golgi)
        xline(meanBinnedGrayLevel_golgi, 'Color', 'r', 'LineWidth', 2);
        saveas(gcf, string(path) + '_golgi_mean.png')
        close


%         data_golgi_filename=string(path_golgi)+'.csv';
%         data_memb_filename=string(path_memb)+'.csv';

        csvwrite(string(path)+'_golgi_data.csv', transpose(data_golgi));
        csvwrite(string(path)+'_memb_data.csv', transpose(data_memb));

end
        % imwrite(memb,baseFileName+'membrane.tif','tiff')
