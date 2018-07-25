function image_processing(input_image, output_image, output_table, scale)
%% import the image and crop the image
%clear
%close all
%disp(scale)
%class(scale)
tic;
%input_image = '/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Liver_cell_nuclei/15_may_2015_IHC/images2run/Image_01.tif';
%input_image = '/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Liver_cell_nuclei/15_may_2015_IHC/images/Image_01_01.tif';
aa = imread(input_image);
%imagesc(aa)

disp('load image...')
toc

%% Crop and select fields
tic;

% crop 60 percents of initial image
scale_crop = 0.6;
[m2, m1, m3] = size(aa);

imwidth = floor(m1*sqrt(scale_crop));
imheight = floor(m2*sqrt(scale_crop));

imxmin = floor(m1/2-imwidth/2);
imymin = floor(m2/2-imheight/2);

rect_image = [imxmin imymin imwidth imheight]; 
bb = imcrop(aa,rect_image);
% imagesc(bb)

% cut croped images into 60 small images 
scale = str2num(scale);
% scale = 0.05;
if(scale>0.6)
    scale = 0.6
end

fprintf('scaling factor is %1.4f\n', scale);

[m2, m1, m3] = size(bb);
%fprintf('size of image is %d and %d \n', m1, m2);

imwidth = floor(m1/6);
imheight = floor(m2/10);
imxymin = zeros(2, 60);

for m = 1:6
    for n = 1:10
        %disp(m)
        %disp(n)
        imxymin(1, (m+(n-1)*6)) = 1+floor(imwidth*(m-1));
        imxymin(2, (m+(n-1)*6)) = 1+floor(imheight*(n-1));
    end
end
   
nb_pieces = floor(scale/0.01);
%rng(1);
rr = randperm(60,nb_pieces);

%bb = imcrop(aa);
%imagesc(bb);

disp('select fields...')
toc

%% to define cutoff of nuclei and cells
% Images around the Clock
nucleus_cutoff = 205; %Red values less than this cutoff in nucleus
cell_cutoff = 170; % Blue values bigger than this cutoff in cytoplasm and nucleus  
nucleus_watershed_dist = 0.3;
cell_watershed_dist = 0.6;
std_scale = 2.0;
min_size_cell = 200;

% E2f KO images
%nucleus_cutoff = 200; %Red values less than this cutoff in nucleus
%nucleus_watershed_dist = 0.2;
%cell_cutoff = 190; % Blue values bigger than this cutoff in cytoplasm and nucleus  
%cell_watershed_dist = 0.4;
%std_scale = 1.5;

%% loop for selected fields
%rr =  (1:3)
cell_index_init = 0;
ftemp = fopen(output_table, 'W');

[mm2, mm1, mm3] = size(aa);
fprintf(ftemp, '%d\t%d\t%d\n', [0 mm1 mm2]); % first line of output is the total pixel size of the whole image

for nb_imp = 1:nb_pieces
    
    disp('.............')
    % nb_imp = 1;
    fprintf('selected field index is %d\n', rr(nb_imp));
    fprintf('initial cell index is %d\n', cell_index_init);
   
    
    rect_image = [imxymin(1,rr(nb_imp)) imxymin(2, rr(nb_imp)) imwidth imheight]; 
    %rect_image = [imxymin(1,1) imxymin(2, 1) imwidth imheight];
    cc = imcrop(bb,rect_image); % small image to analysis
    
    % optimize watershed parameters
    %cc = imread('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Liver_cell_nuclei/15_may_2015_IHC/tests/Image_05_02_crop.tif');
    %imagesc(cc)
    %nucleus_watershed_dist = 0.5;
    %cell_watershed_dist = 0.5;
    %std_scale = 2.0;
    %% Nuclei segmentation
    tic;
    a = cc;
    b1 = bsxfun(@gt, a(:,:,3),a(:,:,2));
    b2 = bsxfun(@gt, a(:,:,3),a(:,:,1));
    b3 = a(:,:,1)<=nucleus_cutoff;

    b = bsxfun(@and, b1, b2);
    b = bsxfun(@and, b, b3);
    
    b = imfill(b,'holes');
    b = imopen(b,strel('disk',8));
    
    %imagesc(b)
    % Watershed for nuclei
    %bw = b;
    %nucleus_watershed_dist = 0.5
    bw = ~bwareaopen(~b, 4);
    D = -bwdist(~bw);
    
    mask = imextendedmin(D,nucleus_watershed_dist);
    %imshowpair(bw,mask,'blend')
    D2 = imimposemin(D,mask);
    Ld2 = watershed(D2);
    bw(Ld2 == 0) = 0;
    %imagesc(bw)

    b = bw; 
    a = imnorm(a);
    b = imnorm(b);
    [labeled_n, test_nb] = bwlabel(b, 8);  % Label each blob so we can make measurements of it
    mn = regionprops(labeled_n, a(:,:,1), {'Area'});
    if(test_nb<20) 
        continue
    end
    
    %imagesc(labeledImage)
    %b = bwmorph(b,'remove');
    %b = repmat(b,[1 1 3]);
    %imagesc((a-b))

    %Cell segementation
    tic;
    a = cc;
    c = a(:,:,3)>cell_cutoff;
    %imagesc(c)

    % Watershed function
    %min_size = min([mn.Area]);
    bw = ~bwareaopen(~c, min_size_cell);
    %imshow(bw2)
    D = -bwdist(~bw);
    mask = imextendedmin(D,cell_watershed_dist);
    %imshowpair(bw,mask,'blend')
    D2 = imimposemin(D,mask);
    Ld2 = watershed(D2);
    %bw3 = bw;
    bw(Ld2 == 0) = 0;
    %imshow(bw3)
    %c = bw;
    %imshow(c)
    c = imfill(bw,'holes');
    c = imopen(c,strel('disk',8));
    %imagesc(c)
    % clearn cell in the border
    c = imclearborder(c);
    %imagesc(c2);

    a = imnorm(a);
    c = imnorm(c);
    [labeled_c, nb_c] = bwlabel(c, 8);     % Label each blob so we can make measurements of it
    mc  = regionprops(labeled_c, a(:,:,1), {'Area'});
    %imagesc(labeledImage)
    disp('segmentation...')
    toc

    %% filter nucleus and cell objects
    tic;
    idx_n = find(([mn.Area]>quantile([mn.Area], 0.01)) & ([mn.Area]<quantile([mn.Area], 0.99)));
    b2 = ismember(labeled_n, idx_n);
    %imagesc(b2)
    labeled_n = bwlabel(b2, 8);  % Label each blob so we can make measurements of it
    %mn = regionprops(labeled_n, a(:,:,1), {'PixelValues','MeanIntensity','Area','Centroid','PixelIdxList','PixelList','Eccentricity'});

    idx_c = find(([mc.Area]>min([mn.Area])) & ([mc.Area]<(50*mean([mn.Area]))) & abs([mc.Area]-mean([mc.Area]))<std_scale*std([mc.Area]));
    c2 = ismember(labeled_c, idx_c);
    [labeled_c, nb_c] = bwlabel(c2, 8);     % Label each blob so we can make measurements of it
    mc  = regionprops(labeled_c, a(:,:,1), {'PixelValues','MeanIntensity','Area','Centroid','PixelIdxList','PixelList','Eccentricity'});
    %c = bwmorph(c,'remove');
    %c = repmat(c,[1 1 3]);
    %imagesc((a-c))
    % check Nuclei and cell segmentations
    %clf
    %imagesc((a-b-c))
    % display segmentation with labeled images
    %imagesc(a-repmat(bwmorph(imnorm(labeled_c>0), 'remove'),[1,1,3])-repmat(bwmorph(imnorm(labeled_n>0), 'remove'),[1,1,3]))
        
    disp('filtering...')
    toc;
    
    %% Find nuclei for each cell object by object
    % output_table = 'test.txt';
    tic;
    index_c = [ ];
    labeled_nf = zeros(size(labeled_n));
    %res = [];
    %disp('initial cell index...')
    %disp(cell_index_init)
    
    for nb = 1:nb_c
        %tic;
        %index_c = [1:2000];
        %idx_c = mc(nb).PixelIdxList;
        %idx_n_all = [mn.PixelIdxList];
        [r,c] = find(labeled_c==nb);
        %toc
    
        %tic;
        %xx = bwselect(labeled_n, c, r);
        %toc
        [bwx, idx] = bwselect(labeled_n, c, r);
        %toc
        %image
        PixelIdx_c = [mc(nb).PixelIdxList];
        %tic;
        if (isempty(idx)<1)
            %clf
            %imagesc(bwx)
            %tic;
            [xL, num] = bwlabel(bwx, 8);
            %toc
            %pause(0.5)
            %mnc  = regionprops(xL, a(:,:,1), {'PixelValues','MeanIntensity','Area','Centroid','PixelIdxList','PixelList','Eccentricity'});
            %tic;
            mnc  = regionprops(xL, a(:,:,1), {'Area','PixelIdxList'});
            %imagesc(labeledImage)
            %toc
            % save results and filter cut nuclei
            %tic;
            for nb_n = 1:num
                if (isempty(setdiff([mnc(nb_n).PixelIdxList], PixelIdx_c)))
                    index_c = [index_c nb]; % save the cell object if nuclei are found inside
                    labeled_nf = bsxfun(@or, labeled_nf, (xL==nb_n)); % save the nucleus objects
                    res = [(nb+cell_index_init) [mc(nb).Area] [mnc(nb_n).Area]];
                
                    % save table in the output
                    %dlmwrite(output_table,res,'-append', 'delimiter','\t');
                    fprintf(ftemp, '%d\t%d\t%d\n',res);
                end
            end
            %toc
        end
        %toc

    end
    
    cell_index_init = cell_index_init + nb_c;
   
    disp('counting ...')
    toc
    %% save the segmented image
    tic;
    %b3 = ismember(labeled_n, index_n);
    %labeled_nf = bwlabel(b3, 8);  % Label each blob so we can make measurements of it
    %mnf = regionprops(labeled_nf, a(:,:,1), {'PixelValues','MeanIntensity','Area','Centroid','PixelIdxList','PixelList','Eccentricity'});
    index_c = unique(index_c);
    c3 = ismember(labeled_c, index_c);
    labeled_cf = bwlabel(c3, 8);     % Label each blob so we can make measurements of it
    %mcf = regionprops(labeled_cf, a(:,:,1), {'PixelValues','MeanIntensity','Area','Centroid','PixelIdxList','PixelList','Eccentricity'});

    % Final segementations for cells adn nuclei
    %clf
    %imagesc(a-repmat(bwmorph(imnorm(labeled_cf>0), 'remove'),[1,1,3])-repmat(bwmorph(imnorm(labeled_nf>0), 'remove'),[1,1,3]))
    cell_nuclei = a-repmat(bwmorph(imnorm(labeled_cf>0), 'remove'),[1,1,3])-repmat(bwmorph(imnorm(labeled_nf>0), 'remove'),[1,1,3]);
    %imagesc(cell_nuclei)
    
    %dlmwrite(output_table,res, 'delimiter','\t');
    imoutput = strcat(output_image, '_part_', num2str(nb_imp), '.png');
    imwrite(uint8(cell_nuclei .* 255),imoutput);

    disp('save image...')
    toc
    
end

fclose(ftemp);


end