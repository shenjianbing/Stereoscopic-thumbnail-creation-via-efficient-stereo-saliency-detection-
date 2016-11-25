clc
clear
addpath( genpath( '.' ) );
options.valScale = 60;
options.salThreshold = 0.2;
options.alpha = 0.002;
options.topRate = 0.01;
left = 1;
right = 2;

% set [w*h] of thumbnail via CPC
options.height = 100;
options.width = 50;

foldername = fileparts( mfilename( 'fullpath' ) );
ImgsFiles{left} = fullfile(foldername, 'images','left');
ImgsFiles{right} = fullfile(foldername, 'images','right');

Imgs{left} = imdir(ImgsFiles{left});
Imgs{right} = imdir(ImgsFiles{right});
nimgs = length( Imgs{left});

for index = 1: nimgs       
        index
        [~, ImgName{left}] = fileparts(Imgs{left}(index).name);
        [~, ImgName{right}] = fileparts(Imgs{right}(index).name);
        if exist(fullfile(ImgsFiles{left}, [ImgName{left} '.png']),'file')
            Img{left} = double(imread(fullfile(ImgsFiles{left}, [ImgName{left} '.png'])));
            Img{right} = double(imread(fullfile(ImgsFiles{right}, [ImgName{right} '.png'])));
        elseif exist(fullfile(ImgsFiles{left}, [ImgName{left} '.jpg']),'file')
            Img{left} = double(imread(fullfile(ImgsFiles{left}, [ImgName{left} '.jpg'])));
            Img{right} = double(imread(fullfile(ImgsFiles{right}, [ImgName{right} '.jpg'])));
        elseif exist(fullfile(ImgsFiles{left}, [ImgName{left} '.bmp']),'file')
            Img{left} = double(imread(fullfile(ImgsFiles{left}, [ImgName{left} '.bmp'])));
            Img{right} = double(imread(fullfile(ImgsFiles{right}, [ImgName{right} '.bmp'])));
        end       
        outfolder = fullfile( foldername, 'results', ImgName{left}(1:end-2));
        if( ~exist( outfolder, 'dir' ) ), mkdir( outfolder ), end;
        
        if( ~exist( fullfile( outfolder, 'regions'), 'dir' ) ), mkdir( fullfile( outfolder, 'regions') ), end;
        
        [ height,width ] = size(Img{left}(:,:,1));
        PixNum = height*width;
        
        %% compute disparity        
        disfile = fullfile( outfolder, 'disparity.mat');
        if( exist( disfile, 'file' ) )
            load( disfile );
        else
            [disparity{left},~,colorflow{left}] = computeDisparity(Img{left},Img{right});
            [disparity{right},~,colorflow{right}] = computeDisparity(Img{right},Img{left});
%             imwrite(colorflow{left},[outfolder '/' 'disparity_left.bmp']);
%             imwrite(colorflow{right},[outfolder '/' 'disparity_right.bmp']);
%             imagesc(disparity{left});
            save( fullfile( outfolder,'disparity.mat'), 'disparity', '-v7.3' );
        end
    
        [dx dy]= gradient(disparity{left});
        magnitude{left} = sqrt( dx.^2 + dy.^2 );
        magnitude{left} = magnitude{left}/max(magnitude{left}(:));
        magnitude{left} = imdilate(magnitude{left},strel('diamond',10));
        
        [dx dy]= gradient(disparity{right});
        magnitude{right} = sqrt( dx.^2 + dy.^2 );
        magnitude{right} = magnitude{right}/max(magnitude{right}(:));
        magnitude{right} = imdilate(magnitude{right},strel('diamond',10));
         
         %% compute superpixel
         Attr=[ height ,width, 1000, 20, PixNum ];
         imgVecR = reshape( Img{left}(:,:,1)', PixNum, 1);
         imgVecG = reshape( Img{left}(:,:,2)', PixNum, 1);
         imgVecB = reshape( Img{left}(:,:,3)', PixNum, 1); 
    
         [ leftLabel, leftSup1, leftSup2, leftSup3, RegionNum{left} ] = SLIC( double(imgVecR), double(imgVecG), double(imgVecB), Attr );
         superpixels{left}.Label = int32(reshape(leftLabel+1,width,height)');
         superpixels{left}.Lab = [leftSup1 leftSup2 leftSup3] ;
         superP{1} = uint32(superpixels{left}.Label);
         [ superpixels{left}.colours, superpixels{left}.centres, superpixels{left}.t ] = getSuperpixelStats(Img(left:left), superP, RegionNum{left} );         
         subImg{left} = superpixel2pixel(double(superpixels{left}.Label),double(superpixels{left}.colours));
         subImg{left} = uint8(reshape(subImg{left},height,width,3));
            
         imgVecR = reshape( Img{right}(:,:,1)', PixNum, 1);
         imgVecG = reshape( Img{right}(:,:,2)', PixNum, 1);
         imgVecB = reshape( Img{right}(:,:,3)', PixNum, 1); 
    
         [ rightLabel, rightSup1, rightSup2, rightSup3, RegionNum{right} ] = SLIC( double(imgVecR), double(imgVecG), double(imgVecB), Attr );
         superpixels{right}.Label = int32(reshape(rightLabel+1,width,height)');
         superpixels{right}.Lab = [rightSup1 rightSup2 rightSup3] ;
         superP{1} = uint32(superpixels{right}.Label);
         [ superpixels{right}.colours, superpixels{right}.centres, superpixels{right}.t ] = getSuperpixelStats(Img(right:right), superP, RegionNum{right} ); 
         subImg{right} = superpixel2pixel(double(superpixels{right}.Label),double(superpixels{right}.colours));
         subImg{right} = uint8(reshape(subImg{right},height,width,3));
         
        %% compute edge correspondence
        edgfile = fullfile( outfolder, 'edg.mat');
        if( exist( edgfile, 'file' ) )
            load( edgfile );
        else
             [gb_thin_CSG{left}] = Gb_CSG(Img{left});
             [gb_thin_CSG{right}] = Gb_CSG(Img{right});
             save( fullfile( outfolder,'edg.mat'), 'gb_thin_CSG', '-v7.3' );
        end
        
        %% compute disparity-edge map and initial saliency via gradient-flow       
         comDiff{left} = computeDE(magnitude{left},gb_thin_CSG{left});        
         [V_Energy1,H_Energy1,V_Energy2,H_Energy2] = energy_map(comDiff{left});            
         Energy{left} = min(min(min(H_Energy1,H_Energy2),V_Energy1),V_Energy2); 
         Energy{left} = Energy{left}/max(Energy{left}(:));

         comDiff{right} = computeDE(magnitude{right},gb_thin_CSG{right});          
         [V_Energy1,H_Energy1,V_Energy2,H_Energy2] = energy_map(comDiff{right});
         Energy{right} = min(min(min(H_Energy1,H_Energy2),V_Energy1),V_Energy2);
         Energy{right} = Energy{right}/max(Energy{right}(:));
        
         %% considering center-bais
        nLabel = double(max(superpixels{left}.Label(:)));
        L{1} = uint32(superpixels{left}.Label);
        S{1} = repmat(Energy{left}(:),[1 3]);
        [ R, ~, ~ ] = getSuperpixelStats(S(1:1),L, nLabel );
        R = double(R(:,1));
        R = (R-min(R))/(max(R)-min(R));
        IniSal_region{left} = R;
        IniSal{left} = double(R(superpixels{left}.Label));
        [IniSal_region{left} IniSal{left} mask] = centerBias(IniSal_region{left},IniSal{left},superpixels{left},options);
        %imwrite(mask,[outfolder '/' ImgName{left} '_mask' '.bmp']);

        nLabel = double(max(superpixels{right}.Label(:)));
        L{1} = uint32(superpixels{right}.Label);
        S{1} = repmat(Energy{right}(:),[1 3]);
        [ R, ~, ~ ] = getSuperpixelStats(S(1:1),L, nLabel );
        R = double(R(:,1));
        R = (R-min(R))/(max(R)-min(R));
        IniSal_region{right} = R;
        IniSal{right} = double(R(superpixels{right}.Label));
        [IniSal_region{right} IniSal{right} mask] = centerBias(IniSal_region{right},IniSal{right},superpixels{right},options);
       
        %imwrite(mask,[outfolder '/' ImgName{right} '_mask' '.bmp']);          
%         imwrite(IniSal{left},[outfolder '/' ImgName{left} '_initial' '.bmp']);
%         imwrite(IniSal{right},[outfolder '/' ImgName{right} '_initial' '.bmp']); 
%         imwrite(comDiff{left},[outfolder '/' ImgName{left} '_comDiff' '.bmp']);
%         imwrite(comDiff{right},[outfolder '/' ImgName{right} '_comDiff' '.bmp']); 
        %% saliency optimization
        Sal = computeUS(IniSal_region,IniSal,superpixels,options);
        imwrite(Sal{left},[outfolder '/' ImgName{left}  '.bmp']);
        imwrite(Sal{right},[outfolder '/' ImgName{right}  '.bmp']);
        
        %% Context-Persistent
         CPC = CPcroping(Sal,Img,options);
        %% Object-Aware
%         OAC = OAcroping(Sal,Img);
       
        
end