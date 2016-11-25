function geoDist = computeBR(superpixels,Sal,geoDist,options)

nLabel = double(max(superpixels.Label(:)));

[height,width] = size(Sal);
PixNum = height*width;
Salmask = Sal > min(options.salThreshold,mean(Sal(:)));
Salmask = imdilate(Salmask,strel('diamond',2));

BTW=bwconncomp(~Salmask,4);
for i=1:BTW.NumObjects
    if numel(BTW.PixelIdxList{i}) < PixNum*0.01
        Salmask(BTW.PixelIdxList{i}(:))=1;
    end
end

BTW=bwconncomp(Salmask,4);
for i=1:BTW.NumObjects
    if numel(BTW.PixelIdxList{i}) < PixNum*0.01
        Salmask(BTW.PixelIdxList{i}(:))=0;
    end
end

if sum(Salmask(:))/PixNum>0.9
    Salmask = UniSal > mean(UniSal(:));
    Salmask = imdilate(Salmask,strel('diamond',5));

    BTW=bwconncomp(~Salmask,4);
    for i=1:BTW.NumObjects
        if numel(BTW.PixelIdxList{i}) < PixNum*0.01
            Salmask(BTW.PixelIdxList{i}(:))=1;
        end
    end

    BTW=bwconncomp(Salmask,4);
    for i=1:BTW.NumObjects
        if numel(BTW.PixelIdxList{i}) < PixNum*0.01
            Salmask(BTW.PixelIdxList{i}(:))=0;
        end
    end
end
 
fd = int32(Salmask).*superpixels.Label;
fd = unique(fd(:));
fd(fd==0) = [];
bdIds = setdiff(unique(superpixels.Label(:)),fd);       
 
Salmask(:) =0;
Salmask(1:end,1)=1;
Salmask(1:end,end)=1;
Salmask(1,1:end)=1;
Salmask(end,1:end)=1;
bd = int32(Salmask).*superpixels.Label;
bd= unique(bd(:));
bd(bd==0) = [];


%colDistM = squareform(pdist(superpixels.Lab));

meanLabCol = colorspace('Lab<-', double(superpixels.colours)/255);
colDistM = GetDistanceMatrix(meanLabCol);
colDistM(bdIds,bdIds)=colDistM(bdIds,bdIds)/5;

[ConSPix conedge]= find_connect_superpixel( superpixels.Label, nLabel, height ,width ); 
ConSPix = full(ConSPix)*2;
ConSPix = ConSPix +eye(size(ConSPix));

[clipVal,geoSigma,neiSigma ]= EstimateDynamicParas(ConSPix,colDistM);
[bgProb, bdCon, bgWeight] = EstimateBgProb(colDistM, ConSPix, double(bd), clipVal, geoSigma);
bgProb = double(bgProb/max(bgProb(:)));
%bgProb(bdIds) = 1;
optwCtr = SaliencyOptimization(ConSPix, double(bdIds), colDistM, neiSigma, bgWeight, geoDist);
geoDist = (optwCtr-min(optwCtr(:))) / (max(optwCtr(:))-min(optwCtr(:))) ; 

