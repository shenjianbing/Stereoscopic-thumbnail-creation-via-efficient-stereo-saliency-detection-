function Sal = computeUS(IniSal_region,IniSal,superpixels,options)
left = 1;
right = 2;
nLabel{left} = double(max(superpixels{left}.Label(:)));
nLabel{right} = double(max(superpixels{right}.Label(:)));
 
UniSal = max(IniSal{left},IniSal{right});
UniSal = (UniSal-min(UniSal(:)))/(max(UniSal(:))-min(UniSal(:)));
[height,width] = size(UniSal);
PixNum = height*width;
Salmask = UniSal > min(options.salThreshold,mean(UniSal(:)));%
Salmask = imdilate(Salmask,strel('diamond',2));


valLAB = [superpixels{left}.Lab;superpixels{right}.Lab];
colDistM = squareform(pdist(valLAB));

ConSPix = zeros(nLabel{left}+nLabel{right},nLabel{left}+nLabel{right}); 
Conedge = [];


[conSPix conedge]= find_connect_superpixel( superpixels{left}.Label, nLabel{left}, height ,width ); 
Conedge = [Conedge;conedge];
[conSPix conedge]= find_connect_superpixel( superpixels{right}.Label, nLabel{right}, height ,width ); 
Conedge = [Conedge;conedge+nLabel{left}];
intralength = size(Conedge,1);
        
[x y] = meshgrid(1:nLabel{left},1:nLabel{right});
conedge = [x(:),y(:)];
connect = sum((superpixels{left}.centres(conedge(:,1),:) - superpixels{right}.centres(conedge(:,2),:)).^2,2 );
cross_po_dis = conedge(connect<800,:);

cross_po_dis(:,2) = cross_po_dis(:,2)+ nLabel{left};
Conedge = [Conedge;cross_po_dis];
               
ConSPix=sparse([Conedge(:,1);Conedge(:,2)],[Conedge(:,2);Conedge(:,1)], ...
         [ones(size(Conedge(:,1)));ones(size(Conedge(:,1)))],nLabel{left}+nLabel{right},nLabel{left}+nLabel{right});
ConSPix = full(ConSPix);
ConSPix = ConSPix +eye(size(ConSPix));
        

fd = superpixels{left}.Label(Salmask);
fd = unique(fd(:));
fd(fd==0) = [];
bd = setdiff(unique(superpixels{left}.Label(:)),fd);       
        
bdIds = [];
bdIds = [bdIds;bd];   
        
fd = superpixels{right}.Label(Salmask);
fd = unique(fd(:));
fd(fd==0) = [];
bd = setdiff(unique(superpixels{right}.Label(:)),fd);
bdIds = [bdIds;bd+nLabel{left}]; 
        
[clipVal,geoSigma,neiSigma ]= EstimateDynamicParas(ConSPix,colDistM);
geoDistall = GeodesicSaliency(ConSPix, double(bdIds), colDistM, clipVal, false,[])';


geoDist{left} = geoDistall(1:nLabel{left});
geoDist{right} = geoDistall(nLabel{left}+1:end);

tmp = sort(geoDist{left}, 'descend');
pos = round(options.topRate * length(tmp));
maxVal = tmp(pos);
geoDist{left} = (geoDist{left}-min(geoDist{left}(:))) / (maxVal-min(geoDist{left}(:))+0.001) ; 
geoDist{left}(geoDist{left} > 1) = 1;
geoDist{left} = geoDist{left}*0.9+IniSal_region{left}*0.1;
geoDist{left} = geoDist{left}/max(geoDist{left}(:));
Sal{left} = double(geoDist{left}(superpixels{left}.Label));

tmp = sort(geoDist{right}, 'descend');
pos = round(options.topRate * length(tmp));
maxVal = tmp(pos);
geoDist{right} = (geoDist{right} -min(geoDist{right}(:)))/ (maxVal-min(geoDist{right}(:))+0.001) ; 
geoDist{right}(geoDist{right} > 1) = 1;
geoDist{right} = geoDist{right}*0.9+IniSal_region{right}*0.1;
geoDist{right} = geoDist{right}/max(geoDist{right}(:));
Sal{right} = double(geoDist{right}(superpixels{right}.Label));

geoDist{left} = computeBR(superpixels{left},Sal{left},geoDist{left},options);
Sal{left} = double(geoDist{left}(superpixels{left}.Label));
geoDist{right} = computeBR(superpixels{right},Sal{right},geoDist{right},options);
Sal{right} = double(geoDist{right}(superpixels{right}.Label));

geoDistall =[geoDist{left};geoDist{right}];

valDistances=sqrt(sum((valLAB(Conedge(:,1),:)-valLAB(Conedge(:,2),:)).^2,2));
valDistances(intralength+1:end)=valDistances(intralength+1:end)/3;   
valDistances=normalize(valDistances);
weights=exp(-options.valScale*valDistances)+ 1e-5;
weights=sparse([Conedge(:,1);Conedge(:,2)],[Conedge(:,2);Conedge(:,1)], ...
    [weights;weights],(nLabel{left}+nLabel{right}),(nLabel{left}+nLabel{right}));
E = sparse(1:(nLabel{left}+nLabel{right}),1:(nLabel{left}+nLabel{right}),ones((nLabel{left}+nLabel{right}),1)); iD = sparse(1:(nLabel{left}+nLabel{right}),1:(nLabel{left}+nLabel{right}),1./sum(weights));
P = iD*weights;
geoDistall = (E-P+10*options.alpha*E)\geoDistall;

geoDist{left} = geoDistall(1:nLabel{left});
geoDist{right} = geoDistall(nLabel{left}+1:end);

tmp = sort(geoDist{left}, 'descend');
pos = round(options.topRate * length(tmp));
maxVal = tmp(pos);
geoDist{left} = (geoDist{left}-min(geoDist{left}(:))) / (maxVal-min(geoDist{left}(:))) ; 
geoDist{left}(geoDist{left} > 1) = 1;

tmp = sort(geoDist{right}, 'descend');
pos = round(options.topRate * length(tmp));
maxVal = tmp(pos);
geoDist{right} = (geoDist{right} -min(geoDist{right}(:)))/ (maxVal-min(geoDist{right}(:))) ; 
geoDist{right}(geoDist{right} > 1) = 1;

Sal{left} = double(geoDist{left}(superpixels{left}.Label));
Sal{right} = double(geoDist{right}(superpixels{right}.Label));
 