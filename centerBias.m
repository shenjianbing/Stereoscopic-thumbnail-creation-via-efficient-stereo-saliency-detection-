function [IniSal_region,IniSal,mask] = centerBias(IniSal_region,IniSal,superpixels,options) 
[height,width] = size(IniSal);
nLabel = double(max(superpixels.Label(:)));
Label = superpixels.Label;
fd = Label(IniSal>mean(IniSal_region(:)));
fd = unique(fd(:));
fd(fd==0) = [];
fd = sort(fd);

bd = Label(IniSal<=mean(IniSal_region(:)));
bd = unique(bd(:));
bd(bd==0) = [];
        
colDistM = squareform(pdist(superpixels.Lab));
colDistMx = colDistM;
colDistMx(bd,:) = colDistMx(bd,:)*5;
colDistMx(:,bd) = colDistMx(:,bd)*5;
[conSPix Conedge]= find_connect_superpixel( Label, nLabel, height ,width );  
ConSPix=sparse([Conedge(:,1);Conedge(:,2)],[Conedge(:,2);Conedge(:,1)], ...
         [ones(size(Conedge(:,1)));ones(size(Conedge(:,1)))],nLabel,nLabel);
ConSPix = full(ConSPix);
ConSPix = ConSPix +eye(size(ConSPix));
       
        
bcenters = superpixels.centres( bd,:);
b_x = repmat(bcenters(:,1)',[nLabel 1]);
b_y = repmat(bcenters(:,2)',[nLabel 1]);        
x = repmat(superpixels.centres(:,1),[1 size(bd,1)]);
y = repmat(superpixels.centres(:,2),[1 size(bd,1)]);        
b_dis = ((x-b_x).^2+(y-b_y).^2).^0.5;
b_dis = min(b_dis');
b_dis = b_dis(:,fd);
b_dis = (b_dis-min(b_dis(:)))/(max(b_dis(:))-min(b_dis(:)));
b_dis = exp(2*b_dis);

fcenters = superpixels.centres( fd,:);
f_x = repmat(fcenters(:,1)',[nLabel 1]);
f_y = repmat(fcenters(:,2)',[nLabel 1]);        
x = repmat(superpixels.centres(:,1),[1 size(fd,1)]);
y = repmat(superpixels.centres(:,2),[1 size(fd,1)]);        
f_dis = ((x-f_x).^2+(y-f_y).^2).^0.5;
f_dis = f_dis(fd,:);
f_dis = sum(f_dis,2)/size(fd,1);
f_dis = (f_dis-min(f_dis(:)))/(max(f_dis(:))-min(f_dis(:)));
f_dis = exp(5*f_dis);

border =  unique([find(fcenters(:,1)<height/4);find(fcenters(:,1)>3*height/4);find(fcenters(:,2)<width/4);find(fcenters(:,2)>3*width/4)]);
        
b_geoDist = GeodesicSaliency(ConSPix, double(bd), colDistM, 0,false,[]);
b_geoDist = b_geoDist / max(b_geoDist(:)); 
b_geoDist(:,bd) = [];
b_geoDist = exp(b_geoDist);
adjcMatrix_virtual = tril(ConSPix, -1);
edgeWeight = colDistMx(adjcMatrix_virtual > 0);
    
f_geoDist = [];
for i = 1:size(fd,1)
    f_geoDist = [f_geoDist;graphshortestpath(sparse(adjcMatrix_virtual), fd(i), 'directed', false, 'Weights', edgeWeight)];             
end
f_geoDist(:,bd) = [];
%geoDistMatrix = graphallshortestpaths(sparse(adjcMatrix), 'directed', false, 'Weights', edgeWeight);
measure = max(f_geoDist).*f_dis'./b_geoDist./b_dis;
if size(border,1)>1
    measure(border) = measure(border)*2;
end
[~,pos] = min(measure);
%[~,pos] = max(dis);
center = fd(pos);

fcenters = superpixels.centres( center,:);
f_x = repmat(fcenters(:,1)',[size(fd,1) 1]);
f_y = repmat(fcenters(:,2)',[size(fd,1) 1]);              
f_disx = (superpixels.centres(fd(:),1) -f_x).^2;
f_disy = (superpixels.centres(fd(:),2) -f_y).^2;
f_disx =  (4*mean(f_disx)).^0.5;
f_disy =  (4*mean(f_disy)).^0.5;
[X,Y] = meshgrid(1:width,1:height); 
mask = (exp(-((Y-fcenters(:,1)).^2/(f_disx^2)+ (X-fcenters(:,2)).^2/(f_disy^2))));
mask = double(mask/max(mask(:)));

IniSal = (0.1+IniSal).*(0.1+mask);
L{1} = uint32(Label);
S{1} = repmat(IniSal(:),[1 3]);
[ R, ~, ~ ] = getSuperpixelStats(S(1:1),L, nLabel );
R = double(R(:,1));

tmp = sort(R, 'descend');
pos = round(options.topRate * length(tmp));
maxVal = tmp(pos);
R = (R-min(R(:))) / (maxVal-min(R(:))) ; 
R(R > 1) = 1;

% R = (R-min(R))/(max(R)-min(R));
IniSal_region = R;
IniSal = IniSal_region(Label);


