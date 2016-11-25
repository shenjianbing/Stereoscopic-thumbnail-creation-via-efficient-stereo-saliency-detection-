function [cropimage] = OAcroping(Sal,Img)
left = 1;
right = 2;
 
UniSal = max(Sal{left},Sal{right});
UniSal = (UniSal-min(UniSal(:)))/(max(UniSal(:))-min(UniSal(:)));
UniSal = imdilate(UniSal,strel('diamond',25));
[height,width] = size(UniSal);
pixelNum = height*width;
[X,Y] = meshgrid(1:width,1:height); 
T = exp(UniSal)-1;
X = X.*T;
Y = Y.*T;

salCenter = [sum(X(:))/sum(T(:)) sum(Y(:))/sum(T(:))];
[X,Y] = meshgrid(1:width,1:height);
X = abs(X-salCenter(1)).*T;
Y = abs(Y-salCenter(2)).*T;
avwidth = 2*sum(X(:))/sum(T(:));
avheight = 2*sum(Y(:))/sum(T(:));
search.x1= max(floor(salCenter(1)-avwidth),1);
search.x2= min(floor(salCenter(1)+avwidth),width);
search.y1= max(floor(salCenter(2)-avheight),1);
search.y2= min(floor(salCenter(2)+avheight),height);
wenergy = 0;
for p=0.5:0.1:1
    windows.height = floor(p*(search.y2-search.y1));
    windows.width  = floor(p*(search.x2-search.x1));
    num = windows.height*windows.width;
    for i=0:search.x2-search.x1-windows.width
        for j=0:search.y2-search.y1-windows.height
            energy = UniSal;
            energy = energy(j+search.y1:j+search.y1+windows.height,i+search.x1:i+search.x1+windows.width);
            energy = sum(energy(:))/(num^0.5);
            if energy>wenergy 
                wenergy = energy;
                crop.x1=i+search.x1;
                crop.x2=i+search.x1+windows.width;
                crop.y1=j+search.y1;
                crop.y2=j+search.y1+windows.height;
            end
        end
    end 
end
cropimage{left} = Img{left}(crop.y1:crop.y2,crop.x1:crop.x2,:);
cropimage{right} = Img{right}(crop.y1:crop.y2,crop.x1:crop.x2,:);

% cropimage2{left} = Img{left};
% cropimage2{right} = Img{right};
% cropmask{left} =zeros(height,width);
% cropmask{left}(crop.y1:crop.y2,crop.x1:crop.x2)=1;
% cropmask{right} =zeros(height,width);
% cropmask{right}(crop.y1:crop.y2,crop.x1:crop.x2)=1;
% 
% mask = zeros(height,width);
% mask(crop.y1:crop.y2,crop.x1,:) = 1;
% mask(crop.y1:crop.y2,crop.x2,:) = 1;
% mask(crop.y1,crop.x1:crop.x2,:) = 1;
% mask(crop.y2,crop.x1:crop.x2,:) = 1;
% mask = imdilate(mask,strel('diamond',5));
% mask = double(repmat(~mask,[1,1,3]));
% cropimage2{left} = cropimage2{left}.*mask;
% cropimage2{right} = cropimage2{right}.*mask;






