function cropimage = CPcroping(Sal,Img,options)
left = 1;
right = 2;
 
UniSal = max(Sal{left},Sal{right});
UniSal = (UniSal-min(UniSal(:)))/(max(UniSal(:))-min(UniSal(:)));
UniSal = imdilate(UniSal,strel('diamond',10));
[height,width] = size(UniSal);

windows.height = height;
windows.width  = floor(options.width*height/options.height);
if  windows.width>width 
    windows.width  = width;
    windows.height = floor(options.height*width/options.width);
    wenergy = 1000000;
    begin = 0;
    for i=0:height-windows.height
        energy = UniSal;
        energy(i+1:i+windows.height,:) = 0;
        energy = sum(energy(:));
        if energy < wenergy
            begin = i;
            wenergy = energy;
        end
    end
    cropimage{left} = Img{left}(begin+1:begin+windows.height,:,:);
    cropimage{right} = Img{right}(begin+1:begin+windows.height,:,:);
else 
    wenergy = 1000000;
    begin = 0;
    for i=0:width-windows.width
        energy = UniSal;
        energy(:,i+1:i+windows.width) = 0;
        energy = sum(energy(:));
        if energy < wenergy
            begin = i;
            wenergy = energy;
        end
    end
    cropimage{left} = Img{left}(:,begin+1:begin+windows.width,:);
    cropimage{right} = Img{right}(:,begin+1:begin+windows.width,:);
end
