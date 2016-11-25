function comDiff = computeDE(magnitude,gb_thin_CSG)
        comDiff = (0.5+magnitude).*(0+gb_thin_CSG);
        comDiff = (comDiff-min(comDiff(:)))/(max(comDiff(:))-min(comDiff(:)));
        comDiff = imdilate(comDiff,strel('diamond',5));

       
