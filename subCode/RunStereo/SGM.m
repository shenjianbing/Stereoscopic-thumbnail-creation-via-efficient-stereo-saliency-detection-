function Disp = SGM(IL, IR)
Folder = 'RunStereo\in\'
i = 1;
pad = 512
f_sgm.maxDisp = 512;
f_sgm.blocksz =3;
f_sgm.outf =  [Folder '\SGM_Disp.png'];
[m,n,c] = size(IL);
ILp = [zeros(m,pad,c) IL];
IRp = [zeros(m,pad,c) IR];

imwrite(ILp, [Folder '\ILp.png']);
imwrite(IRp, [Folder '\IRp.png']);

system(['RunStereo\RunStereo.exe ' [Folder '\ILp.png'] ' '  [Folder '\IRp.png']  ' --max-disparity=' num2str(f_sgm.maxDisp) ' --blocksize=' num2str(f_sgm.blocksz) ' -o ' f_sgm.outf]); 


D = imread(f_sgm.outf);
imwrite(D(:,1:end,:),f_sgm.outf);
Disp = -double(D(:,1:end,:))./3;
