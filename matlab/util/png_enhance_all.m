% enhance all pngs in dir

!ls *.png > pnglist.tmp 
files=importdata('pnglist.tmp');

for k = 1:numel(files)
  
  imgA = imread(files{k});
  imgB = enhance_overlay(imgA,2);
  imwrite(imgB,files{k}); 
  disp(['wrote to file: ' files{k}]);
  
end
