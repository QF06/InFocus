function [imPatch,EPatch, pixelRatio] = userDefineImage(folderName, pixelWH, pixelToMm,Threshold)
%% Stitching images and plot the beam profile, calculate the beam width
%% Read the images
NImage = size(dir([folderName '\*.bmp']),1);
imStitch = [];
switch NImage
  case 5
    for i=NImage-1:-1:2
      image=imread([folderName,'\',num2str(i),'.bmp']);
      image=image(:,:,1);
      imStitch=vertcat(imStitch,image);
    end
  case 1
    image = imread([folderName,'\',num2str(1),'.bmp']);
    image = image(:,:,1);
    imStitch = vertcat(zeros(size(image)),image,zeros(size(image)));
end 
imStitch = double(imStitch);
[ySize, xSize] = size(imStitch);

%% Patch image 
imPatch = zeros(max(xSize,ySize));
for i=(ySize-xSize)/2 :(ySize+xSize)/2-1
  imPatch(:,i) = imStitch(:,i-(ySize-xSize)/2+1);
end
imPatch(imPatch<max(max(imPatch))* Threshold)=0;
imPatch = imresize(imPatch,[pixelWH,pixelWH]);
EPatch = sqrt(imPatch);
pixelRatio = pixelToMm * ySize/pixelWH;
Center = regionprops(imPatch,mat2gray(imPatch),'WeightedCentroid');
XYCentroids = [Center.WeightedCentroid];
YCentroid = nanmean(XYCentroids(1:2:end));
XCentroid = nanmean(XYCentroids(2:2:end));
[Columns,Rows] = size(imPatch);
YShift = round(Columns/2 - YCentroid);
XShift = round(Rows/2 - XCentroid);
imPatch = circshift(imPatch, [XShift YShift]);

%% Phase info
if (contains(folderName,'vortex')||contains(folderName,'Vortex'))
  [YSizeImPatch, XSizeImPatch] = size(imPatch);
  Ph    = zeros(XSizeImPatch,YSizeImPatch);
  Full  = 2*pi;
  for ki = 1:XSizeImPatch
    for kj = 1:XSizeImPatch
          r = pixelRatio*sqrt((ki-XCentroid)^2 + (kj-YCentroid)^2);
          if r<3.5e-3
            Angle = angle((ki-XCentroid)+1i*(kj-YCentroid));
            Ph(ki,kj)=mod(Angle,Full)/Full*2*pi;
          else
            Ph(ki,kj)= 2*pi;
          end
    end
  end
  EPatch = EPatch .* exp(1i*Ph);
  imPatch = abs(EPatch.^2);
else
  imPatch = imPatch;
  EPatch = EPatch;
end
end


