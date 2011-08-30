function [ hist ] = create_hist( x )
%create_hist Create histogram according to the distance data
%   The method is:
%	Step1: For one location, compute structure tensor distances between 
%           each frame in a ???$T$???-long video clip and a reference 
%           frame.
%	Step2: Truncated the range of the distances to [0,M], M is the 95% percentile.
%	Step3: Group the distances to ???$NB$??? bins from the smallest to the largest.
%	Step4: Calculate the probability for each bin.
%	Step5: Repeat Step 1 to Step 4 ???$RT$??? times for another video with
%           same activity pattern at the location, and return the average histogram.
% x
%   The distances between each frame of the video to a reference frame
% hist
%   N*3 matrix containing the histogram
%   Every row is [LowBound UpperBound Percent]

%   Just keep those values larger than 95% percentile
prc95 = prctile(x, 95);
xt = x(x<prc95);
nElm = size(xt,2);
minxt=min(xt);
maxxt=max(xt);
delta=(maxxt-minxt)/(length(xt)-1);

lower=minxt-delta/2;
upper=maxxt+delta/2;
ncell = ceil(sqrt(length(xt)));
space = (upper-lower)/ncell;

hist(1:ncell, 1:3)=0;

y=round( (xt-lower)/(upper-lower)*ncell + 1/2 );
for n=1:nElm
  index=y(n);
  if index >= 1 && index<=ncell
      if index == 1
          hist(index, 1) = 0;
      else
          hist(index, 1) = lower+(index-1)*space;
      end
      hist(index, 2) = lower+(index)*space;
      hist(index, 3) = hist(index, 3)+1;
  end;
end;

figure;
bar(hist(:,3),'hist');

end

