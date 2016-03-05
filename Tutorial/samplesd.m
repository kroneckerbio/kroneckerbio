function [sigma,dsigmady,d2sigmady2] = samplesd(t,yInd,yVal)

sdperc = 0.05;

sigma = max(sdperc*yVal,1);

if sdperc*yVal >= 1
    dsigmady = sdperc + zeros(size(yVal));
else
    dsigmady = zeros(size(yVal));
end

d2sigmady2 = zeros(size(yVal));


end

