function [ bestLambda ] = findBestLambda(im, mosaiced, Omega, iter )
%FINDBESTLAMBDA Summary of this function goes here
%   Detailed explanation goes here

    stepSize = 5;
    error = zeros(1,420);
    parfor k=1:420,
        l = stepSize*k
        abc = demosaic(mosaiced, Omega, l, iter);
        ssd = (im-abc).^2;
        error(k) = sum(abs(ssd(:)));
    end
    
    % show plot
    [~,ind] = min(error);
    plot((1:420)*stepSize, error);

    % find best candidate within retrieved range
    error2 = zeros(1,2*stepSize);
    parfor k=1:2*stepSize,
        l = ind*stepSize-stepSize+k
        abc = demosaic(mosaiced, Omega, l, iter);
        ssd = (im-abc).^2;
        error2(k) = sum(abs(ssd(:)));
    end
    [~,ind2] = min(error2);
    bestLambda = ind*stepSize-stepSize+ind2;
end

