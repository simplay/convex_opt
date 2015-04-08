function [ red, green, blue ] = getMasks(m, n)
%GETMASKS Bayer color filters as n x m matrices.
    % @param m, n are natural numbers > 0. specifying the mask dimensions
    % @return array of red, green, and blue masks.
    
    blue = zeros(m,n);
    blue(1:2:end,1:2:end) = 1;
    red = zeros(m,n);
    red(2:2:end,2:2:end) = 1;
    green = ones(m,n);
    green(1:2:end,1:2:end) = 0;
    green(2:2:end,2:2:end) = 0;

end

