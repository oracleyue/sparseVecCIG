function imshowOm(X, bWidth)
% IMSHOWM convert a matrix into grayscale image and show it.
%
% Inputs:
%    X : matrix
%    borderWidth : integer; border width in pixels

% Copyright [2019] <oracleyue>
% Last modified on 21 Jun 2019

if nargin < 2
    bWidth = 1;
end

set(gcf,'color','white');
imOm = mat2gray(full(abs(X)));
imbOm = addborder(2*imOm, bWidth, 1, 'outer');
imshow(imbOm, 'InitialMagnification','fit');
colormap(1-colormap('gray'));

end  % END of imshowM

function img2 = addBorder(img1, t, c)
% Add borders to image.
% Modified based on "addborder.m" by Eric C. Johnson, 7-Aug-2008.

    [nr1 nc1 d] = size(img1);
    % Add the border thicknesses to the total image size
    nr2 = nr1 + 2*t;
    nc2 = nc1 + 2*t;
    % Create an empty matrix, filled with the border color.
    img2 = repmat(c, [nr2 nc2 1]);
    % Copy IMG1 to the inner portion of the image.
    img2( (t+1):(nr2-t), (t+1):(nc2-t), : ) = img1;

end  % END of addBorder
