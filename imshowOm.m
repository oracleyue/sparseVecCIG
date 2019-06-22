function imshowOm(X, borderWidth)
% IMSHOWM convert a matrix into grayscale image and show it.
%
% Inputs:
%    X : matrix
%    borderWidth : integer; border width in pixels

% Copyright [2019] <oracleyue>
% Last modified on 21 Jun 2019


figure;
set(gcf,'color','white'); borderWidth = 1;
imOm = mat2gray(full(abs(X)));
imbOm = addborder(2*imOm, borderWidth, 1, 'outer');
imshow(imbOm, 'InitialMagnification','fit');
colormap(1-colormap('gray'));

% END of imshowM()