function imshowOm(X, type, dL)
% IMSHOWM convert a matrix into grayscale image and show it.
%
% Inputs:
%    X : matrix
%    type : 'raw' or 'spy'
%    dL : ordered sizes of blocks when X is a partitioned matrix
%
% Note:
%    For convenience, you may use either type or dL as the second
%    argument. For example,
%       imshowOm(X, 'raw')
%       imshowOm(X, [2 3 4 2 1])

% Copyright [2019] <oracleyue>
% Last modified on 28 Jun 2019

if nargin < 2
    type = 'raw';
    dL = [];
end
if nargin == 2  % use "type" or "dL" as the second argument
    if ischar(type)
        dL = [];
    else
        dL = type;
        type = 'raw';
    end
end
assert(any(strcmpi({'raw', 'spy'}, type)), ...
       'The argument "type" must be "raw" or "spy".');

% fill matrix with ones
X = full(abs(X));
if ~isempty(dL)
    d = sum(dL);
    p = length(dL);
    assert(d == size(X, 1), 'The dimension of X fails to match dL!')

    for i = 1:p
        for j = 1:p
            di = dL(i); dj = dL(j);
            iIdx = sum(dL(1:i-1))+1:sum(dL(1:i));
            jIdx = sum(dL(1:j-1))+1:sum(dL(1:j));
            if i == j || sum(sum(X(iIdx, jIdx)))
                X(iIdx, jIdx) = ones(di, dj);
            end
        end
    end
end

switch type
  case 'raw'
    % plot matrix as gray-scale image
    set(gcf,'color','white');
    imOm = mat2gray(X);
    imbOm = addBorder(2*imOm, 1, 1);
    imshow(imbOm, 'InitialMagnification','fit');
    colormap(1-colormap('gray'));
  case 'spy'
    % investiage nonzeros in matrix
    spy(X);
end

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
