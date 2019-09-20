function opt = setOptions(varargin)
% SETOPTIONS generates a valid argument "algOptions" for "calcLambda"
% function. It serves as an interface for easy setting of "algOptions".
%
% opt = setOptions()
% opt = setOptions(Name, Value)
%
% INPUT:
%   Name-Value Pair Arguments:
%   'perm'       :   vector; specifying iteration order
%   'precision'  :   vector or scalar: [tol maxIter], tol, maxIter
%      - tol     :   0 < tol < 1; tolerance to stop iteration
%      - maxIter :   integer > 1; force to stop after maxIter iterations
%   'errorType'  :   char or cell: {tolType, evalType}, tolType, evalType
%      - tolType :   char: 'abs' or 'rel'
%                    choose absolute/relative errors in stopping conditions
%      - evalType:   char: 'val' or 'var'
%                    choose to evaluate convergence of loss functions or variables
%   'initType'   :   char: 'fixed', 'random'
%      - 'fixed' :   use "perm" if provided, or default "perm"
%      - 'random':   use a random permutation in initialization
%
% OUTPUT:
%   opt          :   struct; valid for "algOptions"

% Copyright (c) 2015-2017, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last update on 20 Sep 2019


% default values
defaultPerm = [];
defaultPrecision= [1e-3, 20];
defaultErrorType= {'rel', 'var'};
defaultInitType = 'fixed';
defaultPenalty = 0;

% validating functions
validPerm = @(x) isnumeric(x) || isempty(x);
validPrecision = @(x) isnumeric(x) && isvector(x);
validErrorType = @(x) ischar(x) || iscellstr(x);
validInitType = @(x) any(strcmp(x, {'fixed', 'random'}));
validPenalty = @(x) (x==1) || (x==0);

% setup parser
parser = inputParser;
addParameter(parser, 'perm', defaultPerm, validPerm);
addParameter(parser, 'precision', defaultPrecision, validPrecision);
addParameter(parser, 'errorType', defaultErrorType, validErrorType);
addParameter(parser, 'initType', defaultInitType, validInitType);
addParameter(parser, 'penalty', defaultPenalty, validPenalty);

% parsing
parse(parser, varargin{:});
opt.perm = parser.Results.perm;
opt.precision = parser.Results.precision;
opt.errorType = parser.Results.errorType;
opt.initType = parser.Results.initType;
opt.penalty = parser.Results.penalty;