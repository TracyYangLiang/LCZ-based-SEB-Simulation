function run_all()
% run_all.m (function version)
% Root: D:\Paris Case Study\
% Code: D:\Paris Case Study\Case Simulation Code\
% Results: D:\Paris Case Study\Case Simulation Results\

clc;

% ====== Resolve paths ======
thisFile = mfilename('fullpath');
rootDir  = fileparts(thisFile);  % D:\Paris Case Study
codeDir  = fullfile(rootDir, 'Case Simulation Code');
resDir   = fullfile(rootDir, 'Case Simulation Results');

% ====== Basic checks ======
if ~isfolder(codeDir)
    error('Code folder not found: %s', codeDir);
end
if ~isfolder(resDir)
    warning('Results folder not found. Creating: %s', resDir);
    mkdir(resDir);
end

% ====== Set working directory & MATLAB path ======
cd(codeDir);
addpath(genpath(codeDir));

% Make paths available to the scripts (in base workspace)
assignin('base','rootDir',rootDir);
assignin('base','codeDir',codeDir);
assignin('base','resDir',resDir);

% ====== Script list (execution order) ======
scripts = { ...
    'SEB_Jan_useFunction_imp_veg_soil.m', ...
    'SEB_Feb_to_Dec_useFunction_imp_veg_soil.m', ...
    'Monthly_mean_composite_LST_Jan_Paris.m', ...
    'Monthly_mean_composite_LST_Feb_to_Dec_Paris.m' ...
};

fprintf('=== Running full workflow ===\n');
fprintf('Root dir   : %s\n', rootDir);
fprintf('Code dir   : %s\n', codeDir);
fprintf('Results dir: %s\n\n', resDir);

tAll = tic;

for i = 1:numel(scripts)
    s = scripts{i};

    if ~isfile(fullfile(codeDir, s))
        error('Missing script: %s (expected in %s)', s, codeDir);
    end

    fprintf('--- (%d/%d) Running: %s ---\n', i, numel(scripts), s);
    t = tic;

    try
        % Run in BASE workspace so script-internal "clear" won't kill runner variables
        evalin('base', sprintf('run(''%s'');', s));
    catch ME
        fprintf('\n[ERROR] Script failed: %s\n', s);
        fprintf('Message : %s\n', ME.message);
        if ~isempty(ME.stack)
            fprintf('Location: %s (line %d)\n\n', ME.stack(1).file, ME.stack(1).line);
        end
        rethrow(ME);
    end

    fprintf('Done: %s (%.1f s)\n\n', s, toc(t));
end

fprintf('=== All done. Total time: %.1f s ===\n', toc(tAll));

end

run_all;
%%
clear all, clc
A = load('LST_composite.mat');
B = load('LST_composite_jan.mat');
function A = injectJanIntoCompositeLayer1(A, B, fieldPattern)
% injectJanIntoCompositeLayer1
% Replace A.(field)(:,:,1) with B.(field) for matching fields.
%
% Usage:
%   A = injectJanIntoCompositeLayer1(A, B);
%   A = injectJanIntoCompositeLayer1(A, B, "LST_");   % optional pattern
%
% Inputs:
%   A, B : structs loaded by load()
%   fieldPattern (optional): only process fields containing this pattern
%
% Output:
%   A : updated struct

    if nargin < 3 || isempty(fieldPattern)
        fieldPattern = "LST_";  % default: only touch LST_* fields
    else
        fieldPattern = string(fieldPattern);
    end

    fA = string(fieldnames(A));
    fB = string(fieldnames(B));

    % intersection of field names
    common = intersect(fA, fB);

    % filter by pattern (e.g., "LST_")
    common = common(contains(common, fieldPattern));

    if isempty(common)
        error('No common fields found matching pattern "%s".', fieldPattern);
    end

    for i = 1:numel(common)
        fn = common(i);

        a = A.(fn);
        b = B.(fn);

        % Basic dimensional checks
        if ndims(a) < 3
            error('A.%s is not 3D (size=%s).', fn, mat2str(size(a)));
        end
        if ~ismatrix(b)
            error('B.%s is not 2D (size=%s).', fn, mat2str(size(b)));
        end

        % Size checks
        if size(a,1) ~= size(b,1) || size(a,2) ~= size(b,2)
            error('Size mismatch at %s: A is %s but B is %s.', ...
                fn, mat2str(size(a)), mat2str(size(b)));
        end

        % Inject B into the first layer of A
        a(:,:,1) = b;
        A.(fn) = a;
    end
end

A = injectJanIntoCompositeLayer1(A, B); 
X = A.LST_1200_final;
%if you want to see other data, you can change"A.LST_1200_final" to "A.LST_1300_final"., 
%the data can be changed from 08:00 to 24:00 everydays.
figure('Color','w');
clim = [min(X(:),[],'omitnan') max(X(:),[],'omitnan')];

for m = 1:12
    subplot(3,4,m);
    imagesc(X(:,:,m));
    axis image off;
    title(sprintf(' %d th Months data', m));
    caxis([273,320]); colorbar;
end