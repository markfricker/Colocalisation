function golden = makeColocGolden(fijiExe)
%MAKECOLOCGOLDEN  Regenerate the Fiji Coloc 2 golden values for testColocConformance.
%
%   golden = makeColocGolden()
%   golden = makeColocGolden(fijiExe)
%
% 1. Writes every case of colocConformanceCases as images/<case>_ch1.tif and
%    images/<case>_ch2.tif (uint8).
% 2. Runs coloc2Golden.groovy headless in Fiji, which runs Coloc 2 on each pair
%    with both threshold-regression implementations (Costes = the dialog
%    default, and Bisection) and writes coloc2_golden.csv.
% 3. Saves coloc2_golden.mat: golden.values (table: case, implementation,
%    name, value) and golden.provenance (tool, versions, date, settings,
%    script), which testColocConformance asserts against.
%
% Run only when the cases change or to re-check against a newer Coloc 2;
% the tests read the stored .mat and never call Fiji. fijiExe defaults to
% the Fiji install used for the stored values.

if nargin < 1 || isempty(fijiExe)
    fijiExe = 'C:\Users\dops0035\Documents\Research\Fiji\fiji-windows-x64.exe';
end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));   % colocConformanceCases
imDir = fullfile(here, 'images');
if ~isfolder(imDir)
    mkdir(imDir);
end

C = colocConformanceCases();
for k = 1:numel(C)
    imwrite(C(k).ch1, fullfile(imDir, [C(k).name '_ch1.tif']));
    imwrite(C(k).ch2, fullfile(imDir, [C(k).name '_ch2.tif']));
end

csvFile = fullfile(here, 'coloc2_golden.csv');
if isfile(csvFile)
    delete(csvFile);
end
script = fullfile(here, 'coloc2Golden.groovy');
args = sprintf('inDir=''%s'',outFile=''%s''', strrep(imDir, '\', '/'), strrep(csvFile, '\', '/'));
% The Fiji launcher is a GUI executable: cmd's start /wait blocks until it exits.
cmd = sprintf('start "" /wait "%s" --headless --run "%s" "%s"', fijiExe, script, args);
status = system(cmd);
if status ~= 0 || ~isfile(csvFile)
    error('makeColocGolden:fiji', 'Fiji run failed (status %d); no %s.', status, csvFile);
end

txt = fileread(csvFile);
if contains(txt, sprintf('\nERROR,'))
    error('makeColocGolden:coloc2', 'Coloc 2 script reported an error:\n%s', txt);
end
firstLine = strtok(txt, newline);
% value is read as text (warnings and labels are strings) and parsed into num.
values = readtable(csvFile, 'CommentStyle', '#', 'Delimiter', ',', 'TextType', 'string', ...
    'ReadVariableNames', true, 'VariableNamingRule', 'preserve', 'Format', '%s%s%s%s');
values.Properties.VariableNames{1} = 'caseName';   % 'case' is a MATLAB keyword
for v = {'caseName', 'implementation', 'name', 'value'}
    values.(v{1}) = strip(values.(v{1}), '"');
end
values.num = str2double(values.value);

meta = strtrim(erase(firstLine, '#'));
golden.values = values;
golden.provenance = struct( ...
    'tool', 'Fiji Coloc 2 (sc.fiji.coloc.Coloc_2, Colocalisation_Analysis)', ...
    'versions', meta, ...
    'date', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), ...
    'settings', 'no ROI/mask; PSF 3; 100 Costes randomisations; Manders + Costes test on; both regression implementations', ...
    'script', 'tests/conformance/fixtures/coloc2Golden.groovy', ...
    'images', 'tests/conformance/fixtures/images (from colocConformanceCases)', ...
    'fiji', fijiExe);
save(fullfile(here, 'coloc2_golden.mat'), 'golden');
fprintf('makeColocGolden: %d values from %s\n', height(values), meta);
end
