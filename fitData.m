%% Fit curves (using bootstrapping method)
% by N Mikolajewicz (13.03.19)
clear all; close all;

%% specify inputs and analysis parameters

file = 'liverGrowth.xlsx';          % name of excel file where input data is stored (e.g., 'liverGrowth.xlsx')
sheet = 'DOBDOD';                   % name of excel sheet (e.g., 'Sheet1')
model = 'hyperbolic';               % specify which model to fit to data (linear/hyperbolic) (e.g., 'hyperbolic)
saveResults = true;                 % specify whether results are saved to excel (true/false) (e.g., true)
nBootStraps = 100;                  % specify number of sampling iterations (e.g., 100)
ageOutputUnits = 2;                 % 1 = days, 2 = weeks, 3 = months

%% do not modify

warning('off','all')
fprintf('\nimporting data...\n');
switch model
    case 'hyperbolic'
        fprintf('\ndata will be fit to HYPERBOLIC FUNCTION: \n');
        fprintf('\ty=(a*x)/(b+x) \n');
        fprintf('\twhere \n\t\ta = max value (at plateau) \n\t\tb = x value at half-max of y \n')
    case 'linear'
        fprintf('\ndata will be fit to LINEAR FUNCTION: \n');
        fprintf('\ty=a*x + b: \n');
        fprintf('\twhere \n\t\ta = slope \n\t\tb = intercept \n')      
end

%import data from excel spreadsheet. Refer to template for spreadsheet format.
data = importData(file, sheet);

% check how data input was presented
check(1) = isfield(data,'age');
check(2) = isfield(data,'DOB');
check(3) = isfield(data,'DOD');

if check(1)
    go = 1;
elseif check(2)+check(3) == 2
    go = 2;
    for i = 1:length(data)
        data(i).age = datenum(data(i).DOD) - datenum(data(i).DOB);
        if ageOutputUnits == 2
            data(i).age = data(i).age/7;
        elseif ageOutputUnits == 3
            data(i).age = data(i).age/31;
        end
    end
end

% recode categorical variables to numerical variables
categoricalVar = {'condition'};
[data, codingLegend] = cat2num(data, categoricalVar);
p=numSubplots(length(unique([data.condition])));
uCond = unique([data.condition]);

maxYscale = [];
maxXscale = [];
for i = 1:length(uCond)
    fprintf('\n fitting curves to %s...', codingLegend.condition(i).label);
    MASTER(i).condition = codingLegend.condition(i).label;
    MASTER(i).pData = [data([data.condition]==i)];
    
    x{i} = [MASTER(i).pData.age];
    y{i} = [MASTER(i).pData.outcome];
    
    xx = x{i} ;
    yy = y{i} ;
    
    bootStrapN = 5;
    
    xi = []; yi = [];
    for j = 1:nBootStraps
        
        ind = randsample(length(x{i}), length(x{i}), true);
        subX = x{i}(ind);
        subY = y{i}(ind);
        
        try;
            
            switch model
                case 'hyperbolic'
                    [fitresult{i,j}, gof{i,j}] = hyperbolicFit(subX, subY);
                    coef = coeffvalues(fitresult{i,j});
                    a(i,j) = coef(1);
                    b(i,j) = coef(2);
                    xi(j,:) = [0:max(x{i})];
                    yi(j,:) = (a(i,j).*xi(j,:))./(b(i,j)+xi(j,:));
                case 'linear'
                    [fitresult{i,j}, gof{i,j}] = linearFit(subX, subY);
                    coef = coeffvalues(fitresult{i,j});
                    a(i,j) = coef(1);
                    b(i,j) = coef(2);
                    xi(j,:) = [0:max(x{i})];
                    yi(j,:) = a(i,j).*xi(j,:)+b(i,j);
            end
        catch e; disp(e);
        end
    end
    
    dataStore(i).xi =  xi;
    dataStore(i).yi =  yi;
    
    maxYscale = [maxYscale max(y{i})];
    maxXscale = [maxXscale max(x{i})];
end

maxYscale = max(maxYscale) + 0.1*max(maxYscale);
maxXscale = max(maxXscale);

fprintf('\n plotting results...');

pp=numSubplots(2*length(unique([data.condition])));
offset = length(unique([data.condition]));
minA = []; maxA = [];
minB = []; maxB = [];
colorPal = {'b', 'r', 'm', 'g', 'c', 'k'}; 
for i = 1:length(uCond)
    % plot histrograms (visualize distribution of sampled regression parameters)
    figure(2);
    %     subplot(p(2), p(1), codingLegend.condition(i).code);
    subplot(pp(2), pp(1), 1+(2*(i-1)));
    outIndA{i} = isoutlier(a(i,:));
    aTrim{i} = a(i,~isoutlier(a(i,:)));
    h1{i} = histogram( aTrim{i});
    set(h1{i},'FaceColor','b');
    minA = [minA min(h1{i}.BinLimits)]; maxA = [maxA max(h1{i}.BinLimits)];
    title([codingLegend.condition(i).label ': coefficient a']);
    xlabel('coefficient a'); ylabel('count');
    
    figure(2);
    %     subplot(p(2), p(1), codingLegend.condition(i).code);
    subplot(pp(2), pp(1), (2*(i)));
    outIndB{i} = isoutlier(b(i,:));
    bTrim{i} = b(i,~isoutlier(b(i,:)));
    h2{i} =histogram( bTrim{i});
     minB = [minB min(h2{i}.BinLimits)]; maxB = [maxB max(h2{i}.BinLimits)];
    set(h2{i},'FaceColor','r');
    title([codingLegend.condition(i).label ': coefficient b']);
    xlabel('coefficient b'); ylabel('count');
    
    %remove outliers
    outInd{i} = isoutlier(a(i,:)) | isoutlier(b(i,:));
    subXTrim = dataStore(i).xi(~outInd{i},:);
    subYTrim = dataStore(i).yi(~outInd{i},:);
    figure(3);  subplot(p(2), p(1), codingLegend.condition(i).code); hold on;
    % plot (subXTrim(1,:), subYTrim,'color',[0,0,0]+0.8);
    xiMean = mean(subXTrim,1);  yiMean = mean(subYTrim,1);
    xiSTD = std(subXTrim,1);  yiSTD = std(subYTrim,1);
    
    %confidence and prediction intervals for regression curves
    xiSEMc=  xiSTD/sqrt(length(x{i})-1); yiSEMc = yiSTD/sqrt(length(y{i})-1);
    xiSEMp = xiSTD*sqrt(1+(1/(length(x{i})-2)));    yiSEMp = yiSTD*sqrt(1+(1/(length(y{i})-2)));
    yiCIWidthc = yiSEMc*tinv(1-0.05/2,(length(x{i})-2)); yiCIWidthp = yiSEMp*tinv(1-0.05/2,(length(x{i})-2));
    
    % plot curves
    pl1{i} = plot (xiMean, yiMean,'r');
    pl2{i} = ciplot(yiMean-yiCIWidthc,yiMean+yiCIWidthc, xiMean, 'r',  0.35);
    pl3{i} = ciplot(yiMean-yiCIWidthp,yiMean+yiCIWidthp, xiMean, 'r',  0.15);
    pl4{i} = plot (x{i}, y{i},'ro');
    
    if length(uCond) == i
    legend([pl1{i} pl2{i} pl3{i} pl4{i}], {'fit', '95% CI', '95% PI', 'data'}, 'Location', 'Best'); 
    end
        
    title(codingLegend.condition(i).label)
    xlabel('x'); ylabel('y'); ylim([0 maxYscale]); xlim([0 maxXscale]);
    
    figure(4); hold on;
%     plot (xiMean, yiMean); hold on;   
      
    xPl{i} = xiMean; yPl{i} = yiMean;
    h3{i} = ciplot(yiMean-yiCIWidthc,yiMean+yiCIWidthc, xiMean, colorPal{i},  0.35);
    legendNames{i} = codingLegend.condition(i).label;
    
    % save results to data structure
    results(i).condition = codingLegend.condition(i).label;
    results(i).aMean = mean(aTrim{i});
    results(i).aMedian = median(aTrim{i});
    results(i).aStd = std(aTrim{i});
    results(i).aSem = results(i).aStd/sqrt(length(x{i})-2);
      
    results(i).bMean = mean(bTrim{i});
    results(i).bMedian = median(bTrim{i});
    results(i).bStd = std(bTrim{i});
    results(i).bSem = results(i).bStd/sqrt(length(x{i})-2);
    
    results(i).SampleSize = length(x{i});
end

minA = min(minA); maxA = max(maxA); 
minB = min(minB); maxB = max(maxB); 
for i = 1:length(uCond)
    figure(2);
    subplot(pp(2), pp(1), 1+(2*(i-1)));
    xlim([minA maxA])
    subplot(pp(2), pp(1), (2*(i)));
    xlim([minB maxB])
    
 figure(4); hold on;
  plot(xPl{i}, yPl{i}, colorPal{i});
end
legend([h3{:}],legendNames, 'Location', 'Best');
xlabel('x'); ylabel('y'); 
title('Compare Results (+/- 95% CI)')

% statistical analysis
combo = combnk(1:length(MASTER),2);

for i = 1:size(combo,1)
    statsA(i).condition1 = codingLegend.condition(combo(i,1)).label;
    statsA(i).condition2 = codingLegend.condition(combo(i,2)).label;
    statsA(i).mean1 = results(combo(i,1)).aMean;
    statsA(i).mean2 = results(combo(i,2)).aMean;
    statsA(i).sampleSize1 = results(combo(i,1)).SampleSize;
    statsA(i).sampleSize2 = results(combo(i,2)).SampleSize;
    statsA(i).sem1 = results(combo(i,1)).aSem;
    statsA(i).sem2 = results(combo(i,2)).aSem;
    statsA(i).pooledSEM = sqrt( statsA(i).sem1^2 + statsA(i).sem2^2);
    statsA(i).tscore = abs(statsA(i).mean1 - statsA(i).mean2)/statsA(i).pooledSEM;
    statsA(i).p = 2*(1-(tcdf(statsA(i).tscore,statsA(i).sampleSize1+statsA(i).sampleSize2-4)));
    if statsA(i).p > 1; statsA(i).p = 1; end   

    statsB(i).condition1 = codingLegend.condition(combo(i,1)).label;
    statsB(i).condition2 = codingLegend.condition(combo(i,2)).label;
    statsB(i).mean1 = results(combo(i,1)).bMean;
    statsB(i).mean2 = results(combo(i,2)).bMean;
    statsB(i).sampleSize1 = results(combo(i,1)).SampleSize;
    statsB(i).sampleSize2 = results(combo(i,2)).SampleSize;
    statsB(i).sem1 = results(combo(i,1)).bSem;
    statsB(i).sem2 = results(combo(i,2)).bSem;
    statsB(i).pooledSEM = sqrt( statsB(i).sem1^2 + statsB(i).sem2^2);
    statsB(i).tscore = abs(statsB(i).mean1 - statsB(i).mean2)/statsB(i).pooledSEM;
    statsB(i).p = 2*(1-(tcdf(statsB(i).tscore,statsB(i).sampleSize1+statsB(i).sampleSize2-4)));
    if statsB(i).p > 1; statsB(i).p = 1; end
    
end

if saveResults
    try;
        fprintf('\n saving results...')
        % get current time
        now = datestr(datetime('now'));
        now = strrep(now, ':','-');
        
        % assign file name for results
        file = strrep(file, '.xlsx',' ');
        file = [file 'RESULTS ' now '.xlsx'];
        
        % save results to excel file
        writetable(struct2table(results), file, 'Sheet' ,['Summary']);
        writetable(struct2table(statsA), file, 'Sheet' ,['Statistics_regression coef A']);
        writetable(struct2table(statsB), file, 'Sheet' ,['Statistics_regression coef B']);
        
    catch;
        error(' ')
    end
    % remove first 3 blank sheets in new excel file (Sheet 1, Sheet 2, Sheet 3)
    RemoveSheet123(file);
else
    fprintf('\n results not saved...')
end
fprintf('\n Complete! \n');

function [fitresult, gof] = hyperbolicFit(x, y)


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
% ft = fittype( '(a*(x^n))/((b^n)+(x^n))', 'independent', 'x', 'dependent', 'y' ); % hill function
ft = fittype( '(a*x)/(b+x)', 'independent', 'x', 'dependent', 'y' ); % hyperbolic function
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% opts.StartPoint = [1 2 5];
opts.Lower = [0 0];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

end

function [data, codingLegend] = cat2num(data, categoricalVar);
for i = 1:length(categoricalVar)
    
    % define coding scheme
    a= categorical(cellstr({data.(categoricalVar{i})}));
    [b,GN,~]= grp2idx(a);
    
    % create coding legend
    for j = 1:length(GN)
        codingLegend.(categoricalVar{i})(j).code = j;
        codingLegend.(categoricalVar{i})(j).label = GN{j};
    end
    
    % assign codes to dataset
    for j = 1:length(data)
        data(j).(categoricalVar{i}) = b(j);
    end
    
end
end

function data = importData(file, sheet)

% import spreadsheet
[num,txt,raw]=xlsread(file,sheet);

% extract data and store in data structure
for i = 2:size(raw,1);
    for j = 1:size(raw,2)
        if ~isnan(raw{1,j}); data(i-1).(raw{1,j}) = raw{i,j}; end
    end; end
end


function [p,n]=numSubplots(n)
% function [p,n]=numSubplots(n)
%
% Purpose
% Calculate how many rows and columns of sub-plots are needed to
% neatly display n subplots.
%
% Inputs
% n - the desired number of subplots.
%
% Outputs
% p - a vector length 2 defining the number of rows and number of
%     columns required to show n plots.
% [ n - the current number of subplots. This output is used only by
%       this function for a recursive call.]
%
%
%
% Example: neatly lay out 13 sub-plots
% >> p=numSubplots(13)
% p =
%     3   5
% for i=1:13; subplot(p(1),p(2),i), pcolor(rand(10)), end
%
%
% Rob Campbell - January 2010


while isprime(n) & n>4,
    n=n+1;
end
p=factor(n);
if length(p)==1
    p=[1,p];
    return
end
while length(p)>2
    if length(p)>=4
        p(1)=p(1)*p(end-1);
        p(2)=p(2)*p(end);
        p(end-1:end)=[];
    else
        p(1)=p(1)*p(2);
        p(2)=[];
    end
    p=sort(p);
end
%Reformat if the column/row ratio is too large: we want a roughly
%square design
while p(2)/p(1)>2.5
    N=n+1;
    [p,n]=numSubplots(N); %Recursive!
end
end


function h1 = ciplot(lower,upper,x,colour, transparency);

% ciplot(lower,upper)
% ciplot(lower,upper,x)
% ciplot(lower,upper,x,colour)
%
% Plots a shaded region on a graph between specified lower and upper confidence intervals (L and U).
% l and u must be vectors of the same length.
% Uses the 'fill' function, not 'area'. Therefore multiple shaded plots
% can be overlayed without a problem. Make them transparent for total visibility.
% x data can be specified, otherwise plots against index values.
% colour can be specified (eg 'k'). Defaults to blue.
% Raymond Reynolds 24/11/06
if length(lower)~=length(upper)
    error('lower and upper vectors must be same length')
end
if nargin<4
    colour='b';
end
if nargin<3
    x=1:length(lower);
end
% convert to row vectors so fliplr can work
if find(size(x)==(max(size(x))))<2
    x=x'; end
if find(size(lower)==(max(size(lower))))<2
    lower=lower'; end
if find(size(upper)==(max(size(upper))))<2
    upper=upper'; end
h1 = fill([x fliplr(x)],[upper fliplr(lower)],colour);
set(h1,'EdgeColor','none');
set(h1, 'facealpha', transparency);
end


function [fitresult, gof] = linearFit(xx, yy)
%CREATEFIT(XX,YY)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : xx
%      Y Output: yy
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 10-Mar-2019 17:10:04


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( xx, yy );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );


end


function RemoveSheet123(excelFileName,sheetName)
% RemoveSheet123 - removes the sheets that are automatically added to excel
% file.
% When Matlab writes data to a new Excel file, the Excel software
% automatically creates 3 sheets (the names are depended on the user
% languade). This appears even if the user has defined the sheet name to be
% added.
%
% Usage:
% RemoveSheet123(excelFileName) - remove "sheet1", "sheet2","sheet3" from
% the excel file. excelFileName is a string of the Excel file name.
% RemoveSheet123(excelFileName,sheetName) - enables the user to enter the
% sheet name when the language is other than English.
% sheetName is the default sheet name, without the number.
%
%
%                       Written by Noam Greenboim
%                       www.perigee.co.il
%
%% check input arguments
if nargin < 1 || isempty(excelFileName)
    error('Filename must be specified.');
end
if ~ischar(excelFileName)
    error('Filename must be a string.');
end
try
    excelFileName = validpath(excelFileName);
catch
    error('File not found.');
end
if nargin < 2
    sheetName = 'Sheet'; % EN: Sheet, DE: Tabelle, HE: ?????? , etc. (Lang. dependent)
else
    if ~ischar(sheetName)
        error('Default sheet name must be a string.');
    end
end
%%
% Open Excel file.
objExcel = actxserver('Excel.Application');
objExcel.Workbooks.Open(excelFileName); % Full path is necessary!
% Delete sheets.
try
    % Throws an error if the sheets do not exist.
    objExcel.ActiveWorkbook.Worksheets.Item([sheetName '1']).Delete;
    %     fprintf('\nsheet #1 - deleted.')
    objExcel.ActiveWorkbook.Worksheets.Item([sheetName '2']).Delete;
    %     fprintf('\nsheet #2 - deleted.')
    objExcel.ActiveWorkbook.Worksheets.Item([sheetName '3']).Delete;
    %     fprintf('\nsheet #3 - deleted.\n')
catch
    fprintf('\n')
    O=objExcel.ActiveWorkbook.Worksheets.get;
    if O.Count==1
        error('Can''t delete the last sheet. Excel file must containt at least one sheet.')
    else
        warning('Problem occured. Check excel file.');
    end
end
% Save, close and clean up.
objExcel.ActiveWorkbook.Save;
objExcel.ActiveWorkbook.Close;
objExcel.Quit;
objExcel.delete;
end
function filenameOut = validpath(filename)
% VALIDPATH builds a full path from a partial path specification
%   FILENAME = VALIDPATH(FILENAME) returns a string vector containing full
%   path to a file. FILENAME is string vector containing a partial path
%   ending in a file or directory name. May contain ..\  or ../ or \\. The
%   current directory (pwd) is prepended to create a full path if
%   necessary. On UNIX, when the path starts with a tilde, '~', then the
%   current directory is not prepended.
%
%   See also XLSREAD, XLSWRITE, XLSFINFO.

%   Copyright 1984-2012 The MathWorks, Inc.

%First check for wild cards, since that is not supported.
if strfind(filename, '*') > 0
    error(message('MATLAB:xlsread:Wildcard', filename));
end

% break partial path in to file path parts.
[Directory, file, ext] = fileparts(filename);
if ~isempty(ext)
    filenameOut = getFullName(filename);
else
    extIn = matlab.io.internal.xlsreadSupportedExtensions;
    for i=1:length(extIn)
        try                                                                %#ok<TRYNC>
            filenameOut = getFullName(fullfile(Directory, [file, extIn{i}]));
            return;
        end
    end
    error(message('MATLAB:xlsread:FileDoesNotExist', filename));
end
end
function absolutepath=abspath(partialpath)

% parse partial path into path parts
[pathname, filename, ext] = fileparts(partialpath);
% no path qualification is present in partial path; assume parent is pwd, except
% when path string starts with '~' or is identical to '~'.
if isempty(pathname) && strncmp('~', partialpath, 1)
    Directory = pwd;
elseif isempty(regexp(partialpath,'(.:|\\\\)', 'once')) && ...
        ~strncmp('/', partialpath, 1) && ...
        ~strncmp('~', partialpath, 1);
    % path did not start with any of drive name, UNC path or '~'.
    Directory = [pwd,filesep,pathname];
else
    % path content present in partial path; assume relative to current directory,
    % or absolute.
    Directory = pathname;
end

% construct absolute filename
absolutepath = fullfile(Directory,[filename,ext]);
end
function filename = getFullName(filename)
FileOnPath = which(filename);
if isempty(FileOnPath)
    % construct full path to source file
    filename = abspath(filename);
    if isempty(dir(filename)) && ~isdir(filename)
        % file does not exist. Terminate importation of file.
        error(message('MATLAB:xlsread:FileDoesNotExist', filename));
    end
else
    filename = FileOnPath;
end
end