function outfig = labelplot(fig,varargin)
% Function for labeling plots and making them more readable on documents
%
% Inputs: 
%       fig = figure handle, 
%       varargin = xlabel, ylabel, title, 1 or 0 indicating if a
%                  legend is needed, zlabel (if 3D plot)
% Outputs: 
%       outfig = updated figure handle
% Author: Derrick Choi
%
% Date created: Feb. 6, 2020
% Date last modified: Jan. 30 2022

% Assign inputs
xlab = varargin{1};
ylab = varargin{2};
T = varargin{3};
% Adjust figure position
fig.Position = [400 100 800 600];
% Grab Axes handle from figure handle and start labelling
p = fig.Children;
p.XLabel.Interpreter = 'latex';
p.YLabel.Interpreter = 'latex';
p.XLabel.FontSize = 14;
p.YLabel.FontSize = 14;
p.XLabel.String = xlab;
p.YLabel.String = ylab;
p.Title.Interpreter = 'latex';
p.Title.String = ['\textbf{',T,'}'];
p.Title.FontSize = 16;
p.XAxis.TickLabelInterpreter = 'latex';
p.YAxis.TickLabelInterpreter = 'latex';
p.XAxis.FontSize = 14;
p.YAxis.FontSize = 14;
p.XAxis.MinorTick = 'on';
p.YAxis.MinorTick = 'on';
%add z-label if 3D plot;
if length(varargin) ~= 4
    p.ZLabel.Interpreter = 'latex';
    p.ZLabel.String = varargin{4};
    p.ZLabel.FontSize = 14;
    p.ZAxis.TickLabelInterpreter = 'latex';
    p.ZAxis.FontSize = 14;
    p.ZAxis.MinorTick = 'on';
end
grid on
grid minor

% legend
if varargin{end-1} == 1
lgd = legend(varargin{end});
lgd.Location = 'best';
lgd.FontSize = 12;
lgd.Interpreter = 'latex';
end
% output
outfig = fig;

end