function [low,high] = TorPlot(data,names,sensitivity1,sensitivity2,save,fh)
% TorPlot  A function for generating tornado plots. This program calculates
% a baseline value using the function handle provided, or defaulting to 
% summing the vector 'data'. Each element of 'data' is
% then changed by a certain sensitivity value and the elements in 'data'
% are calculated again, with only the current single element having been
% changed.
%
% Sensitivity plots are useful in assessing importance and weight in
% engineering systems. Preliminary assessments benefit from sensitivity
% analysis, mainly because they highlight the driving factors.
%
%   INPUTS:
%       data = baseline data entered as a vector
%       names = tick mark names for the y axis
%       sensitivity = amount each element is individually changed
%       save = flag determining whether or not to save the plot
%       fh = function handle, ie: @function_name
%
%   OUTPUTS:
%       low = lower sensitivities
%       high = higher sensitivies
%
%   EXAMPLE (Sum by default):
%       
%       data=[0.24,2.55,0.60,0.02,0.09,0.36,0.18];
%       names = {'Distribution';'Crude';'Refinery';'Storage';
%           'Sales Tax';'State Excise Tax';'Federal Excise Tax'};
%       [low, high] = TorPlot(data, names, 0.15, true);
%
%
%   EXAMPLE (With a function):
%       
%       % This creates a matlab file called 'my_fun.m' which
%       % simply returns the cost function related to a mock
%       % manufacturing process based on weight, cost, and transportation
%       fileID=fopen('my_fun.m','w');
%       fprintf(fileID,'function fval=my_fun(InputVector) \r\n');
%       fprintf(fileID,'weight=InputVector(1);');
%       fprintf(fileID,'cost=InputVector(2);');
%       fprintf(fileID,'transportation=InputVector(3);');
%       fprintf(fileID,'fval = 0.1*weight+0.2*cost+0.3*transportation+0.01*weight^2+0.02*cost^2+0.03*transportation^2;');
%       fclose(fileID);
%
%       % Sensitivities are calculated using this function 'my_fun.m'
%       data=[1,2,3];
%       names = {'Weight';'Cost';'Transportation'};
%       [low, high] = TorPlot(data, names, 0.15, true,@my_fun);
%
%       It's important that the function inputs a VECTOR and not
%       individual elements or a matrix, and that it returns a SINGLE value
%
%   Copyright 2013 Richard C.J. McCulloch
%   $Revision.1 $  $Date: 27-Jun-2013 $

if isempty(data) 
    error(message('No data was passed to the function TorPlot')); 
end

if nargin < 4, save=0; % if save isn't specified, don't save
    if nargin < 3, sensitivity1=0.2; % if no sensitivity is given assume 20%
        if nargin < 2
            names = repmat(char(0),length(data),1);
        end
    end
end

% Initialize low and high matricies for speed
Objective_low=zeros(1,length(names));
Objective_high=zeros(1,length(names));
Objective_low_sum=zeros(1,length(names));
Objective_high_sum=zeros(1,length(names));
low=zeros(1,length(names));
high=zeros(1,length(names));

% Calculate the base change due to a single factor at a time
if nargin < 5
    for i=1:length(names)
        Objective_low=data;
        Objective_high=data;
        Objective_low(i)=Objective_low(i)*(1-sensitivity1);
        Objective_high(i)=Objective_high(i)*(1+sensitivity2);
        Objective_low_sum(i)=sum(Objective_low);
        Objective_high_sum(i)=sum(Objective_high);
        low(i)=Objective_low_sum(i);
        high(i)=Objective_high_sum(i);    
        
        % The base value is where the y axis is centered
        Objective_base_value=sum(data);
    end
else
    for i=1:length(names)
        Objective_low=data;
        Objective_high=data;
        Objective_low(i)=Objective_low(i)*(1-sensitivity1);
        Objective_high(i)=Objective_high(i)*(1+sensitivity2);
        Objective_low_sum(i)=fh(Objective_low);
        Objective_high_sum(i)=fh(Objective_high);
        low(i)=Objective_low_sum(i);
        high(i)=Objective_high_sum(i); 
        
                % The base value is where the y axis is centered
        Objective_base_value=fh(data);
    end
end

% Sort the values based on the lower change
% Sort the higher values and the names arrays
%    using the same indices
%[Objective_low_sum,ind]=sort(Objective_low_sum,'descend');
%Objective_high_sum=Objective_high_sum(ind);
%names_Objective=names(ind);
names_Objective=names;
% Create a figure and plot the low and high horizontally
figure
h = barh(Objective_high_sum);
hold on
xmin=min([min(Objective_low_sum),min(Objective_high_sum)]);
xmax=max([max(Objective_low_sum),max(Objective_high_sum)]);
xlim([1.025*xmin 0.975*xmax])
barh(Objective_low_sum,'r')
bh = get(h,'BaseLine');
set(bh,'BaseValue',Objective_base_value);
title('Sensitivities')
if nargin > 1
    set(gca,'yticklabel',names)
    set(gca,'Ytick',[1:length(names)],'YTickLabel',[1:length(names)])
    set(gca,'yticklabel',names_Objective)
end
xlabel('NPV ($ millions)')
if(save)
    saveas(gcf,'Objective.png')
end