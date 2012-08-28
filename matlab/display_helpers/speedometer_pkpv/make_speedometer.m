function [plot_handle] = make_speedometer( vin, vref, fh )
% 
% Draw a speedometer chart 
%
% Assumes you have jcontrol and jfreecharts added to path, 
% for example like:
%
% javaaddpath([pwd '/jfreechart-1.0.14/lib/jcommon-1.0.17.jar'])
% javaaddpath([pwd '/jfreechart-1.0.14/lib/jfreechart-1.0.14.jar'])
%
% returns a plot_handle, 
% which you pass to "set_vin_arrow" and "set_vref_arrow"
%


%% Dialdemo4
import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.GradientPaint;
import java.awt.Point;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.general.DefaultValueDataset;
import org.jfree.chart.plot.dial.DialBackground;
import org.jfree.chart.plot.dial.DialPlot;
import org.jfree.chart.plot.dial.DialCap;
import org.jfree.chart.plot.dial.DialPointer;
import org.jfree.chart.plot.dial.StandardDialFrame;
import org.jfree.chart.plot.dial.StandardDialScale;
import org.jfree.chart.plot.dial.StandardDialRange;
import org.jfree.chart.plot.dial.DialValueIndicator;
import org.jfree.chart.plot.dial.DialTextAnnotation;
import org.jfree.ui.GradientPaintTransformType;
import org.jfree.ui.StandardGradientPaintTransformer;


%% Start

if ~exist('vin','var')
  vin = 0.0;
end
if ~exist('vref','var')
  vref= 70.0;
end

% Create Datasets
this.dataset1 = DefaultValueDataset(  vin );
this.dataset2 = DefaultValueDataset(  vref);
        
plot = DialPlot();                                      % Create DialplotObject
plot.setView(0.0, 0.0, 1.0, 1.0);                       % Setting Viewing Parameters
plot.setDataset(0, this.dataset1);                      % Add dataset to plot, Index in Java starts with zero
plot.setDataset(1, this.dataset2);                      % Add 2. dataset to plot


dialFrame = StandardDialFrame();                        % Create Frame for special Plot
dialFrame.setBackgroundPaint(Color.lightGray);          % Set ForgroundColor of the frame
dialFrame.setForegroundPaint(Color.darkGray);           % Set BackgroundColor of the frame
plot.setDialFrame(dialFrame);                           % Set the Frame to the plot


gp = GradientPaint(Point(), ...  % Create Gradient-Color for DialBackground
           Color(255/255, 255/255, 255/255), Point(), ... 
           Color(100/255, 100/255, 100/255));  

db = DialBackground(gp);    % Set the Color to DialBackground
db.setGradientPaintTransformer(StandardGradientPaintTransformer(GradientPaintTransformType.VERTICAL));      % Set the GradiensPainTransformer to the DialBackground
plot.setBackground(db);             % add the DialBackGround to the plot


annotation1 = DialTextAnnotation(java.lang.String('Speed: km/hr'));   % Create Annotation for Dialplot 
annotation1.setFont(java.awt.Font('Dialog', 1, 14));                  % Set the Font
annotation1.setRadius(0.7);                                           % Set the Radius for Annotation
plot.addLayer(annotation1);                                           % add annotation to plot


dvi = DialValueIndicator(0);                                         % Create DialValueIndicator for showing actual value of dataset1
dvi.setFont(java.awt.Font('Dialog', 0, 10));                         % Set the Font
dvi.setOutlinePaint(Color.darkGray);                                 % Set the OutlineColor
dvi.setRadius(0.60);                                                 % Set the Radius of DVI
dvi.setAngle(-103.0);                                                % set angle of DVI orientation
plot.addLayer(dvi);                                                  % add DVI to plot

dvi2 = DialValueIndicator(1);                                        % Create DialValueIndicator for showing actual value of dataset2
dvi2.setFont(java.awt.Font('Dialog', 0, 10));                        % Set the Font
dvi2.setOutlinePaint(Color.red);                                     % Set the OutlineColor
dvi2.setRadius(0.60);                                                % Set the Radius of DVI
dvi2.setAngle(-77.0);                                                % set angle of DVI orientation
plot.addLayer(dvi2);                                                 % add DVI to plot

scale = StandardDialScale;       % Create the Scale for dataset 1
 % Setting Scale Parameter
 scaleLower = 0;
 scaleUpper = 200; 
scale.setLowerBound( scaleLower );
scale.setUpperBound( scaleUpper );
scale.setStartAngle(-120);
scale.setExtent(-300);
scale.setTickRadius(0.88);
scale.setTickLabelOffset(0.15);
scale.setTickLabelFont(java.awt.Font('Dialog', 0, 14)); 
plot.addScale(0, scale);    % Add scale to plot


scale2 = StandardDialScale;  % Create the Scale for dataset 2
 % Setting Scale Parameter
scale2.setLowerBound(0);
scale2.setUpperBound(100);
scale2.setStartAngle(-120);
scale2.setExtent(-300);
scale2.setTickRadius(0.50);
scale2.setTickLabelOffset(0.15);
scale2.setTickLabelFont(java.awt.Font('Dialog', 0, 10)); 
scale2.setMajorTickPaint(Color.red);
%plot.addScale(1, scale2);     % Add scale to plot
%plot.mapDatasetToScale(1, 1); % map Data o scale


% Create range objects for highlighting ranges
range = StandardDialRange( scaleUpper*0.7, scaleUpper, Color.red);      % "Failure"
range.setInnerRadius(0.53);
range.setOuterRadius(0.56);
plot.addLayer(range);

range2 = StandardDialRange(scaleUpper*0.5, scaleUpper*0.7, Color.orange);  % "Attention"
range2.setInnerRadius(0.53);
range2.setOuterRadius(0.56);
plot.addLayer(range2);

range3 =   StandardDialRange(scaleLower, scaleUpper*0.5, Color.green);% "Easy"
range3.setInnerRadius(0.53);
range3.setOuterRadius(0.56);
plot.addLayer(range3);


needle2 = javaObjectEDT('org.jfree.chart.plot.dial.DialPointer$Pin',1);  %  create needle of pin-type for inner scale and dataset2 
needle2.setRadius(0.55);   % set radius of the inner pin
plot.addLayer(needle2);
 
needle = javaObjectEDT('org.jfree.chart.plot.dial.DialPointer$Pointer',0);   %  create needle of pointer-type for outer scale and dataset1 
plot.addLayer(needle);
 
cap = DialCap();      % Create a Dailcap
cap.setRadius(0.10);  % set the radius of dialcap
plot.setCap(cap);     % add dialcap to plot

plot_handle=plot;     % return this
 
%% Create Chart Area with Panel
chart1 = JFreeChart(plot);
chart1.setTitle('Velameter v0.1');
cp1 =  ChartPanel(chart1);

if ~exist('fh','var')  % New figure
  fh = figure('Units','normalized','position',[0.1,0.1,  0.2,  0.4]);
else
  sfigure(fh); 
end

% ChartPanel with JControl
jp = jcontrol(fh, cp1,'Position',[0.01 0.07 0.98 0.88]);


        
            

% %% Slider Callback for Outer-Needle
% function sh_callback_1(varargin)
% hObject = varargin{1,1};  
% % disp(['Slider moved to ' num2str(get(hObject,'Value'))]);   % diplay stuff in Matlab Command Window
% 
% % Get Handle from java plot object
% plot_cell = get(hObject,'Userdata' );
% plot_h    = plot_cell{1,1};    % handle of plot_object
% 
% % Update DialPlot
% plot_h.setDataset(org.jfree.data.general.DefaultValueDataset(get(hObject,'Value')));   %  change value of dataset1
% 
% 
% %% Slider Callback for Inner-Needle
% function sh_callback_2(varargin)
% hObject = varargin{1,1};  
% %disp(['Slider moved to ' num2str(get(hObject,'Value'))]);   % diplay stuff in Matlab Command Window
% 
% % Get Handle from java plot object
% plot_cell = get(hObject,'Userdata' );
% plot_h    = plot_cell{1,1};    % handle of plot_object
% 
% % Update DialPlot
% plot_h.setDataset(1, org.jfree.data.general.DefaultValueDataset(get(hObject,'Value')));   %  change value of dataset2
% 




% Not used, but if want to add sliders:            
            
% Matlab-Slider for dataset1
% sh = uicontrol(fh,'Style','slider',...
%                 'Max',60,'Min',-40,'Value',10,...
%                 'SliderStep',[0.01 0.01],...
%                 'Units','normalized', ...
%                 'Position',[0.01 0.01 0.90/2  0.03], ...
%                 'UserData', {plot}, ...                       % save the handle of the plot-object to Userdata to change values
%                 'Callback',@sh_callback_1);
%             
% 
%             
% % Matlab-Slider for dataset2
% sh2 = uicontrol(fh,'Style','slider',...
%                 'Max',100,'Min',0,'Value',70,...
%                 'SliderStep',[0.01 0.01],...
%                 'Units','normalized', ...
%                 'Position',[0.95/2 0.01 0.98/2 0.03], ...
%                 'UserData', {plot}, ...                       % save the handle of the plot-object to Userdata to change values
%                 'Callback',@sh_callback_2);
%             
%             
% % Matlab-Text for Sliders
% uicontrol(fh,'Style','Text',...
%                 'Units','normalized', ...
%                 'Position',[0.01 0.04 0.90/2  0.03], ...                
%                 'String', 'Outer Needle:');
%             
% uicontrol(fh,'Style','Text',...
%                 'Units','normalized', ...
%                 'Position',[0.95/2 0.04 0.98/2  0.03], ...                
%                 'String', 'Inner Needle:');
% 
    










