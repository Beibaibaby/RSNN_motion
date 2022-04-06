%computing motion energy of motion stimuli
%By Ruyuan Zhang

%function Constant_stimuli
function [TFenergy, range]=computeEnergy2(duration, speed, contrast, wantfig)
%here, you just need to input duration and speed. you can compute motion
%energy for this stimuli. TFenergy is vector for motion energy and range is
%the vector for "x" axis
% <duration>: milliseconds
% <speed>: deg/sec, default:5.66
% <contrast>: 0~100

clc;

warning('off','MATLAB:dispatcher:InexactMatch');

if notDefined('duration')
    duration = 50;
end
if notDefined('speed')
    speed = 5.66;
end
if notDefined('contrast')
    contrast = 100;
end
if notDefined('wantfig')
    wantfig = 0;
end




%----------------------------------Variables------------------------
frame_rate = 120; % Screen frame rate (hz)
scale_factor = 2; % most important parameter - how many arcmin is one pixel?
background = 126; % background intensity, in gray scale units

%-------------------------------Standard Gebor or Grating----------
angle = 0; % angle of orientation of the Gebor
cycles = 2; % cycles per degree
stimulus_radius = 8*60; % in arcmin
 
if speed >= 0
    direction = 1; %1, left;-1,right
else
    direction = -1;
end

% ********* KEEP DERIVED VARIABLES IN THIS SECTION
%housekeeping stuff
SF = cycles/(stimulus_radius/60);
%SF                              = 1;         %cycles/deg
TF = speed*SF;   %cycles/sec
duration= duration/(1000/frame_rate);   % convert durations to frames
TFstep = (2*pi*TF)/frame_rate;
stimulus_radius  = stimulus_radius/scale_factor;
f=(SF*scale_factor/60)*2*pi;
a=cos(angle)*f; b=sin(angle)*f;
amplitude = background*contrast/100;
yy = round(stimulus_radius);

%make the spatial envelope
[x,y]=meshgrid(-stimulus_radius:stimulus_radius,-yy:yy);
bps = (stimulus_radius)*2+1;
circle=((stimulus_radius)^2-(x.^2+y.^2));
circle = circle>0;
bps= size(circle,1);
R = (sqrt(x.^2 + y.^2) + eps).*circle;
R = R/max(max(R));
cos2D = (cos(R*pi)+1)/2;
circle = (cos2D.*circle);


%---------------------------------create and record motion movie------------------
% calculate the temporal envelope
movie_rect=[0 0 bps bps];

mv_length = duration;
mv_length = round(mv_length/2)*2+1; %ensure that mv_length had odd number of frames
clear xx;
xx=1:mv_length;    	xx=(xx-mean(xx)) ;
time_gauss= exp(- ((xx)/(sqrt(2)*duration)).^2);  % compute the temporal gaussian


Movie=zeros(bps,bps,mv_length);     %3D matrix which record the movie
%make the  movie
motion_step(1) = rand*2*pi;
for i=2:mv_length
    motion_step(i) = motion_step(i-1)-TFstep;
end
for i = 1:mv_length
    %moving_gratting=round(((sin(a*x+b*y+motion_step(i)).*circle*amplitude*time_gauss(i))+background));
    %moving_gratting=round(((sin(a*x+b*y+motion_step(i))*amplitude*time_gauss(i))+background));
    moving_gratting=round(((sin(a*x+b*y+motion_step(i)).*circle*amplitude*1)+background));
    %moving_gratting=round(((sin(a*x+b*y+motion_step(i))*amplitude*1)+background));
    Movie(:,:,i)=moving_gratting;
    %MovieTexture(i)= Screen('MakeTexture',window,moving_gratting);
end 


% convert 
range = 1:(mv_length-1)/2;
range = range/(mv_length/frame_rate);


%--------------------calculate plot the spatial-temporal spectrum--------------------
%xaxis=bps+round((mv_length-1)*SFstep*direction);
slice=(bps+1)/2; 
stimImg1=Movie(slice,:,:); % left/right motion
stimImg2=Movie(:,slice,:); % up/down motion
stimImg1=squeeze(stimImg1)/255;
stimImg2=squeeze(stimImg2)/255;

if wantfig
    figure(1);
    subplot(2,1,1)
    imshow(stimImg1');
    ylabel('time step');
    xlabel('leftrightwidth');

    subplot(2,1,2)
    imshow(stimImg2');
    ylabel('time step');
    xlabel('updownheight');
end
%-------------------do the fourier transform for the spatial---------------
%% true motion axis
F1=fftshift(fft2(stimImg1));
F1=abs(F1)'; % only extract signal amplitude

if wantfig
    figure(2);
    imshow(F1);
    ylabel('temporal frequency');
    xlabel('leftright spatial frequency');
end

TFenergy1 = F1(end-length(range)+1:end, (bps+1)/2+1:end);
TFenergy2 = F1(1:length(range), (bps+1)/2+1:end);

TFenergy1 = sum(TFenergy1(:));
TFenergy2 = sum(TFenergy2(:));

%% othorganal motion axis
F2=fftshift(fft2(stimImg2));
F2=abs(F2)'; % only extract signal amplitude

if wantfig
    figure(3);
    imshow(F2);
    ylabel('temporal frequency');
    xlabel('updown spatial frequency');
end

TFenergy3 = F2(end-length(range)+1:end, (bps+1)/2+1:end);
TFenergy4 = F2(1:length(range), (bps+1)/2+1:end);

TFenergy3 = sum(TFenergy3(:));
TFenergy4 = sum(TFenergy4(:));


%% concatenate
TFenergy = [TFenergy1 TFenergy2 TFenergy3 TFenergy4]';

% figure(3);
% plot(range, TFenergy1, 'r-o','LineWidth',1.5); hold on;
% plot(range, TFenergy2, 'b-o','LineWidth',1.5);
% % % ylim([0 maxPower+200]);
% xlabel('temporal frequency');
% ylabel('amplitude');
%% calculate left/right motion energy

