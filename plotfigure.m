clear all;close all;clc

%%
duration = [20, 40, 80, 160, 250, 350]; % ms
contrast = [1, 4, 8, 16, 32, 52, 100];

nDur = length(duration);
nConst = length(contrast);
%%
motionEnergy=zeros(nDur, nConst, 4); % 2 directions
for i=1:nDur
    for j=1:nConst
        [motionEnergy(i, j, :), ~]=computeEnergy2(duration(i), 5.66, contrast(j));
    end
end

%%
close all;
cpsfigure(1,2);
subplot(1,2,1);
myplot(duration, squeeze(motionEnergy(:, 1, 1))', [], '-ro');
myplot(duration, squeeze(motionEnergy(:, 1, 2))', [], '-bo');
myplot(duration, squeeze(motionEnergy(:, 1, 3))', [], '-ko');
title('Low contrast');
xlabel('Duration (ms)');
ylabel('Energy');

subplot(1,2,2);
myplot(duration, squeeze(motionEnergy(:, end, 1))', [], '-ro');
myplot(duration, squeeze(motionEnergy(:, end, 2))', [], '-bo');
myplot(duration, squeeze(motionEnergy(:, end, 3))', [], '-ko');
title('High contrast');
xlabel('Duration (ms)');
ylabel('Energy');

%%
close all
x1 = squeeze(motionEnergy(:, 1, 1));
x2 = squeeze(motionEnergy(:, 1, 2));
cpsfigure(1,2);
subplot(1,2,1);
myplot(duration, (x1./(x1+x2))', [], '-ro');
myplot(duration, (x2./(x1+x2))', [], '-bo');
title('Low contrast');
xlabel('Duration (ms)');
ylabel('Energy');
set(gca, 'XScale', 'log');


x3 = squeeze(motionEnergy(:, end, 1));
x4 = squeeze(motionEnergy(:, end, 2));

subplot(1,2,2);
myplot(duration, (x3./(x3+x4))', [], '-ro');
myplot(duration, (x4./(x3+x4))', [], '-bo');
title('High contrast');
xlabel('Duration (ms)');
ylabel('Energy');
set(gca, 'XScale', 'log');