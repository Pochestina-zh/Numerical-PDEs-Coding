%% 分数步长法（ADI,预校，LOD）解二维抛物方程

clear;
clc;
%网格参数
m=40; n=40;         %x(y),t方向的网格数
h=1/m; tau=1/n;     %h为x,y空间步长，tau为时间步长


%% 隐格式计算
tic;

toc;
%% 画图

