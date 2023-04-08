%function [x,y]=back_euler_lam(xRange,yInitial,numSteps)
% NOTES
% y'=LAMBDA*(-y+sin(x)) y(0)=0 x0=0, xf=2*pi
% x(i+1)=x(i)+h
% y(i+1)=y(1)+h*LAMBDA*(-y(i+1)+sin(x(i+1)))
% [x,y]=back_euler_lam(xRange,yInitial,numSteps)
% In: xRange=[0,2*pi]
% In: yInitial (y(0)=0)
% In numSteps = number of steps 
% Out: [x y] result of bckwd Euler 
% Victoria Bell March 31st 2023 
clc
clear all
%Define parameter
LAMBDA=1000;
%Independent variable  
x0=0;
xf=2*pi;
xRange=[x0 xf];
numSteps=40;
h=xRange/numSteps;
x=x0:h:xf;
%Initial Condition for y 
y(1)=0;

options=optimset('TolX',1e-06);
i=0;
for j=x0:h:xf-h
    xbe(i+1)=x0+x*h;
    x=xbe(i);
    y(i+1)=y(i)+h*LAMBDA*sin(x)/(1+h*LAMBDA);
end

figure(1)
hold on 
plot(x,y,'bo')
%plot(x,y_exact,'r-')
xlabel('x')
ylabel('y')
title('Backward Euler Method')