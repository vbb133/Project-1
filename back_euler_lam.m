function [x,y]=back_euler_lam(xRange,yInitial,numSteps)
% [x,y]=back_euler_lam(xRange,yInitial,numSteps)
% In: xRange=[0,2*pi]
% In: yInitial (y(0)=0)
% In numSteps = number of steps 
%Out: [x y] result of bckwd Euler 
% Victoria Bell March 31st 2023 
clc
%Define parameter
LAMBDA=10000;
%X range
x0=xRange(1);
x_max=xRange(2);
%step size
S=sign(x_max-x0);
h=S*xRange(2)/numSteps;
%Initial conditions
x=x0;
y=yInitial(:);
xbe=x;
ybe=y.';
xbe=zeros(numSteps,1);
xbe(1) = x;
ybe=zeros(numSteps,size(y,1));
ybe(1,:) = y.';
j=1;

%Implicit Euler
while j<numSteps
  if S*(x+h-x_max)>0
    h=x_max-x;
    xbe(j+1)=x_max;
  else
    xbe(j+1)=x0 +j*h;
  end
  j=j+1;
    xbe(j)=x0+j*h;
    x=xbe(j);
    y=(y+h*LAMBDA*sin(x))/(1+h*LAMBDA);
    ybe(j,:)=y.';
end
x=xbe
y=ybe
abs(y(end)-stiff50_solution(x(end)))
%Input for command:  back_euler_lam([0 2*pi],0,40)

%Notes 
%rhs ODE y'=LAMBDA*(-y+sin(x))
%f=@(x,y) LAMBDA*(-y+sin(x));

%Analytical solution 
%F= @(x,y)LAMBDA/(1+LAMBDA^2)*exp(-LAMBDA*x)+(LAMBDA^2*sin(x)-LAMBDA*cos(x)/(1+LAMBDA^2);

