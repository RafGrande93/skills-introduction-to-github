function [t,u,err]=AB2(fun,t0,T,y0,N,ye)
%  Adams-Bashfort Methods with p=k=2 
% for Cauchy problem y'=f(t,y), y(t0)=y0
h=(T-t0)/N;
t=t0:h:T;
u(1)=y0; f0=feval(fun,t0,y0); 
u1=u(1)+h*f0; % explicit step (Euler's methods)
u(2)=u1; 
for n=1:N-1
    fn=feval(fun,t(n),u(n));
    fn1=feval(fun,t(n+1),u(n+1));  %implicit part, with two steps
    u(n+2)=u(n+1)+(h/2)*(3*fn1-fn);
end
figure(2), subplot(121), plot(t,u,'r.-'),
title('Adams-Bashfort, k=2')
if isempty(ye)
    err=[];
else
    err=abs(ye(t)-u); % Check error
    figure(2),
    tt=linspace(t0,T,200);
    subplot(121), hold on, plot(tt,ye(tt),'b-');
    subplot(122), hold on, semilogy(t,err,'r.-');
end
