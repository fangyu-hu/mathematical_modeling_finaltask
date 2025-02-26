%% 原始数据  
clear;clc
load data1.mat  
global true_s true_i true_r
true_s = p62(:,1);
true_i = p62(:,2);
true_r = p62(:,3);
plot(1:62,true_i,'r*')
legend('I')
hold on
%% 求出参数值并拟合  
lb = [0.002 0.03];  ub = [5.0 0.2];
% fmincom函数  
x0 = [0.1 0.05]
[x, fval] = fmincon(@sir_Obj_fun2,x0,[],[],[],[],lb,ub)

T = 86;
beta = 0.6702;
gamma = 0.0516;
[~,p]=ode45(@(t,p) sir_fun2(t,p,beta,gamma), [1:1:T],[true_s(1) true_i(1) true_r(1)]);
plot(1:T,p(:,2),'r-')
hold on
%% 接种疫苗  
p(86,3)=p(86,3)+p(86,1)*0.75;
p(86,1)=p(86,1)*0.25;
beta = 0.6702;
gamma = 0.0516;
T2 = 146;
[~,p2]=ode45(@(t,p) sir_fun2(t,p,beta,gamma), [86:1:T2],[p(86,1) p(86,2) p(86,3)]);
plot(86:T2,p2(:,2),'b-')

function f = sir_Obj_fun2(x)
    global true_s true_i true_r;
    beta = x(1);
    gamma = x(2);
    [~,p]=ode45(@(t,p) sir_fun2(t,p,beta,gamma), [1:1:62],[true_s(1) true_i(1) true_r(1)]);
    f = sum((true_s - p(:,1)).^2  + (true_i -  p(:,2)).^2  + (true_r - p(:,3)).^2);
end
function dp = sir_fun2(t,p,beta,gamma)
    N =  1324000000;
    dp = zeros(3,1);
    S = p(1);
    I =  p(2);
    dp(1) = -beta * S * I / N;
    dp(2) = beta * S * I / N - gamma * I;
    dp(3) = gamma * I;
end