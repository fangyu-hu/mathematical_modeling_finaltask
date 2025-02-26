%% 源数据  
clear;clc
p24=[1,1,1,2,5,5,5,8,8,12,14,16,25,29,37,40,45,47,49,59,68,78,90,102]';
global true_i
true_i = p24;
plot(1:24,true_i,'r*')
legend('I')
hold on

%% 求参数并拟合  
lb = [0.02];  ub = [0.6];
x0 = [0.1]
[x, fval] = fmincon(@si_Obj_fun,x0,[],[],[],[],lb,ub)

T=24;
beta = 0.2072;
[~,p]=ode45(@(t,p) si_fun(t,p,beta), [1:1:T],true_i(1));
plot(1:T,p,'r-')

function f = si_Obj_fun(x)
    global true_i;
    beta = x;
    [~,p]=ode45(@(t,p) si_fun(t,p,beta), [1:1:24],true_i(1));
    f = sum((true_i -  p).^2);
end

function dp = si_fun(t,p,beta)
    N =  1400000000;
    I =  p;
    dp = beta * I * (1 - I/N);
end