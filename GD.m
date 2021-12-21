close all;
clear all;
clc;
rng(55)
set(0,'defaultfigurecolor','w')
n=50;m=200;
barx = randn(n,1);
A = randn(m,n);
epsilon = randn(m,1)/20;
b= A*barx+epsilon;
num_iter = 10;
initial_guess = randn(n,1);
x = initial_guess;
h_er=[];
h_x=[];
h_delta_xk_with_xbar=[];
h_gamma=[];
gamma = 0.001;
while 1
    gradient = 2*A'*(A*x-b);
    er = norm(A*x-b)^2;
    h_er = [h_er er];
    [gamma, x_test, ~] = goldenSection(0,0.5,@(gamma) (norm ( A*(x - gamma * gradient)-b )^2 )); 
%     gamma=0.0002;
    h_gamma=[h_gamma gamma];
    if  (norm(gamma*gradient) < 1e-8)
        break
    else
        x = x - gamma * gradient;
        h_delta_xk_with_xbar = [h_delta_xk_with_xbar norm(x-barx)];
    end
end
x_pseudo = (A'*A)^(-1)*A'*b;

figure,
set(gcf,'position',[0,300,800,400])
subplot(2,1,1)
semilogy(h_er,'ok-', 'linewidth',    0.1, 'markerfacecolor', [36, 169, 225]/255,"MarkerSize",5)
title("Objective~Function~$ ||Ax-b||^2_2$",'Interpreter','LaTex','FontSize',12)
xlabel(" iteration times",'Interpreter','LaTex')
hl=legend(" $||Ax-b||^2_2$",'Interpreter','LaTex','FontSize',12);
set(hl, 'Box', 'off')

subplot(2,1,2)
semilogy(h_delta_xk_with_xbar,'ok-','linewidth', 1.1, 'markerfacecolor', [29, 191, 151]/255,"MarkerSize",5)
title('The difference from groud truth $||x^{(k)} - \bar{x}||$','Interpreter','LaTex','FontSize',12);
xlabel(" iteration times",'Interpreter','LaTex')
hl=legend("$||x^{(k)} - \bar{x}||$",'Interpreter','LaTex','FontSize',12);
set(hl, 'Box', 'off')


figure,
set(gcf,'position',[300,300,500,400])
stem(barx,'-*','Color',[142 207 201]/255,'LineWidth',4)
hold on
stem(x,'-x','Color',[255 190 122]/255,'LineWidth',2)
hold on
stem(x_pseudo,'-o','Color',[84 134 135]/255,'LineWidth',1)
legend("ground truth $\bar{x}$","$x$ esitimated by GD algorithm","$x$ using pseudo inverse",'Interpreter','LaTex')
title("Result of signal recovery",'Interpreter','LaTex')

figure,
set(gcf,'position',[500,300,800,400])
plot(h_gamma,'ok-','linewidth', 1.1, 'markerfacecolor', [29, 191, 151]/255,"MarkerSize",5)
legend("每次迭代采用的步长")
title("Step size $\gamma$ determined by Golden Section exact line search method",'Interpreter','LaTex','FontSize',12)
