close all;
clear
clc;
set(0,'defaultfigurecolor','w')
rng(55)
n_list=[300, 500, 1000, 2000];m=200;

figure,
for k = 1:4
    n=n_list(k);
    initial_guess_list = [randn(n,1) zeros(n,1) ones(n,1)];
    for j = 1:3
        initial_guess=initial_guess_list(:,j);
        barx = randn(n,1);
        A = randn(m,n);
        epsilon = randn(m,1)/20;
        b= A*barx+epsilon;
        x = initial_guess;
        h_x=[];
        h_delta_xk_with_xbar=[];
        h_gamma=[];
        for i = 1:10000
            index = ceil(rand(1,1)*m); %生成1个1到m的随机整数
            er = norm(A*x-b)^2;

            gamma = 0.0001;
            gradient =  (A(index, :) * x - b(index)) * A(index, :)';
            x = x - gamma * gradient;
            h_delta_xk_with_xbar = [h_delta_xk_with_xbar norm(x-barx)/n];
        end
        x_pseudo =A'*(A*A')^(-1)*b;
        initial_text = ["x~N(0,1)","all 0","all 1"];
        txt = " with initial point is \textbf{"+initial_text(j)+"}, and n=\textbf{"+n+"}";

        subplot(4,3,(k-1)*3+j)
        set(gcf,'position',[300,300,500,400])
        stem(barx,'-*','Color',[142 207 201]/255,'LineWidth',4)
        hold on
        stem(x,'-x','Color',[255 190 122]/255,'LineWidth',2)
        hold on
        stem(x_pseudo,'-o','Color',[84 134 135]/255,'LineWidth',1)
        legend("ground truth $\bar{x}$","$x$ esitimated by GD algorithm","$x$ using pseudo inverse",'Interpreter','LaTex')
        title("Result of signal recovery"+txt,'Interpreter','LaTex')
        hold on
    end 
end

