clear
clc

problem_num = input('Enter problem number: ');
problem_name = "";

R = 1500; % resolution

switch problem_num
    case 1
        problem_name = "Rastrigin";
        disp(problem_name)
        LB = [-5.12 -5.12];  
        UB = [5.12 5.12];
        problem_function = @rastrigin;
    case 2
        problem_name = "Schwefel";
        disp(problem_name)
        LB = [-500 -500];  
        UB = [500 500];
        problem_function = @schwefel;
    case 3
        problem_name = "Easom";
        disp(problem_name)
        LB = [-100 -100];  
        UB = [100 100];
        problem_function = @easom;
    case 4
        problem_name = "Sphere";
        disp(problem_name)
        LB = [-100 -100];  
        UB = [100 100];
        problem_function = @sphere;
    case 5
        problem_name = "Ackley";
        disp(problem_name)
        LB = [-32.768 -32.768];  
        UB = [32.768 32.768];
        problem_function = @ackley;
    case 6
        problem_name = "Shubert";
        disp(problem_name)
        LB = [-10 -10];  
        UB = [10 10];
        problem_function = @shubert;
    case 7
        problem_name = "Langerman";
        disp(problem_name)
        LB = [0 0];  
        UB = [10 10];
        problem_function = @langerman;
    otherwise
        disp(strcat("No problem with number ", num2str(problem_num)))
end

X1=LB(1):(UB(1)-LB(1))/R:UB(1);
X2=LB(2):(UB(2)-LB(2))/R:UB(2);

f = zeros(size(X1));
F = zeros(size(X1));

for j=1:length(X1)
    for i=1:length(X2)
        f(i)=problem_function([X1(j),X2(i)]);
    end
    F(j,:)=f;
end

figure(1)
meshc(X1, X2, F);colorbar;
set(gca,'FontSize',12);
set(gcf, 'Position',  [300, 300, 800, 600])
xlabel('x_2','FontName','Times','FontSize',20,'FontAngle','italic');
set(get(gca,'xlabel'),'rotation',25,'VerticalAlignment','bottom');
ylabel('x_1','FontName','Times','FontSize',20,'FontAngle','italic');
set(get(gca,'ylabel'),'rotation',-25,'VerticalAlignment','bottom');
zlabel('f(X)','FontName','Times','FontSize',20,'FontAngle','italic');
title(problem_name,'FontName','Times','FontSize',24,'FontWeight','bold');

figure(2)
mesh(X1,X2,F);view(0,90);colorbar;set(gca,'FontSize',12);
xlabel('x_2','FontName','Times','FontSize',20,'FontAngle','italic');
ylabel('x_1','FontName','Times','FontSize',20,'FontAngle','italic');
zlabel('f(X)','FontName','Times','FontSize',20,'FontAngle','italic');
title(problem_name,'FontName','Times','FontSize',24,'FontWeight','bold');

