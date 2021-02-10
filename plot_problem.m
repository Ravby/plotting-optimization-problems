function [] = plot_problem(problem, resolution, plot, save)

X1=problem.LB(1):(problem.UB(1)-problem.LB(1))/resolution:problem.UB(1);
X2=problem.LB(2):(problem.UB(2)-problem.LB(2))/resolution:problem.UB(2);

f = zeros(size(X1));
F = zeros(size(X1));

for j=1:length(X1)
    for i=1:length(X2)
        f(i)=problem.evaluate([X1(j),X2(i)]);
    end
    F(j,:)=f;
end

problem3Dfig = figure(1);

if(plot)
    set(problem3Dfig, 'Visible', 'on');
else
    set(problem3Dfig, 'Visible', 'off');
end

meshc(X1, X2, F);colorbar;
view(-135, 30);
set(gca,'FontSize',12);
set(gcf, 'Position',  [300, 50, 1200, 900])
xlabel('x_1','FontName','Helvetica','FontSize',20);
set(get(gca,'xlabel'),'rotation',-20,'VerticalAlignment','middle');
ylabel('x_2','FontName','Helvetica','FontSize',20);
set(get(gca,'ylabel'),'rotation',20,'VerticalAlignment','middle');
set(gca, 'YDir','reverse')
zlabel('f(X)','FontName','Helvetica','FontSize',20);
%title(problem.problem_name,'FontName','Helvetica','FontSize',24,'FontWeight','bold');
if(save)
    saveas(problem3Dfig, 'images\'+problem.problem_name+'3D.png','png');
end

problemHeatmapFig = figure(2);
if(plot)
    set(problemHeatmapFig, 'Visible', 'on');
else
    set(problemHeatmapFig, 'Visible', 'off');
end
mesh(X1, X2, F);view(90, -90);colorbar;set(gca,'FontSize',12);
set(gcf, 'Position',  [300, 50, 1200, 900])
xlabel('x_2','FontName','Helvetica','FontSize',20);
set(get(gca,'xlabel'),'rotation',-90,'VerticalAlignment','middle');
ylabel('x_1','FontName','Helvetica','FontSize',20);
set(get(gca,'ylabel'),'VerticalAlignment','top');
zlabel('f(X)','FontName','Helvetica','FontSize',20);
%title(problem.problem_name,'FontName','Helvetica','FontSize',24,'FontWeight','bold');
if(save)
    saveas(problemHeatmapFig, 'images\'+problem.problem_name+'heatmap.png','png');
end
