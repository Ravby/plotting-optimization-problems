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

% Define Viridis colormap
viridis = [
    0.267004, 0.004874, 0.329415;
    0.282656, 0.099695, 0.422173;
    0.296143, 0.184948, 0.512066;
    0.308428, 0.265254, 0.596338;
    0.319164, 0.341404, 0.674157;
    0.328547, 0.413873, 0.745587;
    0.336496, 0.483614, 0.809956;
    0.343121, 0.550985, 0.865707;
    0.348458, 0.616354, 0.911231;
    0.352538, 0.680057, 0.944781;
    0.355395, 0.742433, 0.964028;
    0.357068, 0.803804, 0.966707;
    0.357588, 0.864492, 0.950075;
    0.356995, 0.924806, 0.911536;
    0.355325, 0.985053, 0.849073;
    0.352615, 1.000000, 0.761464;
    0.348898, 1.000000, 0.647987;
    0.344213, 1.000000, 0.509994;
    0.338594, 1.000000, 0.351695;
    0.332071, 1.000000, 0.179993
];

% Assuming X1, X2, and F are defined properly
% Calculate the size of F
[m, n] = size(F);

% Generate C to match the size of F
C = repmat(viridis, m, n);

meshc(X1, X2, F);
%colormap(viridis); % Set colormap to Viridis
colorbar;

view(-135, 30);
set(gca,'FontSize',12);
set(gcf, 'Position',  [300, 50, 1200, 900])
xlabel('x_1','FontName','Helvetica','FontSize',20);
set(get(gca,'xlabel'),'rotation',-20,'VerticalAlignment','middle');
ylabel('x_2','FontName','Helvetica','FontSize',20);
set(get(gca,'ylabel'),'rotation',20,'VerticalAlignment','middle');
set(gca, 'YDir','reverse')
zlabel('f(X)','FontName','Helvetica','FontSize',20);
title(problem.problem_name,'FontName','Helvetica','FontSize',24,'FontWeight','bold');
if(save)
    exportgraphics(problem3Dfig, 'images\'+problem.problem_name+'3D.png', 'ContentType', 'image', 'BackgroundColor', 'white');
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
title(problem.problem_name,'FontName','Helvetica','FontSize',24,'FontWeight','bold');
if(save)
    exportgraphics(problemHeatmapFig, 'images\'+problem.problem_name+'heatmap.png', 'ContentType', 'image', 'BackgroundColor', 'white');
end
