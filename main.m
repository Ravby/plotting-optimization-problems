clear
clc

resolution = 1500;
plot = false;
save = true;

problems = {HolderTable(), Griewank(), Rastrigin(), Schwefel(), Easom(), Sphere(), Ackley(), Shubert(), Langermann(), ModifiedLangermann()};

for k=1:length(problems)
    problem=problems{k};
    plot_problem(problem, resolution, plot, save);
end