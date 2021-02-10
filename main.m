clear
clc

resolution = 1500;
plot = true;
save = true;

%problems = {Simple2d(), HolderTable(), Griewank(), Rastrigin(), Schwefel(), Easom(), Sphere(), Ackley(), Shubert(), Langermann(), ModifiedLangermann()};
%problems = {CEC2005(1),CEC2005(2),CEC2005(3),CEC2005(4),CEC2005(5),CEC2005(6),CEC2005(7),CEC2005(8),CEC2005(9),CEC2005(10),CEC2005(11),CEC2005(12),CEC2005(13),CEC2005(14),CEC2005(15),CEC2005(16),CEC2005(17),CEC2005(18),CEC2005(19),CEC2005(20),CEC2005(21),CEC2005(22),CEC2005(23),CEC2005(24)};
problems = {SplitDropWave1(),SplitDropWave2()};
for k=1:length(problems)
    problem=problems{k};
    plot_problem(problem, resolution, plot, save);
end