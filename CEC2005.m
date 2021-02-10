classdef CEC2005  < Problem
    properties
        LB;
        UB;
        funNum;
        fhd;
        data;
    end
    methods
        function obj = CEC2005(funNum)
            obj = obj@Problem("CEC 2005 F"+funNum);
            obj.funNum = funNum;
            
            addpath('./CEC2005/');
            
            if funNum==1 obj.LB = [-100 -100]; obj.UB = [100 100]; obj.fhd=str2func('sphere_func'); %[-100,100]
            elseif funNum==2 obj.LB = [-100 -100]; obj.UB = [100 100]; obj.fhd=str2func('schwefel_102');%[-10,10]
            elseif funNum==3 obj.LB = [-100 -100]; obj.UB = [100 100]; obj.fhd=str2func('high_cond_elliptic_rot_func');%[-100,100]
            elseif funNum==4 obj.LB = [-100 -100]; obj.UB = [100 100]; obj.fhd=str2func('schwefel_102_noise_func');%[-100,100]
            elseif funNum==5 obj.LB = [-100 -100]; obj.UB = [100 100]; obj.fhd=str2func('schwefel_206');%[-100,100]
            elseif funNum==6 obj.LB = [-100 -100]; obj.UB = [100 100]; obj.fhd=str2func('rosenbrock_func');%[-100,100]
            elseif funNum==7 obj.LB = [-600 -600]; obj.UB = [600 600]; obj.fhd=str2func('griewank_rot_func');%[-600,600]
            elseif funNum==8 obj.LB = [-32 -32]; obj.UB = [32 32]; obj.fhd=str2func('ackley_rot_func');%[-32,32]
            elseif funNum==9 obj.LB = [-5 -5]; obj.UB = [5 5]; obj.fhd=str2func('rastrigin_func');%[-5,5]
            elseif funNum==10 obj.LB = [-5 -5]; obj.UB = [5 5]; obj.fhd=str2func('rastrigin_rot_func');%[-5,5]
            elseif funNum==11 obj.LB = [-0.5 -0.5]; obj.UB = [0.5 0.5]; obj.fhd=str2func('weierstrass_rot');%[-0.5,0.5]
            elseif funNum==12 obj.LB = [-pi -pi]; obj.UB = [pi pi]; obj.fhd=str2func('schwefel_213');%[-pi,pi]
            elseif funNum==13 obj.LB = [-3 -3]; obj.UB = [1 1]; obj.fhd=str2func('EF8F2_func');%[-3,1]
            elseif funNum==14 obj.LB = [-100 -100]; obj.UB = [100 100]; obj.fhd=str2func('E_ScafferF6_func');%[-100,100]
            elseif funNum==15 obj.LB = [-5 -5]; obj.UB = [5 5]; obj.fhd=str2func('hybrid_func1');%[-5,5]
            elseif funNum==16 obj.LB = [-5 -5]; obj.UB = [5 5]; obj.fhd=str2func('hybrid_rot_func1');%[-5,5]
            elseif funNum==17 obj.LB = [-5 -5]; obj.UB = [5 5]; obj.fhd=str2func('hybrid_rot_func1_noise');%[-5,5]
            elseif funNum==18 obj.LB = [-5 -5]; obj.UB = [5 5]; obj.fhd=str2func('hybrid_rot_func2');%[-5,5]
            elseif funNum==19 obj.LB = [-5 -5]; obj.UB = [5 5]; obj.fhd=str2func('hybrid_rot_func2_narrow');%[-5,5]
            elseif funNum==20 obj.LB = [-5 -5]; obj.UB = [5 5]; obj.fhd=str2func('hybrid_rot_func2_onbound');%[-5,5]        
            elseif funNum==21 obj.LB = [-5 -5]; obj.UB = [5 5]; obj.fhd=str2func('hybrid_rot_func3');%[-5,5]
            elseif funNum==22 obj.LB = [-5 -5]; obj.UB = [5 5]; obj.fhd=str2func('hybrid_rot_func3_highcond');%[-5,5]     
            elseif funNum==23 obj.LB = [-5 -5]; obj.UB = [5 5]; obj.fhd=str2func('hybrid_rot_func3_noncont');%[-5,5]
            elseif funNum==24 obj.LB = [-5 -5]; obj.UB = [5 5]; obj.fhd=str2func('hybrid_rot_func4');%[-5,5]   
            elseif funNum==25 obj.LB = [-5 -5]; obj.UB = [5 5]; obj.fhd=str2func('hybrid_rot_func4');%[-5,5]  
            end
            
            
            obj.data = load('fbias_data');
        end
        function out = evaluate(obj, x)
            out=feval(obj.fhd,x) + obj.data.f_bias(obj.funNum);
        end
    end
end

% benchmark_func.m is the main function for 25 test functions, all minimize
% problems
% e.g. f=benchmark_func(x,func_num)
% x is the variable, f is the function value 
% func_num is the function num,

%       25 TEST FUCNTIONS
% 	    Unimodal Functions (5):
% 1.    Shifted Sphere Function 					                Bounds[-100,100]	f_bias=-450
% 2.	Shifted Schwefel's Problem 1.2				            	Bounds[-100,100]	f_bias=-450
% 3.	Shifted Rotated High Conditioned Elliptic Function			Bounds[-100,100]	f_bias=-450
% 4.	Shifted Schwefel's Problem 1.2 with Noise in Fitness 		Bounds[-100,100]	f_bias=-450
% 5.	Schwefel's  Problem 2.6 with Global Optimum on Bounds		Bounds[-100,100]	f_bias=-310
% 
% 	    Multimodal Functions (20):
% 	    Basic Functions (7):
% 6.	Shifted Rosenbrock's  Function					            Bounds[-100,100]	f_bias=390 
% 7.	Shifted Rotated Griewank's  Function without Bounds	        Initilization Range [0, 600]	f_bias=-180
% 8.	Shifted Rotated Ackley's  Function with Global Optimum on Bounds	Bounds[-32,32]	f_bias=-140
% 9.	Shifted Rastrigin's  Function 					            Bounds[-5,5]	    f_bias=-330
% 10.	Shifted Rotated Rastrigin's  Function 				        Bounds[-5,5]	    f_bias=-330
% 11.	Shifted Rotated Weierstrass Function 				        Bounds[-0.5,0.5]	f_bias=90
% 12.	Schwefel's  Problem 2.13					                Bounds[-100,100]	f_bias=-460 
% 	    Expanded Functions (2):
% 13.	Expanded Extended Griewank's  plus Rosenbrock's  Function (F8F2)	Bounds[-3,1]	f_bias=-130
% 14.	Expanded Rotated Extended Scaffe's  F6 				        Bounds[-100,100]	f_bias=-300
% 	    Hybrid Composition Functions (11):
% 15.	Hybrid Composition Function 1				                Bounds[-5,5]	    f_bias= 120 
% 16.	Rotated Hybrid Composition Function 1				        Bounds[-5,5]	    f_bias= 120
% 17.	Rotated Hybrid Composition Function 1 with Noise in Fitness		Bounds[-5,5]	f_bias= 120
% 18.	Rotated Hybrid Composition Function 2			        	Bounds[-5,5]	    f_bias=10 
% 19.	Rotated Hybrid Composition Function 2 with a Narrow Basin for the Global Optimum	Bounds[-5,5]]	f_bias=10 
% 20.	Rotated Hybrid Composition Function 2 with the Global Optimum on the Bounds		Bounds[-5,5]	f_bias=10
% 21.	Rotated Hybrid Composition Function 3						Bounds[-5,5]    	f_bias=360 
% 22.	Rotated Hybrid Composition Function 3 with High Condition Number Matrix		Bounds[-5,5]	f_bias=360
% 23.	Non-Continuous Rotated Hybrid Composition Function 3		Bounds[-5,5]    	f_bias=360 
% 24.	Rotated Hybrid Composition Function 4				        Bounds[-5,5]	    f_bias=260 
% 25.	Rotated Hybrid Composition Function 4 without Bounds	    Intilization Range[-2,5]	f_bias=260 
%
%J. J. Liang & P. N. Suganthan   2005.Feb 18


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Unimodal%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	1.Shifted Sphere Function 
function fit=sphere_func(x)
persistent fData
[ps,D]=size(x);
if isempty(fData)
    fData = load('sphere_func_data');
    if length(fData.o)>=D
         fData.o=fData.o(1:D);
    else
         fData.o=-100+200*rand(1,D);
    end
end
x=x-repmat(fData.o,ps,1);
fit=sum(x.^2,2);
end

% 	2.Shifted Schwefel's Problem 1.2
function f=schwefel_102(x)
persistent fData
[ps,D]=size(x);
if isempty(fData)
    fData = load('schwefel_102_data');
   if length(fData.o)>=D
         fData.o=fData.o(1:D);
    else
         fData.o=-100+200*rand(1,D);
    end
end

x=x-repmat(fData.o,ps,1);
f=0;
for i=1:D
    f=f+sum(x(:,1:i),2).^2;
end
end

% 	3.Shifted Rotated High Conditioned Elliptic Function
function fit=high_cond_elliptic_rot_func(x)
persistent fData matData
[ps,D]=size(x);
if isempty(fData)
    fData = load('high_cond_elliptic_rot_data');
    if length(fData.o)>=D
         fData.o=fData.o(1:D);
    else
         fData.o=-100+200*rand(1,D);
    end
    c=1;
    if D==2, matData = load('elliptic_M_D2');
    elseif D==10,matData = load('elliptic_M_D10');
    elseif D==30,matData = load('elliptic_M_D30');
    elseif D==50,matData = load('elliptic_M_D50');
    else 
        A=normrnd(0,1,D,D);[matData.M,r]=cGram_Schmidt(A);
    end
end
x=x-repmat(fData.o,ps,1);
x=x*matData.M;
a=1e+6;
fit=0;
for i=1:D
fit=fit+a.^((i-1)/(D-1)).*x(:,i).^2;
end
end

% 	4.Shifted Schwefel's Problem 1.2 with Noise in Fitness 
function f=schwefel_102_noise_func(x)
persistent fData
[ps,D]=size(x);
if isempty(fData)
    fData = load('schwefel_102_data');
    if length(fData.o)>=D
         fData.o=fData.o(1:D);
    else
         fData.o=-100+200*rand(1,D);
    end
end
x=x-repmat(fData.o,ps,1);
f=0;
for i=1:D
    f=f+sum(x(:,1:i),2).^2;
end
f=f.*(1+0.4.*abs(normrnd(0,1,ps,1)));
end

% 	5.Schwefel's Problem 2.6
function f=schwefel_206(x)%after Fletcher and Powell
persistent fData B
[ps,D]=size(x);
if isempty(fData)
    fData = load('schwefel_206_data');
    if length(fData.o)>=D
         fData.A=fData.A(1:D,1:D);fData.o=fData.o(1:D);
    else
         fData.o=-100+200*rand(1,D);
         fData.A=round(-100+2*100.*rand(D,D));
         while det(fData.A)==0
         fData.A=round(-100+2*100.*rand(D,D));
         end
    end
    fData.o(1:ceil(D/4))=-100;fData.o(max(floor(0.75*D),1):D)=100;
    B=fData.A*fData.o';
end
for i=1:ps
f(i,1)=max(abs(fData.A*(x(i,:)')-B));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Multimodal%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	6.Shifted Rosenbrock's Function
function f=rosenbrock_func(x)
persistent fData
[ps,D]=size(x);
if isempty(fData)
    fData = load('rosenbrock_func_data');
    if length(fData.o)>=D
         fData.o=fData.o(1:D);
    else
         fData.o=-90+180*rand(1,D);
    end
end
x=x-repmat(fData.o,ps,1)+1;
f=sum(100.*(x(:,1:D-1).^2-x(:,2:D)).^2+(x(:,1:D-1)-1).^2,2);
end

% 	7.Shifted Rotated Griewank's Function
function f=griewank_rot_func(x)
persistent fData matData
[ps,D]=size(x);
if isempty(fData)
    fData = load('griewank_func_data');
    if length(fData.o)>=D
         fData.o=fData.o(1:D);
    else
         fData.o=-600+0*rand(1,D);
    end
    c=3;
    if D==2,matData = load('griewank_M_D2');
    elseif D==10,matData = load('griewank_M_D10'); 
    elseif D==30,matData = load('griewank_M_D30'); 
    elseif D==50,matData = load('griewank_M_D50'); 
    else 
        matData.M=rot_matrix(D,c);
        matData.M=matData.M.*(1+0.3.*normrnd(0,1,D,D));
    end
    fData.o=fData.o(1:D);
end
x=x-repmat(fData.o,ps,1);
x=x*matData.M;
f=1;
for i=1:D
    f=f.*cos(x(:,i)./sqrt(i));
end
f=sum(x.^2,2)./4000-f+1;
end

% 	8.Shifted Rotated Ackley's Function with Global Optimum on Bounds
function f=ackley_rot_func(x)
persistent fData matData
[ps,D]=size(x);
if isempty(fData)
    fData = load('ackley_func_data');
    if length(fData.o)>=D
         fData.o=fData.o(1:D);
    else
         fData.o=-30+60*rand(1,D);
    end
    fData.o(2.*[1:floor(D/2)]-1)=-32;
    c=100;
    if D==2,matData = load('ackley_M_D2');
    elseif D==10,matData = load('ackley_M_D10');
    elseif D==30,matData = load('ackley_M_D30');
    elseif D==50,matData = load('ackley_M_D50');
    else 
       matData.M=rot_matrix(D,c);
    end
end
x=x-repmat(fData.o,ps,1);
x=x*matData.M;
f=sum(x.^2,2);
f=20-20.*exp(-0.2.*sqrt(f./D))-exp(sum(cos(2.*pi.*x),2)./D)+exp(1);
end

% 	9.Shifted Rastrign's Function
function f=rastrigin_func(x)
persistent fData
[ps,D]=size(x);
if isempty(fData)
    fData = load('rastrigin_func_data');
    if length(fData.o)>=D
         fData.o=fData.o(1:D);
    else
         fData.o=-5+10*rand(1,D);
    end
end
x=x-repmat(fData.o,ps,1);
f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);
end

% 	10.Shifted Rotated Rastrign's Function 
function f=rastrigin_rot_func(x)
persistent fData matData
[ps,D]=size(x);
if isempty(fData)
    fData = load('rastrigin_func_data');
    if length(fData.o)>=D
         fData.o=fData.o(1:D);
    else
         fData.o=-5+10*rand(1,D);
    end
    c=2;
    if D==2,matData = load('rastrigin_M_D2');
    elseif D==10,matData = load('rastrigin_M_D10');
    elseif D==30,matData = load('rastrigin_M_D30');
    elseif D==50,matData = load('rastrigin_M_D50');
    else 
        matData.M=rot_matrix(D,c);
    end
end
x=x-repmat(fData.o,ps,1);
x=x*matData.M;
f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);
end

% 	11.Shifted Rotated Weierstrass Function
function [f]=weierstrass_rot(x)
persistent fData matData
[ps,D]=size(x);
if isempty(fData)
    fData = load('weierstrass_data');
    if length(fData.o)>=D
         fData.o=fData.o(1:D);
    else
         fData.o=-0.5+0.5*rand(1,D);
    end
    c=5;
    if D==2,matData = load('weierstrass_M_D2');
    elseif D==10,matData = load('weierstrass_M_D10');
    elseif D==30,matData = load('weierstrass_M_D30');
    elseif D==50,matData = load('weierstrass_M_D50');
    else 
        matData.M=rot_matrix(D,c);
    end
end
x=x-repmat(fData.o,ps,1);
x=x*matData.M;
x=x+0.5;
a = 0.5;%0<a<1
b = 3;
kmax = 20;
[ps,D]=size(x);

c1(1:kmax+1) = a.^(0:kmax);
c2(1:kmax+1) = 2*pi*b.^(0:kmax);
c=-w(0.5,c1,c2);
f=0;
for i=1:D
f=f+w(x(:,i)',c1,c2);
end
f=f+repmat(c*D,ps,1);

%--------------------------------
end

% 	12.Schwefel's Problem 2.13
function f=schwefel_213(x)%after Fletcher and Powell
persistent fData A 
[ps,D]=size(x);
if isempty(fData)
    fData = load('schwefel_213_data');
    if length(fData.alpha)>=D
        fData.alpha=fData.alpha(1:D);fData.a=fData.a(1:D,1:D);fData.b=fData.b(1:D,1:D);
    else
        fData.alpha=-3+6*rand(1,D);
        fData.a=round(-100+200.*rand(D,D));
        fData.b=round(-100+200.*rand(D,D));
    end
    fData.alpha=repmat(fData.alpha,D,1);
    A=sum(fData.a.*sin(fData.alpha)+fData.b.*cos(fData.alpha),2);
end

for i=1:ps
    xx=repmat(x(i,:),D,1);
    B=sum(fData.a.*sin(xx)+fData.b.*cos(xx),2);
    f(i,1)=sum((A-B).^2,1);
end
end

% 	13. Expanded Extended Griewank's plus Rosenbrock's Function (F8F2)
function fit=EF8F2_func(x)
%-3,1
persistent fData
[ps,D]=size(x);
if isempty(fData)
    fData = load('EF8F2_func_data');
    if length(fData.o)>=D
         fData.o=fData.o(1:D);
    else
         fData.o=-1+1*rand(1,D);
    end
end
x=x-repmat(fData.o,ps,1)+1;
fit=0;
for i=1:(D-1)
    fit=fit+F8F2(x(:,[i,i+1]));
end
    fit=fit+F8F2(x(:,[D,1]));
end

function f=F8F2(x)
f2=100.*(x(:,1).^2-x(:,2)).^2+(1-x(:,1)).^2;
f=1+f2.^2./4000-cos(f2);
end

% ---------------------------------------------------------------  
% 	14. Expanded Rotated Extended Scaffer's F6 	
function f=E_ScafferF6_func(x)
persistent  fData matData
fhd=str2func('ScafferF6');
[ps,D]=size(x);
if isempty(fData)
    fData = load('E_ScafferF6_func_data');
    if length(fData.o)>=D
         fData.o=fData.o(1:D);
    else
         fData.o=-100+200*rand(1,D);
    end

    c=3;
    if D==2,matData = load('E_ScafferF6_M_D2');
    elseif D==10,matData = load('E_ScafferF6_M_D10');
    elseif D==30,matData = load('E_ScafferF6_M_D30');
    elseif D==50,matData = load('E_ScafferF6_M_D50');
    else 
       matData.M=rot_matrix(D,c);
    end
end
x=x-repmat(fData.o,ps,1);
x=x*matData.M;
f=0;
for i=1:(D-1)
    f=f+feval(fhd,(x(:,i:i+1)));
end
    f=f+feval(fhd,x(:,[D,1]));
end

function f=ScafferF6(x)
f=0.5+(sin(sqrt(x(:,1).^2+x(:,2).^2)).^2-0.5)./(1+0.001*(x(:,1).^2+x(:,2).^2)).^2;
end

%---------------------------------------------------
%   15.Hybrid Composition Function 1
function fit=hybrid_func1(x)
persistent fData  fun_num func sigma lamda bias M
if isempty(fData)
    fData = load('hybrid_func1_data'); % saved the predefined optima
    [ps,D]=size(x);
    fun_num=10;
    if length(fData.o(1,:))>=D
         fData.o=fData.o(:,1:D);
    else
         fData.o=-5+10*rand(fun_num,D);
    end
    func.f1=str2func('frastrigin');
    func.f2=str2func('frastrigin');
    func.f3=str2func('fweierstrass');
    func.f4=str2func('fweierstrass');
    func.f5=str2func('fgriewank');
    func.f6=str2func('fgriewank');
    func.f7=str2func('fackley');
    func.f8=str2func('fackley');
    func.f9=str2func('fsphere');
    func.f10=str2func('fsphere');
    bias=((1:fun_num)-1).*100;
    sigma=ones(1,fun_num);
    lamda=[1; 1; 10; 10; 5/60; 5/60; 5/32; 5/32; 5/100; 5/100];
    lamda=repmat(lamda,1,D);
    for i=1:fun_num
        eval(['M.M' int2str(i) '=diag(ones(1,D));']);
    end
end
fit=hybrid_composition_func(x,fun_num,func,fData.o,sigma,lamda,bias,M);
end

%---------------------------------------------------------------------
%   16.Rotated Hybrid Composition Function 1	
function fit=hybrid_rot_func1(x)
persistent  fun_num func fData sigma lamda bias matData
if isempty(fData)
    fData = load('hybrid_func1_data'); % saved the predefined optima
    [ps,D]=size(x);
    fun_num=10;
    if length(fData.o(1,:))>=D
         fData.o=fData.o(:,1:D);
    else
         fData.o=-5+10*rand(fun_num,D);
    end
    func.f1=str2func('frastrigin');
    func.f2=str2func('frastrigin');
    func.f3=str2func('fweierstrass');
    func.f4=str2func('fweierstrass');
    func.f5=str2func('fgriewank');
    func.f6=str2func('fgriewank');
    func.f7=str2func('fackley');
    func.f8=str2func('fackley');
    func.f9=str2func('fsphere');
    func.f10=str2func('fsphere');
    bias=((1:fun_num)-1).*100;
    sigma=ones(1,fun_num); 
    lamda=[1; 1; 10; 10; 5/60; 5/60; 5/32; 5/32; 5/100; 5/100];
    lamda=repmat(lamda,1,D);
    c=[2,2,2,2,2,2,2,2,2,2,2];
    if D==2,matData = load('hybrid_func1_M_D2');
    elseif D==10,matData = load('hybrid_func1_M_D10');
    elseif D==30,matData = load('hybrid_func1_M_D30');
    elseif D==50,matData = load('hybrid_func1_M_D50');
    else 
        for i=1:fun_num
            eval(['M.M' int2str(i) '=rot_matrix(D,c(i));']);
        end
    end
end
fit=hybrid_composition_func(x,fun_num,func,fData.o,sigma,lamda,bias,matData.M);
end
%----------------------------------------------------------------
%   17.	Rotated Hybrid Composition Function 1 with Noise in Fitness	
function fit=hybrid_rot_func1_noise(x)
[ps,D]=size(x);
fit=hybrid_rot_func1(x).*(1+0.2.*abs(normrnd(0,1,ps,1)));
end
%----------------------------------------------------------------
%   18.	Rotated Hybrid Composition Function 2
function fit=hybrid_rot_func2(x)
persistent  fun_num func fData sigma lamda bias matData
if isempty(fData)
    fData = load('hybrid_func2_data'); % saved the predefined optima
    [ps,D]=size(x);
    fun_num=10;
    if length(fData.o(1,:))>=D
         fData.o=fData.o(:,1:D);
    else
         fData.o=-5+10*rand(fun_num,D);
    end
    fData.o(10,:)=0;
    func.f1=str2func('fackley');
    func.f2=str2func('fackley');
    func.f3=str2func('frastrigin');
    func.f4=str2func('frastrigin');
    func.f5=str2func('fsphere');
    func.f6=str2func('fsphere');
    func.f7=str2func('fweierstrass');
    func.f8=str2func('fweierstrass');
    func.f9=str2func('fgriewank');
    func.f10=str2func('fgriewank');
    bias=((1:fun_num)-1).*100;
    sigma=[1 2 1.5 1.5 1 1 1.5 1.5 2 2];
    lamda=[2*5/32; 5/32; 2*1; 1; 2*5/100; 5/100; 2*10; 10; 2*5/60; 5/60];
    lamda=repmat(lamda,1,D);
    c=[2 3 2 3 2 3 20 30 200 300];
    if D==2,matData = load('hybrid_func2_M_D2');
    elseif D==10,matData = load('hybrid_func2_M_D10');
    elseif D==30,matData = load('hybrid_func2_M_D30');
    elseif D==50,matData = load('hybrid_func2_M_D50');
    else 
        for i=1:fun_num
            eval(['M.M' int2str(i) '=rot_matrix(D,c(i));']);
        end
    end
end
fit=hybrid_composition_func(x,fun_num,func,fData.o,sigma,lamda,bias,matData.M);
end
%----------------------------------------------------------------
%   19.	Rotated Hybrid Composition Function 2 with a Narrow Basin for the Global Optimum
function fit=hybrid_rot_func2_narrow(x)
persistent  fun_num func fData sigma lamda bias matData
if isempty(fData)
    fData = load('hybrid_func2_data'); % saved the predefined optima
    [ps,D]=size(x);
    fun_num=10;
    if length(fData.o(1,:))>=D
         fData.o=fData.o(:,1:D);
    else
         fData.o=-5+10*rand(fun_num,D);
    end
    fData.o(10,:)=0;
    func.f1=str2func('fackley');
    func.f2=str2func('fackley');
    func.f3=str2func('frastrigin');
    func.f4=str2func('frastrigin');
    func.f5=str2func('fsphere');
    func.f6=str2func('fsphere');
    func.f7=str2func('fweierstrass');
    func.f8=str2func('fweierstrass');
    func.f9=str2func('fgriewank');
    func.f10=str2func('fgriewank');
    bias=((1:fun_num)-1).*100;
    sigma=[0.1 2 1.5 1.5 1 1 1.5 1.5 2 2];
    lamda=[0.1*5/32; 5/32; 2*1; 1; 2*5/100; 5/100; 2*10; 10; 2*5/60; 5/60];
    lamda=repmat(lamda,1,D);
    c=[2 3 2 3 2 3 20 30 200 300];
    if D==2,matData = load('hybrid_func2_M_D2');
    elseif D==10,matData = load('hybrid_func2_M_D10');
    elseif D==30,matData = load('hybrid_func2_M_D30');
    elseif D==50,matData = load('hybrid_func2_M_D50');
    else 
        for i=1:fun_num
            eval(['M.M' int2str(i) '=rot_matrix(D,c(i));']);
        end
    end
end
fit=hybrid_composition_func(x,fun_num,func,fData.o,sigma,lamda,bias,matData.M);
end
%----------------------------------------------------------------
%   20.	Rotated Hybrid Composition Function 2 with the Global Optimum on the Bounds	
function fit=hybrid_rot_func2_onbound(x)
persistent  fun_num func fData sigma lamda bias matData
if isempty(fData)
    fData = load('hybrid_func2_data'); % saved the predefined optima
    [ps,D]=size(x);
    fun_num=10;
    if length(fData.o(1,:))>=D
         fData.o=fData.o(:,1:D);
    else
         fData.o=-5+10*rand(fun_num,D);
    end
    fData.o(10,:)=0;
    fData.o(1,2.*[1:floor(D/2)])=5;
    func.f1=str2func('fackley');
    func.f2=str2func('fackley');
    func.f3=str2func('frastrigin');
    func.f4=str2func('frastrigin');
    func.f5=str2func('fsphere');
    func.f6=str2func('fsphere');
    func.f7=str2func('fweierstrass');
    func.f8=str2func('fweierstrass');
    func.f9=str2func('fgriewank');
    func.f10=str2func('fgriewank');
    bias=((1:fun_num)-1).*100;
    sigma=[1 2 1.5 1.5 1 1 1.5 1.5 2 2];
    lamda=[2*5/32; 5/32; 2*1; 1; 2*5/100; 5/100; 2*10; 10; 2*5/60; 5/60];
    lamda=repmat(lamda,1,D);
    c=[2 3 2 3 2 3 20 30 200 300];
    if D==2,matData = load('hybrid_func2_M_D2');
    elseif D==10,matData = load('hybrid_func2_M_D10');
    elseif D==30,matData = load('hybrid_func2_M_D30');
    elseif D==50,matData = load('hybrid_func2_M_D50');
    else 
        for i=1:fun_num
            eval(['M.M' int2str(i) '=rot_matrix(D,c(i));']);
        end
    end
end
fit=hybrid_composition_func(x,fun_num,func,fData.o,sigma,lamda,bias,matData.M);
end
%-------------------------------------------------
%    21.Rotated Hybrid Composition Function 3		
function fit=hybrid_rot_func3(x)
persistent  fun_num func fData sigma lamda bias matData
if isempty(fData)
    fData = load('hybrid_func3_data'); % saved the predefined optima, a 10*1000 matrix
    [ps,D]=size(x);
    fun_num=10;
    if length(fData.o(1,:))>=D
         fData.o=fData.o(:,1:D);
    else
         fData.o=-5+10*rand(fun_num,D);
    end
    func.f1=str2func('fE_ScafferF6');
    func.f2=str2func('fE_ScafferF6');
    func.f3=str2func('frastrigin');
    func.f4=str2func('frastrigin');
    func.f5=str2func('fEF8F2');
    func.f6=str2func('fEF8F2');
    func.f7=str2func('fweierstrass');
    func.f8=str2func('fweierstrass');
    func.f9=str2func('fgriewank');
    func.f10=str2func('fgriewank');
    bias=((1:fun_num)-1).*100;
    sigma=[1,1,1,1,1,2,2,2,2,2];
    lamda=[5*5/100; 5/100; 5*1; 1; 5*1; 1; 5*10; 10; 5*5/200; 5/200];
    lamda=repmat(lamda,1,D);
    c=ones(1,D);
    if D==2,matData = load('hybrid_func3_M_D2');
    elseif D==10,matData = load('hybrid_func3_M_D10');
    elseif D==30,matData = load('hybrid_func3_M_D30');
    elseif D==50,matData = load('hybrid_func3_M_D50');
    else 
        for i=1:fun_num
            A=normrnd(0,1,D,D);
            eval(['M.M' int2str(i) '=cGram_Schmidt(A));']);
        end
    end
end
fit=hybrid_composition_func(x,fun_num,func,fData.o,sigma,lamda,bias,matData.M);
end
%-----------------------------------------
%   22.	Rotated Hybrid Composition Function 3 with High Condition Number Matrix
function fit=hybrid_rot_func3_highcond(x)
persistent  fun_num func fData sigma lamda bias matData
if isempty(fData)
    fData = load('hybrid_func3_data'); % saved the predefined optima, a 10*1000 matrix
    [ps,D]=size(x);
    fun_num=10;
    if length(fData.o(1,:))>=D
         fData.o=fData.o(:,1:D);
    else
         fData.o=-5+10*rand(fun_num,D);
    end
    func.f1=str2func('fE_ScafferF6');
    func.f2=str2func('fE_ScafferF6');
    func.f3=str2func('frastrigin');
    func.f4=str2func('frastrigin');
    func.f5=str2func('fEF8F2');
    func.f6=str2func('fEF8F2');
    func.f7=str2func('fweierstrass');
    func.f8=str2func('fweierstrass');
    func.f9=str2func('fgriewank');
    func.f10=str2func('fgriewank');
    bias=((1:fun_num)-1).*100;
    sigma=[1,1,1,1,1,2,2,2,2,2];
    lamda=[5*5/100; 5/100; 5*1; 1; 5*1; 1; 5*10; 10; 5*5/200; 5/200];
    lamda=repmat(lamda,1,D);
    c=[10 20 50 100 200 1000 2000 3000 4000 5000];
    if D==2,matData = load('hybrid_func3_HM_D2');
    elseif D==10,matData = load('hybrid_func3_HM_D10');
    elseif D==30,matData = load('hybrid_func3_HM_D30');
    elseif D==50,matData = load('hybrid_func3_HM_D50');
    else 
        for i=1:fun_num
            eval(['M.M' int2str(i) '=rot_matrix(D,c(i));']);
        end
    end
end
fit=hybrid_composition_func(x,fun_num,func,fData.o,sigma,lamda,bias,matData.M);
end
%-----------------------------------------
%   23.	Non-Continuous Rotated Hybrid Composition Function 3
function fit=hybrid_rot_func3_noncont(x)
persistent  fData 
[ps,D]=size(x);
if isempty(fData)
    fData = load('hybrid_func3_data'); % saved the predefined optima, a 10*1000 matrix
    fData.o=fData.o(1,1:D);
end
fData.o=repmat(fData.o,ps,1);
x=(abs(x-fData.o)<0.5).*x+(abs(x-fData.o)>=0.5).*(round(x.*2)./2);
fit=hybrid_rot_func3(x);
end
%-----------------------------------------
%   24.	Rotated Hybrid Composition Function 4	
function fit=hybrid_rot_func4(x)
persistent  fun_num func fData sigma lamda bias matData
if isempty(fData)
    fData = load('hybrid_func4_data'); % saved the predefined optima, a 10*1000 matrix
    [ps,D]=size(x);
    fun_num=10;
    if length(fData.o(1,:))>=D
         fData.o=fData.o(:,1:D);
    else
         fData.o=-5+10*rand(fun_num,D);
    end
    func.f1=str2func('fweierstrass');
    func.f2=str2func('fE_ScafferF6');
    func.f3=str2func('fEF8F2');
    func.f4=str2func('fackley');
    func.f5=str2func('frastrigin');
    func.f6=str2func('fgriewank');
    func.f7=str2func('fE_ScafferF6_noncont');
    func.f8=str2func('frastrigin_noncont');
    func.f9=str2func('felliptic');
    func.f10=str2func('fsphere_noise');
    bias=((1:fun_num)-1).*100;
    sigma=[2,2,2,2,2,2,2,2,2,2];
    lamda=[10; 5/20; 1; 5/32; 1; 5/100 ; 5/50; 1; 5/100; 5/100; ];
    lamda=repmat(lamda,1,D);
    c=[100 50 30 10 5 5 4 3 2 2];
    if D==2,matData = load('hybrid_func4_M_D2');
    elseif D==10,matData = load('hybrid_func4_M_D10');
    elseif D==30,matData = load('hybrid_func4_M_D30');
    elseif D==50,matData = load('hybrid_func4_M_D50');
    else 
        for i=1:fun_num
            eval(['M.M' int2str(i) '=rot_matrix(D,c(i));']);
        end
    end
end
fit=hybrid_composition_func(x,fun_num,func,fData.o,sigma,lamda,bias,matData.M);
end

%----------------------------------
function fit=hybrid_composition_func(x,fun_num,func,o,sigma,lamda,bias,M)
[ps,D]=size(x);
for i=1:fun_num
    oo=repmat(o(i,:),ps,1);
    weight(:,i)=exp(-sum((x-oo).^2,2)./2./(D*sigma(i)^2));
end

[tmp,tmpid]=sort(weight,2);
for i=1:ps
    weight(i,:)=(weight(i,:)==tmp(i,fun_num)).*weight(i,:)+(weight(i,:)~=tmp(i,fun_num)).*(weight(i,:).*(1-tmp(i,fun_num).^10));
end
weight=weight./repmat(sum(weight,2),1,fun_num);

fit=0;
for i=1:fun_num
    oo=repmat(o(i,:),ps,1);
    eval(['f=feval(func.f' int2str(i) ',((x-oo)./repmat(lamda(i,:),ps,1))*M.M' int2str(i) ');']);
    x1=5*ones(1,D);
    eval(['f1=feval(func.f' int2str(i) ',(x1./lamda(i,:))*M.M' int2str(i) ');']);
    fit1=2000.*f./f1;
    fit=fit+weight(:,i).*(fit1+bias(i));
end
end
%-------------------------------------------------
%basic functions

function f=fsphere(x)
%Please notice there is no use to rotate a sphere function, with rotation
%here just for a similar structure as other functions and easy programming
[ps,D]=size(x);
f=sum(x.^2,2);
end
%--------------------------------
function f=fsphere_noise(x)
[ps,D]=size(x);
f=sum(x.^2,2).*(1+0.1.*normrnd(0,1,ps,1));
end
%--------------------------------
function f=fgriewank(x)
[ps,D]=size(x);
f=1;
for i=1:D
    f=f.*cos(x(:,i)./sqrt(i));
end
f=sum(x.^2,2)./4000-f+1;
end
%--------------------------------
function f=fackley(x)
[ps,D]=size(x);
f=sum(x.^2,2);
f=20-20.*exp(-0.2.*sqrt(f./D))-exp(sum(cos(2.*pi.*x),2)./D)+exp(1);
end
%--------------------------------
function f=frastrigin(x)
[ps,D]=size(x);
f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);
end
%--------------------------------
function f=frastrigin_noncont(x)
[ps,D]=size(x);
x=(abs(x)<0.5).*x+(abs(x)>=0.5).*(round(x.*2)./2);
f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);
end
%--------------------------------
function [f]=fweierstrass(x)
[ps,D]=size(x);
x=x+0.5;
a = 0.5;
b = 3;
kmax = 20;
c1(1:kmax+1) = a.^(0:kmax);
c2(1:kmax+1) = 2*pi*b.^(0:kmax);
f=0;
c=-w(0.5,c1,c2);
for i=1:D
f=f+w(x(:,i)',c1,c2);
end
f=f+c*D;
end

function y = w(x,c1,c2)
y = zeros(length(x),1);
for k = 1:length(x)
	y(k) = sum(c1 .* cos(c2.*x(:,k)));
end
end
%--------------------------------
function f=fE_ScafferF6(x)
fhd=str2func('ScafferF6');
[ps,D]=size(x);

f=0;
for i=1:(D-1)
    f=f+feval(fhd,(x(:,i:i+1)));
end
    f=f+feval(fhd,x(:,[D,1]));
end
%--------------------------------    
function f=fE_ScafferF6_noncont(x)
fhd=str2func('ScafferF6');
[ps,D]=size(x);
x=(abs(x)<0.5).*x+(abs(x)>=0.5).*(round(x.*2)./2);
f=0;
for i=1:(D-1)
    f=f+feval(fhd,(x(:,i:i+1)));
end
    f=f+feval(fhd,x(:,[D,1]));
end    
%------------------------------
function f=fEF8F2(x)
[ps,D]=size(x);
f=0;
for i=1:(D-1)
    f=f+F8F2(x(:,[i,i+1]));
end
    f=f+F8F2(x(:,[D,1]));
end
%--------------------------------
function f=fschwefel_102(x)
[ps,D]=size(x);
f=0;
for i=1:D
    f=f+sum(x(:,1:i),2).^2;
end
end
%--------------------------------
function f=felliptic(x)
[ps,D]=size(x);
a=1e+6;
f=0;
for i=1:D
f=f+a.^((i-1)/(D-1)).*x(:,i).^2;
end
end
%--------------------------------
% classical Gram Schmid 
 function [q,r] = cGram_Schmidt (A)
% computes the QR factorization of $A$ via
% classical Gram Schmid 
% 
 [n,m] = size(A); 
 q = A;    
 for j=1:m
     for i=1:j-1 
         r(i,j) = q(:,j)'*q(:,i);
     end
     for i=1:j-1   
       q(:,j) = q(:,j) -  r(i,j)*q(:,i);
     end
     t =  norm(q(:,j),2 ) ;
     q(:,j) = q(:,j) / t ;
     r(j,j) = t  ;
 end 
 end

function M=rot_matrix(D,c)
A=normrnd(0,1,D,D);
P=cGram_Schmidt(A);
A=normrnd(0,1,D,D);
Q=cGram_Schmidt(A);
u=rand(1,D);
D=c.^((u-min(u))./(max(u)-min(u)));
D=diag(D);
M=P*D*Q;
end