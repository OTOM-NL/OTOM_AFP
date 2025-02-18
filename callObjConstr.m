
function [x, Fval,exitFlag,Output] = callObjConstr

Lb=[-180];
Ub=[-110];


% [x, fval]=ga(@myObjective,1,[] ,[], [], [],Lb,Ub,@ellipseparabola);


A = []; b = []; % No linear inequality constraints
Aeq = []; beq = []; % No linear equality constraints
options = gaoptimset('PlotFcns',@gaplotpareto);

[x,Fval,exitFlag,Output]=gamultiobj(@myObjective,1,A, ...
    b,Aeq,beq,Lb,Ub,options);
fprintf('The number of points on the Pareto front was: %d\n', size(x,1));
fprintf('The number of generations was : %d\n', Output.generations);
plot3(x,Fval(:,1),Fval(:,2),'r*')
s=sprintf('Number of design= %d, Pareto design= %d',Output.funccount,size(x,1));
title(['Multi-obj-optimization, ',s]);
xlabel('Var');
ylabel('Obj1');
zlabel('Obj2');
grid on



% [x fval] = fmincon(@myObjective,x0,[],[],[],[],[],[], ...
%     @ellipseparabola,options);

function f = myObjective(xin)
Temp=thermal_domain_generator_opt(xin);
f=[norm(abs(Temp-300)), std(Temp)];
disp(f);
% f = 3*(xin(1)-xin(2))^4 + 4*(xin(1)+xin(3))^2/(1+sum(xin.^2)) ...
%     + cosh(xin(1)-1) + tanh(xin(2)+xin(3));

function [c,ceq] = ellipseparabola(x)
% c(1) = (x(1)^2)/9 + (x(2)^2)/4 - 1;
% c(2) = x(1)^2 - x(2) - 1;
c=[];
ceq = [];