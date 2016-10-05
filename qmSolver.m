function qmSolver(x, V, m, num)

% x is the distance vector
% V is V(x) the voltage potential (eV) with respect to x.
% m is mass.
% num is number of eigenvalues to evaluate
%%%% By Eric Perez 

%%%%%%%%%%  Constants   %%%%%%%%%%%%
hbar = 6.58211951440e-16;
% e = 1.60217662e-19;
me = 9.10938356e-31;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~nargin
    V = [zeros(10,1) ; 0.2 * ones(20,1); zeros(10,1); ...
        zeros(10,1) ; 0.2 * ones(20,1); zeros(10,1); ...
        zeros(10,1) ; 0.2 * ones(20,1); zeros(10,1); ...
        zeros(10,1) ; 0.2 * ones(20,1); zeros(10,1); ...
        zeros(10,1) ; 0.2 * ones(20,1); zeros(10,1); ...
        zeros(10,1) ; 0.2 * ones(20,1); zeros(10,1); ...
        zeros(10,1) ; 0.2 * ones(20,1); zeros(10,1); ...
        zeros(10,1) ; 0.2 * ones(20,1); zeros(10,1); ...
        zeros(10,1) ; 0.2 * ones(20,1); zeros(10,1); ...
        zeros(10,1) ; 0.2 * ones(20,1); zeros(10,1); ...
        ];
    m = 1*me;
    x = 1:numel(V);
    num = 20;
end

%%%%% Time-Independent Schrodinger Eqn
%%% d^2(u(x)/dx^2 + (2m / hBar) * (E-V(x))*u(x) = 0
%%%%% 
%%%% Will Solve the Hamiltonian to get values of the 
%%%% energy of the system

%%%% H*U(x,t) = i*hbar* d/dt(U(x,t))
%%%% Hamiltonian is equal to:
%%%% -hbar^2/2m * dx^2 + V(x)

% Get stepSize (assume evenly spaced)


%%%% Using Method from
%  https://wiki.physics.udel.edu/phys824/Discretization_of_1D_Hamiltonian

a = x(2) - x(1);
t = (hbar)^2 /(2*m*a^2);

H = zeros(numel(V),numel(V));

for i = 1:numel(V)
    H(i,i) = V(i) + 2*t;
    if i < numel(V)
        H(i,i+1) = -t;
    end
    if i > 1
        H(i,i-1) = -t;
    end
end


[vectors,values] = eig(H);
figure(1)
plot(x,V,'r','LineWidth',3); 
hold on;


% plot eign vectors
for i = 1:num
    u = values(i,i)+vectors(:,i);
%     u = normU(u);
    plot(x,u)
    
    
end
ylabel('Energy (eV)');

hold off;

figure(2)
plot(x,V,'r','LineWidth',3); 
hold on;
for i = 1:num;
    plot(x,ones(1,numel(x))*values(i,i));
    
end
ylabel('Energy (eV')
title('Energy Levels')
hold off;
end

function u = normU(u)
magArray = u.*conj(u);
mag = sum(magArray);

u = u./ (mag^0.5);
end

