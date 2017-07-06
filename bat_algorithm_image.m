
% Usage: bat_algorithm_image_image([20 1000 0.5 0.5]);                 %


% -------------------------------------------------------------------     
% Implimented based on the following papers:                        %
% (Citation details):                                               %
% 1) Yang X.-S., A new metaheuristic bat-inspired algorithm,        %
%    in: Nature Inspired Cooperative Strategies for Optimization    %
%    (NISCO 2010) (Eds. J. R. Gonzalez et al.), Studies in          %
%    Computational Intelligence, Springer, vol. 284, 65-74 (2010).  %
% 2) Yang X.-S., Nature-Inspired Metaheuristic Algorithms,          %
%    Second Edition, Luniver Presss, Frome, UK. (2010).             %
% 3) Yang X.-S. and Gandomi A. H., Bat algorithm: A novel           %
%    approach for global engineering optimization,                  %
%    Engineering Computations, Vol. 29, No. 5, pp. 464-483 (2012).  %
% -------------------------------------------------------------------


% Main programs starts here
%function [best,fmin,psnrc,N_iter]=bat_algorithm_image(para,hi,lamda,gamma,dh,sm)
function [best,fmin,psnrc,N_iter]=bat_algorithm_image(n,N_gen,A,r,hi,lamda,gamma,dh,sm)
% Display help
% help bat_algorithm_image_image.m

if nargin<1
    n = 25;
end


% Default parameters


%if nargin<1,  para=[20 1000 0.5 0.5];  end
%n=para(1);      % Population size, typically 10 to 40
%N_gen=para(2);  % Number of generations
%A=para(3);      % Loudness  (constant or decreasing)
%r=para(4);      % Pulse rate (constant or decreasing)
% This frequency range determines the scalings
n=30;
N=10;
A=0.5;
r=0.5;


% You should change these values if necessary
Qmin=0;         % Frequency minimum
Qmax=2;         % Frequency maximum
% Iteration parameters
N_iter=30;       % Total number of function evaluations
% Dimension of the search variables
d=256;           % Number of dimensions 
% Lower limit/bounds/ a vector
%Lb=-2*ones(1,d);
% Upper limit/bounds/ a vector
%Ub=2*ones(1,d);
% Initializing arrays
Q=zeros(n,1);   % Frequency
v=zeros(n,d);   % Velocities

number_of_solution = 256;
Lb = zeros(1,number_of_solution);
Ub = 255.*ones(1,number_of_solution);


% Initialize the population/solutions
for i=1:n,
  Sol(i,:)=Lb+(Ub-Lb).*rand(1,d);
  Fitness(i)=Fun(Sol(i,:),hi,lamda,gamma,dh,sm);
end
% Find the initial best solution
[fmin,I]=min(Fitness);
best=Sol(I,:);

% ======================================================  %
% Note: As this is a demo, here we did not implement the  %
% reduction of loudness and increase of emission rates.   %
% Interested readers can do some parametric studies       %
% and also implementation various changes of A and r etc  %
% ======================================================  %

% Start the iterations -- Bat Algorithm (essential part)  %
for t=1:N_gen, 
% Loop over all bats/solutions
        for i=1:n,
          Q(i)=Qmin+(Qmin-Qmax)*rand;
          v(i,:)=v(i,:)+(Sol(i,:)-best)*Q(i);
          S(i,:)=Sol(i,:)+v(i,:);
          % Apply simple bounds/limits
          Sol(i,:)=simplebounds(Sol(i,:),Lb,Ub);
          % Pulse rate
          if rand>r
          % The factor 0.001 limits the step sizes of random walks 
              S(i,:)=best+0.001*randn(1,d);
          end

     % Evaluate new solutions
           Fnew=Fun(S(i,:),hi,lamda,gamma,dh,sm);
     % Update if the solution improves, or not too loud
           if (Fnew<=Fitness(i)) & (rand<A) ,
                Sol(i,:)=S(i,:);
                Fitness(i)=Fnew;
           end

          % Update the current best solution
          if Fnew<=fmin,
                best=S(i,:);
                fmin=Fnew;
          end
        end
        N_iter=N_iter+n;
end
% Output/display
disp(['Number of evaluations: ',num2str(N_iter)]);
disp(['Best =',num2str(best),' fmin=',num2str(fmin)]);

% Application of simple limits/bounds
function s=simplebounds(s,Lb,Ub)
  % Apply the lower bound vector
  ns_tmp=s;
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I);
  
  % Apply the upper bound vector 
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  % Update this new move 
  s=ns_tmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective function: your own objective function can be written here
% Note: When you use your own function, please remember to 
%       change limits/bounds Lb and Ub (see lines 52 to 55) 
%       and the number of dimension d (see line 51). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function z=Fun(u)
% Sphere function with fmin=0 at (0,0,...,0)
%z=sum(u.^2);


function z = Fun(u,hi,lamda,gamma,dh,sm)
            dh = dh*transpose(u);
            dh = sum(dh);
            s1 = sum((transpose(u)-hi).^2);
            s2 = sum((transpose(u)-sm).^2);
            s3 = sum(dh.^2);
            z = s1 + lamda*s2 + gamma*s3;

%%%%% ============ end ====================================
