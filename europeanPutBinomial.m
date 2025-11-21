function[]=europeanPutBinomial(S,K,T,r,sigma,baseSteps,iter)
% Input:
% S     = initial asset price
% E     = strike price
% T     = time to maturity
% r     = risk-free interest rate
% sigma = volatility
% M     = initial number of steps in the binomial tree
% k     = number of refinement iterations
%
% Output:
% Returns the approximated value of a European put option using the
% binomial model, the exact Black–Scholes price, and a plot showing
% convergence behavior.

% Binomial method

approx=zeros(1,iter);
x_app=1:iter;

for w=1:iter % At each iteration, increase the number of binomial steps
    steps=baseSteps*w;
    dt = T/steps;
    u=exp(sigma*sqrt(dt));
    d=1/u;
    p=(exp(r*dt)-d)/(u-d);
    
    % Triangular matrix of option values
    optionValues=zeros(steps+1,steps+1);
    % Terminal payoff column
    for i=1:steps+1
        optionValues(i,steps+1)=max(K-(S*u^(steps+1-i)*d^(i-1)),0);
    end

    for j=steps:-1:1
        for i=1:j % Backward induction from second-to-last column to the first
            optionValues(i,j)=exp(-r*dt)*((p*optionValues(i,j+1))+((1-p)*optionValues(i+1,j+1)));
        end
    end
    approx(w)=optionValues(1,1);
end

% Black–Scholes closed-form value
d1 = (log(S/K) + (r + 0.5*sigma^2)*T)/(sigma*sqrt(T));
d2 = d1 - sigma*sqrt(T);
N1 = 0.5*(1+erf(-d1/sqrt(2)));
N2 = 0.5*(1+erf(-d2/sqrt(2)));
BS_value = K*exp(-r*T)*N2 - S*N1;
disp('European put value (Black–Scholes): '), disp(BS_value)
   
% Convergence plot
plot(x_app,approx)
hold on
plot(x_app,BS_value*ones(1,iter))
legend('Binomial approximation','Black–Scholes value','Location','best')