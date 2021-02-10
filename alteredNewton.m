function [Y,X,k] = alteredNewton(f,hessian,gradient,N_max,x0,e1,e2,e3,lowerLimit,upperLimit) 
% Optimization Theory
% Class Assignment -  Newton Algorithm Main File
% 28/11/2020  


k = 0;
xk = x0;

while 1
    
   gradR = gradient(xk(1),xk(2));
   hessR = hessian(xk(1),xk(2)) ;
   
   
   pk = -1*(inv(hessR))*gradR;
             
    
  g = @(sk) ( f(xk(1)+sk*pk(1),xk(2)+sk*pk(2)));
    
  [x1_g,~,~,~,~,~] = goldenSection(g,lowerLimit,upperLimit,1e-9);
    
    sk = x1_g           ;
    xk_old = xk         ;
    xk = xk + sk*pk     ;
    
    k = k +1            ;
    
    
    
     if k>N_max
        fprintf('Iteration exceeds the N_max\n');
        break;
    elseif abs(( f(xk(1),xk(2)) - f(xk_old(1),xk_old(2)) )) < e1
        fprintf('Fx+1 - Fx<e1\n');
        break;
    elseif norm(xk - xk_old)<e2
        fprintf('norm of x+1 - x<e2\n');
        break;
    elseif norm(gradient(xk(1),xk(2)))<e3
        fprintf('norm of Fx+1<e2\n');
        break;
    end
    
    
    
end


    Y = f(xk(1),xk(2));
    X = xk;


end