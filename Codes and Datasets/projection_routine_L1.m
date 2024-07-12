
% Routine per la proiezione della direzione nel cono delle direzioni
% ammissibili

function [lambda_n,f_min, norma, alfa] = projection_routine_L1(Q,q,lambda_value,anti_grad,N,epsilon,C,tol, isin_border,L1)
    if isin_border
         A=zeros(N,N); % Matrice vincoli attivi
         rows_A=0;
            for i=1:N
              if lambda_value(i)<=tol
                   v=zeros(1,N);
                   v(i)=-1;
                   A(i,:)=v;
                   rows_A=rows_A+1;
              end
              if lambda_value(i)>=C-tol
                   v=zeros(1,N);
                   v(i)=1;
                   A(i,:)=v;
                   rows_A=rows_A+1;
              end
            end
           A=A(1:rows_A,:);
           v=A*anti_grad';

           isin_feasible=~any(v>tol); 
           if (isin_feasible && abs(sum(anti_grad(1:N/2))-sum(anti_grad(N/2+1:end)))<=tol)
               d_star=anti_grad;
           else
               grad_lambda=-anti_grad;
               d_star=AP(grad_lambda,lambda_value,N,C,tol); % routine che risolve il  problema CQKnP (sottoproblema 1.5 report)  
           end
    else
        if abs(sum(anti_grad(1:N/2))-sum(anti_grad(N/2+1:end)))<=tol 
            d_star=anti_grad;
        else
           mu=(sum(anti_grad(1:N/2))-sum(anti_grad(N/2+1:end)))/N; % risoluzione P_KKT system
           d_star=-([ones(1,N/2),-ones(1,N/2)]*mu-anti_grad);
        end
    end
   
    norma=norm(d_star);
    
        if norma<=epsilon  
            f_min=1/2*lambda_value*Q*lambda_value'+q*lambda_value'; 
            alfa=0;
            lambda_n=lambda_value;
        else
            u=C*ones(1,N);
            ind = d_star> 0;  
            maxt = min((u(ind)-lambda_value(ind))./d_star(ind));
            ind = d_star<0;  
            maxt = min([maxt, min(-lambda_value(ind)./d_star(ind))]);
            alfa = 1/L1;
            alfa=min(alfa,maxt);
            lambda_n=lambda_value+alfa*d_star; 
            f_min=1/2*lambda_value*Q*lambda_value'+q*lambda_value';
       end
end