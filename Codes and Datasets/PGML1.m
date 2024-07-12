
% Funzione principale algoritmo

function [sol, iter_true] = PGML1(Q,q,epsilon,C,lambda_0,N, max_iter,tol)
    L=eig(Q);
    L1=L(end); % autovalore piu grande
    tic;
    sol=zeros(1,N);
    anti_grad=-((Q*lambda_0')'+q); % -gradiente punto corrente lambda_0 
    isin_border = any(lambda_0<=tol | lambda_0>=C-tol); % isin_border specifica se il punto corrente si trova sul bordo oppure no
    curr_iter=1;
    iter_true=curr_iter;
    disp('----------------------------------------------------');
    fprintf( 'iter\tf(x)\t\t\t||nabla f(x)||\t alpha \n\n' );
    disp('----------------------------------------------------')
    [lambda_n,f_min, norma, alpha]=projection_routine_L1(Q,q,lambda_0,anti_grad,N,epsilon,C,tol, isin_border,L1);
    
    % Riempimento vettore output
    sol(1)=f_min;
    
    fprintf( '%4d\t%1.8e\t\t%1.4e\t %1.4e\n' , curr_iter , f_min , norma , alpha);
     
    while norma>epsilon % condizione terminazione
        if max_iter==curr_iter
            status='Stopped';
            iter_true=curr_iter;
            break
        end
        
       anti_grad=-((Q*lambda_n')'+q); % -gradiente punto corrente lambda_n 
       isin_border = any(lambda_n<= tol | lambda_n >= C-tol);
       [lambda_n,f_min,norma, alpha]=projection_routine_L1(Q,q,lambda_n,anti_grad,N,epsilon,C,tol, isin_border,L1);
       curr_iter=curr_iter+1;
       iter_true=curr_iter;
       sol(curr_iter)=f_min;
       fprintf( '%4d\t%1.8e\t\t%1.4e\t %1.4e\n\n' , curr_iter , f_min , norma , alpha);
    end

    if curr_iter<max_iter
        status='Optimal';
    end

    sol=sol(1,1:curr_iter);
    fprintf( 'stop: %d iter, status = %s, fbest = %1.8e\n,' , curr_iter , status , f_min );
    toc;
   
end


