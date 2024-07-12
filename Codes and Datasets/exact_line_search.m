function [a] = exact_line_search(Q,d_star,grad_lambda,C,N, lambda_value)
   
    u=C*ones(1,N);
    ind = d_star> 0;  
    maxt = min((u(ind)-lambda_value(ind))./d_star(ind));
    ind = d_star<0;  
    maxt = min([maxt, min((-lambda_value(ind))./d_star(ind))]);
    den = d_star*Q*d_star';
    
    if den <= 1e-16  
     a = maxt;     
    else
     a = min( [-( d_star * -grad_lambda' )/ den, maxt ] );
    end
    
end