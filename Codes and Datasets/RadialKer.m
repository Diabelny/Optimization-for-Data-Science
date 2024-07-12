function k= RadialKer(x,y,s)
k=exp(-(norm(x-y))^2/(2*s^2));
end
