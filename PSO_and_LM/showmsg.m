function showmsg(ite_num,maxit,misfit)
if ite_num==maxit
    fprintf(1,'Maximum iteration %d is reached£¡err=%e\n',[maxit;misfit]);
else
    fprintf(1,'Converged£¡Iteration£º%d, err=%e\n',[ite_num,misfit]);
end
