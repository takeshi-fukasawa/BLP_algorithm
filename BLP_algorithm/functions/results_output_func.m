function results=results_output_func(iter_info,s_jt_predict,S_jt_data);
results=NaN(1,5);

TOL_DIST_s_jt=1e-10;
DIST_s_jt=max(abs(s_jt_predict(:)-S_jt_data(:)));

results(1,1)=iter_info.feval;
results(1,2)=iter_info.t_cpu;
results(1,3)=(iter_info.feval<iter_info.ITER_MAX);
results(1,4)=log10(DIST_s_jt);
results(1,5)=(results(1,4)<log10(TOL_DIST_s_jt));

end
