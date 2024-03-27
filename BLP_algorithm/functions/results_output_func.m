function results=results_output_func(iter_info,s_jt_predict,S_jt_data);
results=NaN(1,5);

TOL_DIST_s_jt=1e-10;
DIST_s_jt=max(abs(s_jt_predict(:)-S_jt_data(:)));

id=1;results(1,id)=iter_info.feval;
id=id+1;results(1,id)=round(iter_info.t_cpu,3);
id=id+1;results(1,id)=(iter_info.feval<iter_info.ITER_MAX &...
    iter_info.FLAG_ERROR==0);
id=id+1;results(1,id)=round(log10(DIST_s_jt),4);
id=id+1;results(1,id)=(results(1,4)<log10(TOL_DIST_s_jt));
end
