pro plot_random_walk,file

; read data file
data = float(rd_tfile(file,2))
step = reform(data[0,*])
rmsd = reform(data[1,*])
; calculate diffusion coefficient
tlog = alog(step)
dlog = alog(rmsd)
fit  = poly_fit(tlog,dlog,1)
print,'diffusion coefficient = ',trim(fit[1])
; plot results
p = plot(step,rmsd,layout=[1,2,1],title='RMSD vs Step')
p = plot(tlog,dlog,layout=[1,2,2],/current,title='log-log RMSD vs Step')

end
