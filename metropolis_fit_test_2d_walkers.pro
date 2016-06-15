; sampling test using Metropolis-Hastings random walk with multiple walkers
; uses MCMC hammer method for multiple walkers

pro metropolis_fit_test_2d_walkers

; inputs
; -----------------------------------------
; set monte carlo iterations
mc_iter = 1000L
; set number of walkers
nwalkers = 1000L
; display images
display = 'no'
; -----------------------------------------


; generate sample data
; -----------------------------------------
print,'generating data ...'
; generate x,y grid
dmin = -15
dmax =  15
step = 0.2
n    = (dmax-dmin)/step+1
x = dmin+step*findgen(n)
y = dmin+step*findgen(n)

; set parameters for first 2d gaussian
a1 = 5   ; peak
b1 = 0   ; center 1
c1 = 4   ; width 1
d1 = -5  ; center 2
e1 = 2   ; width 2

; set parameters for second 2d gaussian
a2 = 2   ; peak
b2 = -5  ; center 1
c2 = 2   ; width 1
d2 = 5   ; center 2
e2 = 2   ; width 2

; generate bimodal 2d gaussian
z = fltarr(n,n)
for i = 0, n-1 do begin
  for j = 0, n-1 do begin
    xi = x[i]
    yi = y[j]
    f1 = a1*exp(-((xi-b1)^2)/(2*c1^2))*exp(-((yi-d1)^2)/(2*e1^2))
    f2 = a2*exp(-((xi-b2)^2)/(2*c2^2))*exp(-((yi-d2)^2)/(2*e2^2))
    fi = f1+f2
    z[i,j] = fi
  endfor
endfor
; -----------------------------------------


; begin monte carlo sim
; -----------------------------------------
print,'running monte carlo ...'
; initialization
; ---------------
; initialize walkers
wmin = -15
wmax =  15
xi = randomu(seed,nwalkers)*(wmax-wmin)+wmin
yi = randomu(seed,nwalkers)*(wmax-wmin)+wmin
f1 = a1*exp(-((xi-b1)^2)/(2*c1^2))*exp(-((yi-d1)^2)/(2*e1^2))
f2 = a2*exp(-((xi-b2)^2)/(2*c2^2))*exp(-((yi-d2)^2)/(2*e2^2))
fi = f1+f2
xr = fltarr(mc_iter,nwalkers)
yr = fltarr(mc_iter,nwalkers)
zr = fltarr(mc_iter,nwalkers)
xr[0,*] = xi
yr[0,*] = yi
zr[0,*] = fi
xarr = xi
yarr = yi
zarr = fi

; initialize acceptance
accept = 0L

; generate random numbers for Metropolis
rand = randomu(seed,mc_iter,nwalkers)

; generate variables for MCMC hammer
aham = 2.0
zmin = sqrt(1.0/aham)
zmax = sqrt(aham)
rham = randomu(seed,mc_iter,nwalkers)
zham = (rham*(zmax-zmin)+zmin)^(2.0)
half = fix(0.5*nwalkers)
; ---------------

; begin simulation
for i = 1, mc_iter-1 do begin

  ; loop over sets of half-walkers
  for a = 0, 1 do begin

    ; split walkers into half-groups
    if a eq 0 then begin
      ; randomize walker pairs
      indx = permute(nwalkers)
      w1 = indx[0:half-1]
      w2 = indx[half:*]
    endif
    if a eq 1 then begin
      w1 = indx[half:*]
      w2 = indx[0:half-1]
    endif

    ; loop over half-walkers
    for b = 0, half-1 do begin

      ; get walker k from group 1
      k   = w1[b]
      xk0 = xarr[k]
      yk0 = yarr[k]

      ; get walker j from group 2
      j   = w2[b]
      xj0 = xarr[j]
      yj0 = yarr[j]

      ; compute stretch move
      xk = xj0 + zham[i,k]*(xk0-xj0)
      yk = yj0 + zham[i,k]*(yk0-yj0)

      ; get new values
      f1 = a1*exp(-((xk0-b1)^2)/(2*c1^2))*exp(-((yk0-d1)^2)/(2*e1^2))
      f2 = a2*exp(-((xk0-b2)^2)/(2*c2^2))*exp(-((yk0-d2)^2)/(2*e2^2))
      fold = f1+f2
      f1 = a1*exp(-((xk-b1)^2)/(2*c1^2))*exp(-((yk-d1)^2)/(2*e1^2))
      f2 = a2*exp(-((xk-b2)^2)/(2*c2^2))*exp(-((yk-d2)^2)/(2*e2^2))
      fnew = f1+f2

      ; accept/reject move
      prob = (fnew/fold)*zham[i,k]
      ; accept
      if prob gt 1 then begin
        xarr[k] = xk
        yarr[k] = yk
        zarr[k] = fnew
        accept += 1
      endif else begin
        ; otherwise accept with some probability
        if prob gt rand[i,k] then begin
          xarr[k] = xk
          yarr[k] = yk
          zarr[k] = fnew
          accept += 1
        endif
        ; otherwise reject
      endelse
    endfor
  endfor
  ; store values
  xr[i,*] = xarr
  yr[i,*] = yarr
  zr[i,*] = zarr
endfor
; -----------------------------------------

; display images
if display eq 'yes' then begin
  print,'displaying images...'
  for i = 0, mc_iter-1 do begin
    if i eq 0 then begin
      ; plot surface
      p1 = surface(z,x,y,dimensions=[1200,600],xrange=[dmin,dmax],yrange=[dmin,dmax],layout=[2,1,1])
      p1 = scatterplot3d(xr[i,*],yr[i,*],zr[i,*],symbol='dot',sym_size=2,xrange=[dmin,dmax],yrange=[dmin,dmax],/overplot)
      ; plot contour
      p2 = contour(z,x,y,layout=[2,1,2],xrange=[dmin,dmax],yrange=[dmin,dmax],/current)
      p2 = scatterplot(xr[i,*],yr[i,*],symbol='dot',sym_size=2,xrange=[dmin,dmax],yrange=[dmin,dmax],/overplot)
    endif else begin
      ; update plot
      data1 = p1.data
      data2 = p2.data
      data1[0,*] = xr[i,*]
      data1[1,*] = yr[i,*]
      data1[2,*] = zr[i,*]
      data2[0,*] = xr[i,*]
      data2[1,*] = yr[i,*]
      p1.putdata,data1
      p2.putdata,data2
    endelse
    wait,0.05
  endfor
endif

; save plot
p = contour(z,x,y,layout=[2,1,1],dimensions=[600,300],xrange=[dmin,dmax],yrange=[dmin,dmax],title='Initial',/buffer)
p = scatterplot(xr[0,*],yr[0,*],symbol='dot',sym_size=2,xrange=[dmin,dmax],yrange=[dmin,dmax],/overplot)
p = contour(z,x,y,layout=[2,1,2],xrange=[dmin,dmax],yrange=[dmin,dmax],/current,title='Final')
p = scatterplot(xr[-1,*],yr[-1,*],symbol='dot',sym_size=2,xrange=[dmin,dmax],yrange=[dmin,dmax],/overplot)
p.save,'metropolis_fit_test_2d_walkers.ps'

; print aceeptance
print,'acceptance rate : '+trim(accept/float(nwalkers*(mc_iter-1)))
print,'min/max walkers'
print,'x ',min(xarr),max(xarr)
print,'y ',min(yarr),max(yarr)
print,'mean, median, stddev'
print,'x ',mean(xarr),median(xarr),stddev(xarr)
print,'y ',mean(yarr),median(yarr),stddev(yarr)

end
