; sampling test using Metropolis-Hastings random walk

pro metropolis_fit_test_2d

; inputs
; -----------------------------------------
; set monte carlo iterations
mc_iter = 1000
; set step size
xstep = 4.0
ystep = 4.0
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
; initialize walker
xi = 8
yi = 8
f1 = a1*exp(-((xi-b1)^2)/(2*c1^2))*exp(-((yi-d1)^2)/(2*e1^2))
f2 = a2*exp(-((xi-b2)^2)/(2*c2^2))*exp(-((yi-d2)^2)/(2*e2^2))
fi = f1+f2
xr = fltarr(mc_iter)
yr = fltarr(mc_iter)
zr = fltarr(mc_iter)
xr[0] = xi
yr[0] = yi
zr[0] = fi

; initialize acceptance
accept = 0

; generate random steps
xrand = randomu(seed,mc_iter)*2*xstep-xstep
yrand = randomu(seed,mc_iter)*2*ystep-ystep

; generate random numbers for Metropolis
rand = randomu(seed,mc_iter)
; ---------------

; begin simulation
for i = 1, mc_iter-1 do begin

  xold = xr[i-1]
  yold = yr[i-1]
  f1 = a1*exp(-((xold-b1)^2)/(2*c1^2))*exp(-((yold-d1)^2)/(2*e1^2))
  f2 = a2*exp(-((xold-b2)^2)/(2*c2^2))*exp(-((yold-d2)^2)/(2*e2^2))
  fold = f1+f2

  xnew = xold + xrand[i]
  ynew = yold + yrand[i]
  f1 = a1*exp(-((xnew-b1)^2)/(2*c1^2))*exp(-((ynew-d1)^2)/(2*e1^2))
  f2 = a2*exp(-((xnew-b2)^2)/(2*c2^2))*exp(-((ynew-d2)^2)/(2*e2^2))
  fnew = f1+f2

  ; accept/reject move
  prob = fnew/fold
  ; accept
  if prob gt 1 then begin
    xr[i] = xnew
    yr[i] = ynew
    zr[i] = fnew
    accept += 1
  endif else begin
    ; otherwise accept with some probability
    if prob gt rand[i] then begin
      xr[i] = xnew
      yr[i] = ynew
      zr[i] = fnew
      accept += 1
    endif else begin
      ; otherwise reject
      xr[i] = xold
      yr[i] = yold
      zr[i] = fold
    endelse
  endelse
endfor
; -----------------------------------------

; display images
if display eq 'yes' then begin
  print,'displaying images...'
  for i = 0, mc_iter-1 do begin
    print,i,' ',xr[i],' ',yr[i]
    if i eq 0 then begin
      ; plot surface
      p = surface(z,x,y,dimensions=[1200,600],xrange=[dmin,dmax],yrange=[dmin,dmax],layout=[2,1,1])
      p = plot3d([xr[i]],[yr[i]],[zr[i]],'.',sym_size=2,xrange=[dmin,dmax],yrange=[dmin,dmax],/overplot)
      ; plot contour
      p = contour(z,x,y,layout=[2,1,2],xrange=[dmin,dmax],yrange=[dmin,dmax],/current)
      p = plot([xr[i]],[yr[i]],'.',sym_size=2,xrange=[dmin,dmax],yrange=[dmin,dmax],/overplot)
    endif else begin
      ; plot points
      p = plot([xr[i]],[yr[i]],'.',sym_size=4,xrange=[dmin,dmax],yrange=[dmin,dmax],/overplot)
    endelse
    wait,0.05
  endfor
endif else begin
  print,'making plot ...'
  ; plot surface
  p = surface(z,x,y,dimensions=[800,400],xrange=[dmin,dmax],yrange=[dmin,dmax],layout=[2,1,1])
  p = plot3d(xr,yr,zr,'.',sym_size=2,xrange=[dmin,dmax],yrange=[dmin,dmax],/overplot)
  ; plot contour
  p = contour(z,x,y,layout=[2,1,2],xrange=[dmin,dmax],yrange=[dmin,dmax],/current)
  p = plot(xr,yr,'.',sym_size=2,xrange=[dmin,dmax],yrange=[dmin,dmax],/overplot)
  p = plot([xr[0]],[yr[0]],'or',xrange=[dmin,dmax],yrange=[dmin,dmax],/overplot)
  p.save,'metropolis_fit_test_2d.ps'
endelse

; print aceeptance
print,'acceptance rate : '+trim(accept/float(mc_iter))

end
