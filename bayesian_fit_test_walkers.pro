; gaussian fit model
; ---------------------------------------
function get_gaussian, p,x
a = p[0]  ; peak
b = p[1]  ; center
c = p[2]  ; width
d = p[3]  ; background
fit = a*exp(-((x-b)^2)/(2*c^2))+d
return,fit
end
; ---------------------------------------

; bayesian likelihood function
; ---------------------------------------
function get_ln_likelihood, y,yy,ysig
b = ((y-yy)*(y-yy))/(2.0*ysig*ysig)
ln_likelihood = total(-b,/nan)
return,ln_likelihood
end
; ---------------------------------------

; bayesian prior function
; ---------------------------------------
function get_ln_prior, p,pp,psig
b = ((p-pp)*(p-pp))/(2.0*psig*psig)
ln_prior = total(-b,/nan)
return,ln_prior
end
; ---------------------------------------

; bayesian posterior function
; ---------------------------------------
function get_ln_posterior, p,pp,psig,y,yy,ysig
ln_likelihood = get_ln_likelihood(y,yy,ysig)
ln_prior = get_ln_prior(p,pp,psig)
ln_posterior = ln_likelihood+ln_prior
return,ln_posterior
end
; ---------------------------------------


; fit test using Bayesian Metropolis-Hastings Markov Chain Monte Carlo (MCMC)
; uses MCMC hammer method for multiple walkers
pro bayesian_fit_test_walkers
tic

; inputs
; -----------------------------------------
; set monte carlo iterations
mc_iter = 1000L
; set number of walkers
nwalkers = 1000L
; set initial gaussian fit parameters
a0 = 210  ; peak
b0 = 1.0  ; center
c0 = 1.0  ; width
d0 = 25   ; background
; display parameter plots
pltdisplay = 'no'
; -----------------------------------------


; generate sample data
; -----------------------------------------
; generate x data
n = 101
x = -5+0.1*findgen(n)

; set parameters for gaussian
a = 200  ; peak
b = 0    ; center
c = 1.2  ; width
d = 30   ; background

; generate gaussian with poisson errors
y = get_gaussian([a,b,c,d],x)
y = poidev(y)
ysig = sqrt(y)
; -----------------------------------------


; begin monte carlo sim
; -----------------------------------------
print,'running monte carlo ...'
; initialization
; ---------------
; initialize walkers
p0    = [a0,b0,c0,d0]
npar  = n_elements(p0)
parr  = rebin(transpose(p0),nwalkers,npar)
r0    = randomn(seed,nwalkers,npar,/double)*0.2
parr  = parr+parr*r0

; display parameter plot
if pltdisplay eq 'yes' then begin
  plt1 = plot(parr[*,0],parr[*,3],'.',layout=[2,1,1],title='Background vs Peak',$
      dimensions=[1200,600],xrange=[100,300],yrange=[-30,70])
  plt = plot([a],[d],'or',/overplot)
  plt2 = plot(parr[*,1],parr[*,2],'.',layout=[2,1,2],title='Width vs Center',$
      xrange=[-1,1],yrange=[0,2],/current)
  plt = plot([b],[c],'or',/overplot)
  wait,0.02
endif

; initialize acceptance
accept = 0L

; generate random numbers for Metropolis
rand = randomu(seed,mc_iter,nwalkers)

; set values for prior
pp   = p0
psig = [10,2,2,5]

; generate variables for MCMC hammer
aham = 2.0
zmin = sqrt(1.0/aham)
zmax = sqrt(aham)
rham = randomu(seed,mc_iter,nwalkers)
zham = (rham*(zmax-zmin)+zmin)^(2.0)
half = fix(0.5*nwalkers)

; get initial posterior for all walkers
ln_posterior0 = dblarr(nwalkers)
for k = 0, nwalkers-1 do begin
  pk0 = reform(parr[k,*])
  yy0 = get_gaussian(pk0,x)
  ln_posterior0[k] = get_ln_posterior(pk0,pp,psig,y,yy0,ysig)
endfor
; ---------------

; begin simulation
for i = 1, mc_iter-1 do begin

  ; loop over sets of half-walkers
  for aw = 0, 1 do begin

    ; split walkers into half-groups
    if aw eq 0 then begin
      ; randomize walker pairs
      indx = permute(nwalkers)
      w1 = indx[0:half-1]
      w2 = indx[half:*]
    endif
    if aw eq 1 then begin
      w1 = indx[half:*]
      w2 = indx[0:half-1]
    endif

    ; loop over half-walkers
    for bw = 0, half-1 do begin

      ; get walker k from group 1
      k   = w1[bw]
      pk0 = reform(parr[k,*])

      ; get walker j from group 2
      j   = w2[bw]
      pj0 = reform(parr[j,*])

      ; generate new parameters from stretch move
      pk = pj0 + zham[i,k]*(pk0-pj0)

      ; get new log posterior for walker k
      yy  = get_gaussian(pk,x)
      ln_posterior  = get_ln_posterior(pk,pp,psig,y,yy,ysig)

      ; accept/reject move
      prob = ln_posterior-ln_posterior0[k]+(npar-1)*alog(zham[i,k])
      ; accept
      if prob gt 0 then begin
        parr[k,*] = pk
        ln_posterior0[k] = ln_posterior
        accept += 1
      endif else begin
        ; otherwise accept with some probability
        if prob gt alog(rand[i,k]) then begin
          parr[k,*] = pk
          ln_posterior0[k] = ln_posterior
          accept += 1
        endif else begin
          ; otherwise reject
        endelse
      endelse
    endfor
  endfor
  ; display parameter plot
  if pltdisplay eq 'yes' then begin
    data1 = plt1.data
    data2 = plt2.data
    data1[0,*] = parr[*,0]
    data1[1,*] = parr[*,3]
    data2[0,*] = parr[*,1]
    data2[1,*] = parr[*,2]
    plt1.putdata,data1
    plt2.putdata,data2
    wait,0.02
  endif
endfor
; -----------------------------------------
print, 'making plots ...'

; plot parameter contour
p = plot(parr[*,0],parr[*,3],'.',layout=[2,1,1],title='Background vs Peak',$
    dimensions=[600,300],/buffer)
p = plot([a],[d],'or',/overplot)
p = plot(parr[*,1],parr[*,2],'.',layout=[2,1,2],title='Width vs Center',$
    /current)
p = plot([b],[c],'or',/overplot)
p.save,'bayesian_test_contour_walkers.ps'

; plot parameter distributions
for i = 0, npar-1 do begin
  h = histogram(parr[*,i],nbins=100,locations=l)
  medparam = median(parr[*,i],/double)
  stdparam = stddev(parr[*,i],/double)
  if i eq 0 then begin
    p = plot(l,h,layout=[npar,2,i+1],dimensions=[npar*300,600],title='P1 Histogram',/buffer)
    ax = p.axes
    p = plot([medparam,medparam],[ax[0].yrange[0],ax[0].yrange[1]],'r',/overplot)
    p = plot(lindgen(nwalkers),parr[*,i],layout=[npar,2,i+1+npar],/current,title='P1 Values')
    ax = p.axes
    p = plot([ax[0].xrange[0],ax[1].xrange[1]],[medparam,medparam],'r',/overplot)
  endif else begin
    p = plot(l,h,layout=[npar,2,i+1],/current,title='P'+trim(i+1)+' Histogram')
    ax = p.axes
    p = plot([medparam,medparam],[ax[0].yrange[0],ax[0].yrange[1]],'r',/overplot)
    p = plot(lindgen(nwalkers),parr[*,i],layout=[npar,2,i+1+npar],/current,title='P'+trim(i+1)+' Values')
    ax = p.axes
    p = plot([ax[0].xrange[0],ax[1].xrange[1]],[medparam,medparam],'r',/overplot)
  endelse
endfor
p.save,'bayesian_test_parameters_walkers.ps'

; plot fit to data
medpar = median(parr,dimension=1,/double)
stdpar = stddev(parr,dimension=1,/double)
f = get_gaussian(medpar,x)
p = plot(x,y,xtitle='x',ytitle='y',title='MCMC Gaussian Fit Walkers',/buffer)
p = plot(x,f,'r',/overplot)
p.save,'bayesian_test_fit_walkers.ps'

; print stats
print,'mc iterations : '+trim(mc_iter)
print,'walkers       : '+trim(nwalkers)
acceptrate = (accept/float(nwalkers*(mc_iter-1)))*100
print,'acceptance rate: ',trim(acceptrate,'(i)'),+' %'
print,'input model parameters'
print,'peak       : ',a
print,'center     : ',b
print,'width      : ',c
print,'background : ',d
print,'best-fit model parameters'
print,'peak       : ',medpar[0],' +/- ',stdpar[0]
print,'center     : ',medpar[1],' +/- ',stdpar[1]
print,'width      : ',medpar[2],' +/- ',stdpar[2]
print,'background : ',medpar[3],' +/- ',stdpar[3]
print,'walker min/max on parameters'
print,'peak       : ',min(parr[*,0]),max(parr[*,0])
print,'center     : ',min(parr[*,1]),max(parr[*,1])
print,'width      : ',min(parr[*,2]),max(parr[*,2])
print,'background : ',min(parr[*,3]),max(parr[*,3])

toc
end
