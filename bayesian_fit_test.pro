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
pro bayesian_fit_test

; inputs
; -----------------------------------------
; set monte carlo iterations
mc_iter = 20000L
; set burn-in iterations
mc_burn = 20000L
; set initial gaussian fit parameters
a0 = 220  ; peak
b0 = 1.0  ; center
c0 = 1.0  ; width
d0 = 20   ; background
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
; initialization
; ---------------
; get iterations
nsim = mc_iter+mc_burn

; store initial parameters
p0 = [a0,b0,c0,d0]
npar = n_elements(p0)
params = fltarr(nsim,npar)
avgpar = fltarr(npar)
stdpar = fltarr(npar)
params[0,*] = p0

; initialize acceptance
accept = fltarr(nsim)

; generate random numbers
prop_a = randomn(seed,nsim)*0.80
prop_b = randomn(seed,nsim)*0.01
prop_c = randomn(seed,nsim)*0.01
prop_d = randomn(seed,nsim)*0.80
rand   = randomu(seed,nsim)

; set values for prior
pp   = p0
psig = [10,2,2,5]

; get initial log posterior
yy0 = get_gaussian(p0,x)
ln_posterior0 = get_ln_posterior(p0,pp,psig,y,yy0,ysig)
; ---------------

; begin simulation
for i = 1, nsim-1 do begin

  ; generate new parameters from proposal function
  p = p0+[prop_a[i],prop_b[i],prop_c[i],prop_d[i]]

  ; get new log posterior
  yy = get_gaussian(p,x)
  ln_posterior = get_ln_posterior(p,pp,psig,y,yy,ysig)

  ; accept/reject move
  prob = ln_posterior-ln_posterior0
  ; accept
  if prob gt 0 then begin
    p0 = p
    ln_posterior0 = ln_posterior
    params[i,*] = p
    accept[i] = 1
  endif else begin
    ; otherwise accept with some probability
    if prob gt alog(rand[i]) then begin
      p0 = p
      ln_posterior0 = ln_posterior
      params[i,*] = p
      accept[i] = 1
    endif else begin
    ; otherwise reject
      params[i,*] = p0
      accept[i] = 0
    endelse
  endelse
endfor
; -----------------------------------------

; plot parameter distributions
for i = 0, npar-1 do begin
  h = histogram(params[mc_burn:*,i],nbins=100,locations=l)
  avgparam = mean(params[mc_burn:*,i])
  stdparam = stddev(params[mc_burn:*,i])
  avgpar[i] = avgparam
  stdpar[i] = stdparam
  if i eq 0 then begin
    p = plot(l,h,layout=[npar,2,i+1],dimensions=[npar*300,600],title='P1 Histogram',/buffer)
    ax = p.axes
    p = plot([avgparam,avgparam],[ax[0].yrange[0],ax[0].yrange[1]],'r',/overplot)
    p = plot(lindgen(mc_iter),params[mc_burn:*,i],layout=[npar,2,i+1+npar],/current,title='P1 Values')
    ax = p.axes
    p = plot([ax[0].xrange[0],ax[1].xrange[1]],[avgparam,avgparam],'r',/overplot)
  endif else begin
    p = plot(l,h,layout=[npar,2,i+1],/current,title='P'+trim(i+1)+' Histogram')
    ax = p.axes
    p = plot([avgparam,avgparam],[ax[0].yrange[0],ax[0].yrange[1]],'r',/overplot)
    p = plot(lindgen(mc_iter),params[mc_burn:*,i],layout=[npar,2,i+1+npar],/current,title='P'+trim(i+1)+' Values')
    ax = p.axes
    p = plot([ax[0].xrange[0],ax[1].xrange[1]],[avgparam,avgparam],'r',/overplot)
  endelse
endfor
p.save,'bayesian_test_parameters.ps'

; plot fit to data
f = get_gaussian(avgpar,x)
p = plot(x,y,xtitle='x',ytitle='y',title='MCMC Gaussian Fit',/buffer)
p = plot(x,f,'r',/overplot)
p.save,'bayesian_test_fit.ps'

; print stats
print,'burn-in       : '+trim(mc_burn)
print,'mc iterations : '+trim(mc_iter)
acceptrate = (total(accept[mc_burn:*])/float(mc_iter))*100
print,'acceptance rate: ',trim(acceptrate,'(i)'),+' %'
; print best-fit parameters
print,'input model parameters'
print,'peak       : ',a
print,'center     : ',b
print,'width      : ',c
print,'background : ',d
print,'best-fit model parameters'
print,'peak       : ',avgpar[0],' +/- ',stdpar[0]
print,'center     : ',avgpar[1],' +/- ',stdpar[1]
print,'width      : ',avgpar[2],' +/- ',stdpar[2]
print,'background : ',avgpar[3],' +/- ',stdpar[3]

end
