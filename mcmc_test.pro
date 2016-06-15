; gaussian fit model to be used with 'mpfit'
; ---------------------------------------
function mcmc_test_fit, p,x
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
function mcmc_test_likelihood, y,mu,sig
a = 1.0/(sig*sqrt(2*!pi))
b = (y-mu)
c = 2*sig*sig
likelihood = a*exp(-(b*b)/c)
return,likelihood>1e-30
end
; ---------------------------------------

; bayesian prior function
; ---------------------------------------
function mcmc_test_prior, p,pmin,pmax
  prior = fltarr(n_elements(p))+0.1
  prior[where(p lt pmin or p gt pmax)] = 1e-30
  return,prior
end
; ---------------------------------------


; MCMC test fit
pro mcmc_test

; inputs
; ---------------------------------------------
; set monte carlo iterations
mc_iter = 500000
; set burn-in iterations
mc_burn = 500000
; set limits for uniform prior distribution
y_min = 0
y_max = 2000
; set step size proposal distribution
y_stp = 10
; ---------------------------------------------

; generate sample data
; -----------------------------------------
; generate x data
xmin = -5.0
xmax =  5.0
npts = 10
x = ((xmax-xmin)/npts)*findgen(npts)+xmin

; set parameters for gaussian
a = 1000  ; peak
b = 0     ; center
c = 1.2   ; width
d = 200   ; background

; generate gaussian with poisson errors
y = mcmc_test_fit([a,b,c,d],x)
y = poidev(y)
e = sqrt(y)
; -----------------------------------------

; initialize variables
nsim = mc_iter+mc_burn
nx = n_elements(x)
accept  = fltarr(nsim)
mc_yval = fltarr(nsim,nx)

; initialize variables
mc_y0 = randomu(seed,nx)*(y_max-y_min)+y_min
likelihood0 = mcmc_test_likelihood(y,mc_y0,e)
ln_likelihood0 = total(alog(likelihood0))
prior0 = mcmc_test_prior(mc_y0,y_min,y_max)
ln_prior0 = total(alog(prior0))
ln_posterior0 = ln_likelihood0+ln_prior0

; store initial values
mc_yval[0,*] = mc_y0

; generate random numbers
prop = randomn(seed,nsim,nx)*y_stp
rand = randomu(seed,nsim)

; begin monte carlo simulation
for i = 1, nsim-1 do begin

  ; generate new y
  mc_y = mc_y0+prop[i,*]

  ; get the likelihood
  likelihood = mcmc_test_likelihood(y,mc_y,e)
  ln_likelihood = total(alog(likelihood))

  ; get the prior
  prior = mcmc_test_prior(mc_y,y_min,y_max)
  ln_prior = total(alog(prior))

  ; get the posterior
  ln_posterior = ln_likelihood+ln_prior

  ; accept/reject move
  prob = ln_posterior-ln_posterior0
  ; accept
  if prob gt 0 then begin
    mc_y0 = mc_y
    ln_posterior0 = ln_posterior
    mc_yval[i,*] = mc_y
    accept[i] = 1
  endif else begin
    ; otherwise accept with some probability
    if prob gt alog(rand[i]) then begin
      mc_y0 = mc_y
      ln_posterior0 = ln_posterior
      mc_yval[i,*] = mc_y
      accept[i] = 1
    endif else begin
    ; otherwise reject
      mc_yval[i,*] = mc_y0
      accept[i] = 0
    endelse
  endelse
endfor

; get mean y
avg_y= fltarr(nx)
for i = 0, nx-1 do begin
  avg_y[i] = mean(mc_yval[mc_burn:*,i])
endfor

; plot y
p = plot(x,y,dimensions=[800,400],/stairstep,layout=[2,1,1],title='Actual Data')
ax = p.axes
p = plot(x,avg_y,'r',yrange=[ax[0].yrange[0],ax[0].yrange[1]],$
    /stairstep,layout=[2,1,2],/current,title='MCMC Fit')
p.save,'mcmc_test.ps'

; print stats
print,'burn-in       : '+trim(mc_burn)
print,'mc iterations : '+trim(mc_iter)
acceptrate = (total(accept[mc_burn:*])/float(mc_iter))*100
print,'acceptance rate: ',trim(acceptrate,'(i)'),+' %'
print,'total error : ',mean(abs(y-avg_y)/y)

end
