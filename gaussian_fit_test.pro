; gaussian fit model to be used with 'mpfit'
; ---------------------------------------
function gaussian_fit, p,x
a = p[0]  ; peak
b = p[1]  ; center
c = p[2]  ; width
d = p[3]  ; background
fit = a*exp(-((x-b)^2)/(2*c^2))+d
return,fit
end
; ---------------------------------------

; gaussian deviates to be used with 'mpfit'
; ---------------------------------------
function gaussian_dev, p,x=x,y=y,err=err
fit = gaussian_fit(p,x)
dev = (y-fit)/err
return,dev
end
; ---------------------------------------


; fit test using 'mpfit'
pro gaussian_fit_test

; generate x data
n = 101
x = -5+0.1*findgen(n)

; set parameters for gaussian
a = 200  ; peak
b = 0    ; center
c = 1.2  ; width
d = 30   ; background

; generate gaussian with poisson errors
y = gaussian_fit([a,b,c,d],x)
y = poidev(y)
e = sqrt(y)

; get least squares fit to gaussian
p0 = [210.0,1.,1.0,25.0]  ; initial guess at parameters
fa = {x:x,y:y,e:e}        ; input data and errors
; Note: mpfit calls two functions:
;     1) 'gaussian_dev' to generate deviates
;     2) 'gaussian_fit' to generate model fit
p = mpfit('gaussian_dev',p0,functargs=fa,maxiter=1000,$
           perror=perror,status=status,/quiet)

 ; catch convergence erorrs on fit
 if status lt 0 then print, ' ** failure to converge **'
 if status eq 5 then print, ' ** max iterations reached **'

; print best-fit parameters
print,'input model parameters'
print,'peak       : ',a
print,'center     : ',b
print,'width      : ',c
print,'background : ',d
print,'best-fit model parameters'
print,'peak       : ',p[0],' +/- ',perror[0]
print,'center     : ',p[1],' +/- ',perror[1]
print,'width      : ',p[2],' +/- ',perror[2]
print,'background : ',p[3],' +/- ',perror[3]

; plot data and fit
fit = gaussian_fit(p,x)
plt = plot(x,y,xtitle='x',ytitle='y',title='Least-Squares Gaussian Fit')
plt = plot(x,fit,'r',/overplot)

end
