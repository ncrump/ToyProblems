pro randomwalk2d
; 2d random walk with variable step size

; inputs
; ----------------------
; set sim size
walkers = 1000
steps   = 1000

; video option
video   = 'no'
fps     = 5
xpltlim = [-10,10]
ypltlim = [-10,10]
outfile = 'Plots/Animate/RandomWalk2d.mp4'
; ----------------------

; initialize walkers
x = fltarr(walkers)
y = fltarr(walkers)

; storage arrays
step = indgen(steps)+1
rmsd = []

; start video
if (video eq 'yes') then begin
  p = plot(x,y,'o',xrange=xpltlim,yrange=ypltlim,/buffer)
  t = text(0.76,0.9,'step'+string(0),/normal)
  write_video,outfile,p.copywindow(/antialias),handle=handle,$
              format='mp4',video_fps=fps
endif

; random walk
for i = 1, steps do begin
  ; generate random uniform steps between (-1,1)
  xstp = randomu(seed,walkers)*2 - 1
  ystp = randomu(seed,walkers)*2 - 1
  ; take a step
  x = x + xstp
  y = y + ystp
  ; get rmsd vs step
  msd  = mean(x^2 + y^2)
  rmsd = [rmsd, sqrt(msd)]
  ; save to video
  if (video eq 'yes') then begin
    p = plot(x,y,'o',xrange=xpltlim,yrange=ypltlim,/buffer)
    t = text(0.76,0.9,'step'+string(i),/normal)
    write_video,outfile,p.copywindow(/antialias),handle=handle
  endif
endfor

; end video
if (video eq 'yes') then begin
  write_video,/close,handle=handle
endif

; calculate diffusion coefficient
tlog = alog(step)
dlog = alog(rmsd)
fit  = poly_fit(tlog,dlog,1)
print, 'diffusion coefficient =',fit[1]

; configure plot window
device, decomposed=0                   ; set color mode
!p.background=255                      ; set background to white
!p.color=0                             ; set plot to black
; plot rmsd vs t
p1 = plot(step,rmsd,'-b',$
     title ="Random Walk 2D",$
     xtitle='Steps',ytitle='RMSD',$
     layout=[1,2,1])
p1.font_size=14
; plot log(rmsd) vs log(t)
p2 = plot(tlog,dlog,'-b',$
     xtitle='log(Steps)',ytitle='log(RMSD)',$
     /current,layout=[1,2,2])
p2.font_size=14
; plot best fit line to log-log plot
p3 = plot(tlog,fit[0]+tlog*fit[1],'-r',/overplot)

end
