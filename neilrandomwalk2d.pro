pro neilrandomwalk2d
; 2d random walk with variable step size
; modified with a small displacement towards the origin after each step

; set start parameters
; ----------------------
walkers = 1000
steps   = 1000
; ----------------------

; initialize walkers
x = fltarr(walkers)
y = fltarr(walkers)

; storage arrays
step = indgen(steps)+1
rmsd = []

; random walk
for i = 1, steps do begin
  ; generate random uniform steps between (-1,1)
  xstp = randomu(seed,walkers)*2 - 1
  ystp = randomu(seed,walkers)*2 - 1
  ; take a step
  x = x + xstp
  y = y + ystp
  ; take small step towards origin
  x = x - x*abs(xstp)*0.1
  y = y - y*abs(ystp)*0.1
  ; get rmsd vs step
  msd  = mean(x^2 + y^2)
  rmsd = [rmsd, sqrt(msd)]
endfor

; configure plot window
device, decomposed=0                   ; set color mode
!p.background=255                      ; set background to white
!p.color=0                             ; set plot to black
; plot rmsd vs t
plot,step,rmsd,charsize=2,$
     title ="Neil's Random Walk 2D",$
     xtitle='Steps',ytitle='RMSD'
plt = tvrd(true=1)
write_png, 'Plots/NeilRandomWalkRMSD.png', plt

end
