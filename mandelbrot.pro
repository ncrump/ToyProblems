pro mandelbrot
; generate points in the mandelbrot set
; uses a simple slow method of directly evaluating the equation

; inputs
; ------------------
; domain size
xmin = -1.8
xmax =  0.8
ymin = -1.0
ymax =  1.0

; grid points
xpts = 800
ypts = 600

; max iterations
maxn = 30

; cutoff radius
rcut = 2
; ------------------

; iterate over mandelbrot eq
; z=z^2+c for complex c and z

; initialize
im   = complex(0,1)
z0   = complex(0,0)
magz = fltarr(xpts,ypts)+2*rcut
dx   = (xmax-xmin)/xpts
dy   = (ymax-ymin)/ypts

; evaluate grid points
for i = 0,xpts-1 do begin
  x = xmin + dx*i
  for j = 0,ypts-1 do begin
    y = ymin + dy*j
    c = x + im*y
    z = z0
    k = 0
    ; if abs(z) stays less than cutoff radius then point is in the set
    while (abs(z) lt rcut) and (k lt maxn) do begin
      z = z*z + c
      k += 1
    endwhile
  ; store magnitude of point for plotting
  if (k ge maxn) then begin
    magz[i,j] = z
  endif
  endfor
endfor

; plot image
device, decomposed=0                   ; set color mode
loadct, 1                              ; set color map
window, /free, xsize=xpts, ysize=ypts  ; set window size
tvscl, magz                            ; display scaled plot image
plt = tvrd(true=1)                     ; read display image
write_png, 'Plots/Mandelbrot.png',plt  ; save image

end
