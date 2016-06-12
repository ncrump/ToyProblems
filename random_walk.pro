pro random_walk

; inputs
; ----------------------
; set size
walkers = 1000
steps   = 1000
; show random walk
showwalk = 'no'
; ----------------------

; initialize
dr   = 0.2
wpos = fltarr(walkers,2)
step = indgen(steps)+1
rmsd = fltarr(steps)
thet = randomu(seed,walkers,steps)*2*!pi

if showwalk eq 'yes' then begin
  p = plot(wpos[*,0],wpos[*,1],'bo',xrange=[-10,10],yrange=[-10,10],title='1')
  wait,0.01
endif

; random walk
tic
for i = 0, steps-1 do begin
  msd = 0.0
  for j = 0, walkers-1 do begin
    ; take random step
    wpos[j,0] += dr*cos(thet[j,i])
    wpos[j,1] += dr*sin(thet[j,i])
    xij = wpos[j,0]
    yij = wpos[j,1]
    ; get displacement
    msd += xij*xij + yij*yij
  endfor
  rmsd[i] = sqrt(msd/float(walkers))
  if showwalk eq 'yes' then begin
    data = p.data
    data[0,*] = wpos[*,0]
    data[1,*] = wpos[*,1]
    p.putdata,data
    p.title = trim(i+1)
    wait,0.02
  endif
endfor
toc

; write output to file
openw,1,'random_walk_idl.txt'
for i = 0, steps-1 do printf,1,step[i],rmsd[i],format='(i,1x,f)'
close,1

end
