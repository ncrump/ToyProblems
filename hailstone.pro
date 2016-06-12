pro hailstone
; count sequence length as a function of n for n <= 100

num = indgen(100)
cnt = intarr(100)
for i = 1, 100 do begin
  n = i
  k = 1
  ; compute sequence
  while (n ne 1) do begin
    if (n mod 2 eq 0) then begin
        n  = n/2
        k += 1
    endif else begin
        n  = 3*n+1
        k += 1
    endelse
  endwhile
  cnt[i-1] = k
endfor

; plot results
!p.background=255
!p.color=0
plot,num,cnt,charsize=2,xtitle='Number',ytitle='Hailstone Length'
plt = tvrd(true=1)
write_png, 'Plots/HailstoneLength.png', plt
end
