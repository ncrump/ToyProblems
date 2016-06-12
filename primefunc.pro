function primefunc, x
; function to check if input integer is prime

; check for non-ints and negatives
if (ceil(x) ne x or abs(x) ne x) then begin
    s = 'input must be a positive integer'
endif else begin
  ; check for 0,1,2
  if (x lt 2) then s = string(x)+' is not prime'
  if (x eq 2) then s = string(x)+' is prime'
  ; check for even integers > 2
  if (x gt 2 and x mod 2 eq 0) then s = string(x)+' is not prime'
  ; check for odd integers > 2
  if (x gt 2 and x mod 2 ne 0) then begin
    mx = ceil(sqrt(x))
    p  = 1
    for i = 2, mx do begin
      if (x mod i eq 0) then begin
        p = 0
        break
      endif
    endfor
    if (p eq 1) then s = string(x)+' is prime'
    if (p eq 0) then s = string(x)+' is not prime'
  endif
endelse

return, s

end
