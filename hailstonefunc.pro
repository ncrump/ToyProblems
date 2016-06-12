function hailstonefunc, n
; return the hailstone sequence for input integer

; check input is valid
if (n eq 0 or abs(n) ne n or ceil(n) ne n) then begin
  return, 'input must be a positive integer'
endif else begin
  ; compute sequence
  seq = [n]
  while (n ne 1) do begin
    if (n mod 2 eq 0) then begin
        n   = n/2
        seq = [seq, n]
    endif else begin
        n   = 3*n+1
        seq = [seq, n]
    endelse
  endwhile
endelse

print, 'length of sequence =',n_elements(seq)
print, 'hailstone sequence:'
return, seq
end
