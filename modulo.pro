pro modulo
; find sum of all multiples of 3 or 5 below 1000

s = long(0)
for i = 1, 999 do begin
  if (i mod 3 eq 0 or i mod 5 eq 0) then s = s+i
endfor

print, 'sum =', s
end
