pro worldseries
; monte carlo sim to estimate world series outcome of 7 games
; odds of yanks winning first 3 and socks winning last 4

; set number of runs
runs  = 1000000L

; generate random uniform numbers
rand  = randomu(seed,7*runs)
; set initial variables
cnt   = 0.0
ndx   = 0

; simulate 7 game series
for i = 1, runs do begin
  yanks = 0
  socks = 0
  for j = 1, 3 do begin
    r = rand[ndx]
    if (r le 0.5) then yanks += 1
    if (r gt 0.5) then socks += 1
    ndx   += 1
  endfor
  if (yanks eq 3) then begin
    for j = 1, 4 do begin
      r = rand[ndx]
      if (r le 0.5) then yanks += 1
      if (r gt 0.5) then socks += 1
      ndx   += 1
    endfor
    if (socks eq 4) then cnt += 1
  endif
endfor

print, '# wins:',cnt
print, '% odds:', (cnt/runs)*100

end
