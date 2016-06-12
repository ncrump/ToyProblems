pro tolstoyread
; analyze war&peace text file

; set file location
file = '../data/war&peace.txt'

; read file
text = rd_tfile(file, /compress)

; remove common punctuation
punc  = ['.',',','!','?',':',';','"','-','--']
npunc = n_elements(punc)
for i = 0, npunc-1 do text = str_replace(text, punc[i], ' ')

; split by words per line
text = strsplit(text, ' ', /extract)

; split by words
words = text.toarray(dimension=1)

; make all lower case
words = strlowcase(words)


; analyze word count
totwords = n_elements(words)
momentum = n_elements(where(words eq 'momentum'))


print, '   total word count:', totwords
print, 'momentum word count:', momentum

end
