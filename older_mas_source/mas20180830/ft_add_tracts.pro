;; $Id$
;; This procedure doesn't seem to ever get called

pro ft_add_tracts, files, targetfile

  if (n_elements(files) eq 0) then begin
    ch=1 
    files=strarr(1)
    files[0] = Dialog_PickFile(/read, /must_exist, filter='*.flt',$
		path=Get_Path, Title='Read in 1st *.flt file')
    if (files[0] eq '') then return
    count=1
    while(ch) do begin
      files2=strarr(count+1)
      files2[0:count-1]=reform(files)
      if (count eq 1) then titleword='2nd' else if (count eq 2) then $
	titleword='3rd' else titleword=strtrim(string(count+1),1)+'th'
      files2[count] = Dialog_PickFile(/read, /must_exist,filter='*.flt',$
	path=Get_Path, Title='Read in '+titleword+' *.flt file')
      if (files2[count] eq '') then ch=0 else begin
        files=files2
        files2=0
      endelse
    count=count+1
    endwhile
  endif

  sz_files=size(files)
  num_files=sz_files(1)

  for i=0, num_files-1 do begin
    rtemp=rflt(files[i])
    if (i eq 0) then r=rtemp else r=r+rtemp
  endfor
  
  r = float(r ge 1.)

  if (n_elements(targetfile) eq 0) then begin
    targetfile=Dialog_PickFile(/read, filter='*.flt',$
	path=Get_Path, Title='Enter The Target File Name')
  endif
  openw, 55, targetfile, /SWAP_IF_LITTLE_ENDIAN
  writeu, 55, size(r), r
  free_lun, 55

end
