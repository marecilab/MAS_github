;; $id$
;;

; ------------------------------------------------------------------------------
;
; ft_visualize.pro
;
; The procedure FT_VISUALIZE is the visualization tool for Fiber-Tract Mapping
; Program.  It checks if the file exists and in the right resolution format.
; Skips if the file was not found.  Resizes if it is in wrong resolution factor.
; And creates the .dat file if it is in the wrong format. Then it uses the IDL
; slicer3 software to visualize user-chosen files.
;
; SYNTAX
;	ft_visualize, refs, tractfile, rfx, rfy, rfz
;
; ARGUMENTS
;
;  refs
;	An integer array that has the value 1 if the corresponding file was
;	chosen for visualization and 0 otherwise.
;
;  tractfile
;	Name of the file that the user chose to store the fiber tracts data.
;
;  rfx, rfy, rfz
;	Resolution factors in each dimensions.
;
; Created on October 18, 2001 by E. Ozarslan
;
; ------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; function that does the checking, resizing and type changing (if necessary)
; and returns 1 for success and 0 for failure.

;function resize, name, p_name, rfx, rfy, rfz
;
;  ok=0
;  ;if (findfile(name+'.dat') eq '') then begin
;    if (findfile(name+'.flt') ne '') then begin
;      ;pro_create_3d, name+'.flt', name+'.dat'
;      ok=1
;    endif else begin
;      if (findfile(p_name) ne '') then begin
;	im=rflt(p_name)
;	sz_data=size(im)
;	fdim=sz_data(1)
;	pdim=sz_data(2)
;	sdim=sz_data(3)
;	im2=make_array(rfx*fdim,rfy*pdim,rfz*sdim)
;	if (strmid(p_name,strlen(p_name)-8,8,/reverse_offset) eq 'read.flt') $
;									then im=im(*,*,*,0)
;	im2=rebin(im, rfx*fdim, rfy*pdim, rfz*sdim)
;	openw, 55, name+'.flt', /SWAP_IF_LITTLE_ENDIAN
;	writeu, 55, size(im2), im2
;	free_lun, 55
;	;pro_create_3d, name+'.flt', name+'.dat'
;	ok=1
;      endif
;    endelse
;  ;endif else ok=1
;  if (ok eq 0) then print, 'Visualization file doesnt exist: ', p_name
;  return, ok
;end
;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; produces the name array for the files, and visualizes

pro ft_visualize, b
  compile_opt idl2
  common scan_data



; Adding support for reading from MAS

if strtrim(string((b.rfiles)[0]),2) eq 1 AND ptr_valid(project.dataArray[project.CI].adt) AND ptr_valid(project.dataArray[project.CI].eign_val) AND $
ptr_valid(project.dataArray[project.CI].eign_Vec) AND ptr_valid(project.dataArray[project.CI].frac_Ani) $
	AND ptr_valid(project.dataArray[project.CI].Avg_Dif) THEN BEGIN

p_names=[b.fname+'.trs']

rd_slicer, p_names, aspect_ratio=float(b.resz)/float(b.resx), $
  		rf=[b.rfx, b.rfx, b.rfx], length_limit=float(b.ll)/float(b.ss)


endif else begin

  no_files=5
  p_names=[b.fname+'.trs', 'Fract_Anisotropy.flt', $
		'dir_coherence.flt', 'S0.flt', 'read.flt']


  ;j=b.rfiles[0] ; Not used

  for i=0, no_files-1 do if b.rfiles[i] then vis_files=[vis_files, p_names[i]]

  if n_elements(vis_files) gt 1 then $
  rd_slicer, vis_files[1: n_elements(vis_files)-1], $
  		aspect_ratio=float(b.resz)/float(b.resx), $
  		rf=[b.rfx, b.rfx, b.rfx], length_limit=float(b.ll)/float(b.ss)

endelse

end




pro ft_visualize_newversion
    compile_opt idl2
    common scan_data
    common share_fibers_main
    
    ;; Adding support for reading from MAS

    if ptr_valid(project.dataArray[project.CI].adt)      AND $
       ptr_valid(project.dataArray[project.CI].eign_val) AND $
       ptr_valid(project.dataArray[project.CI].eign_Vec) AND $
       ptr_valid(project.dataArray[project.CI].frac_Ani) AND $
       ptr_valid(project.dataArray[project.CI].Avg_Dif)  THEN BEGIN

        ;;p_names=[fname+'.trs']
        p_names=[fname]
        
        rd_slicer, p_names, aspect_ratio=float(resz)/float(resx), $
          rf=[rfx, rfx, rfx], length_limit=float(ll)/float(ss)
        
        
    endif else begin
        
        no_files=5
        p_names=[b.fname+'.trs', 'Fract_Anisotropy.flt', $
                 'dir_coherence.flt', 'S0.flt', 'read.flt']
        

  ;j=b.rfiles[0] ; Not used

        for i=0, no_files-1 do if b.rfiles[i] then vis_files=[vis_files, p_names[i]]
        
        if n_elements(vis_files) gt 1 then $
          rd_slicer, vis_files[1: n_elements(vis_files)-1], $
          aspect_ratio=float(b.resz)/float(b.resx), $
          rf=[b.rfx, b.rfx, b.rfx], length_limit=float(b.ll)/float(b.ss)
        
    endelse
    
end
