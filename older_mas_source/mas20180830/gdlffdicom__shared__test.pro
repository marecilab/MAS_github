;L+
; LICENSE:
;
; IDL user contributed source code
; Copyright (C) 2006 Robbie Barnett
;
;    This library is free software;
;    you can redistribute it and/or modify it under the
;    terms of the GNU Lesser General Public License as published
;    by the Free Software Foundation; 
;    either version 2.1 of the License, 
;    or (at your option) any later version.
;
;    This library is distributed in the hope that it will
;    be useful, but WITHOUT ANY WARRANTY;
;    without even the implied warranty of MERCHANTABILITY
;    or FITNESS FOR A PARTICULAR PURPOSE. 
;    See the GNU Lesser General Public License for more details.
;
;    You should have received a copy of the GNU Lesser General Public License
;    along with this library; if not, write to the
;    Free Software Foundation, Inc., 59 Temple Place,
;    Suite 330, Boston, MA 02111-1307 USA 
;
; Please send queries to:
; Robbie Barnett
; Nuclear Medicine and Ultrasound
; Westmead Hospital
; +61 2 9845 7223
;L-


; The test data for this routine can be obtained from
; http://www.creatis.insa-lyon.fr/~jpr/PUBLIC/gdcm/gdcmData.tar.gz

pro gdlffdicom__shared__test0

files = file_search('*dcm',COUNT=nfiles)
nprocessed = 0
nidlffdicom = 0
du = obj_new('gdlffdicom__shared')
for i=0l,nfiles-1l do begin
    
    fstat = FILE_INFO(files[i])
 ;   help, fstat, /STRUCTURE
    if (du -> open( files[i], ACCESS_TIME=at, /NO_CATCH, MAX_BYTES=fstat.size)) then begin
        print, "Open " + files[i] + " access time", at
        images = du -> assoc(COUNT=count, ACCESS_TIME=at, TRUE=true, /NO_CATCH)
        help, images
        print, "Associate Acess time", at
        if (count gt 0) then begin
             nprocessed = nprocessed + 1
             window, nprocessed
              if (true) then begin
                 for j=0l,count-1l do  tv, images[*,*,j], TRUE=true
             endif else begin
                 for j=0l,count-1l do tvscl, images[*,*,j], j mod 4
             endelse
        endif else print, "Could not read " + files[i]
        result = du -> NewSOPInstanceUID()
        print, result
        du -> close

                                ; Warning - Causes segmentation fault for some invalid files
    endif
   
endfor
obj_destroy, du
print, "Number of files successfully processed by gdlffdicom__shared: ", nprocessed, ".  Total number of files ", nfiles, FORMAT="(A,I0,A,I0)"
print, "Number of files successfully opened by IDLffDICOM: ", nidlffdicom, ".  Total number of files ", nfiles, FORMAT="(A,I0,A,I0)"
end



pro gdlffdicom__shared__test1
files = file_search('*dcm',COUNT=nfiles)

for i=0l,nfiles-1l do begin
    du = obj_new('gdlffdicom__shared')
    print, "Opening " + files[i]
    fstat = FILE_INFO(files[i])

    if (du -> open( files[i],ACCESS_TIME=at,/INDEX_SEQUENCES, MAX_BYTES=fstat.size, /NO_CATCH)) then begin    
        print, "Dumping"
        du -> dump
        du -> close
    endif
endfor
obj_destroy, du

end


pro gdlffdicom__shared__test, inds, ALL=all

if (keyword_set(all)) then inds = indgen(1)
for i=0,n_elements(inds)-1 do begin
    call_procedure, STRING('gdlffdicom__shared__test',inds[i],FORMAT="(A,I0)")
endfor

end

