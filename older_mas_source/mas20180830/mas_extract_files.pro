;; $Id$
;; Copyright 2003 University of Florida. All Rights Reserved


; Subroutine name: get_dir_slash
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; this function will return the proper slashes for directories depending upon
; which operating system MAS is running on.
; slash[0] is the slash that goes in front of a directory
;           this is useful in unix b/c every thing starts at /
;           in windows this is blank
; slash[1] is the slash to show the difference between directories
;           in unix the slash is / in windows =  \

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

function get_dir_slash
    slash = strarr(2)
    if !d.name eq 'X' then begin
        slash[0] = '/'
        slash[1] = '/'
    endif else begin
        slash[0] = ''
        slash[1] = '\'
    end

    return, slash

end


; Subroutine name: is_file_32_bit
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;this function will return a 1 if the specified lun is stored in 32 bit signed
;integer else this function will return 0 if the specified lun is not stored
;in 32 bit signed integer

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

Function is_file_32_bit, lun, fdim, pdim, sdim, adim

    ;first get the file size in bytes
    sz_file = (FSTAT(lun)).size
    ;convert bytes to bits.
    sz_file *= long(8)

    ;now in order to calculate the proper file size we have to calculate the
    ;dimensions the proper dimensions of the file and if the file
    ;does not equal 32 bit signed integer size (long) then we
    ;will assume that the file is 16 bit signed integer (integer)
    calc_sz_file_32_bit = long(fdim) ;first set file size to fdim
    calc_sz_file_32_bit *= long(pdim) ;then multiply that by pdim
    calc_sz_file_32_bit *= long(sdim) ;then multiply that by sdim
    calc_sz_file_32_bit *= long(adim) ;then multiply that by adim
    calc_sz_file_32_bit *= long(2)    ;then multiply that by 2 for complex
    calc_sz_file_32_bit *= long(32)   ;then multiply that by 32 bits

    if sz_file eq calc_sz_file_32_bit then return, 1
    ;if the program makes it past this point then the file is not
    ;32 bit signed integer
    return, 0
end


; Subroutine name: rd_cmplx_blk_spec
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting


Function rd_cmplx_blk_spec, lun, fdim, word_size=word_size

    COMPILE_OPT IDL2
    
    data_blk = ComplexArr(fdim/2)
    if (word_size eq 32) then begin
        data_in = lonarr(2,fdim/2)
    endif else begin
        data_in = intarr(2,fdim/2)
    endelse
    ReadU, lun, data_in
    data_blk = Complex(Transpose(data_in[0,*]),Transpose(data_in[1,*]))
    Return, data_blk ; return complex data block

End


; Subroutine name: rd_cmplx_data_movie
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Read in single-slice movie 2D data then return the result
; parameters:
;   lun      [in]    file descriptor
;   fdim       [in]   Number of complex integers in freq. dimension
;   pdim       [in]   Number of phase-encode blocks
;   sdim       [in]   should be 1
;   adim        [in]  number of frames in the movie
;   rare_factor   [in]  RARE phase-encoding factor
;   islice       [in]  Vector of slice indices (lowest index is zero)
;   Returns    [out]   a block of complex data
; description:
;   Created (from original rd_cmplx_blk) by Ty Black
;   This function reads standard or RARE encoded data into complex format
;
; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

function rd_cmplx_data_movie,lun,fdim,pdim,sdim,adim,rare_factor,islice
    COMPILE_OPT IDL2

    if is_file_32_bit( lun, fdim, pdim, sdim, adim) then $
        read_file = LonArr(2,fdim, /NOZERO) $
    else $
        read_file = LonArr(2,fdim, /NOZERO)

    p_data_blk = ptr_new(ComplexArr(fdim,pdim,1,adim, /NOZERO))


    for pp=0, pdim-1 do begin
        for aa=0, adim-1 do begin
            ReadU, lun, read_file
            data_ts = Complex( Transpose(read_file[0,*]), Transpose(read_file[1,*]))
            (*p_data_blk)[*,pp,0,aa] = data_ts
        end
    end


    Return, p_data_blk ; return complex data block

end

; Subroutine name: rd_cmplx_data_DTI
; Created by: BT 2009-02
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Read in multiple-slice DTI (DTI_Tensor) 2D data then return the result
; parameters:
;   lun      [in]    file descriptor
;   fdim       [in]   Number of complex integers in freq. dimension
;   pdim       [in]   Number of phase-encode blocks
;   sdim       [in]   Number of slices
;   rare_factor   [in]  RARE phase-encoding factor
;   islice       [in]  Vector of slice indices (lowest index is zero)
;   Returns    [out]   a block of complex data
; description:
;   Created (from original rd_cmplx_blk) by Bill Triplett (*_T2 version)
;   This function reads standard or RARE encoded data into complex format

; Editing Information:

Function rd_cmplx_data_DTI,lun,fdim,pdim,sdim,adim,rare_factor,islice

    COMPILE_OPT IDL2

    read_file = LonArr(2,fdim, /NOZERO)
    
    p_data_blk = ptr_new(ComplexArr(fdim,pdim,sdim,adim, /NOZERO))
    
    a = fdim/8
    check_size = fdim MOD 128   ; Bruker data files are saved in
    check2 = fdim/128           ; blocks of 128; this helps to cut the
    check2 = fix(check2) + 1    ; empty points out
    subdim = pdim/rare_factor
    
                                ;RESP
                                ;open the data in these nested loops phase{slice{echo{freq}}}
    If check_size EQ 0 Then Begin ; if the data fits the blocks exactly

       for pp=0, subdim-1 do begin
          for ee=0, adim-1 do begin
             for ss=0, sdim-1 do begin
                for rare=0, rare_factor-1 do begin
                   ReadU, lun, read_file
                   data_ts = Complex( Transpose(read_file[0,*]), Transpose(read_file[1,*]))
                   (*p_data_blk)[*,pp+subdim*rare,islice[ss],ee] = data_ts
                                ;print, pp, islice[ss] ,ee
                endfor
             endfor
          endfor
       endfor

    endif else begin

       print, 'Data has to be padded b/c it fdim is not divisable by 128 evenly'
       read_file_padding = LonArr(2,check2*128 - fdim) ; this is the protion of the data that is gonig to be thown away
       for pp=0, subdim-1 do begin
          for ee=0, adim-1 do begin
             for ss=0, sdim-1 do begin
                for rare=0, rare_factor-1 do begin
                   ReadU, lun, read_file, read_file_padding
                   data_ts = Complex( Transpose(read_file[0,*]), Transpose(read_file[1,*]))
                   (*p_data_blk)[*,pp+subdim*rare,islice[ss],ee] = data_ts
                                ;print, pp, islice[ss] ,ee
                endfor
             endfor
          endfor
       endfor
       
    endelse
    
    Return, p_data_blk          ; return complex data block
    
End


; Subroutine name: rd_cmplx_data_T2
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Read in multiple-slice T2 2D data then return the result
; parameters:
;   lun      [in]    file descriptor
;   fdim       [in]   Number of complex integers in freq. dimension
;   pdim       [in]   Number of phase-encode blocks
;   sdim       [in]   Number of slices
;   rare_factor   [in]  RARE phase-encoding factor
;   islice       [in]  Vector of slice indices (lowest index is zero)
;   Returns    [out]   a block of complex data
; description:
;   Created (from original rd_cmplx_blk) by Ty Black
;   This function reads standard or RARE encoded data into complex format

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

Function rd_cmplx_data_T2,lun,fdim,pdim,sdim,adim,rare_factor,islice
COMPILE_OPT IDL2

read_file = LonArr(2,fdim, /NOZERO)


p_data_blk = ptr_new(ComplexArr(fdim,pdim,sdim,adim, /NOZERO))

a = fdim/8
check_size = fdim MOD 128       ; Bruker data files are saved in
check2 = fdim/128               ; blocks of 128; this helps to cut the
check2 = fix(check2) + 1        ; empty points out
subdim = pdim/rare_factor

                                ;RESP
                                ;open the data in these nested loops phase{slice{echo{freq}}}
If check_size EQ 0 Then Begin       ; if the data fits the blocks exactly
    for pp=0, subdim-1 do begin
        for ss=0, sdim-1 do begin
            for ee=0, adim-1 do begin
                for rare=0, rare_factor-1 do begin
                    ReadU, lun, read_file
                    data_ts = Complex( Transpose(read_file[0,*]), Transpose(read_file[1,*]))
                    (*p_data_blk)[*,pp+subdim*rare,islice[ss],ee] = data_ts
                                    ;print, pp, islice[ss] ,ee
                end


            end
        end
    end
end else begin
    print, 'Data has to be padded b/c it fdim is not divisable by 128 evenly'
    read_file_padding = LonArr(2,check2*128 - fdim)    ; this is the protion of the data that is gonig to be thown away
    for pp=0, subdim-1 do begin
        for ss=0, sdim-1 do begin
            for ee=0, adim-1 do begin
                for rare=0, rare_factor-1 do begin
                    ReadU, lun, read_file, read_file_padding
                    data_ts = Complex( Transpose(read_file[0,*]), Transpose(read_file[1,*]))
                    (*p_data_blk)[*,pp+subdim*rare,islice[ss],ee] = data_ts
                                    ;print, pp, islice[ss] ,ee
                end
            end
        end
    end

end



Return, p_data_blk              ; return complex data block

End



; Subroutine name: rd_cmplx_data_2d
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Read in multiple-slice 2D data then return the result
; parameters:
;   lun      [in]    file descriptor
;   fdim       [in]   Number of complex integers in freq. dimension
;   pdim       [in]   Number of phase-encode blocks
;   sdim       [in]   Number of slices
;   rare_factor   [in]  RARE phase-encoding factor
;   islice       [in]  Vector of slice indices (lowest index is zero)
;   Returns    [out]   a block of complex data
; description:
;   Created (from original rd_cmplx_blk) by T. H. Mareci, 1 August 2000
;   This function reads standard or RARE encoded data into complex format

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

Function rd_cmplx_data_2d,lun,fdim,pdim,sdim,adim,rare_factor,islice
COMPILE_OPT IDL2

    ;print,'rd_cmplx_data_2d arguments',fdim,pdim,sdim,adim,rare_factor
    ;;return, rd_cmplx_data_DTI(lun,fdim,pdim,sdim,adim,rare_factor,islice)

    if is_file_32_bit( lun, fdim, pdim, sdim, adim) then $
        data_in = LonArr(2,fdim, /NOZERO) $
    else $
        data_in = LonArr(2,fdim, /NOZERO)

    p_data_blk = ptr_new(ComplexArr(fdim,pdim,sdim,adim, /NOZERO))

    for oo = 0, adim-1 do begin
        ; Check data file for size in mod 128
        a = fdim/8
        check_size = fdim MOD 128           ; Bruker data files are saved in
        check2 = fdim/128                   ; blocks of 128; this helps to cut the
        check2 = fix(check2) + 1            ; empty points out

        subdim = pdim/rare_factor

        If check_size EQ 0 Then Begin       ; if the data fits the blocks exactly
            For i=0,subdim-1 Do Begin
               For j=0,sdim-1 Do Begin
                 For k=0,rare_factor-1 Do Begin

                  ReadU, lun, data_in               ; Read complex integer time series
                  ; Convert data_in from integer array to complex vector
                  data_ts = Complex( Transpose(data_in[0,*]), Transpose(data_in[1,*]))
                  (*p_data_blk)[*,i+subdim*k,islice[j],oo] = data_ts
                 EndFor
               EndFor
            EndFor
        EndIf Else Begin                     ; if data does not fit the blocks exactly
            For i=0,subdim-1 Do Begin
               For j=0,sdim-1 Do Begin
                 For k=0,rare_factor-1 Do Begin
                  data_in = LonArr(2,fdim)                  ; data_in(0,*) for real, data_in(1,*) for imag
                  data_in_1 = LonArr(2,check2*128 - fdim)    ; data_in_1(0,*) for real, data_in_1(1,*) for imag
                  ReadU, lun, data_in, data_in_1                   ; Read complex integer time series
                  ; Convert data_in from integer array to complex vector
                  data_ts = Complex( Transpose(data_in[0,*]), Transpose(data_in[1,*]))
                  (*p_data_blk)[*,i+subdim*k,islice[j],oo] = data_ts
                 EndFor
               EndFor
            EndFor
        EndElse

       ;print, oo

    end

Return, p_data_blk ; return complex data block

End


; Subroutine name: read_cmplx_data_3d
; Created by:  Created (from original rd_cmplx_blk) by T. H. Mareci, 1 August 2000
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Read in 3D data then return the result
;
; parameters:
;   lun      [in]    file descriptor
;   fdim       [in]   Number of complex integers in freq. dimension
;   pdim       [in]   Number of phase-encode blocks
;   sdim       [in]   Number of slices
;   rare_factor   [in]  RARE phase-encoding factor
;   islice       [in]  Vector of slice indices (lowest index is zero)
;   Returns    [out]   a block of complex data
; description:
; This function reads standard or RARE encoded data into complex format

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

Function read_cmplx_data_3d,lun,fdim,pdim,sdim,rare_factor,adim
COMPILE_OPT IDL2

;make a nice progress indicator
progressbar = Obj_New('progressbar', Color='red', Text='Extracting 3D Data from file',/NOCANCEL)
progressbar -> Start

if (n_elements(adim) eq 0) then adim = 1L

if is_file_32_bit( lun, fdim, pdim, sdim, 1) then $
        data_in = LonArr(2,fdim, /NOZERO) $ ; data_in(0,*) for real, data_in(1,*) for imaginary
    else $
        data_in = LonArr(2,fdim, /NOZERO)

data_blk = ptr_new(ComplexArr(fdim,pdim,sdim,adim, /NOZERO))

; Check data file for size in mod 128
a = fdim/8
check_size = fdim MOD 128      ; Bruker data files are saved in
check2 = fdim/128             ; blocks of 128; this helps to cut the
check2 = fix(check2) + 1         ; empty points out

subdim = pdim/rare_factor

If check_size EQ 0 Then Begin     ; if the data fits the blocks exactly
  
    For i=0,sdim-1 Do Begin
       For j=0,subdim-1 Do Begin
         For k=0,rare_factor-1 Do Begin
            for v =0,adim-1 do begin
          ReadU, lun, data_in          ; Read complex integer time series
          ; Convert data_in from integer array to complex vector
          data_ts = Complex( Transpose(data_in[0,*]), Transpose(data_in[1,*]))
          (*data_blk)[*,j+subdim*k,i,v] = data_ts
         EndFor
       EndFor
       progressBar -> Update, (float(i)/float(sdim-1))*100.0
    EndFor
  endfor
EndIf Else Begin              ; if data does not fit the blocks exactly

   For i=0,sdim-1 Do Begin
     For j=0,subdim-1 Do Begin
      For k=0,rare_factor-1 Do Begin

       for v=0, adim-1 do begin
  
          data_in_1 = LonArr(2,check2*128 - fdim)     ; data_in_1(0,*) for real, data_in_1(1,*) for imag
          ReadU, lun, data_in, data_in_1          ; Read complex integer time series
          ; Convert data_in from integer array to complex vector
          data_ts = Complex( Transpose(data_in[0,*]), Transpose(data_in[1,*]))
          (*data_blk)[*,j+subdim*k,i,v] = data_ts

       EndFor
     EndFor
     progressBar -> Update, (float(i)/float(sdim-1))*100.0
    EndFor
  endfor
EndElse

progressbar -> Destroy

Return, data_blk ; return complex data block

End


; Subroutine name: increment_dynamics_directory
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; parameters:
;   directoryToincrement  [in]     file pointer string to increment
;                   [out]  a file pointer string
; description:
;   this function takes in a string to a directory and adds 1 to the pre/pst number.
;
; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

function increment_dynamics_directory, directoryToincrement
    COMPILE_OPT IDL2


    if strpos(directoryToincrement,'fid') ne -1 then begin
       ;this is fid file pointer
       ;print, 'directoryToincrement: ', directoryToincrement

       placeToBreak = strpos(directoryToincrement,'pre')
       if placeToBreak eq -1 then $
         placeToBreak = strpos(directoryToincrement,'PRE')
       if placeToBreak eq -1 then $
         placeToBreak = strpos(directoryToincrement,'pst')
       if placeToBreak eq -1 then $
         placeToBreak = strpos(directoryToincrement,'PST')


       part1 = STRMID(directoryToincrement, 0, placeToBreak)
       part2 = STRMID(directoryToincrement, placeToBreak)
       prePstPart = STRMID(part2, 0,3)
       numAndFid = STRMID(part2, 3)
       fidpart = STRMID(numAndFid, strlen(numAndFid)-4)

       NumberLength = strpos(numAndFid,'fid')-1

       ;numberToincrement = fix(STRMID(numAndFid, NumberLength)) + 1
        numberToincrement = fix( STRMID( numAndFid, 0, NumberLength ) ) + 1

           ;numString = string(numberToincrement)
       numString = string(numberToincrement)
       numString = STRTRIM( numString , 2 )
       aatemp = STRJOIN( part1+prePstPart+numString+fidpart, /single)
       ;print, 'numberToincrement:' , numberToincrement
       ;print, aatemp
       return, aatemp
    endif else begin
       ;this is imnd file pointer
       ;print, 'directoryToincrement: ', directoryToincrement

       placeToBreak = strpos(directoryToincrement,'pre')
       if placeToBreak eq -1 then $
         placeToBreak = strpos(directoryToincrement,'PRE')
       if placeToBreak eq -1 then $
         placeToBreak = strpos(directoryToincrement,'pst')
       if placeToBreak eq -1 then $
         placeToBreak = strpos(directoryToincrement,'PST')


       part1 = STRMID(directoryToincrement, 0, placeToBreak)
       part2 = STRMID(directoryToincrement, placeToBreak)
       prePstPart = STRMID(part2, 0,3)
       ;print,'prePstPart ',prePstPart

       numAndimnd = STRMID(part2, 3)
       ;print, 'numAndimnd ',numAndimnd

       imndpart = STRMID(numAndimnd, strlen(numAndimnd)-5)
       ;print, 'imndpart ',imndpart

       NumberLength = strpos(numAndimnd,'imnd')-1
       ;print, 'NumberLength ',NumberLength

       numberToincrement = fix( STRMID( numAndimnd, 0, NumberLength ) ) + 1
       ;print,'numberToincrement ',numberToincrement

       numString = string(numberToincrement)
       numString = STRTRIM( numString , 2 )
       aatemp = STRJOIN( part1+prePstPart+numString+imndpart, /single)
       ;print, 'numberToincrement:' , numberToincrement
       ;print, aatemp
       return, aatemp

    end

end


; Subroutine name: replacePreWithPst
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; parameters:
;   filePointer   [in]  file pointer string to change
;          [out]  returns a string
; description:
;   This function takes in a string and changes the pre to pst in the string.
;
; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

function replacePreWithPst , filePointer
    COMPILE_OPT IDL2

    if strpos(filePointer,'fid') ne -1 then begin
       ;the file pointer contains fid
       ;this test if the windowing enviroment is a WIN system or a X windows system
       if STRCMP( !D.name, 'X') eq 1 then $
         placeToBreak = strpos(filePointer,'pre') $
       else $
         placeToBreak = strpos(filePointer,'PRE')

       part1 = STRMID(filePointer, 0, placeToBreak)
       part2 = STRMID(filePointer, placeToBreak)

       if STRCMP( !D.name, 'X') eq 1 then $
         prePstPart = 'pst' $
       else $
         prePstPart = 'PST'

       numAndFid = STRMID(part2, 3)
       fidpart = STRMID(numAndFid, strlen(numAndFid)-4)


       ;placeToBreak = strpos(directoryToincrement,'pre')
       ;if placeToBreak eq -1 then $
       ; placeToBreak = strpos(directoryToincrement,'PRE')
       ;if placeToBreak eq -1 then $
       ; placeToBreak = strpos(directoryToincrement,'pst')
       ;if placeToBreak eq -1 then $
       ; placeToBreak = strpos(directoryToincrement,'PST')
       ;strput, filePointer ,'pst',placeToBreak
       temp = STRJOIN( part1+prePstPart+'1'+fidpart, /single)
       ;print, temp
       return, temp

    endif else begin
       ;the file pointer does not contain fid so i assume that the oter option is
       ; that this fie contains imnd instead
       ;this test if the windowing enviroment is a WIN system or a X windows system
       if STRCMP( !D.name, 'X') eq 1 then $
         placeToBreak = strpos(filePointer,'pre') $
       else $
         placeToBreak = strpos(filePointer,'PRE')

       part1 = STRMID(filePointer, 0, placeToBreak)
       part2 = STRMID(filePointer, placeToBreak)

       if STRCMP( !D.name, 'X') eq 1 then $
         prePstPart = 'pst' $
       else $
         prePstPart = 'PST'

       numAndFid = STRMID(part2, 3)
       fidpart = STRMID(numAndFid, strlen(numAndFid)-5)


       ;placeToBreak = strpos(directoryToincrement,'pre')
       ;if placeToBreak eq -1 then $
       ; placeToBreak = strpos(directoryToincrement,'PRE')
       ;if placeToBreak eq -1 then $
       ; placeToBreak = strpos(directoryToincrement,'pst')
       ;if placeToBreak eq -1 then $
       ; placeToBreak = strpos(directoryToincrement,'PST')
       ;strput, filePointer ,'pst',placeToBreak
       temp = STRJOIN( part1+prePstPart+'1'+fidpart, /single)
       ;print, temp
       return, temp

    end
end


; Subroutine name: findNumberOrPrePstImages
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; parameters:
;   fid_file_pointer  [in] file pointer string to fid file
;   PRE           [in/out]
;   PST             [in/out]
; description:
;   takes in 3 arguments and returns 2. the number of pre and pst directories
;   is returned in the pre and pst values.
;
; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

Pro findNumberOrPrePstImages ,fid_file_pointer , PRE , PST
    COMPILE_OPT IDL2

    ;writes over the values of PRE and PST to represent the number of scans of each.

    path_in = fid_file_pointer
    done = 0
    ;numberOfPREImages = 0
    ;numberOfPSTImages = 0
    ;print,fid_file_pointer
    ;newDirectory = increment_dynamics_directory( fid_file_pointer )
    ;print,newDirectory

    ;data_in = ComplexArr(fdim,pdim,sdim,1)
    ;temp = ComplexArr(fdim,pdim,sdim,1)
    while done eq 0 do begin

       OpenR, lun, path_in, /Get_LUN, /SWAP_IF_LITTLE_ENDIAN, ERROR = err
       ; If err is nonzero, something happened. Print the error message to
       ; the standard error file (logical unit -2):
       IF (err NE 0) then  begin
         ;once there is an error you must have reached the end of the
         ; preimages
         done = 1
         CONTINUE
       end

       ;so if you can open another preimage then increase the number of preimages
       PRE = PRE + 1
       close, lun
       Free_LUN, lun
       path_in = increment_dynamics_directory ( path_in )

    endwhile

    path_in = replacePreWithPst ( path_in )

    done = 0
    while done eq 0 do begin

       OpenR, lun, path_in, /Get_LUN, /SWAP_IF_LITTLE_ENDIAN, ERROR = err
       ; If err is nonzero, something happened. Print the error message to
       ; the standard error file (logical unit -2):
       IF (err NE 0) then  begin
         ;once there is an error you must have reached the end of the
         ; pstimages
         done = 1
         CONTINUE
       end

       ;so if you can open another preimage then increase the number of preimages
       PST = PST + 1
       close, lun
       Free_LUN, lun
       path_in = increment_dynamics_directory ( path_in )

    endwhile

    ;print, PRE,' ',PST

end


; Subroutine name: extract_dynamics_fids
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; parameters:
;   fid_file_pointer  [in]
;
; description:
;   extracts fids across all the pre and pst folders in a dynamics scan.
;
; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

function extract_dynamics_fids , fid_file_pointer
    COMPILE_OPT IDL2
    COMMON scan_data

    ;get the proper slash depending upon which operation system is in use.
    slash = get_dir_slash()

    fdim = project.imndArray[project.ci].fdim
    pdim = project.imndArray[project.ci].pdim
    sdim = project.imndArray[project.ci].sdim
    adim = project.imndArray[project.ci].adim
    rare_factor = project.imndArray[project.ci].rare
    islice = *project.imndArray[project.ci].slice_scheme_ptr
    slice = project.imndArray[project.ci].slices
    ndims= project.imndArray[project.ci].dimensions
    navg = project.imndArray[project.ci].n_avg
    pdim_shift =  project.imndArray[project.ci].pdim_shift
    dim  =  project.imndArray[project.ci].dimensions
    multi_scan_file_array = *project.imndArray[project.ci].multi_scan_file_array


    n_avg = project.imndArray[project.ci].n_avg
    n_avg_flag =  project.procPramArray[project.ci].n_avg_flag

    ;this will help determine whether or not the imnd file contains certain
    ; words that will help decide whether or not to use little endian
    position = strpos( project.imndArray[project.ci].imnd_File ,'##ORIGIN= Bruker Medizintechnik GMBH' )
    ;print, position
    little_endian_flag = 1
    if position eq -1 then little_endian_flag = 0


    ;if the user selected to NOT divide the image by the number of averages then set averages to 1
    ;dividing by 1 doesn't change the image.
    if n_avg_flag  eq 0 then n_avg = 1


    data_in = ComplexArr(fdim,pdim,sdim,adim)


    FOR ii=0, adim -1 DO BEGIN

;      print, 'extracting'
;      print, string(multi_scan_file_array[ii]+slash[1]+'fid')
;      print, 'placing at position',string(ii)

        if (project.big_endian eq 1 and little_endian_flag eq 0) then begin
            OpenR, lun, multi_scan_file_array[ii]+slash[1]+'fid', /Get_LUN $
              , /SWAP_IF_BIG_ENDIAN ;, ERROR = err
        endif else begin
            OpenR, lun, multi_scan_file_array[ii]+slash[1]+'fid', /Get_LUN $
              , SWAP_IF_LITTLE_ENDIAN=little_endian_flag ; , ERROR = err
        endelse

;;       OpenR, lun, multi_scan_file_array[ii]+slash[1]+'fid', /Get_LUN, SWAP_IF_LITTLE_ENDIAN=little_endian_flag

       data_in[*,*,*,ii] = (*(rd_cmplx_data_2d (lun,fdim,pdim,sdim,1,rare_factor,islice)))/n_avg
       close, lun
       Free_LUN, lun


    END

    return, ptr_new(data_in,/no_copy)
end


; Subroutine name: extract_multiT2_fids
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; parameters:
;   fid_file_pointer  [in]
;
; description:
;   extracts fids across all the sets folders in a dynamics scan.
;
; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

function extract_multiT2_fids
    COMPILE_OPT IDL2

    COMMON scan_data

    fdim = project.imndArray[project.ci].fdim
    pdim = project.imndArray[project.ci].pdim
    sdim = project.imndArray[project.ci].sdim
    adim = project.imndArray[project.ci].adim
    rare_factor = project.imndArray[project.ci].rare
    islice = *project.imndArray[project.ci].slice_scheme_ptr
    slice = project.imndArray[project.ci].slices
    ndims= project.imndArray[project.ci].dimensions
    navg = project.imndArray[project.ci].n_avg
    pdim_shift =  project.imndArray[project.ci].pdim_shift
    dim  =  project.imndArray[project.ci].dimensions
    echoTime = *project.imndArray[project.ci].echo_Time_ptr
    parent_dir_file_listing = *project.imndArray[project.ci].multi_scan_file_array

    slash = get_dir_slash()

    n_avg = project.imndArray[project.ci].n_avg
    n_avg_flag =  project.procPramArray[project.ci].n_avg_flag

    ;this will help determine whether or not the imnd file contains certain
    ; words that will help decide whether or not to use little endian
    position = strpos( project.imndArray[project.ci].imnd_File ,'##ORIGIN= Bruker Medizintechnik GMBH' )
    ;print, position
    little_endian_flag = 1
    if position eq -1 then little_endian_flag = 0

    ;if the user selected to NOT divide the image by the number of averages then set averages to 1
    ;dividing by 1 doesn't change the image.
    if n_avg_flag  eq 0 then n_avg = 1


    data_in = ComplexArr(fdim,pdim,sdim,adim)

    FOR ii=0, adim -1 DO BEGIN
       path_in = parent_dir_file_listing[ii]+slash[1]+'fid'

       if (project.big_endian eq 1 and little_endian_flag eq 0) then begin
           OpenR, lun, path_in, /Get_LUN $
             , /SWAP_IF_BIG_ENDIAN ;, ERROR = err
       endif else begin
           OpenR, lun, path_in, /Get_LUN $
             , SWAP_IF_LITTLE_ENDIAN=little_endian_flag ; , ERROR = err
       endelse

;;       OpenR, lun, path_in, /Get_LUN, SWAP_IF_LITTLE_ENDIAN=little_endian_flag

       data_in[*,*,*,ii] = (*(rd_cmplx_data_2d(lun,fdim,pdim,sdim,1,rare_factor,islice)))/n_avg
       close, lun
       Free_LUN, lun

    END

    return, ptr_new(data_in,/no_copy)
end

; Subroutine name: extract_multiADT_fids
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

function extract_multiADT_fids
    COMPILE_OPT IDL2

    COMMON scan_data

    fdim                    =  project.imndArray[project.ci].fdim
    pdim                    =  project.imndArray[project.ci].pdim
    sdim                    =  project.imndArray[project.ci].sdim
    adim                    =  project.imndArray[project.ci].adim
    rare_factor             =  project.imndArray[project.ci].rare
    islice                  = *project.imndArray[project.ci].slice_scheme_ptr
    slice                   =  project.imndArray[project.ci].slices
    ndims                   =  project.imndArray[project.ci].dimensions
    pdim_shift              =  project.imndArray[project.ci].pdim_shift
    dim                     =  project.imndArray[project.ci].dimensions
    n_bvals                 =  project.imndArray[project.ci].n_bvals
    parent_file_pointer     =  project.imndArray[project.ci].file_Path
    parent_dir_file_listing = *project.imndArray[project.ci].multi_scan_file_array

    ;; this will help determine whether or not the imnd file contains certain
    ;; words that will help decide whether or not to use little endian
    position = strpos( project.imndArray[project.ci].imnd_File ,'##ORIGIN= Bruker Medizintechnik GMBH' )
    little_endian_flag = 1
    if position eq -1 then little_endian_flag = 0

    data_in    = ComplexArr(fdim,pdim,sdim,adim)
    n_scans    = (size(parent_dir_file_listing))[1]
    n_avg      = *(project.imndArray[project.ci].pn_avg)
    n_avg_flag =  project.procPramArray[project.ci].n_avg_flag

    vols_per_file = adim/n_scans
    ;;if the user selected to NOT divide the image by the number of averages then set averages to 1
    ;;dividing by 1 doesn't change the image.
    if n_avg_flag  eq 0 then begin
       n_avg[*] = 1
    end

    progressbar = Obj_New('progressbar', Color='red', Text='Extracting DWI Data from files', /NOCANCEL)
    progressbar -> Start
    index1 = 0

    FOR ii=0, n_scans-1 DO BEGIN

        n_bvals = (*project.imndArray[project.ci].file_adim)[ii]

        update_status_bar, string ('Opening scan ',strtrim(ii,2),'/', strtrim(n_scans-1,2))

        print, 'Opening ',parent_dir_file_listing[ii]

        open_file  = parent_dir_file_listing[ii]+'\fid'
        if !d.name eq 'X' then begin
            open_file  = parent_dir_file_listing[ii]+'/fid'
        endif

;;        OpenR, lun, open_file , /Get_LUN,
;;        SWAP_IF_LITTLE_ENDIAN=little_endian_flag

        if (project.big_endian eq 1 and little_endian_flag eq 0) then begin
            OpenR, lun, open_file, /Get_LUN $
              , /SWAP_IF_BIG_ENDIAN;, ERROR = err
        endif else begin
            OpenR, lun, open_file, /Get_LUN $
              , SWAP_IF_LITTLE_ENDIAN=little_endian_flag; , ERROR = err
        endelse

        if (vols_per_file eq 1) then begin
            index2 = index1+ n_bvals-1
            ;;print, 'writing fid between '+string( index1) +'to '+string(index2)
            ;;store the data at this location and divided by the number of averages for this scan.
            data_in[*,*,*,index1:index2] = (*(rd_cmplx_data_2d(lun,fdim,pdim,sdim,n_bvals,rare_factor,islice)))/n_avg[ii]
            print, 'Dividing by', n_avg[ii]
            close, lun
            Free_LUN, lun
            index1 += n_bvals
        endif else begin
            index2 = index1+ vols_per_file-1
            ;;print, 'writing fid between '+string( index1) +'to '+string(index2)
            ;;store the data at this location and divided by the number of averages for this scan.
            data_in[*,*,*,index1:index2] = (*(rd_cmplx_data_DTI(lun,fdim,pdim,sdim,vols_per_file,rare_factor,islice)))/n_avg[ii]
            print, 'Dividing by', n_avg[ii]
            close, lun
            Free_LUN, lun
            index1 += vols_per_file
        endelse
        progressBar -> Update, (float(ii)/float(n_scans-1))*100.0

    ENDfor

    progressbar -> Destroy

    return, ptr_new(data_in, /no_copy)

end


; Subroutine name: extract_multi_phase_array_fids
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/05.
    ;Fix spelling mistakes and commenting

function extract_multi_phase_array_fids
    COMPILE_OPT IDL2

    COMMON scan_data

    fdim = project.imndArray[project.ci].fdim
    pdim = project.imndArray[project.ci].pdim
    sdim = project.imndArray[project.ci].sdim
    adim = project.imndArray[project.ci].adim
    rare_factor = project.imndArray[project.ci].rare
    islice = *project.imndArray[project.ci].slice_scheme_ptr
    slice = project.imndArray[project.ci].slices
    ndims= project.imndArray[project.ci].dimensions
    navg = project.imndArray[project.ci].n_avg
    pdim_shift =  project.imndArray[project.ci].pdim_shift
    dim  =  project.imndArray[project.ci].dimensions
    n_bvals = project.imndArray[project.ci].n_bvals
    parent_file_pointer = project.imndArray[project.ci].file_Path
    parent_dir_file_listing = *project.imndArray[project.ci].multi_scan_file_array


    n_avg = project.imndArray[project.ci].n_avg
    n_avg_flag =  project.procPramArray[project.ci].n_avg_flag


    ;this will help determine whether or not the imnd file contains certain
    ; words that will help decide whether or not to use little endian
    position = strpos( project.imndArray[project.ci].imnd_File ,'##ORIGIN= Bruker Medizintechnik GMBH' )
    ;print, position
    little_endian_flag = 1
    if position eq -1 then little_endian_flag = 0

    ;if the user selected to NOT divide the image by the number of averages then set averages to 1
    ;dividing by 1 doesn't change the image.
    if n_avg_flag  eq 0 then begin
       n_avg = 1
    end


    n_scans = (size(parent_dir_file_listing))[1]

    data_in = ComplexArr(fdim,pdim,sdim,n_scans)

    progressbar = Obj_New('progressbar', Color='red', Text='Extracting Array Data',/NOCANCEL)
    progressbar -> Start

    ;*************************************************************************************************

    image1 = ComplexArr(fdim,pdim,sdim,n_scans) ;Creating a complex array with four parameters
    image2 = ComplexArr(fdim,pdim,sdim,n_scans) ;Creating a complex array with four parameters
    image3 = ComplexArr(fdim,pdim,sdim,n_scans) ;Creating a complex array with four parameters
    image4 = ComplexArr(fdim,pdim,sdim,n_scans) ;Creating a complex array with four parameters



       open_file = parent_dir_file_listing[0]+'\fid'  ;Seting open_file to the first folder with a fid file

    if !d.name eq 'X' then open_file = parent_dir_file_listing[0]+'/fid'   ;Checking

;    OpenR, lun, open_file , /Get_LUN, SWAP_IF_LITTLE_ENDIAN=little_endian_flag  ;Opening the file

    if (project.big_endian eq 1 and little_endian_flag eq 0) then begin
        OpenR, lun, open_file, /Get_LUN $
          , /SWAP_IF_BIG_ENDIAN ;, ERROR = err
    endif else begin
        OpenR, lun, open_file, /Get_LUN $
          , SWAP_IF_LITTLE_ENDIAN=little_endian_flag ; , ERROR = err
    endelse

    data_in[*,*,*,0] = (*(rd_cmplx_data_2d(lun,fdim,pdim,sdim,1,rare_factor,islice)))/n_avg  ;Use the function rd_cmplx_data_2d
                                                                                             ;to read the data and then transfer
                                                                                             ;it to data_in folder 0 in the fourth parameter
       image1[*,*,*,0] = data_in[*,*,*,0]   ;Transfer the data in data_in to image1


    close, lun     ;Closing the lun
    Free_LUN, lun  ;Deallocates previously-allocated file units

progressBar -> Update, 25

       open_file = parent_dir_file_listing[1]+'\fid'  ;Seting open_file to the second folder with a fid file

    if !d.name eq 'X' then open_file = parent_dir_file_listing[1]+'/fid'   ;checking
    OpenR, lun, open_file , /Get_LUN, SWAP_IF_LITTLE_ENDIAN=little_endian_flag  ;opening the file


    data_in[*,*,*,1] = (*(rd_cmplx_data_2d(lun,fdim,pdim,sdim,1,rare_factor,islice)))/n_avg  ;Use the function rd_cmplx_data_2d
                                                                                             ;to read the data and then transfer
                                                                                             ;it to data_in folder 1 in the fourth parameter

      image2[*,*,*,1] = data_in[*,*,*,1]  ;Transfer the data in data_in to image2

    close, lun     ;Closing the lun
    Free_LUN, lun  ;Deallocates previously-allocated file units

progressBar -> Update, 50

       open_file = parent_dir_file_listing[2]+'\fid'  ;Seting open_file to the third folder with a fid file

    if !d.name eq 'X' then open_file = parent_dir_file_listing[2]+'/fid'   ;checking
    OpenR, lun, open_file , /Get_LUN, SWAP_IF_LITTLE_ENDIAN=little_endian_flag  ;opening the file

    data_in[*,*,*,2] = (*(rd_cmplx_data_2d(lun,fdim,pdim,sdim,1,rare_factor,islice)))/n_avg  ;Use the function rd_cmplx_data_2d
                                                                                             ;to read the data and then transfer
                                                                                             ;it to data_in folder 2 in the fourth parameter

       image3[*,*,*,2] = data_in[*,*,*,2]  ;Transfer the data in data_in to image3

    close, lun     ;Closing the lun
    Free_LUN, lun  ;Deallocates previously-allocated file units


progressBar -> Update, 75

       open_file = parent_dir_file_listing[3]+'\fid'  ;Seting open_file to the fourth folder with a fid file

    if !d.name eq 'X' then open_file = parent_dir_file_listing[3]+'/fid'   ;checking
    OpenR, lun, open_file , /Get_LUN, SWAP_IF_LITTLE_ENDIAN=little_endian_flag  ;opening the file


    data_in[*,*,*,3] = (*(rd_cmplx_data_2d(lun,fdim,pdim,sdim,1,rare_factor,islice)))/n_avg  ;Use the function rd_cmplx_data_2d
                                                                                             ;to read the data and then transfer
                                                                                             ;it to data_in folder 3 in the fourth parameter

      image4[*,*,*,3] = data_in[*,*,*,3]  ;Transfer the data in data_in to image4

   close, lun     ;Closing the lun
   Free_LUN, lun  ;Deallocates previously-allocated file units

progressBar -> Update, 100

;*********************************************************************************************

;    FOR ii=0, n_scans-1 DO BEGIN
;
;        progressBar -> Update, (float(ii)/float(n_scans-1))*100.0
;
;        open_file  = parent_dir_file_listing[ii]+'\fid'
;        if !d.name eq 'X' then open_file  = parent_dir_file_listing[ii/n_bvals]+'/fid'
;
;        OpenR, lun, open_file , /Get_LUN, SWAP_IF_LITTLE_ENDIAN=little_endian_flag
;
;        data_in[*,*,*,ii] = (*(rd_cmplx_data_2d(lun,fdim,pdim,sdim,1,rare_factor,islice)))/n_avg
;       close, lun
;        Free_LUN, lun
;    ENDfor


    progressbar -> Destroy

    return, ptr_new(data_in,/no_copy)



end


; Subroutine name: extract_multi_slice_fids
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/05.
    ;Fix spelling mistakes and commenting


function extract_multi_slice_fids

    COMPILE_OPT IDL2

    COMMON scan_data

    fdim = project.imndArray[project.ci].fdim
    pdim = project.imndArray[project.ci].pdim
    sdim = project.imndArray[project.ci].sdim
    adim = project.imndArray[project.ci].adim
    rare_factor = project.imndArray[project.ci].rare
    islice = *project.imndArray[project.ci].slice_scheme_ptr
    slice = project.imndArray[project.ci].slices
    ndims= project.imndArray[project.ci].dimensions
    navg = project.imndArray[project.ci].n_avg
    pdim_shift =  project.imndArray[project.ci].pdim_shift
    dim  =  project.imndArray[project.ci].dimensions
    n_bvals = project.imndArray[project.ci].n_bvals
    parent_file_pointer = project.imndArray[project.ci].file_Path
    parent_dir_file_listing = *project.imndArray[project.ci].multi_scan_file_array


    n_avg = project.imndArray[project.ci].n_avg
    n_avg_flag =  project.procPramArray[project.ci].n_avg_flag


    ;this will help determine whether or not the imnd file contains certain
    ; words that will help decide whether or not to use little endian
    position = strpos( project.imndArray[project.ci].imnd_File ,'##ORIGIN= Bruker Medizintechnik GMBH' )
    ;print, position
    little_endian_flag = 1
    if position eq -1 then little_endian_flag = 0

    ;if the user selected to NOT divide the image by the number of averages then set averages to 1
    ;dividing by 1 doesn't change the image.
    if n_avg_flag  eq 0 then n_avg = 1

    data_in = ComplexArr(fdim,pdim,sdim,adim)

    progressbar = Obj_New('progressbar', Color='red', Text='Extracting movie Data',/NOCANCEL)
    progressbar -> Start

    FOR ii=0, sdim-1 DO BEGIN

        progressBar -> Update, (float(ii)/float(sdim-1))*100.0

        open_file  = parent_dir_file_listing[ii]+'\fid'
        if !d.name eq 'X' then open_file  = parent_dir_file_listing[ii/n_bvals]+'/fid'

        if (project.big_endian eq 1 and little_endian_flag eq 0) then begin
            OpenR, lun, open_file, /Get_LUN $
              , /SWAP_IF_BIG_ENDIAN;, ERROR = err
        endif else begin
            OpenR, lun, open_file, /Get_LUN $
              , SWAP_IF_LITTLE_ENDIAN=little_endian_flag; , ERROR = err
        endelse

 ;;       OpenR, lun, open_file , /Get_LUN, SWAP_IF_LITTLE_ENDIAN=little_endian_flag

        extracted_data = (*(rd_cmplx_data_movie(lun,fdim,pdim,sdim,adim,rare_factor,islice)))/n_avg
        data_in[*,*,ii,*] = extracted_data
        close, lun
        Free_LUN, lun
    ENDfor



    progressbar -> Destroy

    return, ptr_new(data_in,/no_copy)


end

; Subroutine name: extract_ADT_DICOM
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; This function will extract data from the current indexes multi_scan_file_array and return a pointer
; to that data.
;
; Editing Information:
    ;Edited by HS 2006/10/05.
    ;Fix spelling mistakes and commenting

function extract_ADT_DICOM
    COMPILE_OPT IDL2
    COMMON scan_data
    CI = project.ci
    fdim = project.imndArray[ci].fdim
    pdim = project.imndArray[ci].pdim
    sdim = project.imndArray[ci].sdim
    adim = project.imndArray[ci].adim
    ;in dealing with dimcom images multi_scan_file_array contains the directory for adt images and b_0 only 2 entries.
    multi_scan_file_array = *project.imndArray[ci].multi_scan_file_array


    temp = multi_scan_file_array[0]+'*.dcm'
    adtdicomFileListing = FILE_SEARCH(temp, count=n_adtDicomImages)

    b_value_length = fix(n_adtDicomImages/project.imndArray[ci].n_bvals)

    ;now we can open the data and extract the images and
    ;coalate them into their proper positions.
    p_data = ptr_new(fltArr(fdim, pdim, sdim, adim))


    update_status_bar, 'opening the data'
    progressbar = Obj_New('progressbar', Color='red', Text='Extracting DICOM data from file',/NOCANCEL)
    progressbar -> Start

    for ii = 0 , n_adtDicomImages-1 do  begin

        dicom = get_dicom_header_and_image( adtdicomFileListing[ii] )

        ;this tells is the image nuber is >= b_value_length
        odd_shift =(dicom.image_Number/b_value_length)

        sdim_index   = dicom.image_Number mod sdim
        adim_index    = (((dicom.image_Number - b_value_length * odd_shift)/sdim)*2 + odd_shift)+1

        (*p_data)[*,*,sdim_index,adim_index]= dicom.image

        progressBar -> Update, (float(ii)/float(n_adtDicomImages-1))*100.0

    endfor
    progressbar -> Destroy

    update_status_bar,''


    ;open the b_zero image.

    ;what are the b_zero file names and how many are there
    b_zeroFileListing = FILE_SEARCH(multi_scan_file_array[1]+'*.dcm', count=num_b_zero)

    progressbar = Obj_New('progressbar', Color='red', Text='Extracting DICOM b0 Data from file',/NOCANCEL)
    progressbar -> Start

    ;open all the b_zero files and put them in there right position.
    for ii = 0 , num_b_zero-1 do  begin

        dicom = get_dicom_header_and_image( b_zeroFileListing[ii] )
        sdim_index   = dicom.image_Number mod sdim

        (*p_data)[*,*,sdim_index,0] = dicom.image

        progressBar -> Update, (float(ii)/float(num_b_zero-1))*100.0
    endfor
    progressbar -> Destroy


    return ,p_data
end


;*****************************************************************************************************
;
; NAME:
;   mas_extract_data
;
; PURPOSE:
;   this procedure will extract data directly from the files that a scan was created with.
;
; ARGUMENTS:
;
;
; MODIFICATION HISTORY:
;
; Copyright 2003 University of Florida. All Rights Reserved
; Edited by HS 2006/10/05.
; Fix spelling mistakes and commenting
; HS - 20061028
; Modifying the code to be able to open .raw files in the load state 1
;
; Changed reader for multi-echo 3D data. (Magdoom, 1/13/16)
;*****************************************************************************************************
function mas_extract_data

    COMPILE_OPT IDL2
    COMMON scan_data
    
    CI = project.ci
    fdim = project.imndArray[CI].fdim
    pdim = project.imndArray[CI].pdim
    sdim = project.imndArray[CI].sdim
    adim = project.imndArray[CI].adim
    rare_factor = project.imndArray[CI].rare
    
    slice = project.imndArray[CI].slices
    ndims = project.imndArray[CI].dimensions
    navg =  project.imndArray[CI].n_avg
    dim  =  project.imndArray[CI].dimensions
    
    ; *** Part to read the .raw files.
    path_in = project.imndArray[ci].file_Path
    position = strpos(path_in,'.raw')
    
    if position ne -1 then begin
    
        ; Open the file to get the LUN
        if (project.big_endian eq 1 ) then begin
            OpenR, lun, path_in, /Get_LUN, /SWAP_IF_BIG_ENDIAN
        endif else begin
            OpenR, lun, path_in, /Get_LUN, /SWAP_IF_LITTLE_ENDIAN
        endelse
        
        ; The dimensionality of the data is contained in the first byte of the raw file.
        dimensionality = LonArr(1)
        
        ; Read in the byte number defining the dimensionality of the data
        ReadU, lun, dimensionality
        
        ; Set up the array to determine the dimensionality of the data (up to 5)
        if dimensionality eq 3 then begin
            header_raw = LonArr(6)
            header_temp = LonArr(5)
        endif else begin
            header_raw = LonArr(7)
            header_temp = LonArr(6)
        endelse
        
        header_raw[0] = dimensionality
        
        ; Read in the rest of the header
        ReadU, lun, header_temp
        
        header_raw[1:*] = header_temp[*]
        
        ; Create the complex data array depending in the dimensions
        if  dimensionality eq 3 then begin
            complex_data_in = ComplexArr(header_raw[1],header_raw[2], header_raw[3])
        endif else begin
            complex_data_in = ComplexArr(header_raw[1],header_raw[2], header_raw[3], header_raw[4])
        endelse
        
        ; Finish reading the data.
        ReadU, lun, complex_data_in
        free_lun, lun
        return, ptr_new(complex_data_in)
        
    endif else begin
    
        ; *** This is the original code.
        path_in = project.imndArray[CI].file_Path+'fid'
        if (not file_test(path_in, /read)) then begin
            path_in = project.imndArray[CI].file_Path+'ser'
            if (not file_test(path_in, /read)) then begin
                ;;void = dialog_message(['Unable to read fid or ser file from:', $
                ;;                       path_in], /error, /center)
                ;;return, ptr_new()
            endif
        endif
                                      
        image_type = project.imndArray[CI].image_type
        
        print, 'Image type', image_type
        
        update_status_bar,' Extracting fid data from file(s)'
        
        ;this will help determine whether or not the imnd
        ; file contains certain words that will help decide
        ; whether or not to use little endian
        position = strpos( project.imndArray[CI].imnd_File, '##ORIGIN= Bruker Medizintechnik GMBH' )
        ;print, position
        little_endian_flag = 1
        if position eq -1 then little_endian_flag = 0
        
        ;==========================================================================
        ;==========================================================================
        ;================================== IS THIS REALLY NECESSARY???? ==========
        ;==========================================================================
        ;==========================================================================
        
        ;    ;check to see if there is a fid file or a ser file.
        ;    OpenR, lun, path_in, /Get_LUN , ERROR = err
        ;    IF (err NE 0) then  begin
        ;        path_in = project.imndArray[CI].file_Path+'ser'
        ;        print, 'Opening the ser file.'
        ;        close, /all
        ;        OpenR, lun, path_in, /Get_LUN, ERROR = err
        ;        if (err NE 0) then print, 'There is an error opening the ser and fid file'
        ;    end
        ;    close, /all
        
        ;free_lun, lun ; <- This part stays commented
        
        
        ;if the image_type is 2d or 3d
        if image_type eq 0 or image_type eq 9 then begin
        
            islice = *project.imndArray[CI].slice_scheme_ptr
            n_avg = project.imndArray[project.ci].n_avg
            n_avg_flag =  project.procPramArray[project.ci].n_avg_flag

            ;open the file unit and prepare for
            ;opening
            if (project.big_endian eq 1 and little_endian_flag eq 0) then begin
                OpenR, lun, path_in, /Get_LUN, /SWAP_IF_BIG_ENDIAN, ERROR = err
            endif else begin
                OpenR, lun, path_in, /Get_LUN, SWAP_IF_LITTLE_ENDIAN=little_endian_flag, ERROR = err
            endelse
            
            if (err NE 0) then print, 'There is an error opening the data file.'
            
            ;print, 'opening ',path_in
            ; print, 'little endian flag', little_endian_flag
            ; print, 'dimensions', dim
            ;if the data is 2d then use the 2d open
            if dim eq 2 then begin
                if (strpos(project.imndarray[ci].imnd_file, 'DW_Contrast') ne -1 or $
                    strpos(project.imndarray[ci].imnd_file, 'DW_Tensor') ne -1) then begin
                    
                    p_data =  rd_cmplx_data_DTI(lun,fdim,pdim,sdim,adim,rare_factor,islice)
                    
                endif else begin
                
                    p_data =  rd_cmplx_data_2d(lun,fdim,pdim,sdim,adim,rare_factor,islice)
                    
                endelse
            endif
            
            ;if the data is 3d then use the 3d open
            if dim eq 3 then p_data =  read_cmplx_data_3d(lun,fdim,pdim,sdim,rare_factor,adim)
            print, size(p_data)
            ;if the user selected to divide the image by the number of averages then
            ;divide the signal by the number of averages.
            if n_avg_flag  eq 1 then begin
                print,'dividing by', n_avg
                (*p_data) /=n_avg
            endif
            
            close, lun
            Free_LUN, lun
            HEAP_GC
            mas_redraw_GUI
            update_status_bar,' '
            return, p_data
            
        endif else if  image_type eq 1 then begin
        
            ;open scans in multiple files which are dynamic
            ;extract the fids that are dynamics then return
            p_data = extract_dynamics_fids ( path_in )
            CLOSE, /all
            mas_redraw_GUI
            update_status_bar,' '
            return, p_data
            
        endif else if  image_type eq 2 then begin
        
            update_status_bar, 'Extracting multi T2 fid data from files'
            ;extract fids associated with multi T2 data
            p_data = extract_multiT2_fids ( )
            CLOSE, /all
            mas_redraw_GUI
            update_status_bar,' '
            return, p_data
            
        endif else if  image_type eq 3 then begin
        
            update_status_bar, 'Extracting multi ADT fid data from files'
            ;extract fids associated with multi T2 data
            p_data = extract_multiADT_fids()
            CLOSE, /all
            mas_redraw_GUI
            update_status_bar,' '
            return, p_data
            
        endif else if  image_type eq 4 then begin
        
            update_status_bar, 'Extracting multi Phase Array fid data from files'
            ;extract fids associated with multi T2 data
            p_data = extract_multi_phase_array_fids()
            CLOSE, /all
            mas_redraw_GUI
            update_status_bar,' '
            return, p_data
            
        endif else if image_type eq 7 then begin
        
            update_status_bar, 'Extracting multi ADC fid data from files'
            p_data = extract_multiADT_fids()
            CLOSE, /all
            mas_redraw_GUI
            update_status_bar,' '
            return, p_data
            
        endif else if  image_type eq 11 then begin
        
            update_status_bar, 'Extracting DICOM images'
            p_data = extract_ADT_DICOM()
            CLOSE, /all
            HEAP_GC
            mas_redraw_GUI
            update_status_bar,' '
            return, p_data
            
        endif else if  image_type eq 12 then begin
        
            islice = *project.imndArray[CI].slice_scheme_ptr
            update_status_bar, 'Extracting T2 multi echo single scan'
            print, 'Extracting T2 multi echo single scan'
            ;extract fid associated with T2 data
            ;open the file unit and prepare for
            ;opening
            if (project.big_endian eq 1 and little_endian_flag eq 0) then begin
                OpenR, lun, path_in, /Get_LUN, /SWAP_IF_BIG_ENDIAN;, ERROR = err
            endif else begin
                OpenR, lun, path_in, /Get_LUN, SWAP_IF_LITTLE_ENDIAN=little_endian_flag; , ERROR = err
            endelse
            ;        OpenR, lun, path_in, /Get_LUN,
            ;        SWAP_IF_LITTLE_ENDIAN=little_endian_flag
            if dim eq 3 then p_data = read_cmplx_data_3d(lun,fdim,pdim,sdim,rare_factor,adim) else $
              p_data = rd_cmplx_data_T2(lun,fdim,pdim,sdim,adim,rare_factor,islice)
            
            close, lun
            Free_LUN, lun
            CLOSE, /all
            HEAP_GC
            mas_redraw_GUI
            update_status_bar,' '
            return, p_data
            
        endif else if  image_type eq 13 then begin
        
            islice = *project.imndArray[CI].slice_scheme_ptr
            update_status_bar, 'Extracting single movie'
            
            if (project.big_endian eq 1 and little_endian_flag eq 0) then begin
                OpenR, lun, path_in, /Get_LUN, /SWAP_IF_BIG_ENDIAN ;, ERROR = err
            endif else begin
                OpenR, lun, path_in, /Get_LUN, SWAP_IF_LITTLE_ENDIAN=little_endian_flag; , ERROR = err
            endelse
            ;        OpenR, lun, path_in, /Get_LUN, SWAP_IF_LITTLE_ENDIAN=little_endian_flag
            
            p_data = rd_cmplx_data_movie(lun,fdim,pdim,sdim,adim,rare_factor,islice)
            close, lun
            Free_LUN, lun
            CLOSE, /all
            HEAP_GC
            mas_redraw_GUI
            update_status_bar,' '
            return, p_data
            
        endif else if  image_type eq 14 then begin
        
            update_status_bar, 'Extracting multi movie multi '
            p_data = extract_multi_slice_fids()
            CLOSE, /all
            HEAP_GC
            mas_redraw_GUI
            update_status_bar,' '
            return, p_data
            
        endif else if image_type eq 99 then begin
        
            update_status_bar, 'Extracting onepulse FID'
            if (project.big_endian eq 1 and little_endian_flag eq 0) then begin
                OpenR, lun, path_in, /Get_LUN, /SWAP_IF_BIG_ENDIAN ;, ERROR = err
            endif else begin
                OpenR, lun, path_in, /Get_LUN, SWAP_IF_LITTLE_ENDIAN=little_endian_flag; , ERROR = err
            endelse
            
            ;; acq_size is the number of complex points
            acq_size = project.imndarray[ci].spect_acq_size
            p_data = ptr_new(complexarr(acq_size, adim), /no_copy)
            for a = 0, adim-1 do begin
                ;; function expects number number of data elements 2 x number of complex points
                tmp = rd_cmplx_blk_spec(lun, acq_size*2, word_size=32)
                (*p_data)[*,a] = tmp
            endfor
            close, /all
            free_lun, lun
            mas_redraw_GUI
            update_status_bar, ' '
            return, p_data
            
        endif
        
    endelse ;this endelse is part of the .raw reading
    
    update_status_bar,'  '
    HEAP_GC
    
end


pro mas_extract_files

end
