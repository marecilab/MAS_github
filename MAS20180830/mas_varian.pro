; Routines for reading varian 2d and 3d fid files.
; See: mas_varian_example at the end of this file
; 
; Changes:
;   20100729:
;     - Initial
;   20100806:
;     - Now uses procpar parameter 'seqfil' instead of 'pslabel' to identify
;       pulse sequences. pslabel can be user-defined, but seqfil seems fixed.
;     - Fixed 3D reading based on newly collected data. Now reads ncccn and nccsn.
;       Has not been tested with 3D acquisitons with an arrayed parameter
;   20100810:
;     - Throw an error when trying to read scans with a pairwise navigator
;       echo. This is not supported.
;     - Tested with non-square acquisition matrices and matrices with 
;       non-power-of-two dimensions.
;   20101015:
;     - Added preliminary support for 1d fid reader. Needs to be tested with
;       arrayed spectra (probably does not work), but should work for
;       simple spuls.
;     - Added support for epi with single and (not quite) triple
;       pointwise reference scans, with fract_ky=nv/2 only. Note that
;       reference scans are returned with the data as the first (for
;       single) or first three (for triple) volumes in the 4th data dimension. 
;     - Added a single function (mas_varian_read_fid_data_any) to determine
;       which type from the procpar and call the appropriate reader. There
;       is also a 'recon' keyword which, if present, will contain the fft'd
;       imagery. mas_varian_example has been updated to use this.
;   20101018:
;     - Added a second shift to compensate for n/2 fft shifting in 2nd
;       phase encode dim for 3D acquisitions only.
;   20101124:
;     - Fixed bug in triple reference scan computation in which the final
;       FFT should have been an inverse FFT (it was forward).
;   20110304:
;     - Changlog has been moved to version control
; To Do:
;     - Fix data orientation -- there is a 90 in-plane rotation for 2D and a
;       strange vertical flip of 3D.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;; arg1: pdim  - number of phase-encode steps
;; arg2: etl   - echo train length
;; arg3: kzero - zero offset of kspace
;; returns: a table that indicates how to reorder the traces
;;          to reconstruct kspace
function mas_varian_generate_petable, pdim, etl, kzero

    ;; generate a phase-encode table based on evidence sluthed from
    ;; existing petables. 
    
    dim = pdim/2
    
    top_half = transpose(reform(lindgen(dim)-dim+1,dim/etl,etl))
    bot_half = reverse(transpose(reform(lindgen(dim)+1,dim/etl,etl)))
    
    petable = shift(reverse([ [top_half], [bot_half] ]), kzero-1, 0)
    
    if (kzero gt 2) then begin
        petable[0:kzero-2,*] = reverse(petable[0:kzero-2,*],1)
    endif
    
    return, petable

end

;; arg1: f_hd  - MAS_VARIAN_FILE_HEADER struct as read from a data file
;; returns: a properly sized and typed array for reading.
function mas_varian_make_read_array, f_hd

    if ((size(f_hd, /structure)).type_name eq 'STRUCT') then begin
        if (tag_names(f_hd, /structure_name) ne 'MAS_VARIAN_FILE_HEADER') then begin
            message, "argument should be a structure of type MAS_VARIAN_FILE_HEADER"
        endif
    endif
    
    ;; determine data type of fid contents.
    if ((f_hd.status AND '08'x) eq 0) then begin
        ;; data is integer.
        if ((f_hd.status AND '04'x) eq 0) then begin
            ;; data is 16 bit        
            message, "data is 16-bit integer.", /informational
            tmp = intarr(f_hd.np)
        endif else begin
            ;; data is 32 bit
            message, "data is 32-bit integer", /informational
            tmp = lonarr(f_hd.np)
        endelse
    endif else begin
        ;; data is floating point
        message, "data is 32-bit floating point", /informational
        tmp = fltarr(f_hd.np)
    endelse

    return, tmp

end

;; data (fid) file header information:
;; * nblocks   - the number of data blocks present in the file.
;; * ntraces   - the number of traces in each block.
;; * np        - the number of simple elements (16-bit integers, 32-bit integers,
;;               or 32-bit floating point numbers) in one trace.
;;               It is equal to twice the number of complex data points.
;; * ebytes    - the number of bytes in one element, either 2 (for 16-bit integers
;;               in single precision FIDs) or 4 (for all others).
;; * tbytes    - set to (np*ebytes).
;; * bbytes    - set to (ntraces*tbytes + nbheaders*sizeof(struct datablockhead)).
;;               The size of the datablockhead structure is 28 bytes.
;; * vers_id   - the version identification of present VnmrJ.
;; * nbheaders - the number of block headers per data block.
;; * status    - bits as defined below with their hexadecimal values.
;;
;; 0  S_DATA         0x1    0 = no data, 1 = data
;; 1  S_SPEC         0x2    0 = FID, 1 = spectrum
;; 2  S_32           0x4    *
;; 3  S_FLOAT        0x8    0 = integer, 1 = floating point
;; 4  S_COMPLEX      0x10   0 = real, 1 = complex
;; 5  S_HYPERCOMPLEX 0x20   1 = hypercomplex
;; 6  ----- UNUSED ---- [ per email**: 1 - complex data ] 
;; 7  S_ACQPAR       0x80   0 = not Acqpar, 1 = Acqpar
;; 8  S_SECND        0x100  0 = first FT, 1 = second FT
;; 9  S_TRANSF       0x200  0 = regular, 1 = transposed
;; 10 ----- UNUSED -------------------------------------
;; 11 S_NP           0x800  1 = np dimension is active
;; 12 S_NF           0x1000 1 = nf dimension is active
;; 13 S_NI           0x2000 1 = ni dimension is active
;; 14 S_NI2          0x4000 1 = ni2 dimension is active
;; 15 ----- UNUSED -------------------------------------
;;
;; All other bits must be zero.
;; Bits 0–6: file header and block header status bits (bit 6 is unused):
;; * If S_FLOAT=0, S_32=0 for 16-bit integer, or S_32=1 for 32-bit integer.
;; If S_FLOAT=1, S_32 is ignored.
;; Bits 7–14: file header status bits (bits 10 and 15 are unused):
;; Block headers are defined by the following C specifications:
;; 
;; To make life VERY CONFUSED there is no clarity in the status word 
;; as to the use of bits 4 and 6 .... if one of the two is set 
;; IT'S COMPLEX DATA.
pro mas_varian_file_header__define

    ;; See Vnmrj User Programming Manual, Page 247
    struct = { MAS_VARIAN_FILE_HEADER, $ ;; 32 bytes
               nblocks: long(0), $ ;; number of blocks in file
               ntraces: long(0), $ ;; number of traces per block
               np: long(0), $      ;; number of elements per trace
               ebytes: long(0), $  ;; number of bytes per element
               tbytes: long(0), $  ;; number of bytes per trace
               bbytes: long(0), $  ;; number of bytes per block
               vers_id: 0, $       ;; software version, file_id status bits
               status: 0, $        ;; status of whole file
               nbheaders: long(0)$ ;; number of blk hdrs per blk
             }

end

;; 28 bytes
pro mas_varian_block_header__define

    struct = { MAS_VARIAN_BLOCK_HEADER, $
               scale: 0, $         ;; scaling factor
               status: 0, $        ;; status of data in block
               index: 0, $         ;; block index
               mode: 0, $          ;; mode of data in block
               ctcount: long(0), $ ;; ct valus for FID
               lpval: float(0), $  ;; f2 (2D-f1) left phase in phasefile
               rpval: float(0), $  ;; f2 (2D-f1) right phase in phasefile
               lvl: float(0), $    ;; level drift correction
               tlt: float(0) $     ;; tilt drift correction
              }

end

;; Read a 1D fid file.
;; arg1: lun - the pre-opened file LUN, as from IDL's openr routine
;; arg2: opp - an instance of the mas_varian_procpar
;;
;; returns a pointer to the complex kspace array
function mas_varian_read_fid_data_1d, lun, opp

    f_hd = create_struct(name='MAS_VARIAN_FILE_HEADER')
    b_hd = create_struct(name='MAS_VARIAN_BLOCK_HEADER')

    np = opp->lookup('np', not_found=ntfd)
    
    readu, lun, f_hd
    
    tmp = mas_varian_make_read_array(f_hd)
    
    adim = f_hd.nblocks

    data = complexarr(np/2, adim)

    for a = 0, adim-1 do begin
    
        readu, lun, b_hd
    
        readu, lun, tmp
    
        data[*,a] = complex(tmp[0:np-1:2], -tmp[1:np-1:2])
    
    endfor
    
    return, ptr_new(data)
    
end


;; Read a 2D fid file.
;; arg1: lun - the pre-opened file LUN, as from IDL's openr routine
;; arg2: opp - an instance of the mas_varian_procpar
;;
;; returns a pointer to the complex kspace array
;; Modifications
;; 1)Fixed bug reading mems and mgems code (Magdoom, 10/29/16)

function mas_varian_read_fid_data_2d, lun, opp

    f_hd = create_struct(name='MAS_VARIAN_FILE_HEADER')
    b_hd = create_struct(name='MAS_VARIAN_BLOCK_HEADER')
    
    ;; these should exist for any sequence. 
    names = [ 'seqcon', 'ni', 'np', 'nf', 'arraydim', $
              'ns', 'nD', 'ne', 'nv', 'nv2' ]
    error = bytarr(n_elements(names))
    par = 0
    seqcon     = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    ni         = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    np         = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    nf         = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    arraydim   = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    num_slices = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    num_dims   = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    num_echos  = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    nv         = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    nv2        = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    if (total(error ne 0) ne 0) then begin
        msg = "One or more crucial parameters are missing from the procpar file: "
        msg += strjoin(names[where(error ne 0)], ', ')+'.'
        message, msg
    endif
        
    seqfil = opp->lookup('seqfil')
    navigator = opp->lookup('navigator', not_found=ntfd)
    if (ntfd eq 1 || navigator ne 'y') then begin
        navigator = 'n'
    endif else begin
        nav_type = opp->lookup('nav_type')
        nav_echo = opp->lookup('nav_echo')
        if (nav_type ne 'linear') then begin
            message, /info, "nav_type is '"+nav_type+"', but right now, only linear navigator type is supported."
        endif
    endelse
    
    ;; convert the string into a character array
    seqcon = string(transpose(byte(seqcon)))
    
    data_dims = opp->computeDataDimensions()
    f_dim = data_dims[0] ;; freq. dimension
    p_dim = data_dims[1] ;; phase dimension
    s_dim = data_dims[2] ;; slice (or phase2) dim.
    a_dim = data_dims[3] ;; volume dimension
    data = complexarr(f_dim, p_dim, s_dim, a_dim)
    ;; read in the file header.
    readu, lun, f_hd
    
    tmp = mas_varian_make_read_array(f_hd)

    ;; these will be set in certain circumstances, in particular
    ;; if the sequences is a fsems (rare).
    if (opp->paramExists('etl', value=etl) eq 0) then begin
        etl = 1L ;; echo train length
    endif
    
    if (opp->paramExists('kzero', value=kzero) eq 0) then begin
        kzero = 1 ;; kzero
    endif
    
    ;; In many cases, if petable exists and is not "" then there will
    ;; be also a pelist that contains the petable data. This is not
    ;; always the case, and I have one example where pelist is not
    ;; present even though petable is. 
    ;; See: Vnmrj Imaging Manual, page 74 "Table Ordered Data"
    pelist = lindgen(nv) ;; default pelist shall be sequential
    if (opp->paramExists('petable', value=petable) ne 0 && petable ne "") then begin
        if (strmid(seqfil, 0,3) eq 'fse') then begin
            if (opp->paramExists('pelist', value=pelist) eq 0) then begin
                message, "petable = "+petable+", but there is no pelist. "+$
                         "I'm going to generate one, but be warned that the "+$
                         "phase-encode steps could be incorrect.", /continue
                         ;; note that this makes a sequential pelist for etl=1, kzero=1.
                pelist = (mas_varian_generate_petable(nv, etl, kzero))[*]
            endif else begin
                message, "found pelist for petable = "+petable, /informational
            endelse
            ;; need to massage pelist a little more
            pelist = transpose(reform((long(pelist+abs(min(pelist)))), etl, p_dim/etl))
            pelist = pelist[*]
        endif
    endif else begin
        message, "Using sequential pelist.", /informational
    endelse
    
    ;; There is most certainly a way to optimize these loops to make
    ;; them more generic, but for now I would like to keep them like this
    ;; since it makes explicit the differences between the various
    ;; seqcon codes (which are still a little mysterious.)
    if (seqcon[0] eq 'c' and seqcon[2] eq 'c') then begin
       readu, lun, b_hd
        for p = 0, p_dim-1 do begin
            for s = 0, (s_dim+n_elements(nav_echo))-1 do begin
                for a = 0, a_dim-1 do begin
                    readu, lun, tmp
                    data[*, p, s, a] = complex(tmp[0:np-1:2], tmp[1:np-1:2])
                    if (navigator eq "y") then begin
                       readu, lun, tmp ;; discard the navigator echo
                    endif
                endfor
            endfor
        endfor
    endif else if (seqcon[0] eq 'c' and seqcon[2] eq 's') then begin
          
          for p = 0, p_dim-1 do begin
            readu, lun, b_hd
            for s = 0, (s_dim+n_elements(nav_echo))-1 do begin
              for a = 0, a_dim-1 do begin
                readu, lun, tmp
                data[*, p, s, a] = complex(tmp[0:np-1:2], tmp[1:np-1:2])
                if (navigator eq "y") then begin
                  readu, lun, tmp ;; discard the navigator echo
                endif
              endfor
            endfor
          endfor
      endif else if (seqcon[0] eq 's') then begin
          for a = 0, a_dim-1 do begin
            readu, lun, b_hd
            for p = 0, p_dim-1 do begin
              for s = 0, (s_dim+n_elements(nav_echo))-1 do begin
                readu, lun, tmp
                data[*, p, s, a] = complex(tmp[0:np-1:2], tmp[1:np-1:2])
                if (navigator eq "y") then begin
                  readu, lun, tmp ;; discard the navigator echo
                endif
              endfor
            endfor
          endfor
     endif else if (seqcon[1] eq 'c' and seqcon[2] eq 'c') then begin
        ;; compressed-compressed nccnn
        pe_blocks = p_dim/etl
        
        for a=0L, a_dim-1 do begin ;; the number of volumes
            readu, lun, b_hd
            for p=0L, pe_blocks-1 do begin ;; the number of PE blocks
                for s=0L, s_dim-1 do begin ;; the number of slices
                    echonum = 0
                    for e=0L, etl-1 do begin ;; echo train length
                        readu, lun, tmp
                        echonum++
                        pe_index = p + pe_blocks*e
                        data[*, pelist[pe_index], s, a] = complex(tmp[0:np-1:2], tmp[1:np-1:2])
                        if (navigator eq "y") then begin
                            if (nav_type eq 'linear') then begin
                                readu, lun, tmp ;; discard the navigator echo
                            endif else if (nav_type eq 'pairwise') then begin
                                ;; pairwise nagivator is experimental
                                if (echonum eq etl) then begin
                                    for n = 0, n_elements(nav_echo)-1 do begin
                                        readu, lun, tmp
                                    endfor
                                endif
                            endif
                        endif
                    endfor
                endfor
            endfor
        endfor
    endif else if (seqcon[1] eq 'c' and seqcon[2] eq 's') then begin
        ;; standard ncsnn
        for p=0L, p_dim-1 do begin
            for a=0L, a_dim-1 do begin
                readu, lun, b_hd
                for s=0L, s_dim-1 do begin
                    readu, lun, tmp
                    data[*, p, s, a] = complex(tmp[0:np-1:2], tmp[1:np-1:2])
                    if (navigator eq "y") then begin
                        readu, lun, tmp
                    endif
                endfor
            endfor
        endfor
    endif

    ;; now reorganize slices
    pss = opp->lookup('pss', not_found=not_found)
    if (not_found eq 0) then begin
        reslice = sort(pss)
        data = temporary(data[*,*,reslice,*])
    endif
    
    ;;data = reverse(data, 2, /overwrite)
    return, ptr_new(data, /no_copy)
    
end

;; Read a 3D fid file.
;; arg1: lun - the pre-opened file LUN, as from IDL's openr routine
;; arg2: opp - an instance of the mas_varian_procpar
;; returns a pointer to the complex kspace array
;; Modifications
;; - Fixed bug that could read seqcon - nccsn acquisitions
;; - Now able to read seqcon - cccsn acquisitions (02/19/2016)

function mas_varian_read_fid_data_3d, lun, opp

    f_hd = create_struct(name='MAS_VARIAN_FILE_HEADER')
    b_hd = create_struct(name='MAS_VARIAN_BLOCK_HEADER')

    ;; these should exist for any sequence. 
    names = [ 'seqcon', 'ni', 'np', 'nf', 'arraydim', $
              'ns', 'nD', 'ne', 'nv', 'nv2' ]
    error = bytarr(n_elements(names))
    par = 0
    seqcon     = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    ni         = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    np         = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    nf         = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    arraydim   = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    num_slices = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    num_dims   = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    num_echos  = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    nv         = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    nv2        = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    if (total(error ne 0) ne 0) then begin
        msg = "One or more crucial parameters are missing from the procpar file: "
        msg += strjoin(names[where(error ne 0)], ', ')+'.'
        message, msg
    endif

    seqcon = string(transpose(byte(seqcon)))

    data_dims = opp->computeDataDimensions()
    data = complexarr(data_dims)
    f_dim = data_dims[0] ;; freq. dimension
    p_dim = data_dims[1] ;; phase dimension
    s_dim = data_dims[2] ;; slice (or phase2) dim.
    a_dim = data_dims[3] ;; volume dimension
    if n_elements(data_dims) gt 4 then a_dim2 = data_dims[4]
    
    readu, lun, f_hd
    
    tmp = mas_varian_make_read_array(f_hd)

    ;; these will be set in certain circumstances, in particular
    ;; if the sequences is a fse3d (rare).
    if (opp->paramExists('etl', value=etl) eq 0) then begin
        etl = 1L ;; echo train length
    endif

    if (opp->paramExists('kzero', value=kzero) eq 0) then begin
        kzero = 1 ;; kzero
    endif

    ;; In many cases, if petable exists and is not "" then there will
    ;; be also a pelist that contains the petable data. This is not
    ;; always the case, and I have one example where pelist is not
    ;; present even though petable is. 
    ;; See: Vnmrj Imaging Manual, page 74 "Table Ordered Data"
    pelist = lindgen(nv) ;; default pelist shall be sequential
    seqfil = opp->lookup('seqfil')
    if (opp->paramExists('petable', value=petable) ne 0 && petable ne "") then begin
        if (strmid(seqfil, 0,3) eq 'fse') then begin
            if (opp->paramExists('pelist', value=pelist) eq 0) then begin
                message, "petable = "+petable+", but there is no pelist. "+$
                         "I'm going to generate one, but be warned that the "+$
                         "phase-encode steps could be incorrect.", /continue
                         ;; note that this makes a sequential pelist for etl=1, kzero=1.
                pelist = (mas_varian_generate_petable(nv, etl, kzero))[*]
            endif else begin
                message, "found pelist for petable = "+petable, /informational
            endelse
            ;; need to massage pelist a little more
            pelist = transpose(reform((long(pelist+abs(min(pelist)))), etl, p_dim/etl))
            pelist = pelist[*]
            pelist = float(pelist)/(max(pelist)) * (nv-1)
        endif
    endif else begin
        message, "Using sequential pelist.", /informational
    endelse
    
    ; seqcon : [0] MultiEcho [2] PE1 [3] PE2
    if (seqcon[0] eq 'n' and seqcon[2] eq 'c' and seqcon[3] eq 'c') then begin
        ;; compressed-compressed; n_cc_
        message, "Data is 3D no loop-compressed-compressed.", /informational
        for a=0, a_dim-1 do begin         ;; the number of volumes
            readu, lun, b_hd
            for s=0, s_dim-1 do begin     ;; the number of PE2 lines
                for p=0, p_dim-1 do begin ;; the number of PE1 lines
                    readu, lun, tmp
                    data[*,p,s,a] = complex(tmp[0:np-1:2], tmp[1:np-1:2])
                endfor
            endfor
        endfor
        
    endif else if (seqcon[0] eq 'n' and seqcon[2] eq 'c' and seqcon[3] eq 's') then begin
        ;; compressed; n_cs_
        message, "Data is 3D no loop-compressed-standard.", /informational
        pe_blocks = p_dim/etl
        for s=0L, s_dim-1 do begin ;; the number of PE2 lines
              for a=0L, a_dim-1 do begin ;; the number of volumes
                readu, lun, b_hd
                for p=0L, pe_blocks-1 do begin ;; the number of PE1 blocks
                    echonum = 0
                    for e=0L, etl-1 do begin ;; echo train length
                        readu, lun, tmp
                        echonum++
                        pe_index = p + pe_blocks*e
                        data[*, pelist[pe_index], s, a] = complex(tmp[0:np-1:2], tmp[1:np-1:2])
                    endfor
                endfor
            endfor
        endfor

    endif else if (seqcon[0] eq 'n' and seqcon[2] eq 's' and seqcon[3] eq 'c') then begin
        ;; compressed; n_sc_
        message, "Data is 3D no loop-standard-compressed.", /informational
        for p=0, p_dim-1 do begin             ;; the number of PE1 lines
           for a=0, a_dim-1 do begin           ;; the number of volumes
              readu, lun, b_hd
              for s=0, s_dim-1 do begin           ;; the number of PE2 lines
                readu, lun, tmp
                data[*,p,s,a] = complex(tmp[0:np-1:2], tmp[1:np-1:2])
              endfor
           endfor
         endfor
     endif else if (seqcon[0] eq 'c' and seqcon[2] eq 'c' and seqcon[3] eq 's') then begin
      if n_elements(data_dims) gt 4 then begin 
        ;; compressed; c_cs_
        message, "Data is arrayed 3D multi-echo compressed-compressed-standard.", /informational
        for s=0, s_dim-1 do begin             ;; the number of PE2 lines
          for a2=0, a_dim2-1 do begin
            readu, lun, b_hd                    
            for p=0, p_dim-1 do begin         ;; the number of PE1 lines
              for a=0, a_dim-1 do begin       ;; the number of echoes
                 readu, lun, tmp
                 data[*,p,s,a,a2] = complex(tmp[0:np-1:2], tmp[1:np-1:2])  
               endfor              
            endfor
          endfor
        endfor
       data = reform(transpose(data,[0,1,2,4,3]),f_dim,p_dim,s_dim,a_dim*a_dim2,/overwrite)
      endif else begin
        ;; compressed; c_cs_
        message, "Data is 3D multi-echo compressed-compressed-standard.", /informational
        for s=0, s_dim-1 do begin             ;; the number of PE2 lines
          readu, lun, b_hd
          for p=0, p_dim-1 do begin         ;; the number of PE1 lines
            for a=0, a_dim-1 do begin       ;; the number of echoes
              readu, lun, tmp
              data[*,p,s,a] = complex(tmp[0:np-1:2], tmp[1:np-1:2])
            endfor
          endfor
        endfor
      endelse
    endif else if (seqcon[0] eq 'c' and seqcon[2] eq 's' and seqcon[3] eq 's') then begin
      if n_elements(data_dims) gt 4 then begin
        ;; compressed; c_ss_
        message, "Data is arrayed 3D multi-echo compressed-standard-standard.", /informational
        for s=0, s_dim-1 do begin             ;; the number of PE2 lines
           for p=0, p_dim-1 do begin         ;; the number of PE1 lines
            for a2=0, a_dim2-1 do begin
              readu, lun, b_hd
              for a=0, a_dim-1 do begin       ;; the number of echoes
                readu, lun, tmp
                data[*,p,s,a,a2] = complex(tmp[0:np-1:2], tmp[1:np-1:2])
              endfor
            endfor
          endfor
        endfor
        data = reform(transpose(data,[0,1,2,4,3]),f_dim,p_dim,s_dim,a_dim*a_dim2,/overwrite)
      endif else begin
        ;; compressed; c_ss_
        message, "Data is 3D multi-echo compressed-standard-standard.", /informational
        for s=0, s_dim-1 do begin             ;; the number of PE2 lines
          for p=0, p_dim-1 do begin         ;; the number of PE1 lines
            readu, lun, b_hd
            for a=0, a_dim-1 do begin       ;; the number of echoes
              readu, lun, tmp
              data[*,p,s,a] = complex(tmp[0:np-1:2], tmp[1:np-1:2])
            endfor
          endfor
        endfor
      endelse  
    endif else if (seqcon[0] eq 'c' and seqcon[2] eq 'c' and seqcon[3] eq 'c') then begin
        ;; compressed; c_cc_
        if n_elements(data_dims) gt 4 then begin
          ;; compressed; c_cc_
          ;; Echo and array dimension collapsed into one so the data is still 3D as (array-echo-array-echo....)
          message, "Data is arrayed 3D multi-echo compressed-compressed-compressed.", /informational
          readu, lun, b_hd
          for a2=0, a_dim2-1 do begin
            for s=0, s_dim-1 do begin           ;; the number of PE2 lines
              for p=0, p_dim-1 do begin         ;; the number of PE1 lines
                for a=0, a_dim-1 do begin       ;; the number of echoes
                  readu, lun, tmp
                  data[*,p,s,a,a2] = complex(tmp[0:np-1:2], tmp[1:np-1:2])
                endfor
              endfor
            endfor
          endfor
          data = reform(transpose(data,[0,1,2,4,3]),f_dim,p_dim,s_dim,a_dim*a_dim2,/overwrite)
        endif else begin
          message, "Data is 3D multi-echo compressed-compressed-compressed.", /informational
          readu, lun, b_hd
          for s=0, s_dim-1 do begin           ;; the number of PE2 lines
            for p=0, p_dim-1 do begin         ;; the number of PE1 lines
              for a=0, a_dim-1 do begin       ;; the number of echoes
                readu, lun, tmp
                data[*,p,s,a] = complex(tmp[0:np-1:2], tmp[1:np-1:2])
              endfor
            endfor
          endfor
        endelse
        
    endif else if (seqcon[0] eq 'c' and seqcon[2] eq 's' and seqcon[3] eq 'c') then begin
        ;; compressed; c_sc_
        message, "Data is 3D multi-echo compressed-standard-compressed.", /informational
        for p=0, p_dim-1 do begin             ;; the number of PE1 lines
          readu, lun, b_hd
          for s=0, s_dim-1 do begin           ;; the number of PE2 lines
            for a=0, a_dim-1 do begin         ;; the number of echoes
              readu, lun, tmp
              data[*,p,s,a] = complex(tmp[0:np-1:2], tmp[1:np-1:2])
            endfor
          endfor
        endfor
    endif else begin
        message, "Seqcon not currently supported."
    endelse
    
    return, ptr_new(data, /no_copy)
    
end

;; Read kspace from an EPI experiment.
;; arg1: lun - the pre-opened file LUN, as from IDL's openr routine
;; arg2: opp - an instance of the mas_varian_procpar
;; returns a pointer to the complex kspace array
function mas_varian_read_fid_epi, lun, opp

    f_hd = create_struct(name='MAS_VARIAN_FILE_HEADER')
    b_hd = create_struct(name='MAS_VARIAN_BLOCK_HEADER')

    ;; these should exist for any sequence. 
    names = [ 'seqcon', 'ni', 'np', 'nf', 'arraydim', $
              'ns', 'nD', 'ne', 'nv', 'nv2' ]
    error = bytarr(n_elements(names))
    par = 0
    seqcon     = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    ni         = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    np         = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    nf         = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    arraydim   = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    num_slices = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    num_dims   = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    num_echos  = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    nv         = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    nv2        = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    if (total(error ne 0) ne 0) then begin
        msg = "One or more crucial parameters are missing from the procpar file: "
        msg += strjoin(names[where(error ne 0)], ', ')+'.'
        message, msg
    endif
        
    seqfil = opp->lookup('seqfil')
    navigator = opp->lookup('navigator', not_found=ntfd)
    if (ntfd eq 1 || navigator ne 'y') then begin
        navigator = 'n'
    endif else begin
        nav_type = opp->lookup('nav_type')
        nav_echo = opp->lookup('nav_echo')
        if (nav_type ne 'linear') then begin
            message, "nav_type is '"+nav_type+"', but right now, only linear navigator type is supported."
        endif
    endelse
    
    ;; convert the string into a character array
    seqcon = string(transpose(byte(seqcon)))
    
    data_dims = opp->computeDataDimensions()
    f_dim = data_dims[0] ;; freq. dimension
    p_dim = data_dims[1] ;; phase dimension
    s_dim = data_dims[2] ;; slice (or phase2) dim.
    a_dim = data_dims[3] ;; volume dimension
    data = complexarr(f_dim, p_dim, s_dim, a_dim)
    ;; read in the file header.
    readu, lun, f_hd
    
    tmp = mas_varian_make_read_array(f_hd)

    ;; these will be set in certain circumstances, in particular
    ;; if the sequences is a fsems (rare).
    if (opp->paramExists('etl', value=etl) eq 0) then begin
        etl = 1L ;; echo train length
    endif
    
    if (opp->paramExists('kzero', value=kzero) eq 0) then begin
        kzero = 1 ;; kzero
    endif

    if (opp->paramExists('fract_ky', value=fract_ky) eq 0) then begin
        fract_ky = 0 ;; kzero
    endif

    ;; Number of "shots" for multi-shot.
    if (opp->paramExists('nseg', value=nseg) eq 0) then begin
        nseg = 1 ;; kzero
    endif
    
    ;; In many cases, if petable exists and is not "" then there will
    ;; be also a pelist that contains the petable data. This is not
    ;; always the case, and I have one example where pelist is not
    ;; present even though petable is. 
    ;; See: Vnmrj Imaging Manual, page 74 "Table Ordered Data"
    if (opp->paramExists('petable', value=petable) ne 0 && petable ne "") then begin

        if (opp->paramExists('pe_table', value=pe_table) eq 0) then begin
            message, "petable = "+petable+", but there is no pelist. "+$
                     "I'm going to generate one, but be warned that the "+$
                     "phase-encode steps could be incorrect.", /continue
                     ;; note that this makes a sequential pelist for etl=1, kzero=1.
            pe_table = lindgen(nv)
        endif else begin
            message, "found pe_table for petable = "+petable, /informational
        endelse
        ;; need to massage pelist a little more
        pe_table += abs(min(pe_table))
        pe_table = pe_table[*]
    endif else begin
        message, "Using sequential pe_table.", /informational
        pe_table = lindgen(nv) ;; default pe_table shall be sequential
    endelse
    ;;print, 2*(nv - fract_ky)
    ;; EPI is ncccn 
    if (seqcon[1] eq 'c' and seqcon[2] eq 'c') then begin
        ;; compressed-compressed nccnn
        etl = ((fract_ky eq nv/2) ? nv : nv/2 + fract_ky)/nseg
        pe_blocks = long(nseg) ;; PE blocks will be split into shots
        for a=0L, a_dim-1 do begin ;; the number of volumes
            readu, lun, b_hd
            for p=0L, pe_blocks-1 do begin ;; the number of PE blocks
                for s=0L, s_dim-1 do begin ;; the number of slices
                    for e=0L, etl-1 do begin ;; echo train length
                        readu, lun, tmp
                        pe_index = p + pe_blocks*e
                        if ((e mod 2) ne 0) then begin
                            blk = reverse(complex(tmp[0:np-1:2], tmp[1:np-1:2]))
                        endif else begin
                            blk = complex(tmp[0:np-1:2], tmp[1:np-1:2])
                        endelse
                        data[*, pe_index, s, a] = blk
                        
                        if (navigator eq "y") then begin
                            if (nav_type eq 'linear') then begin
                                readu, lun, tmp ;; discard the navigator echo
                            endif else if (nav_type eq 'pairwise') then begin
                                ;; pairwise not implemented
                            endif
                        endif
                        
                    endfor
                endfor
            endfor
        endfor
    endif
    
    ;; EPI reference scan corrections
    if (opp->paramExists('epiref_type', value=epiref_type)) then begin
        if (epiref_type eq 'single') then begin
            message, "Applying reference scan correction (type=single)", /informational
            ref_data = reform(data[*,*,*,0])
            phase_corr = exp(-complex(0,1)*atan(fft(ref_data, dimension=1), /phase))
            for v = 1, a_dim-1 do begin
                data[*,*,*,v] = fft( fft(data[*,*,*,v], dimension=1)*phase_corr, dimension=1, /inverse)
            endfor
        endif else if (epiref_type eq 'triple') then begin
           ;; triple is not supported yet, even though the code
           ;; exists. Needs to be tested further. 
            if (0) then begin
                message, "WARNING: triple reference scan not supported.", /informational
            endif else begin
                message, "Applying reference scan correction (type=triple)", /informational
                
                ref_indicator = abs(opp->lookup('image'))
                
                ref_data0 =         reform(data[*,*,*,ref_indicator[0]])
                ref_data1 = reverse(reform(data[*,*,*,ref_indicator[1]]), 1)
                ref_data2 = reverse(reform(data[*,*,*,ref_indicator[2]]), 1)

                phase0 = exp(-complex(0,1)*atan(fft(ref_data0, dimension=1), /phase))
                phase2 = exp(-complex(0,1)*atan(fft(ref_data2, dimension=1), /phase))
                
                for a = 3, a_dim-1 do begin
                    data[*,*,*,a] = fft( fft(ref_data1    , dimension=1)*phase2 + $
                                         fft(data[*,*,*,a], dimension=1)*phase0, $
                                         dimension=1, /inverse)
                endfor
            endelse
        endif
    endif
    
    ;; now reorganize slices
    pss = opp->lookup('pss', not_found=not_found)
    if (not_found eq 0) then begin
        reslice = sort(pss)
        data = temporary(data[*,*,reslice,*])
    endif
    
    ;;data = reverse(data, 2, /overwrite)
    return, ptr_new(data, /no_copy)
 

end

;; Read kspace from an EPIP experiment.
;; arg1: lun - the pre-opened file LUN, as from IDL's openr routine
;; arg2: opp - an instance of the mas_varian_procpar
;; returns a pointer to the complex kspace array
;; Modifications (Magdoom,04/18/15)
;; 1) Fixed ghosting issue (Every odd echo reversed instead of every other echo)
;; 2) Incorporated scaled EPI REF phase correction

function mas_varian_read_fid_epip, lun, opp
    f_hd = create_struct(name='MAS_VARIAN_FILE_HEADER')
    b_hd = create_struct(name='MAS_VARIAN_BLOCK_HEADER')

    ;; these should exist for any sequence. 
    names = [ 'seqcon', 'ni', 'np', 'nf', 'arraydim', $
              'ns', 'nD', 'ne', 'nread', 'nphase', 'oversample' ]
    error = bytarr(n_elements(names))
    par = 0
    seqcon     = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    ni         = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    np         = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    nf         = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    arraydim   = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    num_slices = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    num_dims   = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    num_echos  = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    nread      = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    nphase     = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    oversample = opp->lookup(names[par], not_found=ntfd) & error[par++] = ntfd
    nv         = nphase
    
    if (total(error ne 0) ne 0) then begin
        msg = "One or more crucial parameters are missing from the procpar file: "
        msg += strjoin(names[where(error ne 0)], ', ')+'.'
        message, msg
    endif
        
    seqfil = opp->lookup('seqfil')
    navigator = opp->lookup('navigator', not_found=ntfd)
    if (ntfd eq 1 || navigator ne 'y') then begin
        navigator = 'n'
        nav_echo = 1
    endif else begin
        nav_type = opp->lookup('nav_type')
        nav_echo = opp->lookup('nav_echo')
        if (nav_type ne 'linear') then begin
            message, /info, "nav_type is '"+nav_type+"', but right now, only linear navigator type is supported."
        endif
    endelse
    
    ;; convert the string into a character array
    seqcon = string(transpose(byte(seqcon)))
    
    ;data_dims = opp->computeDataDimensions()
    f_dim = oversample*nread/2 ;; freq. dimension
    p_dim = nphase ;; phase dimension
    s_dim = num_slices ;; slice (or phase2) dim.
    a_dim = arraydim ;; volume dimension
    data = complexarr(f_dim, p_dim, s_dim, a_dim)
    ;; read in the file header.
    readu, lun, f_hd
    
    np = nread*oversample
    tmp = mas_varian_make_read_array(f_hd)
    ;; Note we don't want blocks, we want readout-sized array
    ;; so we will truncate what came back up to readout size.
    tmp = tmp[0:np-1]
    
    ;; these will be set in certain circumstances, in particular
    ;; if the sequences is a fsems (rare).
    if (opp->paramExists('etl', value=etl) eq 0) then begin
        etl = 1L ;; echo train length
    endif
    
    if (opp->paramExists('kzero', value=kzero) eq 0) then begin
        kzero = 1 ;; kzero
    endif

    ;; Number of "shots" for multi-shot.
    if (opp->paramExists('nseg', value=nseg) eq 0) then begin
        nseg = 1 ;; kzero
    endif
    
    ;; Load the indicator to tell which images are which
    image_indicator = opp->lookup('image')
    
    ;; Check if the acquistion is compressed
    cseg = opp->lookup('cseg')
    
    ;; Check if the readout is alternated for each shot (only for odd number of shots)
    altread = opp->lookup('altread')
    
    etl = nv/nseg
    pe_blocks = long(nseg) ;; PE blocks will be split into shots    
    if (cseg eq 'n') then order_dim = [pe_blocks-1,s_dim-1] else order_dim = [s_dim-1,pe_blocks-1]
    for a=0L, a_dim-1 do begin ;; the number of volumes
        readu, lun, b_hd
        for i=0L, order_dim[0] do begin ;; the number of slices
           for j=0L, order_dim[1] do begin ;; the number of PE blocks
                if (cseg eq 'n') then begin
                  p = i
                  s = j
                endif else begin
                  p = j
                  s = i
                endelse
                if (nav_echo eq 1) then readu, lun, tmp ;; eat up the navigator                 
                readu, lun, tmp ;; eat the extra read gradient (what is this?)
                for e=0L, etl-1 do begin ;; echo train length
                readu, lun, tmp
                pe_index = pe_blocks*e-nv+p+1
                if (nseg mod 2 ne 0) then begin
                  if (altread eq 'y') then begin                 
                    if (pe_index mod 2 eq 0) then begin
                         blk = reverse(complex(tmp[0:np-1:2], tmp[1:np-1:2]))
                     endif else begin
                         blk = complex(tmp[0:np-1:2], tmp[1:np-1:2])
                     endelse 
                    endif 
                   endif else begin
                     if (e mod 2 eq 0) then begin
                         blk = complex(tmp[0:np-1:2], tmp[1:np-1:2])
                     endif else begin
                         blk = reverse(complex(tmp[0:np-1:2], tmp[1:np-1:2]))
                     endelse
                 endelse
                
             data[*, pe_index, s, a] =  blk  
                                    
;             if (navigator eq "y") then begin ;; don't know if this applies.
;                 if (nav_type eq 'linear') then begin
;                     readu, lun, tmp ;; discard the navigator echo
;                  endif else if (nav_type eq 'pairwise') then begin
;                                ;; pairwise not implemented
;                  endif
;              endif
             endfor ;; etl
         endfor ;; s_dim
      endfor ;; pe_blocks
   endfor ;; a_dim
 
  message, /info, "At EOF: "+strtrim(string(eof(lun)), 2)
    
  ;; EPI reference scan phase corrections
  if opp->paramExists('epi_pc', value=epiref_type) then begin
    if (epiref_type eq 'SINGLE' or epiref_type eq 'POINTWISE') then begin
      message, "Applying reference scan correction (type=single)", /informational
      ref_data = reform(data[*,*,*,0])
      phase_corr = exp(-complex(0,1)*atan(fft(ref_data, dimension=1), /phase))
      for v = 1, a_dim-1 do data[*,*,*,v] = fft(fft(data[*,*,*,v], dimension=1)*phase_corr, dimension=1, /inverse)
      data = shift(data,f_dim/2,0,0,0)
     endif else begin
      message, "Applying reference scan correction (type="+epiref_type+")", /informational
      ;; Prepare the data for correction by FFT in readout direction.
      data = fft(data, dimension=1, /overwrite)    
      for a = 0, a_dim-1 do begin
            
          case image_indicator[a] of 
                
          0: begin ;; Reference scan with no phase encode
             ref_data0 = reform(data[*,*,*,a])
             phase0   = exp(-complex(0,1)*atan(ref_data0, /phase))
             end
         -2: begin ;; Reference scan with no phase encode, inverted RO
             ref_data2 = reform(data[*,*,*,a])
             phase2 = exp(-complex(0,1)*atan(ref_data2, /phase))
             end
          -1: ref_data1 = reform(data[*,*,*,a]) ;; Reference scan with phase encode, inverted RO
 
           1: begin ;; And here is the actual image.
              image1 = reverse(fft(ref_data1*phase2,dimension=2))
              image2 = fft(data[*,*,*,a]*phase0,dimension=2)            
              if a_dim gt 4 and epiref_type eq 'SCALED_TRIPLE' then begin
                 noisefrac = 0.05                   ; Fraction of voxels used to sample noise
                 noiselvl = stddev(abs(image1[f_dim/2-f_dim*noisefrac/2:f_dim/2+f_dim*noisefrac/2,p_dim/2-p_dim*noisefrac/2: p_dim/2+p_dim*noisefrac/2]))
                 eta = abs(image2)/(abs(image1)+abs(image2)*(abs(image1) lt noiselvl))
                 data[*,*,*,a] = (eta*image1 +image2)/2.0 
              endif else data[*,*,*,a] = (image1 +image2)/2.0       
              end                     
            endcase                        
         endfor
            
         data = fft(data, dimension=1, /inverse, /overwrite)
         data = fft(data, dimension=2, /inverse, /overwrite)
         data = shift(data,f_dim/2,0,0,0)
                   
      endelse
  endif
    
  ;; now reorganize slices
  pss = opp->lookup('pss', not_found=not_found)
  if (not_found eq 0) then begin
    reslice = sort(pss)
    data = temporary(data[*,*,reslice,*])
  endif
  return, ptr_new(data, /no_copy)
 
end

;; Read any fid file, type determined from procpar contents
;; arg1: lun - the pre-opened file LUN, as from IDL's openr routine
;; arg2: opp - an instance of the mas_varian_procpar
;; keyword: recon: if present will contain the fftd imagery as a pointer
;; returns a pointer to the complex kspace array
function mas_varian_read_fid_data_any, lun, opp, recon=recon

    if (not opp->paramExists('nD', value=num_dims)) then begin
        if (opp->lookup('seqfil') eq 'spuls') then begin
            num_dims = 1
        endif
    endif else begin
        dims = opp->computeDataDimensions()
    endelse
    
    ;; data will be a pointer to kspace data
    ;; img will be reconstructed imagery
    case num_dims of
       1: begin
          data = mas_varian_read_fid_data_1d(lun, opp)
          if (arg_present(recon)) then begin
              img = shift(fft(data, dimension=1), $
                          n_elements(data)/2)
              recon = ptr_new(img, /no_copy)
          endif
       end
       
       2: begin ;; 2D read and simple reconstruction
          if (opp->lookup('apptype') eq 'imEPI') then begin
              data = mas_varian_read_fid_epi(lun, opp)
          endif else if (opp->lookup('apptype') eq 'im2Depi') then begin
              data = mas_varian_read_fid_epip(lun, opp)
          endif else begin
              data = mas_varian_read_fid_data_2d(lun, opp)
          endelse
          
          if (arg_present(recon)) then begin
              img = fltarr(dims)
              dim1_shift = round(opp->lookup('ppe')/opp->lookup('lpe') * dims[1])
              for a=0, dims[3]-1 do begin
                  for s=0, dims[2]-1 do begin
                      img[*,*,s,a] = shift(abs(fft( (*data)[*,*,s,a] )), $
                                           dims[0]/2, dims[1]/2)
                      img[*,*,s,a] = shift(img[*,*,s,a], 0, -dim1_shift)
                  endfor
              endfor
              recon = ptr_new(img, /no_copy)
          endif
       end

       3: begin ;; 3D read and simple recon.
          data = mas_varian_read_fid_data_3d(lun, opp)
          if (arg_present(recon)) then begin
              img = fltarr(dims)
              dim1_shift = round(opp->lookup('ppe')/opp->lookup('lpe') * dims[1])
              dim2_shift = round(opp->lookup('ppe2')/opp->lookup('lpe2') * dims[2])
              for a=0, dims[3]-1 do begin
                  (*data)[*,*,*,a] = fft( (*data)[*,*,*,a], dimension=3 )
                  for s=0, dims[2]-1 do begin
                      img[*,*,s,a] = shift(abs( fft( (*data)[*,*,s,a] ) ), $
                                           dims[0]/2, dims[1]/2)
                      img[*,*,s,a] = shift(img[*,*,s,a], 0, -dim1_shift)
                  endfor
                  img[*,*,*,a] = shift(img[*,*,*,a], 0, 0, dims[2]/2) 
                  img[*,*,*,a] = shift(img[*,*,*,a], 0, 0, dim2_shift)
              endfor
              recon = ptr_new(img, /no_copy)
          endif
       end

       else: begin ;; not supported
          message, "Data must be 1, 2, or 3 dimensions"
       end
    endcase

    return, data
    
end

function mas_varian_fdf_parse_array, datatype, var_data

    temp_data = strmid(var_data, 1, strlen(var_data)-2)
    temp_data = strjoin(strsplit(temp_data, '"', /extract))
    
    case datatype of
        "char": begin
            array = strsplit(temp_data, ",", /extract)
            
            for i = 0, n_elements(array)-1 do begin
                array[i] = strtrim(array[i], 2)
            endfor
        end
        
        "int": begin
            array = long(strsplit(temp_data, ",", /extract))
         end
         
        "float": begin
            array = float(strsplit(temp_data, ",", /extract))
        end
    endcase
    
    return, array
end


function mas_varian_read_fdf_file, filename, header_only=header_only

    openr, lun, filename, /get_lun
    
    line      = ""
    magic_hdr = ""
    
    readf, lun, magic_hdr
    if (magic_hdr ne "/usr/local/fdf/startup") then begin
        ;; something is wrong with file
    endif
    
    struct = create_struct("magic", magic_hdr)
    read_ok = 1
    while (read_ok eq 1) do begin

        readf, lun, line
                
        test = stregex(line, '^(char|int|float)[ ]+([^ ]+)[ ]+=[ ]+([^;]+);$', $
                       /extract, /subexpr)
        if (n_elements(test) ne 4) then begin
            ;; something's wrong...
        endif
        
        datatype = test[1]
        var_name = test[2]
        var_data = test[3]

        for t = 0, n_tags(struct)-1 do begin
            if ((tag_names(struct))[t] eq strupcase(var_name)) then begin
                var_name = var_name + "_"
                break
            endif
        endfor
        
        case datatype of 
        
            "char": begin
                var_name = strmid(var_name, 1)
                array_pos = strpos(var_name, "[]")
                if (array_pos eq -1) then begin
                    ;; data is scalar
                    struct = create_struct(struct, var_name, strmid(var_data, 1, strlen(var_data)-2))
                endif else begin
                    ;; data is array
                    var_name = strmid(var_name, 0, strlen(var_name)-2)
                    substr = mas_varian_fdf_parse_array(datatype, var_data)
                    struct = create_struct(struct, var_name, substr)
                endelse
             end
             
            "int": begin
                array_pos = strpos(var_name, "[]")
                if (array_pos eq -1) then begin
                    ;; data is scalar
                    struct = create_struct(struct, var_name, long(var_data))
                endif else begin
                    ;; data is array
                    var_name = strmid(var_name, 0, strlen(var_name)-2)
                    substr = mas_varian_fdf_parse_array(datatype, var_data)
                    struct = create_struct(struct, var_name, substr)
                endelse
            end
            
            "float": begin
                array_pos = strpos(var_name, "[]")
                if (array_pos eq -1) then begin
                    ;; data is scalar
                    struct = create_struct(struct, var_name, float(var_data))
                endif else begin
                    ;; data is array
                    var_name = strmid(var_name, 0, strlen(var_name)-2)
                    substr = mas_varian_fdf_parse_array(datatype, var_data)
                    struct = create_struct(struct, var_name, substr)
                endelse
            end
            
            else: read_ok = 0 ;; end of text header
             
        endcase

    endwhile
    
    if (keyword_set(header_only)) then begin
        close, lun 
        free_lun, lun
        return, struct
    endif
    
    point_lun, lun, 0
    nullbyte = byte(0)
    tempbyte = byte(1)
    while (not eof(lun)) do begin
        readu, lun, tempbyte
        if (tempbyte eq nullbyte) then break
    endwhile

    case struct.storage of 
        'float': data_array = fltarr(struct.matrix)
        'integer': begin
            case struct.bits of
                8  : data_array = bytarr(struct.matrix)
                16 : data_array = intarr(struct.matrix)
                32 : data_array = lonarr(struct.matrix)
                64 : data_array = lon64arr(struct.matrix) ;; why?
                else: ;; error
            endcase
        end
    endcase
    
    readu, lun, data_array
    close, lun
    free_lun, lun
    
    data_array = reverse(data_array, 2)
    ori = struct.orientation
    
    struct = create_struct(struct, "img_data", ptr_new(data_array, /no_copy))
    
    return, struct
    
end


;; Usage: mas_varian_example, /path/to/fid/directory
;; Usage: mas_varian_example, dialog_pickfile(/directory)
;; Given an experiment directory (must contain a fid and procpar)
;; reconstruct kspace, and then very simply create images and
;; display the first volume in a window.
;; arg1: fid_dir - a string path to a DIRECTORY containing
;;                 a 'fid' and 'procpar' file.
pro mas_varian_example, fid_dir

    if (0 eq file_test(fid_dir, /directory)) then begin
       message, "Usage: mas_varian_example, /path/to/fid/directory.", $
                /informational
       return
    endif

    ps = path_sep()

    fid_file     = fid_dir+ps+'fid'
    procpar_file = fid_dir+ps+'procpar'

    opp = obj_new('mas_varian_procpar', procpar_file)
    
    ;; NOTE FID FILES ARE BIG ENDIAN!
    openu, lun, fid_file, /get_lun, /swap_if_little_endian, error=error
    
    if (error ne 0) then begin
       message, "Unable to open "+fid_file+". "
       return
    endif

    ;; data will be a pointer to kspace data
    ;; img will be pointer to reconstructed imagery
    data = mas_varian_read_fid_data_any(lun, opp, recon=img)

    close, lun
    free_lun, lun

    ;; display the slices from one volume in a very simple way
    window, /free
    dims = opp->computeDataDimensions()
    for s = 0, dims[2]-1 do begin
       tvscl, (*img)[*,*,s,0], s
    endfor

    ptr_free, data
    ptr_free, img
    obj_destroy, opp

end
