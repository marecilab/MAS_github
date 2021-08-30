;; $Id$
;;

; ------------------------------------------------------------------------------
;
; ft_track.pro
;
; The FT_TRACK procedure calculates the fiber tracts for a data set by solving
; the Frenet Equation using fourth order Runge-Kutta technique and stores the
; trajectory points .trs file.  Bilinear interpolation is used
; to interpolate the vector field.  Tract tracing is terminated with the help
; of a thresholding file.
;
; SYNTAX
;	ft_track, parameters, p_threshold, p_direction, ptrseeds
;
; ARGUMENTS
;  parameters
;	A structure carrying the information of the file path, name of the file
;	that the calculated tracts will be stored in, resolution factors in 3
;	dimensions, step size for tracking, threshold and length limits, and
;	density.
;
;  p_threshold
;	Pointer to the 3d thresholding file.
;
;  p_direction
;	Pointer to the principal vector field.
;
;  ptrseeds
;	The pointer as described in the information in ft_track routine.
;
; Created on October 18, 2001 by E. Ozarslan
; Last modified on September 25, 2003 by E. Ozarslan
;
; ------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; interpolating function

function average, p, r

    compile_opt idl2, hidden

    sum     = fltarr(3)
    gc      = floor(r)
    lengths = fltarr(3, 3)
    weights = fltarr(3, 3, 3)

    for i=0, 2 do begin
        if round(r[i]) eq gc[i] then begin
            lengths[i, 0]=0.5+float(gc[i])-r[i]
            lengths[i, 1]=r[i]+0.5-float(gc[i])
            lengths[i, 2]=0.0
        endif else begin
            lengths[i, 0]=0.0
            lengths[i, 1]=1.5+float(gc[i])-r[i]
            lengths[i, 2]=r[i]-0.5-float(gc[i])
        endelse
    endfor
    
    for kk=0,2 do begin
        for jj=0,2 do begin
            for ii=0,2 do begin
                
                weights[ii, jj, kk] = lengths[0, ii]*lengths[1, jj]*lengths[2, kk]
                aa = gc[0]-1+ii
                bb = gc[1]-1+jj
                cc = gc[2]-1+kk
                
                dproduct=0
                
                for c=0,2 do begin
                    dproduct = dproduct+(*p)[gc[0],gc[1],gc[2],c,0] * (*p)[aa,bb,cc,c,0]
                endfor
                
                if dproduct lt 0 then begin
                    (*p)[aa,bb,cc,*,0] = -(*p)[aa,bb,cc,*,0]
                endif
                
                for component=0,2 do begin
                    sum[component] = sum[component]+weights[ii,jj,kk]*(*p)[aa,bb,cc,component,0]
                endfor
                
            endfor
        endfor
    endfor
    
    return, sum
    
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Seeding function

function seeder, pv, w, sq

    compile_opt idl2, hidden
    
    sz_pv = size(*pv)
    n     = sz_pv[1]                  ; number of voxels

    if w eq 'x' or w eq 'y' or w eq 'z' then begin
        
        case w of
            'x': begin
                ii=1
                jj=2
            end
            'y': begin
                ii=0
                jj=2
            end
            'z': begin
                ii=0
                jj=1
            end
        endcase
        
        imin = min((*pv)[*,ii])
        imax = max((*pv)[*,ii])
        jmin = min((*pv)[*,jj])
        jmax = max((*pv)[*,jj])
        dimi = imax-imin+1
        dimj = jmax-jmin+1
        mask = bytarr(dimi, dimj)

        for i=0,dimi-1 do begin
            for j=0,dimj-1 do begin
                for t=0, n-1 do begin
                    if (*pv)[t,ii] eq imin+i and (*pv)[t,jj] eq jmin+j then $
                      mask[i,j]=1
                endfor
            endfor
        endfor

        lattice = make_array(ceil(dimi/sq)*ceil(dimj/sq),2)
        t = 0

        for i=0,ceil(dimi/sq)-1 do begin
            for j=0,ceil(dimj/sq)-1 do begin
                si = i*sq
                sj = j*sq
                if mask[floor(si),floor(sj)] then begin
                    lattice[t,*]=[imin+si,jmin+sj]
                    t=t+1
                endif
            endfor
        endfor

        seeds=make_array(t,3)
        
        case w of
            'x':begin
                seeds[*,0]=(*pv)[0,0]
                seeds[*,1]=lattice[0:t-1,0]
                seeds[*,2]=lattice[0:t-1,1]
            end
            'y':begin
                seeds[*,1]=(*pv)[0,1]
                seeds[*,0]=lattice[0:t-1,0]
                seeds[*,2]=lattice[0:t-1,1]
            end
            'z':begin
                seeds[*,2]=(*pv)[0,2]
                seeds[*,0]=lattice[0:t-1,0]
                seeds[*,1]=lattice[0:t-1,1]
            end
        endcase

    endif else begin            ; oblique plane
      
        xmin = min((*pv)[*,0])
        ymin = min((*pv)[*,1])
        zmin = min((*pv)[*,2])
        xmax = max((*pv)[*,0])
        ymax = max((*pv)[*,1])
        zmax = max((*pv)[*,2])
        dimx = xmax-xmin+1
        dimy = ymax-ymin+1
        dimz = zmax-zmin+1
        mask = bytarr(dimx, dimy, dimz)
        
        for i=0,dimx-1 do begin
            for j=0,dimy-1 do begin
                for k=0,dimz-1 do for t=0,n-1 do begin
                    if (*pv)[t,0] eq xmin+i and (*pv)[t,1] eq ymin+j and $
                      (*pv)[t,2] eq zmin+k then $
                      mask[i,j,k]=1
                endfor
            endfor
        endfor

        lattice = make_array(ceil(dimx/sq)*ceil(dimy/sq)*ceil(dimz/sq),3)
        t = 0
        for i=0,ceil(dimx/sq)-1 do begin
            for j=0,ceil(dimy/sq)-1 do begin
                for k=0,ceil(dimz/sq)-1 do begin
                    si = i*sq
                    sj = j*sq
                    sk = k*sq
                    if mask[floor(si),floor(sj),floor(sk)] then begin
                        lattice[t,*] = [xmin+si,ymin+sj,zmin+sk]
                        t=t+1
                    endif
                endfor
            endfor
        endfor

        seeds      = make_array(t,3)
        seeds[*,0] = lattice[0:t-1,0]
        seeds[*,1] = lattice[0:t-1,1]
        seeds[*,2] = lattice[0:t-1,2]

    endelse

    ptr_free, pv
    
    return, seeds
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Bookkeeping function

pro ft_journal, b, planes

  compile_opt idl2, hidden

	common share_fibers_main

  if b.tfile eq 0 then tfile='Fract_Anisotropy'
  if b.tfile eq 1 then tfile='dir_coherence'
  if b.tfile eq 2 then tfile='fran_coh'
  if b.dfile then dfile='coherence_dir' else dfile='Eigenvectors'

  sl=1/sqrt(b.sd)
  openw, 5, b.fname+'.txt'
  printf, 5, 'ROOT				:	'+b.path
  printf, 5, 'OUTPUT FILE			:	'+b.fname+'.trs'
  printf, 5, 'DIRECTIONS FILE			:	'+dfile+'.flt'
  printf, 5, 'THRESHOLD FILE			:	'+tfile+'.flt'
  printf, 5, 'THRESHOLD VALUE			:	'+strtrim(string(b.thr),1)
  printf, 5, 'STEP SIZE			:	'+strtrim(string(b.ss),1)
  printf, 5, 'SEED DENSITY			:	'+strtrim(string(b.sd),1)+' per voxels'
  printf, 5, 'DISTANCE BETWEEN SEEDS		:	'+strtrim(string(sl),1)+' voxels'
  printf, 5, ' '
  printf, 5, 'SEED PLANES			:'
  printf, 5, ' '
  printf, 5, '   	Normal Vector (Nx, Ny, Nz)   	'+$
				'	   Center Point (Rx, Ry, Rz)'
  printf, 5, planes
  free_lun, 5
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Tracking procedure

pro ft_track, b, pTHR, pDIR, ptr_seeds, seeds_from=seedsfile

    compile_opt idl2
    
    
    ;;print, b.resx, b.resy, b.resz
    res_fac_x=min([b.resx, b.resy, b.resz])/b.resx
    res_fac_y=min([b.resx, b.resy, b.resz])/b.resy
    res_fac_z=min([b.resx, b.resy, b.resz])/b.resz
    
    offset=5
    
    ;;evecs=rflt(dfile+'.flt')
    sz_data=size(*pDIR)
    ;;if sz_data[0] eq 5 then evecs=evecs[*,*,*,*,0]
    fdim=sz_data[1]
    pdim=sz_data[2]
    sdim=sz_data[3]
    evecs2=fltarr(fdim+2*offset, pdim+2*offset, sdim+2*offset,3)
    evecs2[offset:fdim+offset-1,offset:pdim+offset-1,offset:sdim+offset-1, *]=*pDIR
    p_ev=ptr_new(evecs2, /no_copy)
    fa=fltarr(fdim+2*offset, pdim+2*offset, sdim+2*offset)
    fa[offset:fdim+offset-1,offset:pdim+offset-1,offset:sdim+offset-1]=*pTHR

    ;; maximum # of paces, maximum length of tract is mu*b.ss voxels in each dir.
    mu=2*round(float(max([fdim,pdim,sdim]))/b.ss)

    if n_elements(seedsfile) then begin
        if file_test(seedsfile) eq 0 then return
        restore, seedsfile
        rrx=0
        rry=0
        rrz=0
        start_points=0
        mid_points=0
        matrix=0
        sz_seeds=size(seeds)
        cum_inc=sz_seeds[1]
    endif else begin
        ;; extracting the information pointed by the ptr_seeds
        pa=*ptr_seeds[0]
        planes=*ptr_seeds[1]
        sz_planes=size(*ptr_seeds[1])
        if sz_planes[0] then num_planes=1 else num_planes=sz_planes[2]
        sq=1/sqrt(b.sd)         ; step size
        seeds=0
        cum_inc=0
                                ; seed points are calculated and stored
        for p=0, num_planes-1 do begin
            if ptr_valid(pa[p]) then begin
                inc=0
                if planes[0,p] then w='x' else  $
                  if planes[1,p] then w='y' else $
                  if planes[2,p] then w='z' else w='o'
                seedsp=seeder(pa[p], w, sq)
                sz_seedsp=size(seedsp)
                inc=sz_seedsp[1]
                if cum_inc eq 0 then seeds=seedsp else seeds=[seeds, seedsp]
                cum_inc=cum_inc+inc
            endif
        endfor
        ptr_free, ptr_seeds
        ptr_free, pa
    endelse
    
    sz_seeds=size(seeds)
    
    if sz_seeds[0] ne 0 then begin
        
        x0=seeds[*,0]+offset
        y0=seeds[*,1]+offset
        z0=seeds[*,2]+offset
        trackn=cum_inc
        
        ;; arrays that will store the coordinates for the tracts
        rrx=fltarr(1)
        rry=fltarr(1)
        rrz=fltarr(1)
        
        ;; point counter
        ;;hcount=long(0)
        wbase=widget_base(title='Tracts Calculation in Progress...', $
                          /column, /base_align_center, scr_xsize=450, scr_ysize=80, ypad=10)
        
        wslider=widget_slider(wbase, maximum=trackn, $
                              frame=4, xsize=400)
        widget_control, wslider, set_value=trackn
        WIDGET_CONTROL, wbase, /REALIZE
        widget_control, /hourglass
        
                                ; arrays storing the start points and midpoint locations in the vertex arrays
        start_points=ulon64arr(1)
        mid_points=ulon64arr(1)
        
                                ; main loop
                                ; calculates fiber tracts in both directions for each seed point
        for ttt=0L, trackn-1 do begin
            
            
            rx_forward=fltarr(1)
            ry_forward=fltarr(1)
            rz_forward=fltarr(1)
            deltax=0
            deltay=0
            deltaz=0
            terminationf=0
            xx=x0[ttt]
            yy=y0[ttt]
            zz=z0[ttt]
            evec_forward=(*p_ev)[floor(xx),floor(yy),floor(zz),*]
            rx_forward[0]=xx
            ry_forward[0]=yy
            rz_forward[0]=zz
            hh=0
            
            
            while (terminationf eq 0) do begin
                
                ;; Runge-Kutta
;      k1=average(p_ev,[xx,yy,zz])*b.ss
;      k2=average(p_ev,[xx+k1[0]/2.,yy+k1[1]/2.,zz+k1[2]/2.])*b.ss
;      k3=average(p_ev,[xx+k2[0]/2.,yy+k2[1]/2.,zz+k2[2]/2.])*b.ss
;      k4=average(p_ev,[xx+k3[0],yy+k3[1],zz+k3[2]])*b.ss
;      step=fltarr(3)
;      step=(k1*.5+k2+k3+k4*.5)/3.

;      ; Runge-Kutta new
                
                k1=average(p_ev,[xx                   , yy                   , zz                   ])*b.ss
                k2=average(p_ev,[xx+k1[0]/2.*res_fac_x, yy+k1[1]/2.*res_fac_y, zz+k1[2]/2.*res_fac_z])*b.ss
                k3=average(p_ev,[xx+k2[0]/2.*res_fac_x, yy+k2[1]/2.*res_fac_y, zz+k2[2]/2.*res_fac_z])*b.ss
                k4=average(p_ev,[xx+k3[0]   *res_fac_x, yy+k3[1]   *res_fac_y, zz+k3[2]   *res_fac_z])*b.ss
                
                step=fltarr(3)
                step=(k1*.5+k2+k3+k4*.5)/3.
                step=step*[res_fac_x, res_fac_y, res_fac_z]
                
                if step[0] eq 0. and step[1] eq 0. and step[2] eq 0. then terminationf=1
                
                                ; threshold check
                if (fa[floor(xx), floor(yy), floor(zz)] lt b.thr) then terminationf=1
                
                                ; dot product calculation
                dp=0
                for comp=0,2 do $
                  dp=dp+(*p_ev)[floor(xx),floor(yy),floor(zz),comp]* $
                  (*p_ev)[floor(xx+step[0]),floor(yy+step[1]),floor(zz+step[2]),comp]
                
                                ; direction correction
                if dp lt 0 then $
                  (*p_ev)[floor(xx+step[0]),floor(yy+step[1]),floor(zz+step[2]),*]=$
                  -(*p_ev)[floor(xx+step[0]),floor(yy+step[1]),floor(zz+step[2]),*]
                
                                ; stepping forwards
                xx=xx+step[0]
                yy=yy+step[1]
                zz=zz+step[2]
                
                                ; new coordinates stored
                rx_forward=[rx_forward, xx]
                ry_forward=[ry_forward, yy]
                rz_forward=[rz_forward, zz]
                
                hh=hh+1
                
                                ; loop check
                if hh gt mu-1 then terminationf=1 ; loops
                
            endwhile
            
            
            rx_backward=fltarr(1)
            ry_backward=fltarr(1)
            rz_backward=fltarr(1)
            terminationb=0
            xx=x0[ttt]
            yy=y0[ttt]
            zz=z0[ttt]
            rx_backward[0]=xx
            ry_backward[0]=yy
            rz_backward[0]=zz
            hhh=0
            (*p_ev)[floor(xx),floor(yy),floor(zz),*]=-evec_forward[*]
            
            while (terminationb eq 0) do begin
                
;      ; Runge-Kutta
;      k1=average(p_ev,[xx,yy,zz])*b.ss
;      k2=average(p_ev,[xx+k1[0]/2.,yy+k1[1]/2.,zz+k1[2]/2.])*b.ss
;      k3=average(p_ev,[xx+k2[0]/2.,yy+k2[1]/2.,zz+k2[2]/2.])*b.ss
;      k4=average(p_ev,[xx+k3[0],yy+k3[1],zz+k3[2]])*b.ss
;      step=fltarr(3)
;      step=(k1*.5+k2+k3+k4*.5)/3.


	  ; Runge-Kutta new
                k1=average(p_ev,[xx,yy,zz])*b.ss
                k2=average(p_ev,[xx+k1[0]/2.*res_fac_x,yy+k1[1]/2.*res_fac_y,zz+k1[2]/2.*res_fac_z])*b.ss
                k3=average(p_ev,[xx+k2[0]/2.*res_fac_x,yy+k2[1]/2.*res_fac_y,zz+k2[2]/2.*res_fac_z])*b.ss
                k4=average(p_ev,[xx+k3[0]*res_fac_x,yy+k3[1]*res_fac_y,zz+k3[2]*res_fac_z])*b.ss
                step=fltarr(3)
                step=(k1*.5+k2+k3+k4*.5)/3.
                step=step*[res_fac_x, res_fac_y, res_fac_z]
                


      ; eigenvector check
                if step[0] eq 0. and step[1] eq 0. and step[2] eq 0. then terminationb=1
                
                                ; threshold check
                if fa[floor(xx), floor(yy), floor(zz)] lt b.thr then terminationb=1
                
                                ; dotproduct calculation
                dp=0
                for comp=0,2 do $
                  dp=dp+(*p_ev)[floor(xx),floor(yy),floor(zz),comp]* $
                  (*p_ev)[floor(xx+step[0]),floor(yy+step[1]),floor(zz+step[2]),comp]
                
                                ; direction correction
                if (dp lt 0) then $
                  (*p_ev)[floor(xx+step[0]),floor(yy+step[1]),floor(zz+step[2]),*]=$
                  -(*p_ev)[floor(xx+step[0]),floor(yy+step[1]),floor(zz+step[2]),*]
                
                                ; stepping backwards
                xx=xx+step[0]
                yy=yy+step[1]
                zz=zz+step[2]
                
                                ; new positions stored
                rx_backward=[rx_backward, xx]
                ry_backward=[ry_backward, yy]
                rz_backward=[rz_backward, zz]
                
                hhh=hhh+1
                
                                ; loops check
                if hhh gt mu-1 then terminationb=1 ; loops
                
            endwhile
            
            
            if ttt eq 0 then begin
                rrx=[rx_forward, rx_backward]
                rry=[ry_forward, ry_backward]
                rrz=[rz_forward, rz_backward]
                start_points[0]=0
                mid_points[0]=hh+1
                start_points=[start_points, hh+hhh+2]
            endif else begin
                rrx=[rrx,rx_forward, rx_backward]
                rry=[rry,ry_forward, ry_backward]
                rrz=[rrz,rz_forward, rz_backward]
                mid_points=[mid_points, start_points[n_elements(start_points)-1]+hh+1]
                start_points=[start_points, mid_points[n_elements(mid_points)-1]+hhh+1]
            endelse
            
                                ;hcount=hcount+hhh+hh+1
            
                                ; information written on the terminal for each tract
            print, fix(trackn-fix(1+ttt)),'  ', x0[ttt]-offset,'  ', y0[ttt]-offset, '  ', z0[ttt]-offset
            
            widget_control, wslider, set_value=trackn-ttt-1
            
        endfor
        
        widget_control, wbase, /destroy
        
        rrx=(rrx-offset > 0.) < (float(fdim)-.001)
        rry=(rry-offset > 0.) < (float(pdim)-.001)
        rrz=(rrz-offset > 0.) < (float(sdim)-.001)
        matrix=[fdim,pdim,sdim]
        save, filename=b.fname+'.trs', rrx, rry, rrz, start_points, mid_points, $
          matrix, seeds, planes
        
        ft_journal, b,  planes
        
        
        ft_visualize, b
        
        rename=dialog_message(['Would you like to rename the files?',$
                               'Current '+ 'name is: '+b.fname], /default_no, /question)
        if rename eq 'Yes' then begin
            newname=dialog_pickfile(path='')
            
            if newname ne '' then begin
                file_move, b.fname+'.txt', newname+'.txt', /allow_same, /overwrite
                file_move, b.fname+'.trs', newname+'.trs', /allow_same, /overwrite
            endif
            
        endif
        
    endif else ok=dialog_message('No ROI selected.')
    
    ptr_free, p_ev
    
end



; New ft_track
; Changes the input parameters so that b is no longer sent but rather
; resx, resy, resz.
; Slider is now a progressbar

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; New Tracking procedure

pro ft_track_newversion, ptr_seeds, seeds_from=seedsfile

    common share_fibers_main
    common common_widgets
    compile_opt idl2
    
    print, resx, resy, resz
    res_fac_x=min([resx, resy, resz])/resx
    res_fac_y=min([resx, resy, resz])/resy
    res_fac_z=min([resx, resy, resz])/resz
    
    offset=5

    sz_data=size(*pDIR)

    fdim=sz_data[1]
    pdim=sz_data[2]
    sdim=sz_data[3]
    evecs2=fltarr(fdim+2*offset, pdim+2*offset, sdim+2*offset,3)

    evecs2[offset:fdim+offset-1,offset:pdim+offset-1,offset:sdim+offset-1, *]=*pDIR
    p_ev=ptr_new(evecs2, /no_copy)

    fa=fltarr(fdim+2*offset, pdim+2*offset, sdim+2*offset)
    fa[offset:fdim+offset-1,offset:pdim+offset-1,offset:sdim+offset-1]=*pTHR

    ;; maximum # of paces, maximum length of tract is mu*b.ss voxels in each dir.
    mu=2*round(float(max([fdim,pdim,sdim]))/ss)
    
    if n_elements(seedsfile) then begin

        if file_test(seedsfile) eq 0 then return
        
        restore, seedsfile
        rrx=0
        rry=0
        rrz=0
        start_points=0
        mid_points=0
        matrix=0
        sz_seeds=size(seeds)
        cum_inc=sz_seeds[1]
        
    endif else begin
        ;; extracting the information pointed by the ptr_seeds
        
        pa=*ptr_seeds[0]
        planes=*ptr_seeds[1]
        
        sz_planes=size(*ptr_seeds[1])
        if sz_planes[0] then num_planes=1 else num_planes=sz_planes[2]
        
        ;; Spawning a new window with the ROI information for later reference
        ; HS 20070405
        base_xsize = 250
        base_ysize = 100
        temptext_slicedir = ''
        temptext = 'Initial slice is #'
        temptext2 = 'New slice is #'
        
        wbase = widget_base(title='Selected Slices for ROI', /ALIGN_LEFT , group_leader=WID_BASE_MAIN, $
                            xsize=base_xsize, ysize=base_ysize, XOFFSET=800)

        
        TEXT_file_array =strarr(num_planes+1)

        for p=0, num_planes-1 do begin
            
            if (planes[0,p] eq 0) and (planes[1,p] eq 0) and (planes[2,p] eq 1) then begin
                curslice = round(planes[5,p])
                temptext_slicedir = 'Z'
                
            endif else if (planes[0,p] eq 0) and (planes[1,p] eq 1) and (planes[2,p] eq 0) then begin
                curslice = round(planes[4,p])
                temptext_slicedir = 'Y'
            endif else if (planes[0,p] eq 1) and (planes[1,p] eq 0) and (planes[2,p] eq 0) then begin
                curslice = round(planes[3,p])
                temptext_slicedir = 'X'
            endif
            
            
            temptext = temptext  + STRCOMPRESS(string(curslice), /REMOVE_ALL)
            temptext2= temptext2 + STRCOMPRESS(string(curslice*fix(rfx)), /REMOVE_ALL)
            
            TEXT_file_array[p] = 'In ' + temptext_slicedir + ' plane: ' + temptext + ' - ' + temptext2
            
            temptext_slicedir = ''
            temptext = 'Initial slice is #'
            temptext2 = 'New slice is #'
            
        endfor
        
        TEXT_file_textbox = widget_text(wbase, VALUE=TEXT_file_array, /SCROLL, $
                                        scr_xsize=base_xsize, scr_ysize=base_ysize)
        
        WIDGET_CONTROL, wbase, /REALIZE
        
        
        sq=1/sqrt(sd)           ; step size
        seeds=0
        cum_inc=0
                                ; seed points are calculated and stored
        for p=0, num_planes-1 do begin
            
            if ptr_valid(pa[p]) then begin
                inc=0
                
                if planes[0,p] then w='x' else  $
                  if planes[1,p] then w='y' else $
                  if planes[2,p] then w='z' else w='o'

                seedsp=seeder(pa[p], w, sq)
                sz_seedsp=size(seedsp)
                inc=sz_seedsp[1]
                if cum_inc eq 0 then seeds=seedsp else seeds=[seeds, seedsp]
                cum_inc=cum_inc+inc
            endif
        endfor
        ptr_free, ptr_seeds
        ptr_free, pa
    endelse

    sz_seeds=size(seeds)
    
    if sz_seeds[0] ne 0 then begin
        
        x0=seeds[*,0]+offset
        y0=seeds[*,1]+offset
        z0=seeds[*,2]+offset
        trackn=cum_inc
        
        ;; arrays that will store the coordinates for the tracts
        rrx=fltarr(1)
        rry=fltarr(1)
        rrz=fltarr(1)

        ;; HS 20070123
        ;; Exchanging the old widget slider countdown with the standard progress bar.
        progressbar = Obj_New('progressbar', Color='red', Text='Tracts Calculations in Progress...');,/NOCANCEL)
        progressbar -> Start

        catch, error_status
        
        if (error_status ne 0) then begin
           
           junk = dialog_message('An error occurred that prevented fiber tracking from continuing.')
           progressbar->destroy
           
           catch, /cancel
           return
           
        endif

        ;; arrays storing the start points and midpoint locations in the vertex arrays
        start_points=ulon64arr(1)
        mid_points=ulon64arr(1)

        print, 'Number of Seed Points = '+strcompress(string(trackn-1), /remove_all)
        dp_correct = 1
;        print, 'cos(theta) limit = '+string(angle_limit)+' '+string(acos(angle_limit)*180/!PI)
        show_cos = 0

        ;; main loop
        ;; calculates fiber tracts in both directions for each seed point
        for ttt=0L, trackn-1 do begin

            rx_forward=fltarr(1)
            ry_forward=fltarr(1)
            rz_forward=fltarr(1)
            deltax=0
            deltay=0
            deltaz=0
            terminationf=0
            xx=x0[ttt]
            yy=y0[ttt]
            zz=z0[ttt]
            evec_forward=(*p_ev)[floor(xx),floor(yy),floor(zz),*]
            rx_forward[0]=xx
            ry_forward[0]=yy
            rz_forward[0]=zz
            hh=0

            while (terminationf eq 0) do begin

                ;; Runge-Kutta new
                k1 = average(p_ev,[xx,                    yy,                    zz                   ])*ss
                k2 = average(p_ev,[xx+k1[0]/2.*res_fac_x, yy+k1[1]/2.*res_fac_y, zz+k1[2]/2.*res_fac_z])*ss
                k3 = average(p_ev,[xx+k2[0]/2.*res_fac_x, yy+k2[1]/2.*res_fac_y, zz+k2[2]/2.*res_fac_z])*ss
                k4 = average(p_ev,[xx+k3[0]   *res_fac_x, yy+k3[1]   *res_fac_y, zz+k3[2]   *res_fac_z])*ss
                step = fltarr(3)
                step = (k1 * .5 + k2 + k3 + k4 * .5) / 3.
                step = step*[res_fac_x, res_fac_y, res_fac_z]
                
                if step[0] eq 0. and step[1] eq 0. and step[2] eq 0. then terminationf=1

                ;; threshold check
                if (fa[floor(xx), floor(yy), floor(zz)] lt thr) then terminationf=1

                ;; dot product calculation
                dp = transpose(reform((*p_ev)[floor(xx),floor(yy),floor(zz),*])) # $
                  reform((*p_ev)[floor(xx+step[0]),floor(yy+step[1]),floor(zz+step[2]),*])

                ;; direction correction
                if (dp_correct and dp lt 0) then begin
                    (*p_ev)[floor(xx+step[0]),floor(yy+step[1]),floor(zz+step[2]),*]=$
                      -(*p_ev)[floor(xx+step[0]),floor(yy+step[1]),floor(zz+step[2]),*]
                endif

                mag_1 = sqrt(total( ( (*p_ev)[floor(xx),floor(yy),floor(zz),*] )^2 ))
                mag_2 = sqrt(total( ( (*p_ev)[floor(xx+step[0]),floor(yy+step[1]),floor(zz+step[2]),*] )^2 ))
                cos_theta = dp/(mag_1*mag_2)
                if (abs(cos_theta) lt angle_limit) then terminationf=1

                ;; stepping forwards
                xx=xx+step[0]
                yy=yy+step[1]
                zz=zz+step[2]
                
                ;; new coordinates stored
                rx_forward=[rx_forward, xx]
                ry_forward=[ry_forward, yy]
                rz_forward=[rz_forward, zz]

                hh=hh+1
                
                ;; loop check
                if hh gt mu-1 then terminationf=1 ; loops
                
            endwhile

            
            rx_backward=fltarr(1)
            ry_backward=fltarr(1)
            rz_backward=fltarr(1)
            terminationb=0
            xx=x0[ttt]
            yy=y0[ttt]
            zz=z0[ttt]
            rx_backward[0]=xx
            ry_backward[0]=yy
            rz_backward[0]=zz
            hhh=0
            (*p_ev)[floor(xx),floor(yy),floor(zz),*]=-evec_forward[*]
            
            while (terminationb eq 0) do begin

                ;; Runge-Kutta new
                k1=average(p_ev,[xx,yy,zz])*ss
                k2=average(p_ev,[xx+k1[0]/2.*res_fac_x,yy+k1[1]/2.*res_fac_y,zz+k1[2]/2.*res_fac_z])*ss
                k3=average(p_ev,[xx+k2[0]/2.*res_fac_x,yy+k2[1]/2.*res_fac_y,zz+k2[2]/2.*res_fac_z])*ss
                k4=average(p_ev,[xx+k3[0]*res_fac_x,yy+k3[1]*res_fac_y,zz+k3[2]*res_fac_z])*ss
                step=fltarr(3)
                step=(k1*.5+k2+k3+k4*.5)/3.
                step=step*[res_fac_x, res_fac_y, res_fac_z]
                
                ;; eigenvector check
                if step[0] eq 0. and step[1] eq 0. and step[2] eq 0. then terminationb=1

                ;; threshold check
                if fa[floor(xx), floor(yy), floor(zz)] lt thr then terminationb=1

                ;; dotproduct calculation
                dp = transpose(reform((*p_ev)[floor(xx),floor(yy),floor(zz),*])) # $
                  reform((*p_ev)[floor(xx+step[0]),floor(yy+step[1]),floor(zz+step[2]),*])

                ;; direction correction
                if (dp_correct and dp lt 0) then begin
                    (*p_ev)[floor(xx+step[0]),floor(yy+step[1]),floor(zz+step[2]),*]=$
                      -(*p_ev)[floor(xx+step[0]),floor(yy+step[1]),floor(zz+step[2]),*]
                endif

                mag_1 = sqrt(total( ( (*p_ev)[floor(xx),floor(yy),floor(zz),*] )^2 ))
                mag_2 = sqrt(total( ( (*p_ev)[floor(xx+step[0]),floor(yy+step[1]),floor(zz+step[2]),*] )^2 ))
                cos_theta = dp/(mag_1*mag_2)
                if (abs(cos_theta) lt angle_limit) then terminationb=1

                ;; stepping backwards
                xx=xx+step[0]
                yy=yy+step[1]
                zz=zz+step[2]
                
                ;; new positions stored
                rx_backward=[rx_backward, xx]
                ry_backward=[ry_backward, yy]
                rz_backward=[rz_backward, zz]

                hhh=hhh+1

                ;; loops check
                if hhh gt mu-1 then terminationb=1 ; loops

            endwhile

            if ttt eq 0 then begin
                rrx=[rx_forward, rx_backward]
                rry=[ry_forward, ry_backward]
                rrz=[rz_forward, rz_backward]
                start_points[0]=0
                mid_points[0]=hh+1
                start_points=[start_points, hh+hhh+2]
            endif else begin
                rrx=[rrx,rx_forward, rx_backward]
                rry=[rry,ry_forward, ry_backward]
                rrz=[rrz,rz_forward, rz_backward]
                mid_points=[mid_points, start_points[n_elements(start_points)-1]+hh+1]
                start_points=[start_points, mid_points[n_elements(mid_points)-1]+hhh+1]
            endelse
            
            progressbartemp = float(trackn)

; Explanation of progressbar formula
;	initial_p;fnrogressbar_step is 1 / (number of fibers to track)
;	the number of fibers to track - the current fiber the program is divided by the number of fibers to track
;	this is substracted to 1 (to yield a progress rather than a countdown).
;	then we sum the initial step (initial_progressbar_step) so that it reaches 100% (if not
;	it would stay at 1-initial_progressbar_step

            initial_progressbar_step = 1/progressbartemp
            temp = (1-(progressbartemp-ttt)/progressbartemp)+initial_progressbar_step
            
            progressbar -> Update, temp*100
            
            IF progressBar -> CheckCancel() THEN BEGIN
                progressBar->Destroy
                ptr_free, p_ev
                ptr_free, pDIR
                ptr_free, PS_0
                RETURN
            ENDIF

        endfor
        
        catch, /cancel

        rrx=(rrx-offset > 0.) < (float(fdim)-.001)
        rry=(rry-offset > 0.) < (float(pdim)-.001)
        rrz=(rrz-offset > 0.) < (float(sdim)-.001)
        matrix=[fdim,pdim,sdim]
        
        fname = fname + '.trs'
        
        save, filename=fname, rrx, rry, rrz, start_points, mid_points, $
          matrix, seeds, planes
        
	; Im destroying the progressbar at this point to give it a bit of time to show
	; that it reached 100%.	Doesn't always show it because the save function is fast.
        progressbar -> Destroy

        ft_journal_newversion, planes

;;        ft_visualize_newversion

;;        fibers = make_fibers(fname)
        make_fibers, fname, fiber_data=fibers, transform=invert(transform_matrix)

        view_fibers, fibers, /show_axis, alpha=.5, thick=1, project_dataset=project_dataset, /in_mas

    endif else ok=dialog_message('No ROI selected.', /center)
    
    ptr_free, p_ev
    ptr_free, pDIR
    ptr_free, PS_0
    
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Bookkeeping function

pro ft_journal_newversion, planes

    compile_opt idl2, hidden
    common share_fibers_main
    on_error, 2

    path_to_write = fname
    
    lastslash = STRPOS(path_to_write,'.',REVERSE_SEARCH=1)
    path_to_write = STRMID(fname,0,lastslash)
    
    sl=1/sqrt(sd)
    openw, 5, path_to_write+'.txt'
    printf, 5, 'ROOT				:	'+path
    printf, 5, 'OUTPUT FILE			:	'+fname+'.trs'
    printf, 5, 'DIRECTIONS FILE			:	'+dfile+'.flt'
    printf, 5, 'THRESHOLD FILE			:	'+tfile+'.flt'
    printf, 5, 'THRESHOLD VALUE			:	'+strtrim(string(thr),1)
    printf, 5, 'STEP SIZE			:	'+strtrim(string(ss),1)
    printf, 5, 'SEED DENSITY			:	'+strtrim(string(sd),1)+' per voxels'
    printf, 5, 'DISTANCE BETWEEN SEEDS		:	'+strtrim(string(sl),1)+' voxels'
    printf, 5, ' '
    printf, 5, 'SEED PLANES			:'
    printf, 5, ' '
    printf, 5, '   	Normal Vector (Nx, Ny, Nz)   	'+$
      '	   Center Point (Rx, Ry, Rz)'
    printf, 5, planes
    free_lun, 5
end
