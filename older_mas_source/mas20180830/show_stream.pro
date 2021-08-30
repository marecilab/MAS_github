; $Id: //depot/Release/IDL_81/idl/idldir/examples/doc/objects/show_stream.pro#1 $
;
; Copyright (c) 1999-2011, ITT Visual Information Solutions. All
;       rights reserved.
;+
; NAME:
; SHOW_STREAM
;
; PURPOSE:
; This procedure illustrates the use of PARTICLE_TRACE, STREAMLINE,
;       and VECTOR_FIELD to display 3D vector data.
;
; CATEGORY:
; Graphics.
;
; CALLING SEQUENCE:
;
; SHOW_STREAM, u[, v, w]] [, SEEDS=seeds] $
;           [, /LINES | /TUBES | /VECTOR] $
;           [/ARROWS]
;
;
; INPUTS:
; u:  An array providing the u components of the direction
;               vectors.  If the v argument is not provided, the u argument
;               must be a 3xNXxNYxNZ array specifying the u, v, and
;               w components of the vector data for a volume with
;               dimensions NXxNYxNZ.
;
; OPTIONAL INPUTS:
; v:  An array providing the v components of the direction
;               vectors.  The dimensions of this argument must match
;               the dimensions of the u argument.
; w:  An array providing the w components of the direction
;               vectors.  The dimensions of this argument must match
;               the dimensions of the u argument.
; 
; KEYWORD PARAMETERS:
; SEEDS:  Set this keyword to an array of [x,y,z] locations 
;               ([3,nSeeds]) to be used as the seed points for the 
;               particle trace or as initial positions for the vectors 
;               of the vector field (if the /VECTOR keyword is set).
; LINES:  Set this keyword to a non-zero value to indicate that
;               the particle trace should be rendered as a set of 
;               polylines.  The default is to display the particle
;               trace using ribbon streamlines.
; TUBES:  Set this keyword to a non-zero value to indicate that
;               the particle trace should be rendered as a set of 
;               tubes.  The default is to display the particle
;               trace using ribbon streamlines.
;       VECTOR: Set this keyword to a non-zero value to indicate that
;               a vector-field plot should be displayed.  The default
;               is to display a particle trace using ribbon streamlines.
;       ARROWS: Set this keyword to indicate that arrow heads should
;               be displayed.  This keyword is ignored if the TUBES
;               keyword is set.
;
; MODIFICATION HISTORY:
;   Written by: KB, Sept. 1999.
; July, 2002  Add support for the ARROWS keyword.
;     Update to use XOBJVIEW.
;-

;-----------------------------------------------------------------------
function makeArrowHeadSymbol, SIZE=size
    revolutionType = 6
    shape = [[-1,-0.5,0],[0,0,0]]
    nFacets = 16
    rotVec = [1,0,0]
    MESH_OBJ, revolutionType, verts, conn, shape, $
        P1=nFacets, P3=rotVec 

    oPoly = OBJ_NEW('IDLgrPolygon', verts, POLY=conn)
    oSymbol = OBJ_NEW('IDLgrSymbol', oPoly, SIZE=size)

    return, oSymbol
end

;-----------------------------------------------------------------------
function countPolys, inconn
    nElements = N_ELEMENTS(inconn)
    nPolys = 0
    i=0
    while (i lt nElements) do begin
        nVerts = inconn[i]
        if (nVerts eq 0) then continue
        if (nVerts eq -1) then break
        i = i + nVerts + 1
        nPolys = nPolys + 1
    endwhile
   
    return, nPolys
end

;-----------------------------------------------------------------------
pro SHOW_STREAM, u, v, w, SEEDS = seeds, LINES = lines, TUBES = tubes, $
                VECTOR = vector, ARROWS = arrows, anisotropy = anisotropy, image = image, zpos = zpos

    if(N_ELEMENTS(u) ne 0) then begin ;data supplied, check it
        nudims = SIZE(u, /N_DIMENSIONS)
        udims = SIZE(u, /DIMENSIONS)
        if(nudims eq 4) then begin
            if(udims[0] ne 3) then MESSAGE, 'Input must be 3D array of' + $
                ' 3-vectors or 3 3D arrays of scalars.'
            data = u            
            nx = udims[1] &  ny = udims[2] & nz = udims[3]
        end else if(nudims eq 3) then begin
            if((N_ELEMENTS(v) ne 0) and (N_ELEMENTS(w) ne 0)) then begin   
                nx = udims[0] &  ny = udims[1] & nz = udims[2]
            end else MESSAGE, 'Input must be 3D array of' + $
                ' 3-vectors or 3 3D arrays of scalars.'
            data = FLTARR(3,nx, ny, nz)
            data[0, *, *, *] = u
            data[1, *, *, *] = v
            data[2, *, *, *] = w        
        end
    end else begin ;data not supplied, compute helical flow test data.
        nx = 15 & ny = 15 & nz = 15
        data = FLTARR(3,nx, ny, nz)
        b = .034*7 & c = .14*7
        for i=0,nx-1 do $
            for j=0,ny-1 do $
            for k=0,nz-1 do begin
            x = i-nx/2
            y = j-ny/2
            z = k-nz/2
            data[0,i,j,k] = FLOAT(-b*y)
            data[1,i,j,k] = FLOAT(b*x)
            data[2,i,j,k] = FLOAT(c)
        end
    end

    if(N_ELEMENTS(seeds) eq 0) then begin ;Compute seed points.
        ;Set seed spacing for a rake at the bottom of the grid.
        xstep=LONG(nx/4) & ystep=LONG(ny/4)
        if(not KEYWORD_SET(vector)) then zstep=nz else zstep = 1

        nseeds = 3*LONG((nx*ny*nz)/(xstep*ystep*zstep))
        seeds = FLTARR(nseeds)
        iseed=0L
        for i=0,nx-1 do $
            for j=0,ny-1 do $
            for k=0,nz-1 do begin
            if( ((k mod zstep) eq 0) and ((i mod xstep) eq 0) and $
                ((j mod ystep) eq 0) and (iseed lt (nseeds-2)) ) then begin
                seeds[iseed] = FLOAT(i)
                seeds[iseed+1] = FLOAT(j)
                seeds[iseed+2] = FLOAT(k)
                iseed = iseed+3
            end
        end
    end

    maxIterations=100
    stepSize=.5
    width=.5 ;ribbon thickness

    ;Create streamlines graphic.
    oModel = OBJ_NEW('IDLgrModel')
    if not KEYWORD_SET(anisotropy) then anisotropy = [1,1,1]
    if(not KEYWORD_SET(vector)) then begin ;Streamlines/ribbons/tubes
        PARTICLE_TRACE,data,seeds,outverts,outconn,outnormals, $
            MAX_ITERATIONS=maxIterations, MAX_STEPSIZE=stepSize,  $
            INTEGRATION=0,ANISOTROPY=anisotropy, SEED_NORMAL=[0, 0, 1]

        if (outconn[0] eq -1l) then $
            MESSAGE, 'No particle trace.'

        maxDim = MAX(outverts[0,*]) > $
            MAX(outverts[1,*]) > $
            MAX(outverts[2,*])

        ; If requested, prepare arrow head labels.
        if (KEYWORD_SET(arrows) and not KEYWORD_SET(tubes)) then begin 
            symScale = (KEYWORD_SET(lines) ? 0.03 : 0.07)
            oSymbol = makeArrowHeadSymbol(SIZE=maxDim*symScale)
            nPolys = countPolys(outconn)
            ; Two arrow heads per particle trace path.
            lblPolys = LINDGEN(nPolys*2) / 2
            oLblSymbols = REPLICATE(oSymbol, nPolys*2)
            evens = LINDGEN(nPolys) * 2
            odds = evens+1
            lblOffsets = FLTARR(nPolys*2)
            lblOffsets[evens] = 0.5
            lblOffsets[odds] = 1.0
        endif

        if(KEYWORD_SET(lines)) then begin ;lines
            oStreamlines = OBJ_NEW('IDLgrPolyline',outverts, $
                POLYLINES=outconn, $
                LABEL_OBJECTS=oLblSymbols, $
                LABEL_POLYLINES=lblPolys, $
                LABEL_OFFSETS=lblOffsets, $
                /LABEL_USE_VERTEX_COLOR, $
                /LABEL_NOGAPS)
            oModel->Add, oStreamlines 
            title = 'Particle Trace'
        endif else begin ;ribbons/tubes
            ; If arrows are to be added, keep a copy of the particle trace 
            ; vertex and connectivity.
            if (KEYWORD_SET(arrows) and not KEYWORD_SET(tubes)) then begin 
                averts = outverts
                aconn = outconn
            endif

            if(KEYWORD_SET(tubes)) then $ ;square profile for stream-tubes.
                profile = [[-1,-1],[-1,1],[1,1],[1,-1],[-1,-1]]
            nverts = N_ELEMENTS(outverts)/3
            STREAMLINE, TEMPORARY(outverts),TEMPORARY(outconn), $
                outnormals*width,outverts,outconn, PROFILE=profile,ANISOTROPY=anisotropy, rgb_table =33, auto_range = [-1,1]

            oStreamlines = OBJ_NEW('IDLgrPolygon',outverts, POLYGONS=outconn, $
                                  SHADING = 1)
            oModel->Add, oStreamlines 

            if (KEYWORD_SET(arrows) and not KEYWORD_SET(tubes)) then begin 
                oArrows = OBJ_NEW('IDLgrPolyline',averts, $
                    POLYLINES=aconn, $
                    LINESTYLE=6, $
                    LABEL_OBJECTS=oLblSymbols, $
                    LABEL_POLYLINES=lblPolys, $
                    LABEL_OFFSETS=lblOffsets, $
                    /LABEL_USE_VERTEX_COLOR, $
                    /LABEL_NOGAPS)
                oModel->Add,oArrows
            endif
            title = 'Streamline'
        end
    end else begin ;Hedgehog vector plot
        VECTOR_FIELD,data,outverts,outconn, VERTICES=seeds,ANISOTROPY=anisotropy, scale = 0.1

        if (KEYWORD_SET(arrows)) then begin 
            maxDim = MAX(outverts[0,*]) > $
                MAX(outverts[1,*]) > $
                MAX(outverts[2,*])
            oSymbol = makeArrowHeadSymbol(SIZE=maxDim*0.03)
            nPolys = countPolys(outconn)
            lblPolys = LINDGEN(nPolys)
            oLblSymbols = REPLICATE(oSymbol, nPolys)
            lblOffsets = REPLICATE(1.0, nPolys)
        endif

        oStreamlines=OBJ_NEW('IDLgrPolyline',outverts,POLYLINES=outconn, $
            COLOR=[255,255,0], $
            LABEL_OBJECTS=oLblSymbols, $
            LABEL_POLYLINES=lblPolys, $
            LABEL_OFFSETS=lblOffsets, $
            /LABEL_USE_VERTEX_COLOR, $
            /LABEL_NOGAPS)

        oModel->Add, oStreamlines
        title = 'Vector Field'
    end

    ;Compute velocity magnitude
    magdata = SQRT(data[0,*, *]^2 + data[1,*, *]^2 + data[2,*, *]^2)
    ;Interpolate velocity magnitude at streamline vertices, and
    ;use values to color streamlines.
    vertX =  REFORM(outverts[0,*],N_ELEMENTS(outverts)/3)
    vertY =  REFORM(outverts[1,*],N_ELEMENTS(outverts)/3)
    vertZ =  REFORM(outverts[2,*],N_ELEMENTS(outverts)/3)
    vertcolors = BYTSCL(INTERPOLATE(magdata,vertX, vertY, vertZ))
    oPalette = OBJ_NEW('IDLgrPalette')
    oPalette->LOADCT, 33
    oStreamlines->SetProperty, PALETTE = oPalette, VERT_COLORS = vertcolors
    if (OBJ_VALID(oArrows)) then begin
        oArrows->GetProperty, DATA=averts
        vertX =  REFORM(averts[0,*],N_ELEMENTS(averts)/3)
        vertY =  REFORM(averts[1,*],N_ELEMENTS(averts)/3)
        vertZ =  REFORM(averts[2,*],N_ELEMENTS(averts)/3)
        vertcolors = BYTSCL(INTERPOLATE(magdata,vertX, vertY, vertZ))
        oArrows->SetProperty, PALETTE = oPalette, VERT_COLORS = vertcolors
    endif
    
;    if keyword_set(image) then begin
;      oShader = Obj_New('IDLgrShader')
;      oImage = OBJ_NEW('IDLgrImage',TRANSFORM_MODE = 1,DEPTH_TEST_DISABLE =2,SHADER=oShader)
;      oImage->SetProperty, DATA = BYTSCL(image)
;      oImage->SetProperty, LOCATION = [0,0,zpos]
;      oStreamlines->SetProperty,SHADER = oShader 
;      oStreamlines->SetProperty,SHADING = 1          
;    endif
        
    ; Apply standard initial rotation.
    oModel->Rotate, [1,0,0], -90
    oModel->Rotate, [0,1,0], 30
    oModel->Rotate, [1,0,0], 30

    XOBJVIEW, oModel, SCALE=1.0, TITLE=title, /BLOCK

    OBJ_DESTROY, [oModel, oPalette]
    if (OBJ_VALID(oArrows)) then $
        OBJ_DESTROY, oArrows
    if (OBJ_VALID(oSymbol)) then begin
        oSymbol->GetProperty, DATA=oPoly
        OBJ_DESTROY, [oPoly, oSymbol]
    endif
;    if (OBJ_VALID(oImage)) then $
;        OBJ_DESTROY, oImage
end















