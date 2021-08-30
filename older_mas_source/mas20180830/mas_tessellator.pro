;; $Id$
;;

; Subroutine name: mas_tessellator_make_icos
; Created by: BT, 2008-07
; Calling Information:
;
; vertlist  - OUTPUT - set this to a named variable to receive the
;                      vertex list
; polylist  - OUTPUT - set this to a named variable to receive the
;                      polygon adjency list
; level     - INPUT  - integer [0,3] specifying the order of the tess.
; show      - INPUT  - set this keyword to display an object graphics
;                      representation of the tess.

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; creates an nth order tesselation of the icosahedron, where n=0,1,2,3

; Editing Information:

pro mas_tessellator_make_icos, $
   vertlist=vertlist, $
   polylist=polylist, $
   connectivity=connectivity, $
   voronoi=voronoi, $
   vvertices=vvertices, $
   level=level, $
   hemisphere=hemisphere, $
   show=show

    phi = double(0.5 + 0.5*sqrt(5))
    norm_factor = double(sqrt(1.0 + 4.0*(cos(!pi/5))^2))

    ;; icos. vertices
    vertlist = double([ [ 0, phi, 1], $
                        [ 0,-phi, 1], $
                        [ 0, phi,-1], $
                        [ 0,-phi,-1], $
                        [ 1, 0, phi], $
                        [-1, 0, phi], $
                        [ 1, 0,-phi], $
                        [-1, 0,-phi], $
                        [ phi, 1, 0], $
                        [-phi, 1, 0], $
                        [ phi,-1, 0], $
                        [-phi,-1, 0] ])
    
    if not keyword_set(level) then level = 0 else begin
        if (level gt 5) then level = 5
    endelse

    nv = (size(vertlist, /dimensions))[1]

    ;; normalize
    vertlist /= norm_factor
    connectivity = 0
    
    for lev = 1, level do begin

        vertlist_new = vertlist    
        norms = fltarr(nv) + 99.0

        ;; for each vertex, find its nearest neighbors
        for v = 0, nv-1 do begin
            for ov = 0, nv-1 do begin
                
                if (ov ne v) then begin
                    ;; distance between this vert and all others
                    nm = reform(sqrt(total( (vertlist[*,v] - vertlist[*,ov])^2 ) ))
                    if (nm gt 1e-6) then norms[ov] = nm; else norms[ov] = 99.0
                endif else norms[ov] = 99.0

            endfor
            ;; known values for adjency
            case lev of 
                1: min_n = 1.5
                2: min_n = 0.7
                3: min_n = 0.4
                4: min_n = 0.2
                5: min_n = 0.1
                else:
            endcase

            close_points = where(norms lt min_n and norms gt 1e-5)
            sames = 0
            ;; for each nearest neighbor, bisect the line segment
            ;; joining it to this current vertex and add a new point there.
            for cp = 0, n_elements(close_points)-1 do begin
                
                new_pt = (vertlist[*,v] + vertlist[*,close_points[cp]])/2
                new_pt = new_pt/sqrt(total(new_pt^2))
                
                ;; it may have already been added.
                c1  = abs(new_pt[0] - reform(vertlist_new[0,*]))
                c1 += abs(new_pt[1] - reform(vertlist_new[1,*]))
                c1 += abs(new_pt[2] - reform(vertlist_new[2,*]))

                j = where(c1 lt 1e-6, ct)
                
                if (ct eq 0) then begin
                   vertlist_new = [ [vertlist_new], [new_pt] ]
                endif
                
            endfor
        endfor
        
        vertlist = vertlist_new
        nv = (size(vertlist, /dimensions))[1]

    endfor
    
    if (keyword_set(hemisphere)) then begin
        tmp = temporary(vertlist)
        tmp = tmp[*,where(tmp[2,*] ge 0)]
        n_tmp = n_elements(tmp)/3
        dupes = bytarr(n_tmp)
        for i = 0L, n_tmp-1 do begin 
            vert = tmp[*,i]
            for j = i+1, n_tmp-1 do begin 
                other = tmp[*,j]
                if (vert[0] eq -other[0] and vert[1] eq -other[1] and vert[2] eq -other[2]) then begin
                    dupes[j] = 1
                endif
            endfor
        endfor
        vertlist = tmp[*,where(dupes eq 0)]
        nv = n_elements(vertlist)/3
    endif
    
    qhull, vertlist, tetra, /delaunay;, vdiagram=voronoi, vvertices=vvertices
    polylist = tetra_surface(vertlist, tetra)
    npolys = n_elements(polylist)/4
    adj = bytarr(n_elements(vertlist)/3, n_elements(vertlist)/3)
    for i=0L, npolys-1 do begin
        pnum = long(i*4)
        poly = polylist[pnum+1:pnum+3]
        adj[poly[0],poly[1]] = 1;
        adj[poly[0],poly[2]] = 1;
        adj[poly[1],poly[2]] = 1;
    endfor
;; this part will export the vertices and a connectivity
;; array that represents the adjacency list for the vertices
;    adj += transpose(adj)
;    adj[where(adj ne 0)] = 1
;    openw, lun, "~/Desktop/tess.dat", /get_lun
;    writeu, lun, size(vertlist, /dimensions)
;    writeu, lun, float(vertlist)
;    connectivity = lonarr(9, n_elements(vertlist)/3)
;    max_conn = 0L
;    for p=0, n_elements(vertlist)/3-1 do begin
;        tmp = where(adj[p,*] ne 0)
;        if (n_elements(tmp) gt max_conn) then max_conn = n_elements(tmp)
;        connectivity[0,p] = n_elements(tmp)
;        connectivity[1:n_elements(tmp),p] = tmp
;        writeu, lun, n_elements(tmp)
;        writeu, lun, tmp
;    endfor
;    close, lun
;    connectivity = connectivity[0:max_conn,*]
    
    if (keyword_set(show)) then begin
        print, "# Vertices: "+string(nv, format='(1I5)')
        o  = obj_new('idlgrpolygon', vertlist, polygons=polylist, style=2, thick=1, color=[128,128,128])
        o1 = obj_new('idlgrpolygon', vertlist, polygons=polylist, style=0, thick=3,$
                     color=[255,0,0])
        xobjview, [o,o1], renderer=1, /block
        obj_destroy, [o,o1]
    endif
    
end

; Subroutine name: mas_tessellator_make_sphere
; Created by: BT, 2008-07
; Calling Information:
;
; vdim      - INPUT  - the size of the sphere
; vertlist  - OUTPUT - set this to a named variable to receive the
;                      vertex list
; polylist  - OUTPUT - set this to a named variable to receive the
;                      polygon adjency list
; level     - INPUT  - integer [0,3] specifying the order of the tess.
; show      - INPUT  - set this keyword to display an object graphics
;                      representation of the tess.
; normalize - INPUT  - set this keyword to have the vertices returned
;                      in normalized coordinates
; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Provides vertex and polygon list for points on a sphere.

; Editing Information:

pro mas_tessellator_make_sphere, $
     vdim=vdim, $
     decimate_pct=decimate_pct, $
     vertlist=vertlist, $
     polylist=polylist, $
     show=show, $
     normalize=normalize

    volume = dblarr(vdim,vdim,vdim)
  
    for x=0, vdim-1 do begin
       for y=0, vdim-1 do begin
          for z=0, vdim-1 do begin
             volume[x,y,z] = -total( (vdim/2.-.5-[x,y,z])^2 )
          endfor
       endfor
    endfor
    
    scale3, xrange=[0, vdim], yrange=[0, vdim], zrange=[0, vdim]    
    st=-(vdim/4.)^2*3.25
    shade_volume, volume, st, vertlist, polylist, low=1

    if keyword_set(decimate_pct) then begin

       result = mesh_decimate(vertlist, polylist, new_plist, $
                              vertices=new_vlist, $
                              percent_vertices=decimate_pct)
    
       polylist = temporary(new_plist)
       vertlist = temporary(new_vlist)

    endif

    nv = (size(vertlist, /dimensions))[1]
    shift = replicate(-(vdim/2.0), nv)
    vertlist += transpose([ [shift],[shift],[shift] ])

    if (keyword_set(normalize)) then begin
        for i = 0,nv-1 do begin
            vertlist[*,i] /= sqrt(total(vertlist[*,i]^2))
        endfor
    endif

    if (keyword_set(show)) then begin
        omodel = obj_new('idlgrmodel', depth_test_disable=2)
        o = obj_new('idlgrpolygon', vertlist, polygons=polylist, style=1, thick=1)
        omodel->add, o
        xobjview, omodel, renderer=1, /block
        obj_destroy, o
        obj_destroy, omodel
    endif

end

; Subroutine name: mas_tessellator_sa_calc
; Created by: BT, 2009-03
; Calling Information:
;
;   VERT:   a list of vertices that make up the triangles of the
;           tessellation
;
;   POLY:   the polygon connectivity list 
;  
;   AREA:   set to a named variable that will contain the surface
;           area. the surface area is the sum of the area of the
;           triangles
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; computes surface area of a triangulation. for example
; only. use IDL's built-in MESH_SURFACEAREA for speed.
;
; Editing Information:

pro mas_tessellator_sa_calc, vert=vert, poly=poly, area=area

    p = 0     ;; current triangle
    area = 0. ;; total area
    a = 0.    ;; array to hold individual tri's areas

    while (p lt n_elements(poly)-1) do begin

       sum = 0.

       n_pts = poly[p]
       poly_verts = vert[*, poly[p+1:p+n_pts]]
       print, poly_verts & print, '--------------------'
       mat = [ [reform(poly_verts[0,*])], $
               [reform(poly_verts[1,*])], $
               [1.,1.,1.] ]
       sum += determ(mat)^2

       mat = [ [reform(poly_verts[1,*])], $
               [reform(poly_verts[2,*])], $
               [1.,1.,1.] ]
       sum += determ(mat)^2

       mat = [ [reform(poly_verts[0,*])], $
               [reform(poly_verts[2,*])], $
               [1.,1.,1.] ]
       sum += determ(mat)^2
       
       a = (p eq 0) ? 0.5*sqrt(sum) : [a, 0.5*sqrt(sum)]

       p = p+n_pts+1

    endwhile

    area = total(a)

    print, "Num Tris.: "+string(n_elements(a))
    print, "Mean.....: "+string(mean(a))
    print, "Stdev....: "+string(stdev(a))
    print, "Min......: "+string(min(a))  
    print, "Max......: "+string(max(a))  
    print, "Area.....: "+string(area)

;; uncomment to have areas of individual triangles plotted
;;    window, /free, xsize=1100, ysize=350
;;    plot, a, yrange=[0,0.02], xstyle=1, linestyle=0,  title="Area of Triangles"
    
end

pro mas_tessellator_voronoi, index, $
                             vertlist=vertlist,$
                             polylist=polylist,$
                             vor_area=vor_area,$
                             vor_vertices=vor_vertices

    nverts = n_elements(vertlist)/3
    
    testv = vertlist[*,index]

    testv = transpose([ [replicate(testv[0], nverts)], $
                        [replicate(testv[1], nverts)], $
                        [replicate(testv[2], nverts)] ])
    near = sqrt( total((testv-vertlist)^2, 1) )
    near[where(abs(near) lt 1e-6)] = 9999.0
    keep = (sort(near))[0:4]
    vor_verts = vertlist[*,keep]
    vor_near  = near[keep]

    testv = testv[*,0:4]
    vor_poly = (testv + vor_verts)/2.0
        
    vor_area = mean(vor_near)^2 * 1.720477401
    vor_vertices = temporary(vor_poly)

end

pro mas_tessellator_voronoi_test, level=level

    mas_tessellator_make_icos, level=level, vertlist=vlist, polylist=plist, /show
    
    vor_verts = 0
    va_sum = 0
    for v = 0, n_elements(vlist)/3 - 1 do begin
        mas_tessellator_voronoi, v, vertlist=vlist, polylist=plist, vor_area=va, vor_vertices=vv
        ;print, va
        va_sum += va
        vor_verts = (v eq 0) ? vv : [[vor_verts], [vv]]
    endfor
    print, va_sum

end
