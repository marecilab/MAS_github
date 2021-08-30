; Procedure to generate an rank-2 tensor ellipsoid taken from IDL developer website
; http://www.harrisgeospatial.com/Learn/Blogs/Blog-Details/TabId/2716/ArtMID/10198/ArticleID/16044/Generating-an-Ellipsoid-using-IDL-8-Graphics.aspx
; Modifications (Magdoom)
; 1) Fixed a bug when dealing with unequal x,y major axis lengths as mentioned in the website
; 2) Added x and y axis titles
; 3) Included optional title text for the graphics window included
; 4) Added optional keyword to set the axes range
; 5) Added an keyword to incorporate the eigenvectors of the tensor

function get_rad, theta, phi, eval, evec
  
  ;use the formula for the ellipsoid to
  ;determine the radius for various values
  ;of theta and phi.
  x = cos(theta)*cos(phi)
  y = sin(theta)*cos(phi)
  z = sin(phi)
  N = size(z,/DIMENSIONS)
  radius = make_array(N, /float)
  for i = 0,N[0]-1 do begin
    for j = 0,N[0]-1 do begin
      R = [[x[i,j]],[y[i,j]],[z[i,j]]]
      Rp = evec^2##R^2
      radius[i,j] = real_part(1/sqrt((1/eval)##Rp))
    endfor
  endfor

  RETURN, radius
end

pro dj_ellipsoid_ng, eval, evec, window_title_txt = window_title_txt
   
   N = 100;
  ;Create a grid of theta values
  ;These are the longitude of the of the
  ;ellipsoid and go from 0 to 360 degrees
  theta_1d = (findgen(N+1)*360/N)*!DTOR
  theta = CONGRID(TRANSPOSE(theta_1d),N+1,N+1)
  theta = REFORM(TRANSPOSE(theta))

  ;Create a grid of phi values. These
  ;are the latitude of the ellipsoid and go
  ;from -90 to 90 degrees
  phi_1d = ((findgen(N+1)-N/2)*180/N)*!DTOR
  phi = congrid(TRANSPOSE(phi_1d),N+1,N+1)
    
  ;Use "get_rad" function to determine the
  ;radius values at each point of the theta
  ;and phi grid
  radius = get_rad(theta,phi,eval,evec)
  
  ;Create a mesh using the radius values
  ;This outputs the vertices (vert) and polygons
  ;that will be used to generate the plot
  mesh_obj,4,vert,pol,radius
 
  ;Determine the biggest value and create a scale from it
  scale = (findgen(10)-5)*max(radius)/5

  ;Plot the scale and axis with no data. Set CLIP=0, to
  ;prevent edges of ellipsoid from getting cut off.
  p_scale = plot3d(scale,scale,scale,/NODATA,CLIP=0, ASPECT_Z=1, ASPECT_RATIO=1, xtitle = 'X', ytitle = 'Y', ztitle = 'Z')
  if keyword_set(window_title_txt) then p_scale.window_title = window_title_txt
  
  ;Use the POLYGON to plot the mesh. Use the
  ;vertices and polygons output from the MESH_OBJ
  ;to fill out the DATA argument and CONNECTIVITY keyword.
  ;Use bright green as the color.
  p = POLYGON(vert,CONNECTIVITY=pol,/data,CLIP=0,fill_color=[255,0,0], LINESTYLE = 6)
 
end
