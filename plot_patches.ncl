load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" ;copyatt

; Plots patches, read from ascii, on Lambert Conformal map projection

begin ; PRESETS
  fname1 = "lvl01_trace.txt"
  fname2 = "lvl02_trace.txt"
  fname3 = "lvl03_trace.txt"
  
  ncheader = "/home/adam/data/FLX_MERGED/yearly/nexc_d02_VEI7_200803-200902.nc"

  header1 = "lvl01_x0y0_x1y1.txt"
  header2 = "lvl02_x0y0_x1y1.txt"
  header3 = "lvl03_x0y0_x1y1.txt"
  
  outname="leadbetter_vertical_t1"

  ; The coordinate used to convert pixels to lat-lon:
  xmin_in_lon = todouble(-50.)
  ymin_in_lat = todouble(30.)

  xmax_in_lon = todouble(50.)
  ymax_in_lat = todouble(80.)

  ; Output type
  ;type = "x11"
  ;type = "pdf"
  ;type "ps"
  type = "eps"
  
  ; Map projection settings as in WRF grid (Lambert Conformal is a very special case!)
  res = True
  res@mpProjection = "LambertConformal"
  res@mpLambertMeridianF  = -3.0
  res@mpLambertParallel1F = 50.0
  res@mpLambertParallel2F = 70.0
  ; More map projection under misc plot settings (should be automatic)
end   ;PRESETS

; Function to convert xy-coordinate (on a lat-lon grid) to lat-lon coordinates
function xy2ll(xy,x0,y0,x1,y1,lon0,lat0,lon1,lat1)
local x2lon,y2lat ;variable list
begin
  ; Get factors to convert xy-distance to ll-distance:
  printVarSummary(lon1)
  printVarSummary(x1)
  dlondx = (lon1-lon0)/(x1-x0)  ; distance
  dlatdy = (lat1-lat0)/(y1-y0)  ;
  ; Convert coordinate to relative (lon,lat)=(0,0) and to lat-lon units
  ll = xy
  ll(:,0) = (xy(:,0)-x0)*dlondx+lon0
  ll(:,1) = (xy(:,1)-y0)*dlatdy+lat0
  return ll
end

begin ; load data and convert to lat-lon coordinates
  ncol = 2

  ; Read data from files/headers:
  c1xy = readAsciiTable(fname1,ncol,"double",0)
  c1minmax = readAsciiTable(header1,ncol,"double",0)
  c1xmin = c1minmax(0,0)
  c1xmax = c1minmax(1,0)
  c1ymin = c1minmax(0,1)
  c1ymax = c1minmax(1,1)

  c2xy = readAsciiTable(fname2,ncol,"double",0)
  c2minmax = readAsciiTable(header2,ncol,"double",0)
  c2xmin = c2minmax(0,0)
  c2xmax = c2minmax(1,0)
  c2ymin = c2minmax(0,1)
  c2ymax = c2minmax(1,1)

  c3xy = readAsciiTable(fname3,ncol,"double",0)
  c3minmax = readAsciiTable(header3,ncol,"double",0)
  c3xmin = c3minmax(0,0)
  c3xmax = c3minmax(1,0)
  c3ymin = c3minmax(0,1)
  c3ymax = c3minmax(1,1)
end

begin ; Convert xy-coord to ll-coord:
  c1ll = xy2ll(c1xy,c1xmin,c1ymin,c1xmax,c1ymax,xmin_in_lon,ymin_in_lat,xmax_in_lon,ymax_in_lat)
  c2ll = xy2ll(c2xy,c2xmin,c2ymin,c2xmax,c2ymax,xmin_in_lon,ymin_in_lat,xmax_in_lon,ymax_in_lat)
  c3ll = xy2ll(c3xy,c3xmin,c3ymin,c3xmax,c3ymax,xmin_in_lon,ymin_in_lat,xmax_in_lon,ymax_in_lat)
end

begin ; Read lat2d,lon2d from netcdf file
  f    = addfile(ncheader,"r")
  lat2d = f->XLAT  ; Native latitude (2D)
  lon2d = f->XLONG ; Native longitude (2D)
  nxy   = dimsizes(lat2d)
  nlat  = nxy(0)
  nlon  = nxy(1)
end

begin ; Create work station and define colors
  wks  = gsn_open_wks (type, outname)
  gsn_define_colormap(wks,(/"white","black","tan","LightBlue","Blue"/))
end


begin ; Polygon/Polyline resources:
  pgres                  = True
  pgres@gsEdgesOn        = True
  pgres@gsLineColor      = "black"
  pgres@gsLineThicknessF = 2.0
  pgres@gsFillColor      = "tan"
end

begin ; misc plot resources:
  res@gsnDraw                     = False
  res@gsnFrame                    = False

  ; Coast-lines
  res@mpDataBaseVersion           = "Ncarg4_1"
  res@mpDataBaseVersion           = "MediumRes"          ; Finer coast-lines
  ;res@mpGeophysicalLineColor      = "black"
  res@mpGeophysicalLineThicknessF = 0.5
  ;res@mpFillColors                = "black"
  ;res@mpFillOn                    = True                ; turn on map fill
  ;res@mpAreaMaskingOn             = True
  ;res@mpLandFillColor             = "tan"
  res@mpOutlineOn                 = True
  res@mpOutlineDrawOrder          = "PostDraw"
  res@mpFillOn                    = False

  ; Grid and ticks
  res@pmTickMarkDisplayMode       = "Always"             ; turn on tick marks
  res@mpGridAndLimbOn             = True
  res@mpGridSpacingF              = 5.0
  res@mpGridLineThicknessF        = 0.1
  res@mpGridLineColor             = "black"


  ;res@tiMainString                = var2d@description
  ;res@tiMainString                = "Regional Domain on Stereographic Proj"

  ;res@gsnAddCyclic                = False                ; regional data

  ; Zoom in on world map (a bit ugly, could be done better, e.g. manually)
  res@mpLimitMode               = "Corners"           ; method to zoom
  res@mpLeftCornerLatF          = lat2d(0,0)-(lat2d(0,0)-lat2d(0,0))
  res@mpLeftCornerLonF          = lon2d(0,0)-(lon2d(0,0)-lon2d(0,0))
  res@mpRightCornerLatF         = 2.*lat2d(nlat-1,nlon-1)-lat2d(nlat-1,nlon-1)
  res@mpRightCornerLonF         = 2.*lon2d(nlat-1,nlon-1)-lon2d(nlat-1,nlon-1)

  ;res@tfDoNDCOverlay              = True                 ; don't transform 
  res@tfDoNDCOverlay              = False ; Yes, please!
end ;misc plot settings

begin ;create plot
  ;setvalues wks            
  ;  "wkForegroundColor" : (/0.,0.,0./)  
  ;  "wkBackgroundColor" : (/1.,1.,1./) 
  ;end setvalues 
  map = gsn_csm_map(wks,res)

  dum1 = gsn_add_polygon(wks,map,c1ll(:,0),c1ll(:,1),pgres)
  dum2 = gsn_add_polygon(wks,map,c2ll(:,0),c2ll(:,1),pgres)
  dum3 = gsn_add_polygon(wks,map,c3ll(:,0),c3ll(:,1),pgres)

  draw(map)
  frame(wks)
end 
