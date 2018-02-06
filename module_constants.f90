module module_constants
      integer, parameter :: nclg = 4               ! no of bulk cloud layers: total,high,middle,low
      real, dimension(nclg) :: pt             ! pressure (hPa) at top of bulk cloud layers
      real, dimension(nclg) :: pb             ! pressure (hPa) at top of bulk cloud layers
      data pt / 50.0,  50.0, 440.0, 680.0 /
      data pb /1200.0, 440.0, 680.0,1200.0/
      logical, parameter :: doflip = .true.        ! do vertical flip because input level is surface->top
      
end module module_constants
