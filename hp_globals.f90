! Empirical Mode Decomposition for HEALPix maps

module hp_globals

implicit none

character(len=80)  :: fin, arg, header(43)
integer, parameter :: RING = 1, NEST = 2	! ordering literals
integer, parameter :: CART = 1, HLPX = 2	! vector components

end module
