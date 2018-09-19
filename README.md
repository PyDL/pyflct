# pyflct
07 December 2017 - Jiajia Liu, SP2RC, University of Sheffield

## Update 19 September 2018
Fisher & Welsch updated FLCT to 1.0.6, while the original FLCT version I used was 1.0.4.
A new function of "bias correction" has been implemented in the new version by Fisher & Welsch.

A preliminary test of mine shows that:
* Velocity field resulted from 1.0.6 has **no difference** with that from 1.0.4
* If the "bias correction" is turned on, **exactly the same swirls** will be 
detected by ASDA. But the rotating and expanding/shrinking speed of the detected 
swirls would in average **increase by around 30%**.
* More test will be made in the future.

## Discription: 
 Python wrapper for FLCT code written in C from Fisher & Welsch
 2008. You can download the original C code from the following link:
 http://cgem.ssl.berkeley.edu/cgi-bin/cgem/FLCT/home

 Before a proper run of this program, you need first install
 the FLCT libraries. Extract the downloaded C source files, go
 to the fold. Then check source/README-install.txt
 and Makefile to find out how to install the FLCT libraries
 properly.

## Inputs:
         data1 - image data at time T
         data2 - image data at time T+dT
         deltat - time difference between two images
         deltas - pixel size of each image, units of velocities will be
                  based on deltas/deltat
         sigma - pixel width of the Gausian filter. Default is 10 according
                 to [Louis et al. 2015, Solar Phys. 290, 1135], which found
                 the optimal value of sigma for their simulated intensity
                 map ranges from 10 to 15. If sigma is set to 0, only the
                 overall correlation between two input images will be
                 calculated.
         infile - name of file which will be generated storing data1 and data2
         outfile - name of file which stores the FLCT velocity field
         thresh - Do not compute the velocity at a given pixel if the average
                  absolute value between the 2 images at that location is less
                  than thresh.
         kr - from 0 to 1. Perform gaussian, low pass filtering on the
              sub-images that are used to construct the cross-correlation
              function. The value of kr is expressed in units of the maximum
              wavenumber (Nyquist frequency) in each direction of the
              sub-image. This option is most useful when the images contain
              significant amounts of uncorrelated, pixel-to-pixel noise-like
              structure. Empirically, values of kr in the range of 0.2 to 0.5
              seem to be most useful, with lower values resulting in stronger
              filtering.
         skip - if is not None, the program will only compute the velocity
                every N pixels in both the x and y direction, where N is
                the value of skip.
         xoff, yoff - only valid when skip is not None. xoff and yoff
                      are the offset in x and y direction where the compuation
                      starts.
         interp - only valid when skip is not None. If true, interpolation
                  at pixels skipped will be done.
         pc - If set, then the input images are assumed to be in Plate Carree
              coordinates (uniformly spaced in longitude and latitude). This
              is useful when e.g. input images are SHARP magnetic field data.
         bc - If set, bias correction will be turned on (new feature in 1.0.6)
         latmin, latmax - minimum and maximum latitude. Only valid if pc is
                          set. In units of radian.
         quiet - If set, no non-error output will be shown.

## Outputs:
         a tupe in form of (vx, vy, vm)
         where, vx is the velocity field in x direction;
                vy is the velocity field in y direction;
                vm is the velocity field mask. Pixels with vm value of 0 have
                    not been included when calculating the velocity field.
                    Pixels with vm value of 0.5 have interpolated velocity
                    field.

