
!------------------------------------------------------------------------!
! Purpose:                                                              -!
!  Produce VTK files for 2D discretized NDE, xy, rz and rth geometries  -!  
!  Writes : Points, Cells, Cell_Types, Cell_Data                        -!
!------------------------------------------------------------------------!
!! Record of revisions:                                                 -!
!   Date       Programmer     Description of change                     -!
!   ====       ==========     =====================                     -!
! 14/11/25     T. Charlton    Implemented 2D matrix                     -!
! 15/11/25     T. Charlton    Implemented rth matrix                    -!
! 17/11/25     T. Charlton    Implemented rz matrix                     -! 
!------------------------------------------------------------------------! 

module m_outpVTK
   implicit none
   contains
   subroutine outpVTK(phi, nx, ny, dx, dy)
      integer, intent(in) :: nx, ny
      real(8), intent(in) :: phi(:), dx, dy 
      integer :: ii, jj, kk, n0, n1, n2, n3, count
      real(8) :: x, y, z
      real(8), allocatable ::  phi_v(:,:)

      allocate(phi_v(0:nx,0:ny))
      !-----------------------------------------
      ! Compute vertex-centred phi_v by averaging
      !-----------------------------------------
      do jj = 0, ny
         do ii = 0, nx
            phi_v(ii,jj) = 0.0
            count = 0

            ! cell (ii,jj)
            if (ii < nx .and. jj < ny) then
               kk = jj*nx + ii + 1
               phi_v(ii,jj) = phi_v(ii,jj) + phi(kk)
               count = count + 1
            end if

            ! cell (ii-1,jj)
            if (ii > 0 .and. jj < ny) then
               kk = jj*nx + (ii-1) + 1
               phi_v(ii,jj) = phi_v(ii,jj) + phi(kk)
               count = count + 1
            end if

            ! cell (ii,jj-1)
            if (ii < nx .and. jj > 0) then
               kk = (jj-1)*nx + ii + 1
               phi_v(ii,jj) = phi_v(ii,jj) + phi(kk)
               count = count + 1
            end if

            ! cell (ii-1,jj-1)
            if (ii > 0 .and. jj > 0) then
               kk = (jj-1)*nx + (ii-1) + 1
               phi_v(ii,jj) = phi_v(ii,jj) + phi(kk)
               count = count + 1
            end if

            if (count > 0) phi_v(ii,jj) = phi_v(ii,jj)/count

         end do
      end do
        
      open(unit=11, file="flux_points.vtk", status="replace", action="write")

      ! Write header
      write(11,'(A)') '# vtk DataFile Version 2.0'
      write(11,'(A)') 'cell-centred output'
      write(11,'(A)') 'ASCII'
      write(11,'(A)') 'DATASET UNSTRUCTURED_GRID'
      write(11,'(A,I0,A)') 'POINTS ', (nx+1)*(ny+1), ' double'

      ! Loop over 2D domain nodes
      do jj = 0, ny
         y = jj*dy        ! y-direction
         do ii = 0, nx    ! x-direction
            x = ii*dx
            z = nx*dx*phi_v(ii,jj)/phi_v(nx/2,ny/2)
            write(11,'(F9.5 F9.5 ES16.8)') x, y, z
         end do
      end do

      write(11,'(A)') ''
      write(11,'(A,I0,A,I0)') 'CELLS ', nx*ny,' ', 5*nx*ny

      do jj = 0, (ny-1)
         do ii = 0, (nx-1)
         ! point indices for quad cell
         ! node numbering: row-major
            n0 = jj*(nx+1) + ii
            n1 = n0 + 1
            n2 = n1 + (nx+1)
            n3 = n0 + (nx+1)
            write(11,'(I6,4I6)') 4, n0, n1, n2, n3
         end do
      end do

      write(11,'(A)') ''
      write(11,'(A,I0,A,I0)') 'CELL_TYPES ', nx*ny

      do jj = 1, nx*ny
         write(11,'(I5)') 9    
      end do
        
      write(11,'(A)') ''
      write(11,'(A,I0,A,I0)') 'CELL_DATA ', nx*ny
      write(11,'(A)') 'SCALARS phi double'
      write(11,'(A)') 'LOOKUP_TABLE default'
      do jj = 1, ny
         do ii = 1, nx
            kk = (jj-1)*nx + ii
            write(11,'(ES16.8)') phi(kk)
         end do
      end do

      close(11)
   end subroutine outpVTK

   subroutine outpVTK_xyz(phi, nx, ny, nz, dx, dy, dz)
      integer, intent(in) :: nx, ny, nz
      real(8), intent(in) :: phi(:), dx, dy, dz
      integer :: ii, jj, kk, n0, n1, n2, n3, n4, n5, n6, n7
        
      open(unit=11, file="flux_points.vtk", status="replace", action="write")

      ! Write header
      write(11,'(A)') '# vtk DataFile Version 2.0'
      write(11,'(A)') 'cell-centred output'
      write(11,'(A)') 'ASCII'
      write(11,'(A)') 'DATASET UNSTRUCTURED_GRID'
      write(11,'(A,I0,A)') 'POINTS ', (nx+1)*(ny+1)*(nz+1), ' double'

      ! Loop over 2D domain nodes
      do kk = 0, nz
         do jj = 0, ny 
            do ii = 0, nx
               write(11,'(F9.5 F9.5 ES16.8)') ii*dx, jj*dy, kk*dz
            end do
         end do
      end do

      write(11,'(A)') ''
      write(11,'(A,I0,A,I0)') 'CELLS ', nx*ny*nz,' ', 9*nx*ny*nz

      do kk = 0, (nz-1)
         do jj = 0, (ny-1)
            do ii = 0, (nx-1)
               n0 = kk*(ny+1)*(nx+1) + jj*(nx+1) + ii
               n1 = n0 + 1
               n2 = n1 + (nx+1)
               n3 = n0 + (nx+1)
               n4 = n0 + (nx+1)*(ny+1)
               n5 = n4 + 1
               n7 = n4 + (nx+1)
               n6 = n7 + 1

               write(11,'(I4,8I8)') 8, n0, n1, n2, n3, n4, n5, n6, n7
            end do
         end do
      end do

      write(11,'(A)') ''
      write(11,'(A,I0,A,I0)') 'CELL_TYPES ', nx*ny*nz
      do ii = 1, nx*ny*nz
         write(11,'(I5)') 12    
      end do
        
      write(11,'(A)') ''
      write(11,'(A,I0,A,I0)') 'CELL_DATA ', nx*ny*nz
      write(11,'(A)') 'SCALARS phi double'
      write(11,'(A)') 'LOOKUP_TABLE default'
      do kk = 1, nz
         do jj = 1, ny
            do ii = 1, nx
               write(11,'(ES16.8)') phi((kk-1)*nx*ny + (jj-1)*nx + ii)
            end do
         end do
      end do

      close(11)
   end subroutine outpVTK_xyz

   subroutine outpVTK_rz(phi, nr, nz, dr, dz)
      integer, intent(in) :: nr, nz
      real(8), intent(in) :: phi(:), dr, dz 
      integer :: ii, jj, kk, n0, n1, n2, n3
      real(8) :: r, z
        
      open(unit=11, file="flux_points.vtk", status="replace", action="write")

      ! Write header
      write(11,'(A)') '# vtk DataFile Version 2.0'
      write(11,'(A)') 'cell-centred output'
      write(11,'(A)') 'ASCII'
      write(11,'(A)') 'DATASET UNSTRUCTURED_GRID'
      write(11,'(A,I0,A)') 'POINTS ', (nr+1)*(nz+1), ' double'

      ! Loop over 2D domain nodes
      do jj = 0, nz
         z = jj*dz        ! y-direction
         do ii = 0, nr    ! x-direction
            r = ii*dr
            write(11,'(F9.5 F9.5 ES16.8)') r, z, 0.0
         end do
      end do

      write(11,'(A)') ''
      write(11,'(A,I0,A,I0)') 'CELLS ', nr*nz,' ', 5*nr*nz

      do jj = 0, (nz-1)
         do ii = 0, (nr-1)
         ! point indices for quad cell
         ! node numbering: row-major
            n0 = jj*(nr+1) + ii
            n1 = n0 + 1
            n2 = n1 + (nr+1)
            n3 = n0 + (nr+1)
            write(11,'(I6,4I6)') 4, n0, n1, n2, n3
         end do
      end do

      write(11,'(A)') ''
      write(11,'(A,I0,A,I0)') 'CELL_TYPES ', nr*nz

      do jj = 1, nr*nz
         write(11,'(I5)') 9    
      end do
        
      write(11,'(A)') ''
      write(11,'(A,I0,A,I0)') 'CELL_DATA ', nr*nz
      write(11,'(A)') 'SCALARS phi double'
      write(11,'(A)') 'LOOKUP_TABLE default'
      do jj = 1, nz
         do ii = 1, nr
            kk = (jj-1)*nr + ii
            write(11,'(ES16.8)') phi(kk)
         end do
      end do

      close(11)
   end subroutine outpVTK_rz

   subroutine outpVTK_rth(phi, nr, nth, dr, dth)
      implicit none
      integer, intent(in) :: nr, nth
      real(8), intent(in) :: phi(:), dr, dth
      integer :: ir, it, kk, cell_id
      real(8) :: r, theta
      integer :: np, nc, total_indices

      open(unit=11, file="flux_cylinder.vtk", status="replace", action="write")

      !-----------------------------------
      ! Header
      !-----------------------------------
      write(11,'(A)') '# vtk DataFile Version 2.0'
      write(11,'(A)') 'cell-centred output'
      write(11,'(A)') 'ASCII'
      write(11,'(A)') 'DATASET UNSTRUCTURED_GRID'

      !-----------------------------------
      ! Points
      !-----------------------------------
      
      np = 1 + nr * nth   ! 1 center point + nr radial layers * nth angular points
      write(11,'(A,I0,A)') 'POINTS ', np, ' double'

      ! Central point
      write(11,'(3F10.6)') 0.0d0, 0.0d0, 0.0d0

      ! Outer points
      do ir = 1, nr
         r = ir * dr
         do it = 0, nth-1
            theta = it * dth
            write(11,'(3F10.6)') r*cos(theta), r*sin(theta), 0.0
         end do
      end do

      !-----------------------------------
      ! Cells
      !-----------------------------------
      nc = (nth) + (nr-1)*nth   ! nth triangles at center + remaining quads
      total_indices = nth*4 + (nr-1)*nth*5   ! 3+1 for triangles, 4+1 for quads
      
      write(11,'(A)') ''
      write(11,'(A,I0,A,I0)') 'CELLS ', nc, ' ', total_indices

      !----- Triangles at center -----
      do it = 1, nth-1
         write(11,'(I6,3I6)') 3, 0, it, it+1
      end do
      ! Last triangle wraps around
      write(11,'(I6,3I6)') 3, 0, nth, 1

      !----- Quads for remaining radial layers -----
      do ir = 1, nr-1
         do it = 1, nth
            kk = 1 + (ir-1)*nth + it -1
            if (it < nth) then
               write(11,'(I6,4I6)') 4, kk, kk+nth, kk+nth+1, kk+1
            else
               ! wrap around last point in angular direction
               write(11,'(I6,4I6)') 4, kk, kk+nth, kk+1, kk-nth+1
            end if
         end do
      end do

      !-----------------------------------
      ! Cell types
      !-----------------------------------
      write(11,'(A,I0)') 'CELL_TYPES ', nc
      ! Triangles
      do it = 1, nth
         write(11,'(I5)') 5   ! VTK_TRIANGLE
      end do
      ! Quads
      do ir = 1, (nr-1)*nth
         write(11,'(I5)') 9   ! VTK_QUAD
         end do

         !-----------------------------------
         ! Cell data
         !-----------------------------------
         write(11,'(A,I0)') 'CELL_DATA ', nc
         write(11,'(A)') 'SCALARS phi double'
         write(11,'(A)') 'LOOKUP_TABLE default'

      do cell_id = 1, nc
         kk = (cell_id-1) + 1  ! simple mapping; adjust if needed for your phi array
         write(11,'(E16.8)') phi(kk)
      end do

      close(11)
   end subroutine outpVTK_rth

end module m_outpVTK