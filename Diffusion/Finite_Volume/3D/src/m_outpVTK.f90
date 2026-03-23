
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
   use mesh_types
   implicit none
   contains
   subroutine outpVTK_xyz_multi(filename, phi, mesh, G, use_adjoint)
      type(MeshGrid),intent(in) :: mesh
      logical, intent(in) :: use_adjoint
      character(len=*), intent(in) :: filename
      real(8), intent(in) :: phi(:,:)
      integer, intent(in) :: G

      integer :: nx, ny, nz
      real(8) :: dx, dy, dz
      character(len=256) :: filename_final
      integer :: ii, jj, kk, n0, n1, n2, n3, n4, n5, n6, n7, N, gg

      dx = mesh%dx
      dy = mesh%dy
      dz = mesh%dz

      nx = mesh%nx
      ny = mesh%ny
      nz = mesh%nz

      N = mesh%N

      filename_final = trim(filename)
      if (use_adjoint) then
         filename_final = trim(filename) // "_adjoint"
      endif

      open(unit=11, file="../" // trim(filename_final) // ".vtk", status="replace", action="write")


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
               write(11,'(F9.5, F9.5, ES16.8)') ii*dx, jj*dy, kk*dz
            end do
         end do
      end do

      write(11,'(A)') ''
      write(11,'(A,I0,A,I0)') 'CELLS ', N,' ', 9*N

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
      write(11,'(A,I0,A,I0)') 'CELL_TYPES ', N
      do ii = 1, N
         write(11,'(I5)') 12    
      end do
        
      write(11,'(A)') ''
      write(11,'(A,I0,A,I0)') 'CELL_DATA ', N

      do gg = 1,G
         write(11,'(A,I0,A)') 'SCALARS phi_g', gg,' double'
         write(11,'(A)') 'LOOKUP_TABLE default'
         do kk = 1, nz
            do jj = 1, ny
               do ii = 1, nx
                  write(11,'(ES16.8)') phi(gg, (kk-1)*nx*ny + (jj-1)*nx + ii)
               end do
            end do
         end do
      end do

      close(11)
   end subroutine outpVTK_xyz_multi

   subroutine outpvtk_current(filename, mesh, Jx_g, Jy_g, Jz_g)
      use Mesh_types, only: MeshGrid
      implicit none

      character(len=*), intent(in) :: filename
      type(MeshGrid),   intent(in) :: mesh

      real(8), intent(in) :: Jx_g(:)   ! size N = nx*ny*nz
      real(8), intent(in) :: Jy_g(:)
      real(8), intent(in) :: Jz_g(:)

      integer :: iounit
      integer :: ii, jj, kk, row

      ! Open file
      open(newunit=iounit, file="../"//filename, status='replace', action='write', form='formatted')

      ! VTK header
      write(iounit,'(A)')'# vtk DataFile Version 3.0'
      write(iounit,'(A)')'Current vector field'
      write(iounit,'(A)')'ASCII'
      write(iounit,'(A)')'DATASET STRUCTURED_POINTS'
      write(iounit,'(A,3I8)')        'DIMENSIONS', mesh%nx, mesh%ny, mesh%nz
      write(iounit,'(A,3F16.8)')     'ORIGIN', 0.0d0, 0.0d0, 0.0d0
      write(iounit,'(A,3F16.8)')     'SPACING', mesh%dx, mesh%dy, mesh%dz
      write(iounit,'(A,I8)')         'POINT_DATA', mesh%N

      ! Treat each cell as a "point" for visualization.
      write(iounit,*) 'VECTORS J double'

      do kk = 1, mesh%nz
         do jj = 1, mesh%ny
            do ii = 1, mesh%nx
            row = ii + (jj-1)*mesh%nx + (kk-1)*mesh%nx*mesh%ny
            write(iounit,'(3E20.12)') Jx_g(row), Jy_g(row), Jz_g(row)
            end do
         end do
      end do

      close(iounit)
   end subroutine outpvtk_current

end module m_outpVTK