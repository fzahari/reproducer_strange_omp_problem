    program test
    use mod_grid_storage, ONLY: reset_leb_grid_stack, push_last_leb_grid
    use mod_dft_molgrid, ONLY: init_slices, find_neighbours, get_sorted_lebedev_pts, split_grid
    USE mod_nosp_basis, ONLY: split_sp_basis, nosp_basis
    USE mod_dft_fuzzycell, ONLY: comp_basis_mxdists, dft_fc_blk
    USE prec, ONLY: fp
    USE params, ONLY: dft_bfc_algo
    implicit none
    integer :: i, nat, npt, ncntr, ngrids, npoints, maxpts, nrad
    real(KIND=fp) :: rij(1), const1, rad, prunerads(10),txyz(3*302),twght(302)
    real(KIND=fp) :: xdat(302,1,1),ydat(302,1,1),zdat(302,1,1),wght(302,1,1)
    real(KIND=fp) :: atmxvec(1,1), atmyvec(1,1), atmzvec(1,1), wtab(1,1,24*302)
    real(KIND=fp) :: ptrad(24), wtrad(24)
    real(KIND=fp), parameter :: PI=3.141592653589793238E+00_fp
    integer :: iangn(1,2,1)

    dft_bfc_algo=1

    const1=23.0258_fp
    rij=0.0_fp
    atmxvec=0.0_fp
    atmyvec=0.0_fp
    atmzvec=0.0_fp
    wtab=0.0_fp

    nat=1
    npt=7248
    ncntr=1
    ngrids=1
    npoints=110
    maxpts=302
    nrad=24

    rad=0.58819622124651449_fp

    prunerads=0.0_fp
    prunerads(1)=1.0E30_fp

!    call split_sp_basis
!    call comp_basis_mxdists(nosp_basis,const1)

    call radpt(ptrad,wtrad,nrad)

    call init_slices(500*NAT,NAT,NPT)
    call find_neighbours(rij,NAT)

    call reset_leb_grid_stack

    call get_sorted_lebedev_pts(txyz,twght,npoints,maxpts)
    call push_last_leb_grid

    call split_grid(ncntr,ngrids,rad,prunerads)    
    call dft_fc_blk(atmxvec,atmyvec,atmzvec,rij,nat,wtab)

    write(*,*) "Success"

    end program test 
