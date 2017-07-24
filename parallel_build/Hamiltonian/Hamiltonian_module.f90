module Hamiltonian_module
	use perturbation
	use accuracy
	use mpi_mod
	use Timer
implicit none

subroutine construct_hamiltonian(jrot)

    integer(ik),intent(in) :: jrot ! rotational quantum number
    !
    integer(ik)  :: alloc,alloc_p,dimen,dimen_row,dimen_maxrow
    integer(hik) :: matsize
    integer(ik) :: irow,isym,rlevel,k_i,k_t,istart,jstart,iend,jelem,k_j
    type(PTcoeffsT) :: smat(sym%Nrepresen)
    integer(ik) :: kblock(sym%Nrepresen,0:jrot,2)
    !
    real(rk),allocatable :: mat_t(:,:,:),a(:,:)
    integer(ik) :: ielem,icase,Nterms(sym%Nrepresen),i_irr(sym%Nrepresen),nroots
    integer(ik),allocatable :: ijterm(:,:),k_row(:,:),bterm(:,:)
    !
    real(rk) :: zpe
    integer  :: slevel,dimen_s,max_dim,iterm,jterm,total_roots,icontr
    !
    integer(ik) :: iunit,unitO,unitC,rec_len,irec_len,chkptIO
    integer(ik) :: ncontr,maxcontr,maxcontr0
    character(len=cl)   :: task
    character(len=4)   :: jchar
    character(len=cl)  :: unitfname,filename,statusf,symchar
    logical :: only_store = .false.
    logical :: no_diagonalization = .false.	

    !Determine the size of Htotal:
    !
    dimen = PT%Maxsymcoeffs ! PT%Maxcontrcoeffs
    !
    if (job%verbose>=2) then 
        write (out,"(/'Size of the contracted matrix  = ',i7)") dimen
    endif

    allocate (ijterm(dimen,sym%Nrepresen),stat=alloc)
    call ArrayStart('PThamiltonian_contract:a',alloc,size(ijterm),kind(ijterm))

    ! Count and distribute the rows in the symmetrized representaion: 
    !
    Nterms = 0 
    !
    do irow = 1,dimen
      !
      do isym = 1,sym%Nrepresen
        i_irr(isym) = PT%irr(isym)%N(irow) 
      enddo 
      !
      !jterm = 0 
      !
      ijterm(irow,:) = Nterms(:)
      !
      Nterms(:) = Nterms(:) + i_irr(:)
      !
    enddo
    !
    if (job%verbose>=3) write(out,"(/'Calculating the symmetrized mat. elements...'/)")
    !
    k_t = 0 
    kblock(:,:,1) = 1 ; kblock(:,:,2) = 0
    do irow = 1,dimen
      !
      iterm = PT%contractive_space(0,irow)
      k_i = PT%rot_index(iterm,1)%k
      !
      if (k_i==k_t) then 
        !
        do isym = 1,sym%Nrepresen
           !
           kblock(isym,k_i,2) = ijterm(irow,isym)+PT%irr(isym)%N(irow)
           !
        enddo
        !
      else
        !
        do isym = 1,sym%Nrepresen
           !
           k_t = k_i
           kblock(isym,k_i,1) = ijterm(irow,isym) + 1
           kblock(isym,k_i,2) = ijterm(irow,isym)+PT%irr(isym)%N(irow)
           !
        enddo
        !
      endif 
      !
    enddo

    ! forget about the kblock structure in case of the TD symmetry or if the rot. symmetry is based on the euler angles transformation properties
    !
    if (job%vib_rot_contr.or.(job%rotsym_do.and.any( trim( trove%symmetry )==(/'TD','TD(M)','C2V(M)','C2V','C2V','C3V(M)','C3V','D3H(M)','D3H','C2H(M)','C2H','G4(M)','G4','G4(EM)','D2H(M)','D2H'/) ) ) ) then 
      !
      do isym = 1,sym%Nrepresen
        !
        dimen_s = PT%Max_sym_levels(isym)
        !
        kblock(isym,:,1) = 1
        kblock(isym,:,2) = dimen_s
        !
      enddo
    endif


    task = 'top'
    call PTrestore_rot_kinetic_matrix_elements(jrot,task,iunit,dimen,ncontr,maxcontr)
    
      	
end subroutine


end module Hamiltonian_module
