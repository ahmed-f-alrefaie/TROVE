program pbuild
	use accuracy
	use fields
	use perturbation
	use MPI_module

implicit none

      integer(ik) :: NPTorder,Natoms,Nmodes,Npolyads
      integer(ik) :: j
      !type(FLbasissetT),pointer  :: basisset(:)     ! Basis set specifications: range and type

      !
      !type(FLbasissetT)          :: rotbasis            ! Rotational basis set specifications
      !
      ! Begin with constants intitialization
      !
      call TimerStart('TROVE')
      !
      call accuracyInitialize
      
      call mpi_mod_start(.true.)
      !
      if (job%verbose>=4) then 
        write(out,"('spacing around 1.0 is ',d18.8)") spacing(1.0d0)
      endif 

      call FLReadInput(NPTorder,Npolyads,Natoms,Nmodes,j)
      call FLsetMolecule
      call FLinitilize_Kinetic
      call FLinitilize_Potential
      call FLinit_External_field
      call FLbsetInit
      call PTinit(NPTorder,Nmodes,Npolyads)	
      if (trim(job%IOcontr_action)=='READ') then
        !
        call PTcheck_point_contracted_space('READ') 
      else 
      		stop "This is only for building, anything else is not supported"
      endif
      
      call PT_conctracted_rotational_bset(j)
      !
      call PTsymmetrization(j)
      
      
      if (job%contrci_me_fast) then 
        !
        call PTstore_contr_matelem(j)
        !
        call PTcontracted_matelem_class_fast(j) 
        !
      else
        !
        call PTcontracted_matelem_class(j)
        !
      endif
      
      !This is what I'm going to replace      
      call PThamiltonian_contract(j)
      
      
      call TimerStop('TROVE')


end program pbuild
