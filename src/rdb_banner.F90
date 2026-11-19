!! This file contains the cute banner we will print when a simulation is run
module rdb_banner
   use pic_logger, only: logger => global_logger
   use pic_knowledge, only: get_knowledge
   implicit none
   private
   public :: print_banner
contains

   subroutine print_banner(sim_type)
    !! Print a cute banner to the console when the simulation starts.  The
    !! subtitle line adapts to the regime: "ocean circulation model" for
    !! `sim_type == "ocean"`, the 2D/3D coastal line otherwise (the default).
      character(len=*), intent(in), optional :: sim_type
      character(len=54) :: subtitle   ! 54 = interior box width (the ═ border)

      subtitle = "       2D/3D Coastal Hydrodynamics Model"
      if (present(sim_type)) then
         if (trim(sim_type) == "ocean") subtitle = "            Global Circulation Model"
      end if

      call logger%info("")
      call logger%info("  ╔══════════════════════════════════════════════════════╗")
      call logger%info("  ║                                                      ║")
      call logger%info("  ║     ██████╗  ██████╗ ██╗   ██╗███╗   ██╗██████╗      ║")
      call logger%info("  ║     ██╔══██╗██╔═══██╗██║   ██║████╗  ██║██╔══██╗     ║")
      call logger%info("  ║     ██████╔╝██║   ██║██║   ██║██╔██╗ ██║██║  ██║     ║")
      call logger%info("  ║     ██╔══██╗██║   ██║██║   ██║██║╚██╗██║██║  ██║     ║")
      call logger%info("  ║     ██║  ██║╚██████╔╝╚██████╔╝██║ ╚████║██████╔╝     ║")
      call logger%info("  ║     ╚═╝  ╚═╝ ╚═════╝  ╚═════╝ ╚═╝  ╚═══╝╚═════╝      ║")
      call logger%info("  ║      █████╗ ██████╗  ██████╗ ██╗   ██╗████████╗      ║")
      call logger%info("  ║     ██╔══██╗██╔══██╗██╔═══██╗██║   ██║╚══██╔══╝      ║")
      call logger%info("  ║     ███████║██████╔╝██║   ██║██║   ██║   ██║         ║")
      call logger%info("  ║     ██╔══██║██╔══██╗██║   ██║██║   ██║   ██║         ║")
      call logger%info("  ║     ██║  ██║██████╔╝╚██████╔╝╚██████╔╝   ██║         ║")
      call logger%info("  ║     ╚═╝  ╚═╝╚═════╝  ╚═════╝  ╚═════╝    ╚═╝         ║")
      call logger%info("  ║                                                      ║")
      call logger%info("  ║        A GPU accelerated, Fortran native             ║")
      call logger%info("  ║"//subtitle//"║")
      call logger%info("  ║                                                      ║")
      call logger%info("  ╚══════════════════════════════════════════════════════╝")
      call logger%info("")
      call logger%info("  Contributors:")
      call logger%info("    Jorge Luis Galvez Vallejo")
      call logger%info("    Albert Sietze Thie")
      ! add your name here after you've contributed!
      call logger%info("")
      call get_knowledge()
   end subroutine print_banner

end module rdb_banner
