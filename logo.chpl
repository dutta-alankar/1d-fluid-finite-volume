module logo {

    proc print_logo (): void {
        // **CHHA** *Chapel-based Hydrodynamics for HPC in Astrophysics*
        writeln("""
       (  )   (   )  )
        ) (   )  (  (
        ( )  (    ) )
       ________________
      <     CHHA       > ___
      |  Chapel-based  |/ _ \
      | Hydrodynamics  | | | |
      |    for HPC     |_|_| |
      |      in        |\___/
   ___|  Astrophysics  |___
  /    \______________/    \
  \________________________/ 
    """);
    }
}