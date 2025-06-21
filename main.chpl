use Help;
use structures;

proc main(args: [] string) {  
    var w: wall;
    const npoints: int(64)  = 5;
    const nghosts:  int(64)  = 2;
    const xmin:    real(64) = 0.0;
    const xmax:    real(64) = 5.0;
    var Grid = new owned grid(xmin, xmax, npoints, nghosts);
    writeln("cell centers: ");
    for i in Grid.cells_tot.domain {
      write(Grid.cells_tot[i].center);
      write(" ");
    }
    writeln();
    writeln("left wall positions: ");
    for i in Grid.cells_tot.domain {
      write(Grid.cells_tot[i].wall_left.position);
      write(" ");
    }
    writeln();
    writeln("right wall positions: ");
    for i in Grid.cells_tot.domain {
      write(Grid.cells_tot[i].wall_right.position);
      write(" ");
    }
    writeln();
}
