use Help;
use structures;

proc main(args: [] string) {  
    
    const npoints: int(64)  = 5;
    const nghosts:  int(64)  = 2;
    const xmin:    real(64) = 0.0;
    const xmax:    real(64) = 5.0;
    var Grid = new owned grid(xmin, xmax, npoints, nghosts);
    // Grid.create_grid();
}
