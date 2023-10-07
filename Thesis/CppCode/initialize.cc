void initialize() {

  int t, x, y, site, dir, next_t, next_x, next_y, next_site, i;
#if dim==4 
  int z, next_z;
#endif
#ifdef __F4RootLattice__
  int halfInt, next_halfInt;
#endif

  void (*norm_sun)(dc *v);
  
#if Ncol==2 
  norm_sun=norm_su2; 
#endif

#if Ncol==3 
  norm_sun=norm_su3;
#endif

#if Ncol==4
  norm_sun=norm_su4;
#endif

#if Ncol==5
  norm_sun=norm_su5;
#endif

#if Ncol==6
  norm_sun=norm_su6;
#endif

#if Ncol==7
  norm_sun=norm_su7;
#endif

#if Ncol==8
  norm_sun=norm_su8;
#endif

#if Ncol==9
  norm_sun=norm_su9;
#endif

#if Ncol==10
  norm_sun=norm_su10;
#endif

#ifdef __wanna_Wilson_loops__
  if (nx<=ny) {
    max_Wilson_loop_size=nx/2;
  }
  else {
    max_Wilson_loop_size=ny/2;
  }
#if dim==4  
  if ((nz/2)<max_Wilson_loop_size) {
    max_Wilson_loop_size=nz/2;
  }
#endif
  if ((first_mu==0) && ((nt/2)<max_Wilson_loop_size)) {
    max_Wilson_loop_size=nt/2;
  }
#endif

#ifdef __wanna_Jarzynski_SF__

//   if ((nx!=ny)
// #if dim==4 
//       || (nx!=nz)
// #endif
//   ) {
//     printf("In the present version of the code, the Jarzynski's theorem routines require the spatial sizes of the lattice to be equal\n");
//     exit(0);
//   }
  
#ifdef __wanna_multilevel__
  printf("In the present version of the code, it is not possible to use the multilevel in combination with Jarzynski's theorem\n");
  exit(0);
#endif
  
#ifdef __wanna_improvement__
  printf("In the present version of the code, the Jarzynski's theorem routines do not support the improved gauge action\n");
  exit(0);
#endif

#endif


// Conventions for the site labelling:
//
//   site =          z +
//                y*nz +
//             x*ny*nz +
//          t*nx*ny*nz
//  
// Note that, in particular, on a 10^4 lattice, the site of index
// txyz has coordinates (t,x,y,z).
  
/* For F4 lattice (t,x,y,z) can be either all integer, or all semi-integer.
 * For the first case:
 *   site =      0+2*z +
 *              y*2*nz +
 *           x*ny*2*nz +
 *        t*nx*ny*2*nz
 * For the latter (t,x,y,z)->(t-1/2,x-1/2,y-1/2,z-1/2):
 *   site =      1+2*z +
 *              y*2*nz +
 *           x*ny*2*nz +
 *        t*nx*ny*2*nz
 */
  
#ifdef __MYDEBUGNEIGHBOURS__
  std::cout << "Site\tDir\tneigh+\tneigh-\n";
#endif
  
  for (t=0;t<nt;t++)
  for (x=0;x<nx;x++)
  for (y=0;y<ny;y++)
#if dim==4  
  for (z=0;z<nz;z++)
#endif
#ifdef __F4RootLattice__
  for (halfInt=0;halfInt<2;halfInt++)
#endif
  {
#ifndef __F4RootLattice__
    site=y+ny*(x+nx*t);
  #if dim==4
    site*=nz;
    site+=z;
  #endif
#else
    site=(z+nz*(y+ny*(x+nx*t)))*2+halfInt;
#endif
    for (dir=0;dir<maxDir;dir++) {
#ifdef __wanna_Jarzynski_SF__
      if ((t==(nt-1)) || ((t==0) && (dir!=0))) {
// Links to be held fixed: spatial links at t=0
// and all links at t=nt-1:
        locked_link[dir*nsites+site]=true;
      }
      else {
        locked_link[dir*nsites+site]=false;
      };
#endif
#ifdef __wanna_multilevel__
      if ((t%slab_size==0) && (dir!=0)) {
        locked_link[dir*nsites+site]=true;
      }
      else {
        locked_link[dir*nsites+site]=false;
      };
#ifdef __wanna_improvement__
      if ((t+1)%slab_size==0) {
        locked_link[dir*nsites+site]=true;
      }
#endif
#endif
      next_t=t;
      next_x=x;
      next_y=y;
#if dim==4  
      next_z=z;
#endif
#ifdef __F4RootLattice__
      next_halfInt=halfInt;
#endif
      switch (dir) {
        case 0: {
          next_t=(t+1)%nt;
#ifndef __F4RootLattice__
          next_site=next_y+ny*(next_x+nx*next_t);
#if dim==4  
          next_site*=nz;
          next_site+=next_z;
#endif
#else
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
#endif
          neighbor_plus[dir*nsites+site]=next_site;
          next_t=(t+nt-1)%nt;
#ifndef __F4RootLattice__
          next_site=next_y+ny*(next_x+nx*next_t);
#if dim==4  
          next_site*=nz;
          next_site+=next_z;
#endif
#else
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
#endif
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
        case 1: {
          next_x=(x+1)%nx;
#ifndef __F4RootLattice__
          next_site=next_y+ny*(next_x+nx*next_t);
#if dim==4  
          next_site*=nz;
          next_site+=next_z;
#endif
#else
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
#endif
          neighbor_plus[dir*nsites+site]=next_site;
          next_x=(x+nx-1)%nx;
#ifndef __F4RootLattice__
          next_site=next_y+ny*(next_x+nx*next_t);
#if dim==4  
          next_site*=nz;
          next_site+=next_z;
#endif
#else
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
#endif
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
        case 2: {
          next_y=(y+1)%ny;
#ifndef __F4RootLattice__
          next_site=next_y+ny*(next_x+nx*next_t);
#if dim==4  
          next_site*=nz;
          next_site+=next_z;
#endif
#else
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
#endif
          neighbor_plus[dir*nsites+site]=next_site;
          next_y=(y+ny-1)%ny;
#ifndef __F4RootLattice__
          next_site=next_y+ny*(next_x+nx*next_t);
#if dim==4  
          next_site*=nz;
          next_site+=next_z;
#endif
#else
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
#endif
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
#if dim==4  
        case 3: {
          next_z=(z+1)%nz;
#ifndef __F4RootLattice__
          next_site=next_z+nz*(next_y+ny*(next_x+nx*next_t));
#else
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
#endif
          neighbor_plus[dir*nsites+site]=next_site;
          next_z=(z+nz-1)%nz;
#ifndef __F4RootLattice__
          next_site=next_z+nz*(next_y+ny*(next_x+nx*next_t));
#else
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
#endif
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
#endif
#ifdef __F4RootLattice__
        case 4: { // (+0.5 +0.5 +0.5 +0.5)
          if(halfInt){
            next_halfInt=0;
            next_t=(t+1)%nt;
            next_x=(x+1)%nx;
            next_y=(y+1)%ny;
            next_z=(z+1)%nz;
          }
          else next_halfInt=1;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_plus[dir*nsites+site]=next_site;
          next_halfInt=halfInt;
          next_t=t;
          next_x=x;
          next_y=y;
          next_z=z;
          if(halfInt) next_halfInt=0;
          else{
            next_halfInt=1;
            next_t=(t+nt-1)%nt;
            next_x=(x+nx-1)%nx;
            next_y=(y+ny-1)%ny;
            next_z=(z+nz-1)%nz;
          }
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
        
        case 5: { // (+0.5 +0.5 +0.5 -0.5)
          if(halfInt){
            next_halfInt=0;
            next_t=(t+1)%nt;
            next_x=(x+1)%nx;
            next_y=(y+1)%ny;
          }
          else{
            next_halfInt=1;
            next_z=(z+nz-1)%nz;
          }
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_plus[dir*nsites+site]=next_site;
          next_halfInt=halfInt;
          next_t=t;
          next_x=x;
          next_y=y;
          next_z=z;
          if(halfInt){
            next_halfInt=0;
            next_z=(z+1)%nz;
          }
          else{
            next_halfInt=1;
            next_t=(t+nt-1)%nt;
            next_x=(x+nx-1)%nx;
            next_y=(y+ny-1)%ny;
          }
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
        
        case 6: { // (+0.5 +0.5 -0.5 +0.5)
          if(halfInt){
            next_halfInt=0;
            next_t=(t+1)%nt;
            next_x=(x+1)%nx;
            next_z=(z+1)%nz;
          }
          else{
            next_halfInt=1;
            next_y=(y+ny-1)%ny;
          }
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_plus[dir*nsites+site]=next_site;
          next_halfInt=halfInt;
          next_t=t;
          next_x=x;
          next_y=y;
          next_z=z;
          if(halfInt){
            next_halfInt=0;
            next_y=(y+1)%ny;
          }
          else{
            next_halfInt=1;
            next_t=(t+nt-1)%nt;
            next_x=(x+nx-1)%nx;
            next_z=(z+nz-1)%nz;
          }
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
        
        case 7: { // (+0.5 -0.5 +0.5 +0.5)
          if(halfInt){
            next_halfInt=0;
            next_t=(t+1)%nt;
            next_y=(y+1)%ny;
            next_z=(z+1)%nz;
          }
          else{
            next_halfInt=1;
            next_x=(x+nx-1)%nx;
          }
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_plus[dir*nsites+site]=next_site;
          next_halfInt=halfInt;
          next_t=t;
          next_x=x;
          next_y=y;
          next_z=z;
          if(halfInt){
            next_halfInt=0;
            next_x=(x+1)%nx;
          }
          else{
            next_halfInt=1;
            next_t=(t+nt-1)%nt;
            next_y=(y+ny-1)%ny;
            next_z=(z+nz-1)%nz;
          }
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
        
        case 8: { // (-0.5 +0.5 +0.5 +0.5)
          if(halfInt){
            next_halfInt=0;
            next_x=(x+1)%nx;
            next_y=(y+1)%ny;
            next_z=(z+1)%nz;
          }
          else{
            next_halfInt=1;
            next_t=(t+nt-1)%nt;
          }
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_plus[dir*nsites+site]=next_site;
          next_halfInt=halfInt;
          next_t=t;
          next_x=x;
          next_y=y;
          next_z=z;
          if(halfInt){
            next_halfInt=0;
            next_t=(t+1)%nt;
          }
          else{
            next_halfInt=1;
            next_x=(x+nx-1)%nx;
            next_y=(y+ny-1)%ny;
            next_z=(z+nz-1)%nz;
          }
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
        
        case 9: { // (+0.5 +0.5 -0.5 -0.5)
          if(halfInt){
            next_halfInt=0;
            next_t=(t+1)%nt;
            next_x=(x+1)%nx;
          }
          else{
            next_halfInt=1;
            next_y=(y+ny-1)%ny;
            next_z=(z+nz-1)%nz;
          }
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_plus[dir*nsites+site]=next_site;
          next_halfInt=halfInt;
          next_t=t;
          next_x=x;
          next_y=y;
          next_z=z;
          if(halfInt){
            next_halfInt=0;
            next_y=(y+1)%ny;
            next_z=(z+1)%nz;
          }
          else{
            next_halfInt=1;
            next_t=(t+nt-1)%nt;
            next_x=(x+nx-1)%nx;
          }
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
        
        case 10: { // (+0.5 -0.5 -0.5 +0.5)
          if(halfInt){
            next_halfInt=0;
            next_t=(t+1)%nt;
            next_z=(z+1)%nz;
          }
          else{
            next_halfInt=1;
            next_x=(x+nx-1)%nx;
            next_y=(y+ny-1)%ny;
          }
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_plus[dir*nsites+site]=next_site;
          next_halfInt=halfInt;
          next_t=t;
          next_x=x;
          next_y=y;
          next_z=z;
          if(halfInt){
            next_halfInt=0;
            next_x=(x+1)%nx;
            next_y=(y+1)%ny;
          }
          else{
            next_halfInt=1;
            next_t=(t+nt-1)%nt;
            next_z=(z+nz-1)%nz;
          }
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
        
        case 11: { // (+0.5 -0.5 +0.5 -0.5)
          if(halfInt){
            next_halfInt=0;
            next_t=(t+1)%nt;
            next_y=(y+1)%ny;
          }
          else{
            next_halfInt=1;
            next_x=(x+nx-1)%nx;
            next_z=(z+nz-1)%nz;
          }
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_plus[dir*nsites+site]=next_site;
          next_halfInt=halfInt;
          next_t=t;
          next_x=x;
          next_y=y;
          next_z=z;
          if(halfInt){
            next_halfInt=0;
            next_x=(x+1)%nx;
            next_z=(z+1)%nz;
          }
          else{
            next_halfInt=1;
            next_t=(t+nt-1)%nt;
            next_y=(y+ny-1)%ny;
          }
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
        
        case 12: { // (+1 +1 0 0)
          next_t=(t+1)%nt;
          next_x=(x+1)%nx;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_plus[dir*nsites+site]=next_site;
          next_t=(t+nt-1)%nt;
          next_x=(x+nx-1)%nx;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
        
        case 13: { // (+1 0 +1 0)
          next_t=(t+1)%nt;
          next_y=(y+1)%ny;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_plus[dir*nsites+site]=next_site;
          next_t=(t+nt-1)%nt;
          next_y=(y+ny-1)%ny;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
        
        case 14: { // (+1 0 0 +1)
          next_t=(t+1)%nt;
          next_z=(z+1)%nz;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_plus[dir*nsites+site]=next_site;
          next_t=(t+nt-1)%nt;
          next_z=(z+nz-1)%nz;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
        
        case 15: { // (0 +1 0 +1)
          next_x=(x+1)%nx;
          next_z=(z+1)%nz;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_plus[dir*nsites+site]=next_site;
          next_x=(x+nx-1)%nx;
          next_z=(z+nz-1)%nz;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
        
        case 16: { // (0 +1 +1 0)
          next_x=(x+1)%nx;
          next_y=(y+1)%ny;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_plus[dir*nsites+site]=next_site;
          next_x=(x+nx-1)%nx;
          next_y=(y+ny-1)%ny;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
        
        case 17: { // (0 0 +1 +1)
          next_y=(y+1)%ny;
          next_z=(z+1)%nz;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_plus[dir*nsites+site]=next_site;
          next_y=(y+ny-1)%ny;
          next_z=(z+nz-1)%nz;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
        
        case 18: { // (+1 -1 0 0)
          next_t=(t+1)%nt;
          next_x=(x+nx-1)%nx;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_plus[dir*nsites+site]=next_site;
          next_t=(t+nt-1)%nt;
          next_x=(x+1)%nx;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
        
        case 19: { // (+1 0 -1 0)
          next_t=(t+1)%nt;
          next_y=(y+ny-1)%ny;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_plus[dir*nsites+site]=next_site;
          next_t=(t+nt-1)%nt;
          next_y=(y+1)%ny;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
        
        case 20: { // (+1 0 0 -1)
          next_t=(t+1)%nt;
          next_z=(z+nz-1)%nz;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_plus[dir*nsites+site]=next_site;
          next_t=(t+nt-1)%nt;
          next_z=(z+1)%nz;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
        
        case 21: { // (0 +1 0 -1)
          next_x=(x+1)%nx;
          next_z=(z+nz-1)%nz;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_plus[dir*nsites+site]=next_site;
          next_x=(x+nx-1)%nx;
          next_z=(z+1)%nz;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
        
        case 22: { // (0 +1 -1 0)
          next_x=(x+1)%nx;
          next_y=(y+ny-1)%ny;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_plus[dir*nsites+site]=next_site;
          next_x=(x+nx-1)%nx;
          next_y=(y+1)%ny;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
        
        case 23: { // (0 0 +1 -1)
          next_y=(y+1)%ny;
          next_z=(z+nz-1)%nz;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_plus[dir*nsites+site]=next_site;
          next_y=(y+ny-1)%ny;
          next_z=(z+1)%nz;
          next_site=next_halfInt+2*(next_z+nz*(next_y+ny*(next_x+nx*next_t)));
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
#endif
        default: {
          break;
        }
      }
#ifdef __MYDEBUGNEIGHBOURS__
      std::cout << site << "\t" << dir << "\t" << neighbor_plus[dir*nsites+site];
      std::cout << "\t" << neighbor_minus[dir*nsites+site] << "\n";
#endif
      if (type_of_start==0) {
        for (i=0;i<Ncolsquare;i++) {
          v[i]=dc(1.-2.*random_real(xx),1.-2.*random_real(xx));
        }
        norm_sun(v);
        for (i=0;i<Ncolsquare;i++) {
          ufield[(dir*nsites+site)*Ncolsquare+i]=v[i];
        }
      }

      if (type_of_start==1) {
        for (i=0;i<Ncolsquare;i++) {
          ufield[(dir*nsites+site)*Ncolsquare+i]=dc(0.,0.);
        }
        for (i=0;i<Ncolsquare;i+=Ncol_plus_one) {
          ufield[(dir*nsites+site)*Ncolsquare+i]=dc(1.,0.);
        }
      }
/*
//       if (type_of_start==3) {
// 	if ( ( (site!=0) && (site!=ny) ) || (dir!=0) ) {
//           for (i=0;i<Ncolsquare;i++) {
//             ufield[(dir*nsites+site)*Ncolsquare+i]=dc(0.,0.);
//           }
//           for (i=0;i<Ncolsquare;i+=Ncol_plus_one) {
//             ufield[(dir*nsites+site)*Ncolsquare+i]=dc(1.,0.);
//           }
//         }
//         else {
//           for (i=0;i<Ncolsquare;i++) {
//             ufield[(dir*nsites+site)*Ncolsquare+i]=dc(0.,0.);
//           }
//           for (i=0;i<Ncolsquare;i+=Ncol_plus_one) {
//             ufield[(dir*nsites+site)*Ncolsquare+i]=dc(1.,0.);
//           }
//           ufield[(dir*nsites+site)*Ncolsquare+0]=dc(0.,0.);
//           ufield[(dir*nsites+site)*Ncolsquare+1]=dc(1.,0.);
//           ufield[(dir*nsites+site)*Ncolsquare+Ncol_plus_one]=dc(0.,0.);
//           ufield[(dir*nsites+site)*Ncolsquare+Ncol]=dc(-1.,0.);
// 	}
//       }
*/
    }
  }
  
#ifdef __MYDEBUGSTAPLE__
  // site=(y+ny*(x+nx*t))*2*nz+2*z+halfInt;
  ufield[(0*nsites+0)*Ncolsquare+0] = dc(0., 0.);
  ufield[(0*nsites+0)*Ncolsquare+1] = dc(0., 0.);
  ufield[(0*nsites+0)*Ncolsquare+2] = dc(0., 0.);
  ufield[(0*nsites+0)*Ncolsquare+3] = dc(0., 0.);
  /*
  ufield[(0*nsites+0)*Ncolsquare+0] = dc(1.0/sqrt(2), 0.);
  ufield[(0*nsites+0)*Ncolsquare+1] = dc(1.0/sqrt(2), 0.);
  ufield[(0*nsites+0)*Ncolsquare+2] = dc(-1.0/sqrt(2), 0.);
  ufield[(0*nsites+0)*Ncolsquare+3] = dc(1.0/sqrt(2), 0.);
  * */
  /*
  int site1=(2+ny*(2+nx*2))*2*nz+2*2+1;
  ufield[(0*nsites+site1)*Ncolsquare+0] = dc(0.5, 0.);
  ufield[(0*nsites+site1)*Ncolsquare+1] = dc(0., 0.5*sqrt(3));
  ufield[(0*nsites+site1)*Ncolsquare+2] = dc(0., 0.5*sqrt(3));
  ufield[(0*nsites+site1)*Ncolsquare+3] = dc(0.5, 0.);
  */
#endif
  
#ifdef __F4RootLattice__ // initializing vector with all possible staples in each direction
  for(dir=0; dir<maxDir; dir++) switch(dir){
    case 0:{
      stapleDirs[dir*nStaples+ 0]= +4; stapleDirs[dir*nStaples+ 1]= -8;
      stapleDirs[dir*nStaples+ 2]= +5; stapleDirs[dir*nStaples+ 3]=+10;
      stapleDirs[dir*nStaples+ 4]= +6; stapleDirs[dir*nStaples+ 5]=+11;
      stapleDirs[dir*nStaples+ 6]= +7; stapleDirs[dir*nStaples+ 7]= +9;
      /*
      stapleDirs[dir*nStaples+ 8]= +1; stapleDirs[dir*nStaples+ 9]=-12;
      stapleDirs[dir*nStaples+10]= +2; stapleDirs[dir*nStaples+11]=-13;
      stapleDirs[dir*nStaples+12]= +3; stapleDirs[dir*nStaples+13]=-14;
      stapleDirs[dir*nStaples+14]= -1; stapleDirs[dir*nStaples+15]=-18;
      stapleDirs[dir*nStaples+16]= -2; stapleDirs[dir*nStaples+17]=-19;
      stapleDirs[dir*nStaples+18]= -3; stapleDirs[dir*nStaples+19]=-20;
      */
      break;
    }
    case 1:{
      stapleDirs[dir*nStaples+ 0]= +4; stapleDirs[dir*nStaples+ 1]= -7;
      stapleDirs[dir*nStaples+ 2]= +5; stapleDirs[dir*nStaples+ 3]=-11;
      stapleDirs[dir*nStaples+ 4]= +6; stapleDirs[dir*nStaples+ 5]=-10;
      stapleDirs[dir*nStaples+ 6]= +8; stapleDirs[dir*nStaples+ 7]= +9;
      /*
      stapleDirs[dir*nStaples+ 8]=+24; stapleDirs[dir*nStaples+ 9]=-12;
      stapleDirs[dir*nStaples+10]= +2; stapleDirs[dir*nStaples+11]=-16;
      stapleDirs[dir*nStaples+12]= +3; stapleDirs[dir*nStaples+13]=-15;
      stapleDirs[dir*nStaples+14]=-24; stapleDirs[dir*nStaples+15]=+18;
      stapleDirs[dir*nStaples+16]= -2; stapleDirs[dir*nStaples+17]=-22;
      stapleDirs[dir*nStaples+18]= -3; stapleDirs[dir*nStaples+19]=-21;
      */
      break;
    }
    case 2:{
      stapleDirs[dir*nStaples+ 0]= +4; stapleDirs[dir*nStaples+ 1]= -6;
      stapleDirs[dir*nStaples+ 2]= +5; stapleDirs[dir*nStaples+ 3]= -9;
      stapleDirs[dir*nStaples+ 4]= +7; stapleDirs[dir*nStaples+ 5]=-10;
      stapleDirs[dir*nStaples+ 6]= +8; stapleDirs[dir*nStaples+ 7]=+11;
      /*
      stapleDirs[dir*nStaples+ 8]=+24; stapleDirs[dir*nStaples+ 9]=-13;
      stapleDirs[dir*nStaples+10]= +1; stapleDirs[dir*nStaples+11]=-16;
      stapleDirs[dir*nStaples+12]= +3; stapleDirs[dir*nStaples+13]=-17;
      stapleDirs[dir*nStaples+14]=-24; stapleDirs[dir*nStaples+15]=+19;
      stapleDirs[dir*nStaples+16]= -1; stapleDirs[dir*nStaples+17]=+22;
      stapleDirs[dir*nStaples+18]= -3; stapleDirs[dir*nStaples+19]=-23;
      */
      break;
    }
    case 3:{
      stapleDirs[dir*nStaples+ 0]= +4; stapleDirs[dir*nStaples+ 1]= -5;
      stapleDirs[dir*nStaples+ 2]= +6; stapleDirs[dir*nStaples+ 3]= -9;
      stapleDirs[dir*nStaples+ 4]= +7; stapleDirs[dir*nStaples+ 5]=-11;
      stapleDirs[dir*nStaples+ 6]= +8; stapleDirs[dir*nStaples+ 7]=+10;
      /*
      stapleDirs[dir*nStaples+ 8]=+24; stapleDirs[dir*nStaples+ 9]=-14;
      stapleDirs[dir*nStaples+10]= +1; stapleDirs[dir*nStaples+11]=-15;
      stapleDirs[dir*nStaples+12]= +2; stapleDirs[dir*nStaples+13]=-17;
      stapleDirs[dir*nStaples+14]=-24; stapleDirs[dir*nStaples+15]=+20;
      stapleDirs[dir*nStaples+16]= -1; stapleDirs[dir*nStaples+17]=+21;
      stapleDirs[dir*nStaples+18]= -2; stapleDirs[dir*nStaples+19]=+23;
      */
      break;
    }
    case 4:{
      stapleDirs[dir*nStaples+ 0]=+24; stapleDirs[dir*nStaples+ 1]= +8;
      stapleDirs[dir*nStaples+ 2]= +1; stapleDirs[dir*nStaples+ 3]= +7;
      stapleDirs[dir*nStaples+ 4]= +2; stapleDirs[dir*nStaples+ 5]= +6;
      stapleDirs[dir*nStaples+ 6]= +3; stapleDirs[dir*nStaples+ 7]= +5;
      /*
      stapleDirs[dir*nStaples+ 8]= +9; stapleDirs[dir*nStaples+ 9]=-12;
      stapleDirs[dir*nStaples+10]=+10; stapleDirs[dir*nStaples+11]=-14;
      stapleDirs[dir*nStaples+12]=+11; stapleDirs[dir*nStaples+13]=-13;
      stapleDirs[dir*nStaples+14]= -9; stapleDirs[dir*nStaples+15]=-17;
      stapleDirs[dir*nStaples+16]=-10; stapleDirs[dir*nStaples+17]=-16;
      stapleDirs[dir*nStaples+18]=-11; stapleDirs[dir*nStaples+19]=-15;
      */
      break;
    }
    case 5:{
      stapleDirs[dir*nStaples+ 0]=+24; stapleDirs[dir*nStaples+ 1]=-10;
      stapleDirs[dir*nStaples+ 2]= +1; stapleDirs[dir*nStaples+ 3]=+11;
      stapleDirs[dir*nStaples+ 4]= +2; stapleDirs[dir*nStaples+ 5]= +9;
      stapleDirs[dir*nStaples+ 6]= -3; stapleDirs[dir*nStaples+ 7]= +4;
      /*
      stapleDirs[dir*nStaples+ 8]= +6; stapleDirs[dir*nStaples+ 9]=-12;
      stapleDirs[dir*nStaples+10]= +7; stapleDirs[dir*nStaples+11]=-13;
      stapleDirs[dir*nStaples+12]= +8; stapleDirs[dir*nStaples+13]=-16;
      stapleDirs[dir*nStaples+14]= -6; stapleDirs[dir*nStaples+15]=-23;
      stapleDirs[dir*nStaples+16]= -7; stapleDirs[dir*nStaples+17]=-21;
      stapleDirs[dir*nStaples+18]= -8; stapleDirs[dir*nStaples+19]=-20;
      */
      break;
    }
    case 6:{
      stapleDirs[dir*nStaples+ 0]=+24; stapleDirs[dir*nStaples+ 1]=-11;
      stapleDirs[dir*nStaples+ 2]= +1; stapleDirs[dir*nStaples+ 3]=+10;
      stapleDirs[dir*nStaples+ 4]= -2; stapleDirs[dir*nStaples+ 5]= +4;
      stapleDirs[dir*nStaples+ 6]= +3; stapleDirs[dir*nStaples+ 7]= +9;
      /*
      stapleDirs[dir*nStaples+ 8]= +5; stapleDirs[dir*nStaples+ 9]=-12;
      stapleDirs[dir*nStaples+10]= +7; stapleDirs[dir*nStaples+11]=-14;
      stapleDirs[dir*nStaples+12]= +8; stapleDirs[dir*nStaples+13]=-15;
      stapleDirs[dir*nStaples+14]= -5; stapleDirs[dir*nStaples+15]=+23;
      stapleDirs[dir*nStaples+16]= -7; stapleDirs[dir*nStaples+17]=-22;
      stapleDirs[dir*nStaples+18]= -8; stapleDirs[dir*nStaples+19]=-19;
      */
      break;
    }
    case 7:{
      stapleDirs[dir*nStaples+ 0]=+24; stapleDirs[dir*nStaples+ 1]= -9;
      stapleDirs[dir*nStaples+ 2]= -1; stapleDirs[dir*nStaples+ 3]= +4;
      stapleDirs[dir*nStaples+ 4]= +2; stapleDirs[dir*nStaples+ 5]=+10;
      stapleDirs[dir*nStaples+ 6]= +3; stapleDirs[dir*nStaples+ 7]=+11;
      /*
      stapleDirs[dir*nStaples+ 8]= +5; stapleDirs[dir*nStaples+ 9]=-13;
      stapleDirs[dir*nStaples+10]= +6; stapleDirs[dir*nStaples+11]=-14;
      stapleDirs[dir*nStaples+12]= +8; stapleDirs[dir*nStaples+13]=-17;
      stapleDirs[dir*nStaples+14]= -5; stapleDirs[dir*nStaples+15]=+21;
      stapleDirs[dir*nStaples+16]= -6; stapleDirs[dir*nStaples+17]=+22;
      stapleDirs[dir*nStaples+18]= -8; stapleDirs[dir*nStaples+19]=-18;
      */
      break;
    }
    case 8:{
      stapleDirs[dir*nStaples+ 0]=-24; stapleDirs[dir*nStaples+ 1]= +4;
      stapleDirs[dir*nStaples+ 2]= +1; stapleDirs[dir*nStaples+ 3]= -9;
      stapleDirs[dir*nStaples+ 4]= +2; stapleDirs[dir*nStaples+ 5]=-11;
      stapleDirs[dir*nStaples+ 6]= +3; stapleDirs[dir*nStaples+ 7]=-10;
      /*
      stapleDirs[dir*nStaples+ 8]= +5; stapleDirs[dir*nStaples+ 9]=-16;
      stapleDirs[dir*nStaples+10]= +6; stapleDirs[dir*nStaples+11]=-15;
      stapleDirs[dir*nStaples+12]= +7; stapleDirs[dir*nStaples+13]=-17;
      stapleDirs[dir*nStaples+14]= -5; stapleDirs[dir*nStaples+15]=+20;
      stapleDirs[dir*nStaples+16]= -6; stapleDirs[dir*nStaples+17]=+19;
      stapleDirs[dir*nStaples+18]= -7; stapleDirs[dir*nStaples+19]=+18;
      */
      break;
    }
    case 9:{
      stapleDirs[dir*nStaples+ 0]=+24; stapleDirs[dir*nStaples+ 1]= -7;
      stapleDirs[dir*nStaples+ 2]= +1; stapleDirs[dir*nStaples+ 3]= -8;
      stapleDirs[dir*nStaples+ 4]= -2; stapleDirs[dir*nStaples+ 5]= +5;
      stapleDirs[dir*nStaples+ 6]= -3; stapleDirs[dir*nStaples+ 7]= +6;
      /*
      stapleDirs[dir*nStaples+ 8]= +4; stapleDirs[dir*nStaples+ 9]=-12;
      stapleDirs[dir*nStaples+10]=+10; stapleDirs[dir*nStaples+11]=-19;
      stapleDirs[dir*nStaples+12]=+11; stapleDirs[dir*nStaples+13]=-20;
      stapleDirs[dir*nStaples+14]= -4; stapleDirs[dir*nStaples+15]=+17;
      stapleDirs[dir*nStaples+16]=-10; stapleDirs[dir*nStaples+17]=-21;
      stapleDirs[dir*nStaples+18]=-11; stapleDirs[dir*nStaples+19]=-22;
      */
      break;
    }
    case 10:{
      stapleDirs[dir*nStaples+ 0]=+24; stapleDirs[dir*nStaples+ 1]= -5;
      stapleDirs[dir*nStaples+ 2]= -1; stapleDirs[dir*nStaples+ 3]= +6;
      stapleDirs[dir*nStaples+ 4]= -2; stapleDirs[dir*nStaples+ 5]= +7;
      stapleDirs[dir*nStaples+ 6]= +3; stapleDirs[dir*nStaples+ 7]= -8;
      /*
      stapleDirs[dir*nStaples+ 8]= +4; stapleDirs[dir*nStaples+ 9]=-14;
      stapleDirs[dir*nStaples+10]= +9; stapleDirs[dir*nStaples+11]=-19;
      stapleDirs[dir*nStaples+12]=+11; stapleDirs[dir*nStaples+13]=-18;
      stapleDirs[dir*nStaples+14]= -4; stapleDirs[dir*nStaples+15]=+16;
      stapleDirs[dir*nStaples+16]= -9; stapleDirs[dir*nStaples+17]=+21;
      stapleDirs[dir*nStaples+18]=-11; stapleDirs[dir*nStaples+19]=+23;
      */
      break;
    }
    case 11:{
      stapleDirs[dir*nStaples+ 0]=+24; stapleDirs[dir*nStaples+ 1]= -6;
      stapleDirs[dir*nStaples+ 2]= -1; stapleDirs[dir*nStaples+ 3]= +5;
      stapleDirs[dir*nStaples+ 4]= +2; stapleDirs[dir*nStaples+ 5]= -8;
      stapleDirs[dir*nStaples+ 6]= -3; stapleDirs[dir*nStaples+ 7]= +7;
      /*
      stapleDirs[dir*nStaples+ 8]= +4; stapleDirs[dir*nStaples+ 9]=-13;
      stapleDirs[dir*nStaples+10]= +9; stapleDirs[dir*nStaples+11]=-20;
      stapleDirs[dir*nStaples+12]=+10; stapleDirs[dir*nStaples+13]=-18;
      stapleDirs[dir*nStaples+14]= -4; stapleDirs[dir*nStaples+15]=+15;
      stapleDirs[dir*nStaples+16]= -9; stapleDirs[dir*nStaples+17]=+22;
      stapleDirs[dir*nStaples+18]=-10; stapleDirs[dir*nStaples+19]=-23;
      */
      break;
    }
    /*
    case 12:{
      stapleDirs[dir*nStaples+ 0]=-24; stapleDirs[dir*nStaples+ 1]= -1;
      stapleDirs[dir*nStaples+ 2]= -4; stapleDirs[dir*nStaples+ 3]= -9;
      stapleDirs[dir*nStaples+ 4]= -5; stapleDirs[dir*nStaples+ 5]= -6;
      stapleDirs[dir*nStaples+ 6]=-13; stapleDirs[dir*nStaples+ 7]=-22;
      stapleDirs[dir*nStaples+ 8]=-14; stapleDirs[dir*nStaples+ 9]=-21;
      stapleDirs[dir*nStaples+10]=-15; stapleDirs[dir*nStaples+11]=-20;
      stapleDirs[dir*nStaples+12]=-16; stapleDirs[dir*nStaples+13]=-19;
      for(int i=14;i<nStaples;i++) stapleDirs[dir*nStaples+i]=0;
      break;
    }
    case 13:{
      stapleDirs[dir*nStaples+ 0]=-24; stapleDirs[dir*nStaples+ 1]= -2;
      stapleDirs[dir*nStaples+ 2]= -4; stapleDirs[dir*nStaples+ 3]=-11;
      stapleDirs[dir*nStaples+ 4]= -5; stapleDirs[dir*nStaples+ 5]= -7;
      stapleDirs[dir*nStaples+ 6]=-12; stapleDirs[dir*nStaples+ 7]=+22;
      stapleDirs[dir*nStaples+ 8]=-14; stapleDirs[dir*nStaples+ 9]=-23;
      stapleDirs[dir*nStaples+10]=-16; stapleDirs[dir*nStaples+11]=-18;
      stapleDirs[dir*nStaples+12]=-17; stapleDirs[dir*nStaples+13]=-20;
      for(int i=14;i<nStaples;i++) stapleDirs[dir*nStaples+i]=0;
      break;
    }
    case 14:{
      stapleDirs[dir*nStaples+ 0]=-24; stapleDirs[dir*nStaples+ 1]= -3;
      stapleDirs[dir*nStaples+ 2]= -4; stapleDirs[dir*nStaples+ 3]=-10;
      stapleDirs[dir*nStaples+ 4]= -6; stapleDirs[dir*nStaples+ 5]= -7;
      stapleDirs[dir*nStaples+ 6]=-12; stapleDirs[dir*nStaples+ 7]=+21;
      stapleDirs[dir*nStaples+ 8]=-13; stapleDirs[dir*nStaples+ 9]=+23;
      stapleDirs[dir*nStaples+10]=-15; stapleDirs[dir*nStaples+11]=-18;
      stapleDirs[dir*nStaples+12]=-17; stapleDirs[dir*nStaples+13]=-19;
      for(int i=14;i<nStaples;i++) stapleDirs[dir*nStaples+i]=0;
      break;
    }
    case 15:{
      stapleDirs[dir*nStaples+ 0]= -1; stapleDirs[dir*nStaples+ 1]= -3;
      stapleDirs[dir*nStaples+ 2]= -4; stapleDirs[dir*nStaples+ 3]=+11;
      stapleDirs[dir*nStaples+ 4]= -6; stapleDirs[dir*nStaples+ 5]= -8;
      stapleDirs[dir*nStaples+ 6]=-12; stapleDirs[dir*nStaples+ 7]=+20;
      stapleDirs[dir*nStaples+ 8]=-14; stapleDirs[dir*nStaples+ 9]=+18;
      stapleDirs[dir*nStaples+10]=-16; stapleDirs[dir*nStaples+11]=+23;
      stapleDirs[dir*nStaples+12]=-17; stapleDirs[dir*nStaples+13]=-22;
      for(int i=14;i<nStaples;i++) stapleDirs[dir*nStaples+i]=0;
      break;
    }
    case 16:{
      stapleDirs[dir*nStaples+ 0]= -1; stapleDirs[dir*nStaples+ 1]= -2;
      stapleDirs[dir*nStaples+ 2]= -4; stapleDirs[dir*nStaples+ 3]=+10;
      stapleDirs[dir*nStaples+ 4]= -5; stapleDirs[dir*nStaples+ 5]= -8;
      stapleDirs[dir*nStaples+ 6]=-12; stapleDirs[dir*nStaples+ 7]=+19;
      stapleDirs[dir*nStaples+ 8]=-13; stapleDirs[dir*nStaples+ 9]=+18;
      stapleDirs[dir*nStaples+10]=-15; stapleDirs[dir*nStaples+11]=-23;
      stapleDirs[dir*nStaples+12]=-17; stapleDirs[dir*nStaples+13]=-21;
      for(int i=14;i<nStaples;i++) stapleDirs[dir*nStaples+i]=0;
      break;
    }
    case 17:{
      stapleDirs[dir*nStaples+ 0]= -2; stapleDirs[dir*nStaples+ 1]= -3;
      stapleDirs[dir*nStaples+ 2]= -4; stapleDirs[dir*nStaples+ 3]= +9;
      stapleDirs[dir*nStaples+ 4]= -7; stapleDirs[dir*nStaples+ 5]= -8;
      stapleDirs[dir*nStaples+ 6]=-13; stapleDirs[dir*nStaples+ 7]=+20;
      stapleDirs[dir*nStaples+ 8]=-14; stapleDirs[dir*nStaples+ 9]=+19;
      stapleDirs[dir*nStaples+10]=-15; stapleDirs[dir*nStaples+11]=+22;
      stapleDirs[dir*nStaples+12]=-16; stapleDirs[dir*nStaples+13]=+21;
      for(int i=14;i<nStaples;i++) stapleDirs[dir*nStaples+i]=0;
      break;
    }
    case 18:{
      stapleDirs[dir*nStaples+ 0]=-24; stapleDirs[dir*nStaples+ 1]= +1;
      stapleDirs[dir*nStaples+ 2]= -7; stapleDirs[dir*nStaples+ 3]= +8;
      stapleDirs[dir*nStaples+ 4]=-10; stapleDirs[dir*nStaples+ 5]=-11;
      stapleDirs[dir*nStaples+ 6]=-13; stapleDirs[dir*nStaples+ 7]=+16;
      stapleDirs[dir*nStaples+ 8]=-14; stapleDirs[dir*nStaples+ 9]=+15;
      stapleDirs[dir*nStaples+10]=-19; stapleDirs[dir*nStaples+11]=+22;
      stapleDirs[dir*nStaples+12]=-20; stapleDirs[dir*nStaples+13]=+21;
      for(int i=14;i<nStaples;i++) stapleDirs[dir*nStaples+i]=0;
      break;
    }
    case 19:{
      stapleDirs[dir*nStaples+ 0]=-24; stapleDirs[dir*nStaples+ 1]= +2;
      stapleDirs[dir*nStaples+ 2]= -6; stapleDirs[dir*nStaples+ 3]= +8;
      stapleDirs[dir*nStaples+ 4]= -9; stapleDirs[dir*nStaples+ 5]=-10;
      stapleDirs[dir*nStaples+ 6]=-12; stapleDirs[dir*nStaples+ 7]=+16;
      stapleDirs[dir*nStaples+ 8]=-14; stapleDirs[dir*nStaples+ 9]=+17;
      stapleDirs[dir*nStaples+10]=-18; stapleDirs[dir*nStaples+11]=-22;
      stapleDirs[dir*nStaples+12]=-20; stapleDirs[dir*nStaples+13]=+23;
      for(int i=14;i<nStaples;i++) stapleDirs[dir*nStaples+i]=0;
      break;
    }
    case 20:{
      stapleDirs[dir*nStaples+ 0]=-24; stapleDirs[dir*nStaples+ 1]= +3;
      stapleDirs[dir*nStaples+ 2]= -5; stapleDirs[dir*nStaples+ 3]= +8;
      stapleDirs[dir*nStaples+ 4]= -9; stapleDirs[dir*nStaples+ 5]=-11;
      stapleDirs[dir*nStaples+ 6]=-12; stapleDirs[dir*nStaples+ 7]=+15;
      stapleDirs[dir*nStaples+ 8]=-13; stapleDirs[dir*nStaples+ 9]=+17;
      stapleDirs[dir*nStaples+10]=-18; stapleDirs[dir*nStaples+11]=-21;
      stapleDirs[dir*nStaples+12]=-19; stapleDirs[dir*nStaples+13]=-23;
      for(int i=14;i<nStaples;i++) stapleDirs[dir*nStaples+i]=0;
      break;
    }
    case 21:{
      stapleDirs[dir*nStaples+ 0]= -1; stapleDirs[dir*nStaples+ 1]= +3;
      stapleDirs[dir*nStaples+ 2]= -5; stapleDirs[dir*nStaples+ 3]= +7;
      stapleDirs[dir*nStaples+ 4]= -9; stapleDirs[dir*nStaples+ 5]=+10;
      stapleDirs[dir*nStaples+ 6]=-12; stapleDirs[dir*nStaples+ 7]=+14;
      stapleDirs[dir*nStaples+ 8]=-16; stapleDirs[dir*nStaples+ 9]=+17;
      stapleDirs[dir*nStaples+10]=+18; stapleDirs[dir*nStaples+11]=-20;
      stapleDirs[dir*nStaples+12]=-22; stapleDirs[dir*nStaples+13]=-23;
      for(int i=14;i<nStaples;i++) stapleDirs[dir*nStaples+i]=0;
      break;
    }
    case 22:{
      stapleDirs[dir*nStaples+ 0]= -1; stapleDirs[dir*nStaples+ 1]= +2;
      stapleDirs[dir*nStaples+ 2]= -6; stapleDirs[dir*nStaples+ 3]= +7;
      stapleDirs[dir*nStaples+ 4]= -9; stapleDirs[dir*nStaples+ 5]=+11;
      stapleDirs[dir*nStaples+ 6]=-12; stapleDirs[dir*nStaples+ 7]=+13;
      stapleDirs[dir*nStaples+ 8]=-15; stapleDirs[dir*nStaples+ 9]=+17;
      stapleDirs[dir*nStaples+10]=+18; stapleDirs[dir*nStaples+11]=-19;
      stapleDirs[dir*nStaples+12]=-21; stapleDirs[dir*nStaples+13]=+23;
      for(int i=14;i<nStaples;i++) stapleDirs[dir*nStaples+i]=0;
      break;
    }
    case 23:{
      stapleDirs[dir*nStaples+ 0]= -2; stapleDirs[dir*nStaples+ 1]= +3;
      stapleDirs[dir*nStaples+ 2]= -5; stapleDirs[dir*nStaples+ 3]= +6;
      stapleDirs[dir*nStaples+ 4]=+10; stapleDirs[dir*nStaples+ 5]=-11;
      stapleDirs[dir*nStaples+ 6]=-13; stapleDirs[dir*nStaples+ 7]=+14;
      stapleDirs[dir*nStaples+ 8]=+15; stapleDirs[dir*nStaples+ 9]=-16;
      stapleDirs[dir*nStaples+10]=+19; stapleDirs[dir*nStaples+11]=-20;
      stapleDirs[dir*nStaples+12]=-21; stapleDirs[dir*nStaples+13]=+22;
      for(int i=14;i<nStaples;i++) stapleDirs[dir*nStaples+i]=0;
      break;
    }
    */
    default:{
      break;
    }
  }
#endif
  
#ifdef __wanna_Jarzynski_SF__

  int ind;
  double one_over_nt_minus_one=1./(nt-1);
  for (i=0; i<Ncolsquare; i++) {
    C0[i]=dc(0.,0.);
    C1[i]=dc(0.,0.);
  }
  for (i=0; i<Ncol; i++) {
    phi0[i]=largephi0[i]*division_by_L;
    newphi0[i]=newlargephi0[i]*division_by_L;
    phi1[i]=largephi1[i]*division_by_L;
    newphi1[i]=newlargephi1[i]*division_by_L;
    C0[i*(Ncol_plus_one)]=exp(dc(0.,phi0[i]));
    C1[i*(Ncol_plus_one)]=exp(dc(0.,phi1[i]));
  }

  for (x=0;x<nx;x++)
  for (y=0;y<ny;y++)
#if dim==4  
  for (z=0;z<nz;z++)
#endif
  {
    site=(y+ny*(x+nx*(nt-1)));
#if dim==4
    site*=nz;
    site+=z;
#endif
    for (dir=0;dir<dim;dir++) {
        ind=(dir*nsites+site)*Ncolsquare;
// Set all links in t=nt-1 slice to C1:
        for (i=0;i<Ncolsquare;i++) {
          ufield[ind+i]=C1[i];
        }
    }
    site=y+ny*x;
#if dim==4
    site*=nz;
    site+=z;
#endif
    for (int dir=1;dir<dim;dir++) {
        ind=(dir*nsites+site)*Ncolsquare;
// Set spatial links in t=0 slice to C0:
        for (i=0;i<Ncolsquare;i++) {
          ufield[ind+i]=C0[i];
        }
// Set spatial links in intermediate-t slices to matrices interpolating between C0 and C1:
        for (int t=1;t<nt-1;t++) {
          for (i=0;i<Ncolsquare;i++) {
            ufield[ind+i+t*spatial_volume*Ncolsquare]=dc(0.,0.);
          }
          for (i=0;i<Ncolsquare;i+=Ncol+1) {
            ufield[ind+i+t*spatial_volume*Ncolsquare]=
              exp(dc(0.,one_over_nt_minus_one*((nt-1-t)*phi0[i]+t*phi1[i])));;
          }
        }
    }
  }

#endif

  if (type_of_start==2) {
    read_last_conf();
//     char inputfilename[maximum_filename_length];
//     FILE *inputfile;
//     size_t fread_result;
// #if dim==3
//     sprintf(inputfilename,"saved_conf_Ncol_%d_nt_%d_nx_%d_ny_%d_beta_%12.10lf_thread_%d.dat", Ncol, nt, nx, ny, beta, thread_number);
// #endif
// #if dim==4
//     sprintf(inputfilename,"saved_conf_Ncol_%d_nt_%d_nx_%d_ny_%d_nz_%d_beta_%12.10lf_thread_%d.dat", Ncol, nt, nx, ny, nz, beta, thread_number);
// #endif
//     inputfile=fopen(inputfilename,"r");
//     fread_result=fread( ufield, sizeof(ufield), 1, inputfile);
//     fclose(inputfile);
  }

}
