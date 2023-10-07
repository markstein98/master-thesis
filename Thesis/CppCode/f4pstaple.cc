void f4pstaple(dc *stot, int site, int mu, int first_direction, int use_smeared) {

// Computes the triangular staples from (site) to (site+a mu) in the positive 
// directions (starting from first_direction) and stores their sum into stot.
    
  int i, j, nextSite, vertex, dir1, dir2, aux_ind1, aux_ind2;
  bool dagger1, dagger2;
  
  nextSite=neighbor_plus[mu*nsites+site];

  for (i=0;i<Ncolsquare;i++) {
    stot[i]=dc(0.,0.);
  }
  
  for(j=0; j<nStaples; j+=2){
    dir1 = stapleDirs[mu*nStaples+j];
    dir2 = stapleDirs[mu*nStaples+j+1];
    if(dir1==0 || dir2==0) break;
    if(dir1>0){
      dagger1 = false;
      dir1 = dir1%24; // because direction 0 was indicated as 24
      vertex = neighbor_plus[dir1*nsites+site];
      aux_ind1 = (dir1*nsites+site)*Ncolsquare;
    }
    else{
      dagger1 = true;
      dir1 = -dir1%24;
      vertex = neighbor_minus[dir1*nsites+site];
      aux_ind1 = (dir1*nsites+vertex)*Ncolsquare;
    }
    if(dir2>0){
      dagger2 = false;
      dir2 = dir2%24;
      aux_ind2 = (dir2*nsites+vertex)*Ncolsquare;
#ifdef __DebugStaple__
      if(nextSite!=neighbor_plus[dir2*nsites+vertex]){
        fprintf(stderr, "[ERROR]: pStaple does not close: site=%d, mu=%d, dir1=", site, mu);
        if(dagger1) fprintf(stderr, "-");
        fprintf(stderr, "%d, dir2=-%d\n", dir1, dir2);
        fprintf(stderr, "         Link points to %d, staple ends in %d\n", nextSite, neighbor_plus[dir2*nsites+vertex]);
        exit(1);
      }
#endif
    }
    else{
      dagger2 = true;
      dir2 = -dir2%24;
      aux_ind2 = (dir2*nsites+nextSite)*Ncolsquare;
#ifdef __DebugStaple__
      if(nextSite!=neighbor_minus[dir2*nsites+vertex]){
        fprintf(stderr, "[ERROR]: pStaple does not close: site=%d, mu=%d, dir1=", site, mu);
        if(dagger1) fprintf(stderr, "-");
        fprintf(stderr, "%d, dir2=-%d\n", dir1, dir2);
        fprintf(stderr, "         Link points to %d, staple ends in %d\n", nextSite, neighbor_minus[dir2*nsites+vertex]);
        exit(1);
      }
#endif
    }
    for(i=0; i<Ncolsquare; i++){
      u1[i] = ufield[aux_ind1+i];
      u2[i] = ufield[aux_ind2+i];
    }
    
    if(dagger1){
      if(dagger2){ // u1^dag * u2^dag = (u2 * u1)^dag
        mult_C_equals_AB(u3,u2,u1);
        hermitian_conjugate(product,u3);
      }
      else{ // u1^dag * u2
        mult_C_equals_AdaggerB(product,u1,u2);
      }
    }
    else{
      if(dagger2){ // u1 * u2^dag
        mult_C_equals_ABdagger(product,u1,u2);
      }
      else{ // u1 * u2
        mult_C_equals_AB(product,u1,u2);
      }
    }
    
    for(i=0; i<Ncolsquare; i++) stot[i] += product[i];
  }
}
