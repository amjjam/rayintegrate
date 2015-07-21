/*********************************************************************
 * This program tests rayintegrate by creating a simple model which  *
 * is just a sphere of radius 5 with a hollow center of radius 1.    *
 * The sphere has uniform brightness except that the +X portion is   *
 * twice as bright as the -X portion.                                *
 * They are imaged from coordinate (10, -60, +10) with direction     *
 * vector (x) toward the origin, y vector (0,0,1), z vector (-1,0,0) *
 * I think that strictly speaking this is a problem because I think  *
 * y and z should be perpendicular to x which they are not in this   *
 * case. But it will just lead to a minor distortion of the image    *
 *********************************************************************/

#include <stdio.h>

#include "../../include/rayintegrate.H"

class MODEL: public RAYINTEGRATEMODEL{
public:
  MODEL(){};
  ~MODEL(){};
  double get(std::vector<double>);
private:
};

double MODEL::get(std::vector<double> p){
  double r=0;
  for(int i=0;i<3;i++)
    r+=p[i]*p[i];
  r=sqrt(r);

  if(r<=6&&r>=1){
    if(p[0]>0){
      if(p[2]<0)
	return 3;
      return 2;
    }
    return 1;
  }
  
  // Return 0 if not inside one of the regions
  return 0;
}



int main(int argc, char *argv[]){
  MODEL m;
  RAYINTEGRATE r(&m);

  std::vector<double> pos(3);
  pos[0]=10; pos[1]=-60; pos[2]=10;
  std::vector<double> x(3);
  x[0]=-10; x[1]=60; x[2]=-10;
  std::vector<double> y(3);
  y[0]=0; y[1]=0; y[2]=1;
  std::vector<double> z(3);
  z[0]=-1; z[1]=0; z[2]=0;
  double startdistance=30,stopdistance=90,delta=0.1;
  double fovy=20,fovz=20;
  int ny=150,nz=150;

  std::vector<std::vector<double> > img=
    r.image(pos,x,y,z,startdistance,stopdistance,delta,fovy,fovz,ny,nz);

  FILE *fp=fopen("../dat/test.dat","w");
  fwrite(&ny,sizeof(int),1,fp);
  fwrite(&nz,sizeof(int),1,fp);
  for(int i=0;i<ny;i++)
    fwrite(&img[i][0],sizeof(double),nz,fp);
  fclose(fp);

  return 0;
}
