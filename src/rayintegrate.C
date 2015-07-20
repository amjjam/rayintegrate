#include "../include/rayintegrate.H"

RAYINTEGRATE::RAYINTEGRATE(MODEL *m){
  RAYINTEGRATE::m=m;
}


RAYINTEGRATE::~RAYINTEGRATE(){

}


/*====================================================================
  double integrate(std::vector<double> start, std::vector<double>
  direction, double startdistance, double stopdistance, double delta)
  
  integrate the model m along the line specified by start and
  direction, at spatial steps delta
  ===================================================================*/
double RAYINTEGRATE::integrate(std::vector<double> start, 
			       std::vector<double> direction,
			       double startdistance, 
			       double stopdistance,
			       double delta){
  int i;
  double r;
  
  verify(p);
  verify(direction);

  // Position
  std::vector<double> p=start;
  
  // Step vector
  std::vector<double> d=direction;
  double f;
  for(i=0;i<3;i++)
    f+=direction[i]*direction[i];
  f=sqrt(f)*delta;
  for(i=0;i<3;i++)
    d[i]*=f;
  
  for(r=0;distance(start,p)<stopdistance;addvector(p,d)){
    if(distance(start,p)>=startdistance)
      r+=m->(p);
  }
  
  return r;
}

/*============================================================================
  vector<vector<double> > imagestd::vector<double> pos,
  std::vector<double> x, std::vector<double> y, std::vector<double> z,
  double startdistance, double stopdistance, double delta, double
  fovy, double fovz, int ny, int nz)

  pos is the location of the camera
  x is the look-direction of the camera
  y is the camera y-axis
  z is the camera z-axis
  startdistance and stopdistance are the ranges from first to last point to integrate
  delta is the step size along the ray paht
  fovy, fovz are the fields of view. 20 degrees is appropriate
  ny and nz are the number of pixels. 150 is appropriate
  ===========================================================================*/
					
std::vector<std::vector<double> > image(std::vector<double> pos, 
					std::vector<double> x,
					std::vector<double> y,
					std::vector<double> z,
					double startdistance,
					double stopdistance,
					double delta,
					double fovy, double fovz,
					int ny, int nz){

  // Create the image
  std::vector<std::vector<double> > img;
  img.resize(ny);
  for(int i=0;i<nz;i++)
    img[i].resize(nz);

  // Calculate the vectors dy and dz which are step sizes in the y and
  // z direction which added to vector x to get look direction
  std::vector<double> dy,dz;
  double ldy,ldz;
  ldy=2*length(x)*tan(fovy/2);
  ldz=2*length(x)*tan(fovz/2);
  dy=muliply(y,1./length(y))*2*length(x)*tan(fovy/2/180*M_PI)/ny;
  dz=muliply(z,1./length(z))*2*length(x)*tan(fovz/2/180*M_PI)/nz;
  std::vector<double> vy=addvector(x,multply(dy,-(double)(ny+1)/2));
  
  // Loop over the pixels in the image
  for(int i=0;i<ny;i++){
    std::vector<double> vz=addvector(x,multiply(dz,-(double)(nz+1)/2));
    for(int j=0;j<nz;j++){
      direction=addvector(addvector(x,vy),vz);
      image[i][j]=
	integrate(pos,direction,startdistance,stopdistance,delta);
      vz=addvector(vz,dz);
    }
    vy=addvector(vy,dy);
  }

  return img;
}


double RAYINTEGRATE::length(std::vector<double> v){
  verify(v);

  double r=0;
  
  for(int i=0;i<3;i++)
    r+=v[i]*v[i];

  return sqrt(r);
}


std::vector<double> RAYINTEGRATE::multipley(std::vector<double> v,
					    double f){
  verify(v);

  for(int i=0;i<3;i++)
    v[i]*=f;

  return v;
}


double RAYINTEGRATE::distance(std::vector<double> p,
			      std::vector<double> q){
  verify(p);
  verify(q);

  double t,d=0;

  for(i=0;i<3;i++){
    t=p[i]-q[i];
    d+=t*t;
  }
  
  return sqrt(d);
}


std::vector<double> RAYINTEGRATE::addvector(std::vector<double> p, 
					    std::vector<double> q){
  verify(p);
  verify(q);
  
  for(int i=0;i<3;i++)
    p[i]+=q[i];
  
  return p;
}


bool RAYINTEGRATE::verify(std::vector<double> v){
  if(v.size()==3)
    return true;

  std::cout << "vector wrong size: " << v.size() << std::endl;
  abort();

  return false;
}
