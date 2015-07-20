/*********************************************************************
 * This simply does a integration along a ray from a starting        *
 * location to a finite distance.                                    *
 *********************************************************************/

class RAYINTEGRATEMODEL{
public:
  RAYINTEGRATEMODEL();
  ~RAYINTEGRATEMODEL();
  virtual double get(std::vector<double>)=0;
private:

};


class RAYINTEGRATE{
public:
  RAYINTEGRATE(RAYINTEGRATEMODEL *m);
  ~RAYINTEGRATE();
  double integrate(std::vector<double> start, 
		   std::vector<double> direction,
		   double startdistance, double stopdistance,
		   double delta);
  std::vector<std::vector<double> > image(std::vector<double> pos, 
					  std::vector<double> x,
					  std::vector<double> y,
					  std::vector<double> z,
					  double startdistance,
					  double stopdistance,
					  double delta,
					  double fovy, double fovz,
					  int ny, int nz);
private:
  RAYINTEGRATEMODEL *m;
  
  double length(std::vector<double> v);
  std::vector<double> multiply(std::vector<double> v, double f);
  double distance(std::vector<double> p, std::vector<double> q);
  std::vector<double> addvector(std::vector<double> p, std::vector<double> q);
  bool verify(std::vector<double> v);
};