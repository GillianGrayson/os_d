#ifndef __CONFIG__
#define __CONFIG__

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif


struct ConfigParam
{
	int    N; // number of particles, system size = N+1
	double J; // hopping constant
	double E0; // on-site potential
	double U; // on-site interaction
	double g; // gamma;
	double A0;
	double w;
	double T;
	double h;
	int hasDriving;
	int driving_type;
	int avg_purity;
	int multiplicators;
	int deep_dump;
	int NSTEP;
	int N_T;
	int CalcEig;
	ConfigParam(
		int    _N = 10,
		double _J = 1,
		double _E0 = 0,
		double _U = 1,
		double _g = 1,
		double _A0 = 0.0,
		double _w = 1.0,
		int _NSTEP = 1000,
		int _N_T = 100,
		int _hasDriving = 0,
		int _driving_type = 0,
		int _avg_purity = 0,
		int _multiplicators = 0,
		int _deep_dump = 0,
		int _CalcEig = 1
		)
	{
		N = _N; // number of particles, system size = N+1
		J = _J; // hopping constant
		E0 = _E0; // on-site potential
		U = _U; // on-site interaction
		g = _g / _N; // gamma;
		A0 = _A0;
		w = _w;
		T = 2.0 * PI / w;
		NSTEP = _NSTEP;
		h = T / NSTEP;
		hasDriving = _hasDriving;
		driving_type = _driving_type;
		avg_purity = _avg_purity;
		multiplicators = _multiplicators;
		deep_dump = _deep_dump;
		N_T = _N_T;
		CalcEig = _CalcEig;
	}
};

#endif