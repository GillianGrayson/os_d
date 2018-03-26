#include "newdel.h"

void DimerNewDelBehaviour::init_sizes(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp; 
	MainData * md = ad->md;

	int num_threads = rp->num_threads;
	omp_set_num_threads(num_threads);

	int N = int(cp->params.find("N")->second);
	double T = 2.0 * PI / double(cp->params.find("dimer_drv_freq")->second);
	
	md->sys_size = N + 1;
	md->num_diss = 1;
	md->num_ham_qj = 2;
	md->T = T;
}

void JCSNewDelBehaviour::init_sizes(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int num_threads = rp->num_threads;
	omp_set_num_threads(num_threads);

	int N = int(cp->params.find("N")->second);
	double T = double(cp->params.find("jcs_prm_alpha")->second);

	md->sys_size = N;
	md->num_diss = 1;
	md->num_ham_qj = 2;
	md->T = T;
}

void DimerNewDelBehaviour::init_hamiltonians(AllData * ad) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int N = int(cp->params.find("N")->second);
	int sys_size = md->sys_size;

	double dimer_prm_E = double(cp->params.find("dimer_prm_E")->second);
	double dimer_prm_U = double(cp->params.find("dimer_prm_U")->second);
	dimer_prm_U = dimer_prm_U / double(md->sys_size - 1);
	double dimer_prm_J = double(cp->params.find("dimer_prm_J")->second);

	md->hamiltonian = new double[md->sys_size * md->sys_size];
	md->hamiltonian_drv = new double[md->sys_size * md->sys_size];

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			md->hamiltonian[index] = 0.0;
			md->hamiltonian_drv[index] = 0.0;
		}
	}

	for (int st_id = 0; st_id < md->sys_size; st_id++)
	{
		int index = st_id * md->sys_size + st_id;

		md->hamiltonian[index] = 2.0 * dimer_prm_U * double(st_id * (st_id - 1) + (sys_size - (st_id + 1)) * (sys_size - (st_id + 1) - 1));
		md->hamiltonian_drv[index] = double((sys_size - (st_id + 1)) - st_id);
	}

	for (int st_id = 0; st_id < (md->sys_size - 1); st_id++)
	{
		int down_left = (st_id + 1) * md->sys_size + st_id;
		int up_right = st_id * md->sys_size + (st_id + 1);

		md->hamiltonian[down_left] -= dimer_prm_J * sqrt(double((md->sys_size - (st_id + 1)) * (st_id + 1)));
		md->hamiltonian[up_right] -= dimer_prm_J * sqrt(double((st_id + 1) * (md->sys_size - (st_id + 1))));
	}
}

void JCSNewDelBehaviour::init_hamiltonians(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	double alpha = double(cp->params.find("jcs_prm_alpha")->second);

	int sys_size = md->sys_size;

	double * a_std = new double[md->sys_size * md->sys_size];
	double * a_dag = new double[md->sys_size * md->sys_size];

	double * mult_tmp_1 = new double[md->sys_size * md->sys_size];
	double * mult_tmp_2 = new double[md->sys_size * md->sys_size];
	double * mult_tmp_3 = new double[md->sys_size * md->sys_size];

	md->hamiltonian = new double[md->sys_size * md->sys_size];
	md->hamiltonian_drv = new double[md->sys_size * md->sys_size];
	md->special = new MKL_Complex16[md->sys_size * md->sys_size];

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			a_std[index] = 0.0;
			a_dag[index] = 0.0;

			mult_tmp_1[index] = 0.0;
			mult_tmp_2[index] = 0.0;
			mult_tmp_2[index] = 0.0;

			md->hamiltonian[index] = 0.0;
			md->hamiltonian_drv[index] = 0.0;
		}
	}

	for (int st_id = 0; st_id < sys_size - 1; st_id++)
	{
		int index_std = st_id * md->sys_size + (st_id + 1);
		int index_dag = (st_id + 1) * md->sys_size + st_id;
		a_std[index_std] = sqrt(double(st_id + 1));
		a_dag[index_dag] = sqrt(double(st_id + 1));
	}

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, sys_size, sys_size, sys_size, 1.0, a_dag, sys_size, a_dag, sys_size, 0.0, mult_tmp_1, sys_size);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, sys_size, sys_size, sys_size, 1.0, mult_tmp_1, sys_size, a_std, sys_size, 0.0, mult_tmp_2, sys_size);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, sys_size, sys_size, sys_size, 1.0, mult_tmp_2, sys_size, a_std, sys_size, 0.0, mult_tmp_3, sys_size);

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			md->hamiltonian[index] = 0.5 / pow(alpha, 3) * mult_tmp_3[index];
			md->hamiltonian_drv[index] = (a_dag[index] - a_std[index]);

			md->special[index].real = a_std[index];
			md->special[index].imag = 0.0;
		}
	}

	delete[] a_std;
	delete[] a_dag;

	delete[] mult_tmp_1;
	delete[] mult_tmp_2;
	delete[] mult_tmp_3;
}

void DimerNewDelBehaviour::init_dissipators(AllData * ad) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int N = int(cp->params.find("N")->second);
	int sys_size = md->sys_size;

	md->dissipators = new MKL_Complex16 *[md->num_diss];
	for (int diss_id = 0; diss_id < md->num_diss; diss_id++)
	{
		md->dissipators[diss_id] = new MKL_Complex16[md->sys_size * md->sys_size];
	}

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			for (int diss_id = 0; diss_id < md->num_diss; diss_id++)
			{
				md->dissipators[diss_id][index].real = 0.0;
				md->dissipators[diss_id][index].imag = 0.0;
			}
		}
	}

	for (int st_id = 0; st_id < md->sys_size; st_id++)
	{
		int index = st_id * md->sys_size + st_id;

		md->dissipators[0][index].real = double(st_id - (sys_size - (st_id + 1)));
		md->dissipators[0][index].imag = 0.0;
	}

	for (int st_id = 0; st_id < (md->sys_size - 1); st_id++)
	{
		int down_left = (st_id + 1) * md->sys_size + st_id;
		int up_right = st_id * md->sys_size + (st_id + 1);

		md->dissipators[0][down_left].real -= sqrt(double((md->sys_size - (st_id + 1)) * (st_id + 1)));
		md->dissipators[0][up_right].real += sqrt(double((st_id + 1) * (md->sys_size - (st_id + 1))));
	}
}

void JCSNewDelBehaviour::init_dissipators(AllData * ad) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int sys_size = md->sys_size;

	md->dissipators = new MKL_Complex16 *[md->num_diss];
	for (int diss_id = 0; diss_id < md->num_diss; diss_id++)
	{
		md->dissipators[diss_id] = new MKL_Complex16[md->sys_size * md->sys_size];
	}

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			for (int diss_id = 0; diss_id < md->num_diss; diss_id++)
			{
				md->dissipators[diss_id][index].real = 0.0;
				md->dissipators[diss_id][index].imag = 0.0;
			}
		}
	}

	for (int st_id = 0; st_id < sys_size - 1; st_id++)
	{
		int index= st_id * md->sys_size + (st_id + 1);

		md->dissipators[0][index].real = sqrt(double(st_id + 1));
		md->dissipators[0][index].imag = 0.0;
	}
}

void DimerNewDelBehaviour::init_hamiltonians_qj(AllData * ad) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	double dimer_prm_E = double(cp->params.find("dimer_prm_E")->second);
	double dimer_drv_ampl = double(cp->params.find("dimer_drv_ampl")->second);

	double diss_gamma = double(cp->params.find("diss_gamma")->second);
	diss_gamma = diss_gamma / double(md->sys_size - 1);

	MKL_Complex16 diss_gamma_cmplx = { diss_gamma, 0.0 };
	MKL_Complex16 zero_cmplx = { 0.0, 0.0 };

	md->hamiltonians_qj = new MKL_Complex16*[md->num_ham_qj];
	for (int qj_ham_id = 0; qj_ham_id < md->num_ham_qj; qj_ham_id++)
	{
		md->hamiltonians_qj[qj_ham_id] = new MKL_Complex16[md->sys_size * md->sys_size];
	}

	MKL_Complex16 * diss_part = new MKL_Complex16[md->sys_size * md->sys_size];
	MKL_Complex16 * hamitlonian_part = new MKL_Complex16[md->sys_size * md->sys_size];
	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			diss_part[index].real = 0.0;
			diss_part[index].imag = 0.0;

			hamitlonian_part[index].real = md->hamiltonian[index];
			hamitlonian_part[index].imag = 0.0;

			for (int qj_ham_id = 0; qj_ham_id < md->num_ham_qj; qj_ham_id++)
			{
				md->hamiltonians_qj[qj_ham_id][index].real = 0.0;
				md->hamiltonians_qj[qj_ham_id][index].imag = 0.0;
			}
		}
	}

	cblas_zgemm(
		CblasRowMajor,
		CblasConjTrans,
		CblasNoTrans,
		md->sys_size,
		md->sys_size,
		md->sys_size,
		&diss_gamma_cmplx,
		md->dissipators[0],
		md->sys_size,
		md->dissipators[0],
		md->sys_size,
		&zero_cmplx,
		diss_part,
		md->sys_size
	);

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;
			hamitlonian_part[index].real += 0.5 * diss_part[index].imag;
			hamitlonian_part[index].imag -= 0.5 * diss_part[index].real;
		}
	}

	md->non_drv_part = new MKL_Complex16[md->sys_size * md->sys_size];
	md->drv_part = new MKL_Complex16[md->sys_size * md->sys_size];
	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;
			md->non_drv_part[index].real = hamitlonian_part[index].real;
			md->non_drv_part[index].imag = hamitlonian_part[index].imag;

			md->drv_part[index].real = md->hamiltonian_drv[index];
			md->drv_part[index].imag = 0.0;
		}
	}

	double E_0 = dimer_prm_E + dimer_drv_ampl;
	double E_1 = dimer_prm_E - dimer_drv_ampl;

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			md->hamiltonians_qj[0][index].real += (hamitlonian_part[index].imag + E_0 * 0.0);
			md->hamiltonians_qj[0][index].imag -= (hamitlonian_part[index].real + E_0 * md->hamiltonian_drv[index]);

			md->hamiltonians_qj[1][index].real += (hamitlonian_part[index].imag + E_1 * 0.0);
			md->hamiltonians_qj[1][index].imag -= (hamitlonian_part[index].real + E_1 * md->hamiltonian_drv[index]);
		}
	}

	delete[] diss_part;
	delete[] hamitlonian_part;
}

void JCSNewDelBehaviour::init_hamiltonians_qj(AllData * ad) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	double jcs_drv_ampl = double(cp->params.find("jcs_drv_ampl")->second);
	double jcs_prm_alpha = double(cp->params.find("jcs_prm_alpha")->second);

	double diss_gamma = 1.0 / (4.0 * jcs_prm_alpha);

	MKL_Complex16 diss_gamma_cmplx = { diss_gamma, 0.0 };
	MKL_Complex16 zero_cmplx = { 0.0, 0.0 };

	md->hamiltonians_qj = new MKL_Complex16*[md->num_ham_qj];
	for (int qj_ham_id = 0; qj_ham_id < md->num_ham_qj; qj_ham_id++)
	{
		md->hamiltonians_qj[qj_ham_id] = new MKL_Complex16[md->sys_size * md->sys_size];
	}

	MKL_Complex16 * diss_part = new MKL_Complex16[md->sys_size * md->sys_size];
	MKL_Complex16 * hamitlonian_part = new MKL_Complex16[md->sys_size * md->sys_size];
	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			diss_part[index].real = 0.0;
			diss_part[index].imag = 0.0;

			hamitlonian_part[index].real = md->hamiltonian[index];
			hamitlonian_part[index].imag = 0.0;

			for (int qj_ham_id = 0; qj_ham_id < md->num_ham_qj; qj_ham_id++)
			{
				md->hamiltonians_qj[qj_ham_id][index].real = 0.0;
				md->hamiltonians_qj[qj_ham_id][index].imag = 0.0;
			}
		}
	}

	cblas_zgemm(
		CblasRowMajor,
		CblasConjTrans,
		CblasNoTrans,
		md->sys_size,
		md->sys_size,
		md->sys_size,
		&diss_gamma_cmplx,
		md->dissipators[0],
		md->sys_size,
		md->dissipators[0],
		md->sys_size,
		&zero_cmplx,
		diss_part,
		md->sys_size
	);

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;
			hamitlonian_part[index].real += 0.5 * diss_part[index].imag;
			hamitlonian_part[index].imag -= 0.5 * diss_part[index].real;
		}
	}

	md->non_drv_part = new MKL_Complex16[md->sys_size * md->sys_size];
	md->drv_part = new MKL_Complex16[md->sys_size * md->sys_size];
	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;
			md->non_drv_part[index].real = hamitlonian_part[index].real;
			md->non_drv_part[index].imag = hamitlonian_part[index].imag;

			md->drv_part[index].real = 0.0;
			md->drv_part[index].imag = md->hamiltonian_drv[index];
		}
	}

	double E_0 = jcs_drv_ampl;
	double E_1 = 0.0;

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			md->hamiltonians_qj[0][index].real += (hamitlonian_part[index].imag + E_0 * md->hamiltonian_drv[index]);
			md->hamiltonians_qj[0][index].imag -= (hamitlonian_part[index].real + E_0 * 0.0);

			md->hamiltonians_qj[1][index].real += (hamitlonian_part[index].imag + E_1 * md->hamiltonian_drv[index]);
			md->hamiltonians_qj[1][index].imag -= (hamitlonian_part[index].real + E_1 * 0.0);
		}
	}

	delete[] diss_part;
	delete[] hamitlonian_part;
}

void DimerNewDelBehaviour::free_hamiltonians(AllData * ad) const
{
	MainData * md = ad->md;

	delete[] md->hamiltonian;
	delete[] md->hamiltonian_drv;
}

void JCSNewDelBehaviour::free_hamiltonians(AllData * ad) const
{
	MainData * md = ad->md;

	delete[] md->hamiltonian;
	delete[] md->hamiltonian_drv;
	delete[] md->special;
}

void DimerNewDelBehaviour::free_dissipators(AllData * ad) const
{
	MainData * md = ad->md;

	for (int diss_id = 0; diss_id < md->num_diss; diss_id++)
	{
		delete[] md->dissipators[diss_id];
	}
	delete[] md->dissipators;
}

void JCSNewDelBehaviour::free_dissipators(AllData * ad) const
{
	MainData * md = ad->md;

	for (int diss_id = 0; diss_id < md->num_diss; diss_id++)
	{
		delete[] md->dissipators[diss_id];
	}
	delete[] md->dissipators;
}

void DimerNewDelBehaviour::free_hamiltonians_qj(AllData * ad) const
{
	MainData * md = ad->md;

	delete[] md->non_drv_part;
	delete[] md->drv_part;

	for (int qj_ham_id = 0; qj_ham_id < md->num_ham_qj; qj_ham_id++)
	{
		delete[] md->hamiltonians_qj[qj_ham_id];
	}
	delete[] md->hamiltonians_qj;
}

void JCSNewDelBehaviour::free_hamiltonians_qj(AllData * ad) const
{
	MainData * md = ad->md;

	delete[] md->non_drv_part;
	delete[] md->drv_part;

	for (int qj_ham_id = 0; qj_ham_id < md->num_ham_qj; qj_ham_id++)
	{
		delete[] md->hamiltonians_qj[qj_ham_id];
	}
	delete[] md->hamiltonians_qj;
}
