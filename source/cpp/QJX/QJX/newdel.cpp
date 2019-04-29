#include "newdel.h"

void DimerNewDelBehaviour::init_sizes(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp; 
	MainData * md = ad->md;

	cout << "max_num_threads: " << omp_get_max_threads() << endl;
#pragma omp parallel
	{
		int t_id = omp_get_thread_num();
		if (t_id == 0)
		{
			cout << "num_threads_before_init: " << omp_get_num_threads() << endl;
		}
	}
	int num_threads = rp->num_threads;
	cout << "num_threads_target: " << num_threads << endl;
	omp_set_num_threads(num_threads);
#pragma omp parallel
	{
		int t_id = omp_get_thread_num();
		if (t_id == 0)
		{
			cout << "num_threads_after_init: " << omp_get_num_threads() << endl;
		}
	}

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

	cout << "max_num_threads: " << omp_get_max_threads() << endl;
#pragma omp parallel
	{
		int t_id = omp_get_thread_num();
		if (t_id == 0)
		{
			cout << "num_threads_before_init: " << omp_get_num_threads() << endl;
		}
	}
	int num_threads = rp->num_threads;
	cout << "num_threads_target: " << num_threads << endl;
	omp_set_num_threads(num_threads);
#pragma omp parallel
	{
		int t_id = omp_get_thread_num();
		if (t_id == 0)
		{
			cout << "num_threads_after_init: " << omp_get_num_threads() << endl;
		}
	}

	int N = int(cp->params.find("N")->second);
	double T = double(cp->params.find("jcs_prm_alpha")->second);

	md->sys_size = N;
	md->num_diss = 1;
	md->num_ham_qj = 2;
	md->T = T;
}

void PSNewDelBehaviour::init_sizes(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	cout << "max_num_threads: " << omp_get_max_threads() << endl;
#pragma omp parallel
	{
		int t_id = omp_get_thread_num();
		if (t_id == 0)
		{
			cout << "num_threads_before_init: " << omp_get_num_threads() << endl;
		}
	}
	int num_threads = rp->num_threads;
	cout << "num_threads_target: " << num_threads << endl;
	omp_set_num_threads(num_threads);
#pragma omp parallel
	{
		int t_id = omp_get_thread_num();
		if (t_id == 0)
		{
			cout << "num_threads_after_init: " << omp_get_num_threads() << endl;
		}
	}

	int num_spins = int(cp->params.find("ps_num_spins")->second);
	int ps_num_photons_states = int(cp->params.find("ps_num_photons_states")->second);
	// real number of photons is ps_num_photons = ps_num_photons_states - 1

	double T = double(cp->params.find("ps_prm_alpha")->second);

	md->sys_size = ps_num_photons_states * std::round(std::pow(2, num_spins));

	int diss_type = int(cp->params.find("diss_type")->second);
	if (diss_type == 0)
	{
		md->num_diss = 1;
	}
	else if(diss_type == 1)
	{
		md->num_diss = 1 + num_spins;
	}
	
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
			mult_tmp_3[index] = 0.0;

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

void PSNewDelBehaviour::init_hamiltonians(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	double alpha = double(cp->params.find("ps_prm_alpha")->second);
	double g = double(cp->params.find("ps_prm_g")->second);
	double d = double(cp->params.find("ps_prm_d")->second);

	int sys_size = md->sys_size;

	int num_spins = int(cp->params.find("ps_num_spins")->second);
	int s_num_states = std::round(std::pow(2, num_spins));
	int p_num_states = int(cp->params.find("ps_num_photons_states")->second);

	Eigen::MatrixXd s_identity = Eigen::MatrixXd::Identity(s_num_states, s_num_states);
	Eigen::MatrixXd p_identity = Eigen::MatrixXd::Identity(p_num_states, p_num_states);

	Eigen::MatrixXd p_a_std = Eigen::MatrixXd::Zero(p_num_states, p_num_states);
	Eigen::MatrixXd p_a_dag = Eigen::MatrixXd::Zero(p_num_states, p_num_states);
	for (int st_id = 0; st_id < p_num_states - 1; st_id++)
	{
		p_a_std(st_id, st_id + 1) = sqrt(double(st_id + 1));
		p_a_dag(st_id + 1, st_id) = sqrt(double(st_id + 1));
	}

	Eigen::MatrixXd p_hamiltonian = 0.5 / pow(alpha, 3) * p_a_dag * p_a_dag * p_a_std * p_a_std;
	Eigen::MatrixXd p_hamiltonian_drv = p_a_dag - p_a_std;
	Eigen::MatrixXd p_special = p_a_std;

	//if (rp->is_debug)
	//{
	//	string file_name = "p_hamiltonian.txt";
	//	const Eigen::IOFormat common_fmt(Eigen::FullPrecision, 0, "", "\n", "", "", "", "");
	//	std::ofstream f1(file_name.c_str());
	//	f1 << p_hamiltonian.format(common_fmt);

	//	file_name = "p_hamiltonian_drv.txt";
	//	std::ofstream f2(file_name.c_str());
	//	f2 << p_hamiltonian_drv.format(common_fmt);
	//}

	Eigen::Matrix2d sigma_z;
	sigma_z(0, 0) = 1.0;
	sigma_z(0, 1) = 0.0;
	sigma_z(1, 0) = 0.0;
	sigma_z(1, 1) = -1.0;

	Eigen::Matrix2d sigma_minus;
	sigma_minus(0, 0) = 0.0;
	sigma_minus(0, 1) = 0.0;
	sigma_minus(1, 0) = 1.0;
	sigma_minus(1, 1) = 0.0;

	Eigen::Matrix2d sigma_plus;
	sigma_plus(0, 0) = 0.0;
	sigma_plus(0, 1) = 1.0;
	sigma_plus(1, 0) = 0.0;
	sigma_plus(1, 1) = 0.0;

	Eigen::Matrix2d sigma_0;
	sigma_0(0, 0) = 1.0;
	sigma_0(0, 1) = 0.0;
	sigma_0(1, 0) = 0.0;
	sigma_0(1, 1) = 1.0;

	std::vector<Eigen::MatrixXd> s_sigmas_z;
	std::vector<Eigen::MatrixXd> s_sigmas_minus;
	std::vector<Eigen::MatrixXd> s_sigmas_plus;
	for (int s_id = 0; s_id < num_spins; s_id++)
	{
		Eigen::MatrixXd s_sigma_z(2, 2);
		Eigen::MatrixXd s_sigma_plus(2, 2);
		Eigen::MatrixXd s_sigma_minus(2, 2);
		if (s_id == 0)
		{
			s_sigma_z(0, 0) = 1.0;
			s_sigma_z(0, 1) = 0.0;
			s_sigma_z(1, 0) = 0.0;
			s_sigma_z(1, 1) = -1.0;

			sigma_minus(0, 0) = 0.0;
			sigma_minus(0, 1) = 1.0;
			sigma_minus(1, 0) = 0.0;
			sigma_minus(1, 1) = 0.0;

			sigma_plus(0, 0) = 0.0;
			sigma_plus(0, 1) = 0.0;
			sigma_plus(1, 0) = 1.0;
			sigma_plus(1, 1) = 0.0;
		}
		else
		{
			s_sigma_z(0, 0) = 1.0;
			s_sigma_z(0, 1) = 0.0;
			s_sigma_z(1, 0) = 0.0;
			s_sigma_z(1, 1) = 1.0;

			s_sigma_minus(0, 0) = 1.0;
			s_sigma_minus(0, 1) = 0.0;
			s_sigma_minus(1, 0) = 0.0;
			s_sigma_minus(1, 1) = 1.0;

			s_sigma_plus(0, 0) = 1.0;
			s_sigma_plus(0, 1) = 0.0;
			s_sigma_plus(1, 0) = 0.0;
			s_sigma_plus(1, 1) = 1.0;
		}

		for (int it_id = 1; it_id < num_spins; it_id++)
		{
			if (it_id == s_id)
			{
				s_sigma_z = Eigen::kroneckerProduct(s_sigma_z, sigma_z).eval();
				s_sigma_minus = Eigen::kroneckerProduct(s_sigma_minus, sigma_minus).eval();
				s_sigma_plus = Eigen::kroneckerProduct(s_sigma_plus, sigma_plus).eval();
			}
			else
			{
				s_sigma_z = Eigen::kroneckerProduct(s_sigma_z, sigma_0).eval();
				s_sigma_minus = Eigen::kroneckerProduct(s_sigma_minus, sigma_0).eval();
				s_sigma_plus = Eigen::kroneckerProduct(s_sigma_plus, sigma_0).eval();
			}
		}

		s_sigmas_z.push_back(s_sigma_z);
		s_sigmas_minus.push_back(s_sigma_minus);
		s_sigmas_plus.push_back(s_sigma_plus);
	}

	Eigen::MatrixXd s_J_z = s_sigmas_z[0];
	Eigen::MatrixXd s_J_minus = s_sigmas_minus[0];
	Eigen::MatrixXd s_J_plus = s_sigmas_plus[0];
	for (int s_id = 1; s_id < num_spins; s_id++)
	{
		s_J_z += s_sigmas_z[s_id];
		s_J_minus += s_sigmas_minus[s_id];
		s_J_plus += s_sigmas_plus[s_id];
	}
	Eigen::MatrixXd s_hamiltonian = d * 0.5 * s_J_z;

	Eigen::MatrixXd p_hamiltonian_kron = Eigen::kroneckerProduct(s_identity, p_hamiltonian);

	Eigen::MatrixXd s_hamiltonian_kron = Eigen::kroneckerProduct(s_hamiltonian, p_identity);

	Eigen::MatrixXd p_hamiltonian_drv_kron = Eigen::kroneckerProduct(s_identity, p_hamiltonian_drv);

	Eigen::MatrixXd a_dag_J_minus = Eigen::kroneckerProduct(s_identity, p_a_dag) * Eigen::kroneckerProduct(s_J_minus, p_identity);
	Eigen::MatrixXd J_plus_a_std = Eigen::kroneckerProduct(s_J_plus, p_identity) * Eigen::kroneckerProduct(s_identity, p_a_dag);
	Eigen::MatrixXd s_p_hamiltonian = g * 0.5 * (a_dag_J_minus - J_plus_a_std);

	Eigen::MatrixXd hamiltonian = s_hamiltonian_kron + s_p_hamiltonian + p_hamiltonian_kron;
	Eigen::MatrixXd hamiltonian_drv = p_hamiltonian_drv_kron;
	Eigen::MatrixXd special = Eigen::kroneckerProduct(s_identity, p_special);

	md->hamiltonian = new double[md->sys_size * md->sys_size];
	md->hamiltonian_drv = new double[md->sys_size * md->sys_size];
	md->special = new MKL_Complex16[md->sys_size * md->sys_size];

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			md->hamiltonian[index] = hamiltonian(st_id_1, st_id_2);

			md->hamiltonian_drv[index] = hamiltonian_drv(st_id_1, st_id_2);

			md->special[index].real = special(st_id_1, st_id_2);
			md->special[index].imag = 0.0;
		}
	}
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

void PSNewDelBehaviour::init_dissipators(AllData * ad) const
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

	int num_spins = int(cp->params.find("ps_num_spins")->second);
	int s_num_states = std::round(std::pow(2, num_spins));
	int p_num_states = int(cp->params.find("ps_num_photons_states")->second);

	Eigen::MatrixXd s_identity = Eigen::MatrixXd::Identity(s_num_states, s_num_states);
	Eigen::MatrixXd p_identity = Eigen::MatrixXd::Identity(p_num_states, p_num_states);

	Eigen::MatrixXd p_a_std = Eigen::MatrixXd::Zero(p_num_states, p_num_states);
	Eigen::MatrixXd p_a_dag = Eigen::MatrixXd::Zero(p_num_states, p_num_states);
	for (int st_id = 0; st_id < p_num_states - 1; st_id++)
	{
		p_a_std(st_id, st_id + 1) = sqrt(double(st_id + 1));
		p_a_dag(st_id + 1, st_id) = sqrt(double(st_id + 1));
	}

	Eigen::MatrixXd dissipator_cav = Eigen::kroneckerProduct(s_identity, p_a_std);

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			md->dissipators[0][index].real = dissipator_cav(st_id_1, st_id_2);
			md->dissipators[0][index].imag = 0.0;
		}
	}

	int diss_type = int(cp->params.find("diss_type")->second);

	if (diss_type == 1)
	{
		
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

void PSNewDelBehaviour::init_hamiltonians_qj(AllData * ad) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	double ps_drv_ampl = double(cp->params.find("ps_drv_ampl")->second);
	double ps_prm_alpha = double(cp->params.find("ps_prm_alpha")->second);

	double diss_gamma_cav = 1.0 / (4.0 * ps_prm_alpha);

	MKL_Complex16 diss_gamma_cav_cmplx = {diss_gamma_cav, 0.0};
	MKL_Complex16 zero_cmplx = {0.0, 0.0};

	md->hamiltonians_qj = new MKL_Complex16*[md->num_ham_qj];
	for (int qj_ham_id = 0; qj_ham_id < md->num_ham_qj; qj_ham_id++)
	{
		md->hamiltonians_qj[qj_ham_id] = new MKL_Complex16[md->sys_size * md->sys_size];
	}

	MKL_Complex16 * diss_part_cav = new MKL_Complex16[md->sys_size * md->sys_size];
	MKL_Complex16 * hamitlonian_part = new MKL_Complex16[md->sys_size * md->sys_size];
	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			diss_part_cav[index].real = 0.0;
			diss_part_cav[index].imag = 0.0;

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
		&diss_gamma_cav_cmplx,
		md->dissipators[0],
		md->sys_size,
		md->dissipators[0],
		md->sys_size,
		&zero_cmplx,
		diss_part_cav,
		md->sys_size
	);

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;
			hamitlonian_part[index].real += 0.5 * diss_part_cav[index].imag;
			hamitlonian_part[index].imag -= 0.5 * diss_part_cav[index].real;
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

	double E_0 = ps_drv_ampl;
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

	delete[] diss_part_cav;
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

void PSNewDelBehaviour::free_hamiltonians(AllData * ad) const
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

void PSNewDelBehaviour::free_dissipators(AllData * ad) const
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

void PSNewDelBehaviour::free_hamiltonians_qj(AllData * ad) const
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
