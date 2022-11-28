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

void DimerSyncNewDelBehaviour::init_sizes(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

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

	int N = int(cp->params.find("dimersync_N")->second);
	double T1 = 2.0 * PI / double(cp->params.find("dimersync_drv_freq_1")->second);
	double T2 = 2.0 * PI / double(cp->params.find("dimersync_drv_freq_2")->second);
	cout << "T1 = " << T1 << endl;
	cout << "T2 = " << T2 << endl;
	double T = T1;
	if (T2 < T1)
	{
		T = T2;
	}
	cout << "T = " << T << endl;

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
	cout << "sys_size: " << md->sys_size << endl;

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

void MBLNewDelBehaviour::init_sizes(AllData * ad) const
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

	md->mbl_Nc = int(cp->params.find("mbl_Nc")->second); // Number of cells
	if (md->mbl_Nc % 2 != 0)
	{
		stringstream msg;
		msg << "Nc must divide by 2 without any remainder" << endl;
		Error(msg.str());
	}

	md->mbl_Np = md->mbl_Nc / 2;
	md->mbl_Ns = n_choose_k(md->mbl_Nc, md->mbl_Np);

	int num_global_states = pow(2, md->mbl_Nc);

	md->mbl_adjacement = new int[num_global_states];
	md->mbl_x_to_id = new int[num_global_states];
	md->mbl_id_to_x = new int[md->mbl_Ns];

	int state_id = 0;
	for (int global_state_id = 0; global_state_id < num_global_states; global_state_id++)
	{
		if ((bit_count(global_state_id) == 2) && (bit_count(global_state_id & (global_state_id << 1)) == 1))
		{
			md->mbl_adjacement[global_state_id] = 1;
		}
		else
		{
			md->mbl_adjacement[global_state_id] = 0;
		}

		if (bit_count(global_state_id) == md->mbl_Np)
		{
			md->mbl_x_to_id[global_state_id] = state_id + 1;
			md->mbl_id_to_x[state_id] = global_state_id;
			state_id++;
		}
		else
		{
			md->mbl_x_to_id[global_state_id] = 0;
		}
	}

	if (md->mbl_Ns != state_id)
	{
		stringstream msg;
		msg << "wrong num states calculation" << endl;
		Error(msg.str());
	}

	md->sys_size = md->mbl_Ns;
	int diss_type = int(cp->params.find("diss_type")->second);
	if (diss_type == 0)
	{
		md->num_diss = md->mbl_Nc;
	}
	else if (diss_type == 1)
	{
		md->num_diss = md->mbl_Nc - 1;
	}
	
	md->num_ham_qj = 1;

	double T = double(cp->params.find("mbl_T")->second);
	md->T = T;
}

void LndHamNewDelBehaviour::init_sizes(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

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

	md->sys_size = int(cp->params.find("lndham_N")->second);
	md->num_diss = md->sys_size * md->sys_size - 1;

	md->num_ham_qj = 1;

	double T = double(cp->params.find("lndham_T")->second);
	md->T = T;
}

void IntegrableNewDelBehaviour::init_sizes(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

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
	int N = int(cp->params.find("integrable_N")->second);
	md->sys_size = std::round(std::pow(2, N));
	md->num_diss = N;

	md->num_ham_qj = 1;

	double T = double(cp->params.find("integrable_T")->second);
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

	md->hamiltonian = new MKL_Complex16[md->sys_size * md->sys_size];
	md->hamiltonian_drv = new MKL_Complex16[md->sys_size * md->sys_size];

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			md->hamiltonian[index].real = 0.0;
			md->hamiltonian[index].imag = 0.0;

			md->hamiltonian_drv[index].real = 0.0;
			md->hamiltonian_drv[index].imag = 0.0;
		}
	}

	for (int st_id = 0; st_id < md->sys_size; st_id++)
	{
		int index = st_id * md->sys_size + st_id;

		md->hamiltonian[index].real = 2.0 * dimer_prm_U * double(st_id * (st_id - 1) + (sys_size - (st_id + 1)) * (sys_size - (st_id + 1) - 1));
		md->hamiltonian_drv[index].real = double((sys_size - (st_id + 1)) - st_id);
	}

	for (int st_id = 0; st_id < (md->sys_size - 1); st_id++)
	{
		int down_left = (st_id + 1) * md->sys_size + st_id;
		int up_right = st_id * md->sys_size + (st_id + 1);

		md->hamiltonian[down_left].real -= dimer_prm_J * sqrt(double((md->sys_size - (st_id + 1)) * (st_id + 1)));
		md->hamiltonian[up_right].real -= dimer_prm_J * sqrt(double((st_id + 1) * (md->sys_size - (st_id + 1))));
	}

	init_random_obs(ad);
}

void DimerSyncNewDelBehaviour::init_hamiltonians(AllData* ad) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	int N = int(cp->params.find("dimersync_N")->second);
	int sys_size = md->sys_size;

	double prm_E = double(cp->params.find("dimersync_prm_E")->second);
	double prm_U = double(cp->params.find("dimersync_prm_U")->second);
	prm_U = prm_U / double(md->sys_size - 1);
	double prm_J = double(cp->params.find("dimersync_prm_J")->second);

	md->hamiltonian = new MKL_Complex16[md->sys_size * md->sys_size];
	md->hamiltonian_drv = new MKL_Complex16[md->sys_size * md->sys_size];
	md->hamiltonian_drv_2 = new MKL_Complex16[md->sys_size * md->sys_size];

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			md->hamiltonian[index].real = 0.0;
			md->hamiltonian[index].imag = 0.0;

			md->hamiltonian_drv[index].real = 0.0;
			md->hamiltonian_drv[index].imag = 0.0;


			md->hamiltonian_drv_2[index].real = 0.0;
			md->hamiltonian_drv_2[index].imag = 0.0;
		}
	}

	for (int st_id = 0; st_id < md->sys_size; st_id++)
	{
		int index = st_id * md->sys_size + st_id;

		md->hamiltonian[index].real = 2.0 * prm_U * double(st_id * (st_id - 1) + (sys_size - (st_id + 1)) * (sys_size - (st_id + 1) - 1));
		md->hamiltonian_drv[index].real = double((sys_size - (st_id + 1)) - st_id);
	}

	for (int st_id = 0; st_id < (md->sys_size - 1); st_id++)
	{
		int down_left = (st_id + 1) * md->sys_size + st_id;
		int up_right = st_id * md->sys_size + (st_id + 1);

		md->hamiltonian[down_left].real -= prm_J * sqrt(double((md->sys_size - (st_id + 1)) * (st_id + 1)));
		md->hamiltonian[up_right].real -= prm_J * sqrt(double((st_id + 1) * (md->sys_size - (st_id + 1))));
	}

	for (int st_id = 0; st_id < (md->sys_size - 1); st_id++)
	{
		int down_left = (st_id + 1) * md->sys_size + st_id;
		int up_right = st_id * md->sys_size + (st_id + 1);

		md->hamiltonian_drv_2[down_left].real -= sqrt(double((md->sys_size - (st_id + 1)) * (st_id + 1)));
		md->hamiltonian_drv_2[up_right].real -= sqrt(double((st_id + 1) * (md->sys_size - (st_id + 1))));
	}

	init_random_obs(ad);
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

	md->hamiltonian = new MKL_Complex16[md->sys_size * md->sys_size];
	md->hamiltonian_drv = new MKL_Complex16[md->sys_size * md->sys_size];
	md->special = new MKL_Complex16[md->sys_size * md->sys_size];
	md->special_2 = new MKL_Complex16[md->sys_size * md->sys_size];
	md->special_3 = new MKL_Complex16[md->sys_size * md->sys_size];

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

			md->hamiltonian[index].real = 0.0;
			md->hamiltonian[index].imag = 0.0;

			md->hamiltonian_drv[index].real = 0.0;
			md->hamiltonian_drv[index].imag = 0.0;
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

			md->hamiltonian[index].real = 0.5 / pow(alpha, 3) * mult_tmp_3[index];
			md->hamiltonian_drv[index].real = (a_dag[index] - a_std[index]);

			md->special[index].real = a_std[index];
			md->special[index].imag = 0.0;

			md->special_2[index].real = a_std[index];
			md->special_2[index].imag = 0.0;

			md->special_3[index].real = a_std[index];
			md->special_3[index].imag = 0.0;
		}
	}

	delete[] a_std;
	delete[] a_dag;

	delete[] mult_tmp_1;
	delete[] mult_tmp_2;
	delete[] mult_tmp_3;

	init_random_obs(ad);
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
	Eigen::MatrixXd p_special_4 = p_a_dag * p_a_std;

	if (rp->is_debug)
	{
		string file_name = "p_hamiltonian.txt";
		const Eigen::IOFormat common_fmt(Eigen::FullPrecision, 0, "", "\n", "", "", "", "");
		std::ofstream f1(file_name.c_str());
		f1 << p_hamiltonian.format(common_fmt);

		file_name = "p_hamiltonian_drv.txt";
		std::ofstream f2(file_name.c_str());
		f2 << p_hamiltonian_drv.format(common_fmt);
	}

	Eigen::Matrix2d sigma_z(2, 2);
	sigma_z(0, 0) = 1.0;
	sigma_z(0, 1) = 0.0;
	sigma_z(1, 0) = 0.0;
	sigma_z(1, 1) = -1.0;

	Eigen::Matrix2d sigma_minus(2, 2);
	sigma_minus(0, 0) = 0.0;
	sigma_minus(0, 1) = 0.0;
	sigma_minus(1, 0) = 1.0;
	sigma_minus(1, 1) = 0.0;

	Eigen::Matrix2d sigma_plus(2, 2);
	sigma_plus(0, 0) = 0.0;
	sigma_plus(0, 1) = 1.0;
	sigma_plus(1, 0) = 0.0;
	sigma_plus(1, 1) = 0.0;

	Eigen::Matrix2d sigma_0(2, 2);
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
		Eigen::MatrixXd s_sigma_minus(2, 2);
		Eigen::MatrixXd s_sigma_plus(2, 2);
		if (s_id == 0)
		{
			s_sigma_z(0, 0) = 1.0;
			s_sigma_z(0, 1) = 0.0;
			s_sigma_z(1, 0) = 0.0;
			s_sigma_z(1, 1) = -1.0;

			s_sigma_minus(0, 0) = 0.0;
			s_sigma_minus(0, 1) = 0.0;
			s_sigma_minus(1, 0) = 1.0;
			s_sigma_minus(1, 1) = 0.0;

			s_sigma_plus(0, 0) = 0.0;
			s_sigma_plus(0, 1) = 1.0;
			s_sigma_plus(1, 0) = 0.0;
			s_sigma_plus(1, 1) = 0.0;
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

	s_J_z *= 0.5;
	s_J_minus *= 0.5;
	s_J_plus *= 0.5;

	if (rp->is_debug)
	{
		string file_name = "s_J_z.txt";
		const Eigen::IOFormat common_fmt(Eigen::FullPrecision, 0, "", "\n", "", "", "", "");
		std::ofstream f1(file_name.c_str());
		f1 << s_J_z.format(common_fmt);

		file_name = "s_J_minus.txt";
		std::ofstream f2(file_name.c_str());
		f2 << s_J_minus.format(common_fmt);

		file_name = "s_J_plus.txt";
		std::ofstream f3(file_name.c_str());
		f3 << s_J_plus.format(common_fmt);
	}

	Eigen::MatrixXd s_hamiltonian = d * 0.5 * s_J_z;

	Eigen::MatrixXd p_hamiltonian_kron = Eigen::kroneckerProduct(s_identity, p_hamiltonian);

	Eigen::MatrixXd s_hamiltonian_kron = Eigen::kroneckerProduct(s_hamiltonian, p_identity);

	Eigen::MatrixXd p_hamiltonian_drv_kron = Eigen::kroneckerProduct(s_identity, p_hamiltonian_drv);

	Eigen::MatrixXd a_dag_J_minus = Eigen::kroneckerProduct(s_identity, p_a_dag) * Eigen::kroneckerProduct(s_J_minus, p_identity);
	Eigen::MatrixXd J_plus_a_std = Eigen::kroneckerProduct(s_J_plus, p_identity) * Eigen::kroneckerProduct(s_identity, p_a_std);
	Eigen::MatrixXd s_p_hamiltonian = g * 0.5 * (a_dag_J_minus + J_plus_a_std);

	Eigen::MatrixXd hamiltonian = s_hamiltonian_kron + s_p_hamiltonian + p_hamiltonian_kron;
	Eigen::MatrixXd hamiltonian_drv = p_hamiltonian_drv_kron;
	Eigen::MatrixXd special = Eigen::kroneckerProduct(s_identity, p_special);

	Eigen::MatrixXd special_2 = Eigen::kroneckerProduct(s_J_z, p_identity);
	Eigen::MatrixXd special_3 = Eigen::kroneckerProduct(s_J_plus, p_identity);
	Eigen::MatrixXd special_4 = Eigen::kroneckerProduct(s_identity, p_special_4);

	md->hamiltonian = new MKL_Complex16[md->sys_size * md->sys_size];
	md->hamiltonian_drv = new MKL_Complex16[md->sys_size * md->sys_size];
	md->special = new MKL_Complex16[md->sys_size * md->sys_size];
	md->special_2 = new MKL_Complex16[md->sys_size * md->sys_size];
	md->special_3 = new MKL_Complex16[md->sys_size * md->sys_size];
	md->special_4 = new MKL_Complex16[md->sys_size * md->sys_size];

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			md->hamiltonian[index].real = hamiltonian(st_id_1, st_id_2);
			md->hamiltonian[index].imag = 0.0;

			md->hamiltonian_drv[index].real = hamiltonian_drv(st_id_1, st_id_2);
			md->hamiltonian_drv[index].imag = 0.0;

			md->special[index].real = special(st_id_1, st_id_2);
			md->special[index].imag = 0.0;

			md->special_2[index].real = special_2(st_id_1, st_id_2);
			md->special_2[index].imag = 0.0;

			md->special_3[index].real = special_3(st_id_1, st_id_2);
			md->special_3[index].imag = 0.0;

			md->special_4[index].real = special_4(st_id_1, st_id_2);
			md->special_4[index].imag = 0.0;
		}
	}

	init_random_obs(ad);
}

void MBLNewDelBehaviour::init_hamiltonians(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	md->hamiltonian = new MKL_Complex16[md->sys_size * md->sys_size];
	md->special = new MKL_Complex16[md->sys_size * md->sys_size];

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			md->hamiltonian[index].real = 0.0;
			md->hamiltonian[index].imag = 0.0;

			md->special[index].real = 0.0;
			md->special[index].imag = 0.0;
		}
	}

	// ========= disorder data ==========
	double mbl_prm_W = double(cp->params.find("mbl_prm_W")->second);
	int mbl_seed = double(cp->params.find("mbl_seed")->second);
	int mbl_mns = double(cp->params.find("mbl_mns")->second);

	double* energies = new double[md->mbl_Nc];

	VSLStreamStatePtr stream;
	vslNewStream(&stream, VSL_BRNG_MCG31, 77778888);
	vslLeapfrogStream(stream, mbl_seed, mbl_mns);
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, md->mbl_Nc, energies, -1.0, 1.0);

	string energies_fn = "energies" + cp->fn_suffix;
	cout << "save energies to file:" << endl << energies_fn << endl << endl;
	save_double_data(energies_fn, energies, md->mbl_Nc, 16, false);

	for (int state_id = 0; state_id < md->mbl_Ns; state_id++)
	{
		vector<int> vb = convert_int_to_vector_of_bits(md->mbl_id_to_x[state_id], md->mbl_Nc);
		double sum = 0.0;
		for (int cell_id = 0; cell_id < md->mbl_Nc; cell_id++)
		{
			sum += double(vb[cell_id]) * energies[cell_id];
		}
		sum *= 1.0 * mbl_prm_W;

		md->hamiltonian[state_id * md->mbl_Ns + state_id].real += sum;
	}

	delete[] energies;
	// =================================

	// ============= interaction =============
	double mbl_prm_U = double(cp->params.find("mbl_prm_U")->second);
	for (int state_id = 0; state_id < md->mbl_Ns; state_id++)
	{
		md->hamiltonian[state_id * md->mbl_Ns + state_id].real += mbl_prm_U * bit_count(md->mbl_id_to_x[state_id] & (md->mbl_id_to_x[state_id] << 1));
	}
	// =======================================

	// ============ hopping ===============
	double mbl_prm_J = double(cp->params.find("mbl_prm_J")->second);
	for (int state_id_1 = 0; state_id_1 < md->mbl_Ns; state_id_1++)
	{
		for (int state_id_2 = 0; state_id_2 < md->mbl_Ns; state_id_2++)
		{
			md->hamiltonian[state_id_1 * md->mbl_Ns + state_id_2].real -= mbl_prm_J * md->mbl_adjacement[md->mbl_id_to_x[state_id_1] ^ md->mbl_id_to_x[state_id_2]];
		}
	}
	// ====================================

	// ========== init special matrix ==========
	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			md->special[index].real = md->hamiltonian[index].real;
			md->special[index].imag = 0.0;
		}
	}
	// =========================================

	init_random_obs(ad);
}

void LndHamNewDelBehaviour::init_hamiltonians(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	md->hamiltonian = new MKL_Complex16[md->sys_size * md->sys_size];
	md->special = new MKL_Complex16[md->sys_size * md->sys_size];

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			md->hamiltonian[index].real = 0.0;
			md->hamiltonian[index].imag = 0.0;

			md->special[index].real = 0.0;
			md->special[index].imag = 0.0;
		}
	}

	// ========= disorder data ==========
	int lndham_seed = int(cp->params.find("lndham_seed")->second);
	int lndham_num_seeds = int(cp->params.find("lndham_num_seeds")->second);
	double lndham_alpha = double(cp->params.find("lndham_alpha")->second);

	VSLStreamStatePtr stream;
	vslNewStream(&stream, VSL_BRNG_MCG31, 13371337);
	vslLeapfrogStream(stream, lndham_seed, lndham_num_seeds);
	double* disorder = new double[md->sys_size * md->sys_size];
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, md->sys_size * md->sys_size, disorder, 0.0, 1.0);
	Eigen::MatrixXcd X(md->sys_size, md->sys_size);
	for (auto st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (auto st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			auto index = st_id_1 * md->sys_size + st_id_2;
			X(st_id_1, st_id_2) = std::complex<double>(disorder[index], disorder[index]);
		}
	}
	delete[] disorder;

	Eigen::MatrixXcd y = 0.5 * (X + X.adjoint());
	Eigen::MatrixXcd y2 = y * y;
	const auto y2_tr = y2.trace();
	y = y / std::sqrt(y2_tr.real());


	// ====================================

	// ========== init special matrix ==========
	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			complex<double> val = lndham_alpha / std::sqrt(static_cast<double>(md->sys_size)) * y(st_id_1, st_id_2);
			md->hamiltonian[index].real = val.real();
			md->hamiltonian[index].imag = val.imag();

			md->special[index].real = val.real();
			md->special[index].imag = val.imag();
		}
	}
	// =========================================

	init_random_obs(ad);
}

void IntegrableNewDelBehaviour::init_hamiltonians(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	md->hamiltonian = new MKL_Complex16[md->sys_size * md->sys_size];
	md->special = new MKL_Complex16[md->sys_size * md->sys_size];

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			md->hamiltonian[index].real = 0.0;
			md->hamiltonian[index].imag = 0.0;

			md->special[index].real = 0.0;
			md->special[index].imag = 0.0;
		}
	}

	// ========= disorder data ==========
	int integrable_seed = int(cp->params.find("integrable_seed")->second);
	int integrable_num_seeds = int(cp->params.find("integrable_num_seeds")->second);
	double integrable_tau = double(cp->params.find("integrable_tau")->second);
	double integrable_k = double(cp->params.find("integrable_k")->second);

	VSLStreamStatePtr stream;
	vslNewStream(&stream, VSL_BRNG_MCG31, 13371337);
	vslLeapfrogStream(stream, integrable_seed, integrable_num_seeds);
	double* disorder = new double[md->sys_size * md->sys_size];
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, md->sys_size * md->sys_size, disorder, 0.0, 1.0);
	Eigen::MatrixXcd X(md->sys_size, md->sys_size);
	for (auto st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (auto st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			auto index = st_id_1 * md->sys_size + st_id_2;
			X(st_id_1, st_id_2) = std::complex<double>(disorder[index], disorder[index]);
		}
	}
	delete[] disorder;

	Eigen::MatrixXcd y = 0.5 * (X + X.adjoint());
	Eigen::MatrixXcd y2 = y * y;
	const auto y2_tr = y2.trace();
	y = y / std::sqrt(y2_tr.real());

	// ========== init special matrix ==========
	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;
			complex<double> val = y(st_id_1, st_id_2) / std::sqrt(static_cast<double>(md->sys_size));
			md->special[index].real = val.real();
			md->special[index].imag = val.imag();
		}
	}
	// =========================================

	init_random_obs(ad);
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

void DimerSyncNewDelBehaviour::init_dissipators(AllData* ad) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	int N = int(cp->params.find("dimersync_N")->second);
	int sys_size = md->sys_size;

	md->dissipators = new MKL_Complex16 * [md->num_diss];
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
		Eigen::Matrix2d sigma_minus(2, 2);
		sigma_minus(0, 0) = 0.0;
		sigma_minus(0, 1) = 0.0;
		sigma_minus(1, 0) = 1.0;
		sigma_minus(1, 1) = 0.0;

		Eigen::Matrix2d sigma_0(2, 2);
		sigma_0(0, 0) = 1.0;
		sigma_0(0, 1) = 0.0;
		sigma_0(1, 0) = 0.0;
		sigma_0(1, 1) = 1.0;

		for (int s_id = 0; s_id < num_spins; s_id++)
		{
			Eigen::MatrixXd s_sigma_minus(2, 2);
			if (s_id == 0)
			{
				s_sigma_minus(0, 0) = 0.0;
				s_sigma_minus(0, 1) = 0.0;
				s_sigma_minus(1, 0) = 1.0;
				s_sigma_minus(1, 1) = 0.0;
			}
			else
			{
				s_sigma_minus(0, 0) = 1.0;
				s_sigma_minus(0, 1) = 0.0;
				s_sigma_minus(1, 0) = 0.0;
				s_sigma_minus(1, 1) = 1.0;
			}

			for (int it_id = 1; it_id < num_spins; it_id++)
			{
				if (it_id == s_id)
				{
					s_sigma_minus = Eigen::kroneckerProduct(s_sigma_minus, sigma_minus).eval();
				}
				else
				{
					s_sigma_minus = Eigen::kroneckerProduct(s_sigma_minus, sigma_0).eval();
				}
			}

			Eigen::MatrixXd s_sigma_minus_kron = Eigen::kroneckerProduct(s_sigma_minus, p_identity);

			for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
			{
				for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
				{
					int index = st_id_1 * md->sys_size + st_id_2;

					md->dissipators[s_id + 1][index].real = s_sigma_minus_kron(st_id_1, st_id_2);
					md->dissipators[s_id + 1][index].imag = 0.0;
				}
			}
		}
	}
}

void MBLNewDelBehaviour::init_dissipators(AllData * ad) const
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

	int diss_type = int(cp->params.find("diss_type")->second);
	double diss_phase = double(cp->params.find("diss_phase")->second);

	if (diss_type == 0)
	{
		for (int diss_id = 0; diss_id < md->num_diss; diss_id++)
		{
			for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
			{
				int index = st_id_1 * md->sys_size + st_id_1;
				md->dissipators[diss_id][index].real = double(bit_at(md->mbl_id_to_x[st_id_1], diss_id));
			}
		}
	}
	else if (diss_type == 1)
	{
		for (int diss_id = 0; diss_id < md->num_diss; diss_id++)
		{
			for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
			{
				int index = st_id_1 * md->sys_size + st_id_1;
				md->dissipators[diss_id][index].real = double(bit_at(md->mbl_id_to_x[st_id_1], diss_id)) - double(bit_at(md->mbl_id_to_x[st_id_1], diss_id + 1));
			}

			int tmp = 0;
			for (int state_id_1 = 0; state_id_1 < md->mbl_Ns; state_id_1++)
			{
				int row_id = state_id_1;

				tmp = bit_at(md->mbl_id_to_x[state_id_1], diss_id) - bit_at(md->mbl_id_to_x[state_id_1], diss_id + 1);

				int col_id = 0;

				if (tmp == 0)
				{
					col_id = state_id_1;
				}
				else
				{
					for (int state_id_2 = 0; state_id_2 < md->mbl_Ns; state_id_2++)
					{
						if (md->mbl_adjacement[md->mbl_id_to_x[state_id_1] ^ md->mbl_id_to_x[state_id_2]])
						{
							vector<int> adjacency_bits = convert_int_to_vector_of_bits(md->mbl_id_to_x[state_id_1] ^ md->mbl_id_to_x[state_id_2], md->mbl_Nc);
							vector<int> hop;
							for (int cell_id = 0; cell_id < md->mbl_Nc; cell_id++)
							{
								if (adjacency_bits[cell_id])
								{
									hop.push_back(cell_id);
								}
							}

							for (int ad_cell_id = 0; ad_cell_id < hop.size(); ad_cell_id++)
							{
								hop[ad_cell_id] = (md->mbl_Nc - 1) - hop[ad_cell_id];
							}

							if (hop[1] == diss_id)
							{
								if (bit_at(md->mbl_id_to_x[state_id_1], diss_id))
								{
									col_id = state_id_2;
									int index = row_id * md->sys_size + col_id;

									md->dissipators[diss_id][index].real = -cos(-diss_phase);
									md->dissipators[diss_id][index].imag = -sin(-diss_phase);
								}
								else
								{
									col_id = state_id_2;
									int index = row_id * md->sys_size + col_id;

									md->dissipators[diss_id][index].real = cos(diss_phase);
									md->dissipators[diss_id][index].imag = sin(diss_phase);
								}
							}

							adjacency_bits.clear();
							hop.clear();
						}
					}
				}
			}
		}
	}
}

void LndHamNewDelBehaviour::init_dissipators(AllData* ad) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	std::vector<sp_mtx> f_basis;

	cout << "F-basis generation started" << endl;
	std::vector<triplet> vec_triplets;
	vec_triplets.reserve(md->sys_size);
	for (int st_id = 0; st_id < md->sys_size; st_id++)
	{
		vec_triplets.push_back(triplet(st_id, st_id, std::complex<double>(1.0, 0.0)));
	}
	sp_mtx sp_eye(md->sys_size, md->sys_size);
	sp_eye.setFromTriplets(vec_triplets.begin(), vec_triplets.end());

	sp_mtx tmp = sp_eye / std::sqrt(double(md->sys_size));
	f_basis.push_back(tmp);

	double sqrt_2 = std::sqrt(2.0);

	for (auto i = 0; i < md->sys_size; i++)
	{
		for (auto j = i + 1; j < md->sys_size; j++)
		{
			std::vector<triplet> vec_triplets(2);
			vec_triplets[0] = triplet(i, j, std::complex<double>(1.0 / sqrt_2, 0.0));
			vec_triplets[1] = triplet(j, i, std::complex<double>(1.0 / sqrt_2, 0.0));
			tmp = sp_mtx(md->sys_size, md->sys_size);
			tmp.setFromTriplets(vec_triplets.begin(), vec_triplets.end());
			f_basis.push_back(tmp);

			vec_triplets[0] = triplet(i, j, std::complex<double>(0.0, -1.0 / sqrt_2));
			vec_triplets[1] = triplet(j, i, std::complex<double>(0.0, 1.0 / sqrt_2));
			tmp = sp_mtx(md->sys_size, md->sys_size);
			tmp.setFromTriplets(vec_triplets.begin(), vec_triplets.end());
			f_basis.push_back(tmp);
		}
	}

	for (auto i = 0; i < md->sys_size - 1; i++)
	{
		std::vector<triplet> vec_triplets(i + 2);
		for (auto j = 0; j < i + 1; j++)
		{
			vec_triplets[j] = triplet(j, j, std::complex<double>(1.0 / std::sqrt(double((i + 1) * (i + 2))), 0.0));
		}
		vec_triplets[i + 1] = triplet(i + 1, i + 1, std::complex<double>(-double(i + 1) / std::sqrt(double((i + 1) * (i + 2))), 0.0));
		tmp = sp_mtx(md->sys_size, md->sys_size);
		tmp.setFromTriplets(vec_triplets.begin(), vec_triplets.end());
		f_basis.push_back(tmp);
	}
	cout << "F-basis generation complete" << endl;

	int lndham_seed = int(cp->params.find("lndham_seed")->second);
	int lndham_num_seeds = int(cp->params.find("lndham_num_seeds")->second);

	int M = md->sys_size * md->sys_size - 1;

	cout << "G generation started" << endl;
	VSLStreamStatePtr stream;
	vslNewStream(&stream, VSL_BRNG_MCG31, 77778888);
	vslLeapfrogStream(stream, lndham_seed, lndham_num_seeds);
	double* disorder_real = new double[M * M];
	double* disorder_imag = new double[M * M];
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, M * M, disorder_real, 0.0, 1.0);
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, M * M, disorder_imag, 0.0, 1.0);

	Eigen::MatrixXcd X(M, M);
	for (auto st_id_1 = 0; st_id_1 < M; st_id_1++)
	{
		for (auto st_id_2 = 0; st_id_2 < M; st_id_2++)
		{
			auto index = st_id_1 * M + st_id_2;
			X(st_id_1, st_id_2) = std::complex<double>(0.5 * disorder_real[index], 0.5 * disorder_imag[index]);
		}
	}
	delete[] disorder_real;
	delete[] disorder_imag;

	Eigen::MatrixXcd G = X * X.adjoint();
	const auto trace = G.trace();
	G = double(md->sys_size) * G / trace.real();
	cout << "G generation complete" << endl;

	cout << "G eigen started" << endl;
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es;
	es.compute(G, true);
	cout << "G eigen complete" << endl;

	auto evals_tmp = es.eigenvalues();
	std::vector<std::complex<double>> G_evals(evals_tmp.data(), evals_tmp.data() + evals_tmp.rows() * evals_tmp.cols());
	auto evecs = es.eigenvectors();

	cout << "dissipators generation started" << endl;
	for (int k1 = 0; k1 < M; k1++)
	{
		sp_mtx diss = sp_mtx(md->sys_size, md->sys_size);
		for (int k2 = 0; k2 < M; k2++)
		{
			diss += evecs(k2, k1) * f_basis[k2 + 1];
		}
		diss *= std::sqrt(evals_tmp[k1]);

		md->dissipators_eigen.push_back(diss);
	}
	cout << "dissipators generation complete" << endl;
}

void IntegrableNewDelBehaviour::init_dissipators(AllData* ad) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	int integrable_N = int(cp->params.find("integrable_N")->second);
	int integrable_seed = int(cp->params.find("integrable_seed")->second);
	int integrable_num_seeds = int(cp->params.find("integrable_num_seeds")->second);
	double integrable_tau = double(cp->params.find("integrable_tau")->second);
	double integrable_k = double(cp->params.find("integrable_k")->second);

	int sys_size = md->sys_size;
	md->dissipators = new MKL_Complex16 * [md->num_diss];
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

	Eigen::MatrixXcd con = Eigen::MatrixXcd::Zero(4, 4);
	con(0, 0) = integrable_tau;
	con(3, 3) = integrable_k;
	con(1, 2) = 1;
	con(2, 1) = 1;

	Eigen::MatrixXcd eye = Eigen::MatrixXcd::Identity(2, 2);

	// ======================================================================================
	//                                  Left Dissipator (1, 2)
	// ======================================================================================
	Eigen::MatrixXcd left_diss = con;
	for (int spin_id = 2; spin_id <= integrable_N - 1; spin_id++)
	{
		left_diss = Eigen::kroneckerProduct(left_diss, eye).eval();
	}
	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			md->dissipators[0][index].real = left_diss(st_id_1, st_id_2).real();
			md->dissipators[0][index].imag = left_diss(st_id_1, st_id_2).imag();
		}
	}

	// ======================================================================================
	//                               Right Dissipator (N-1, N)
	// ======================================================================================
	Eigen::MatrixXcd right_diss = eye;
	for (int spin_id = 1; spin_id <= integrable_N - 3; spin_id++)
	{
		right_diss = Eigen::kroneckerProduct(right_diss, eye).eval();
	}
	right_diss = Eigen::kroneckerProduct(right_diss, con).eval();
	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			md->dissipators[integrable_N - 2][index].real = right_diss(st_id_1, st_id_2).real();
			md->dissipators[integrable_N - 2][index].imag = right_diss(st_id_1, st_id_2).imag();
		}
	}

	// ======================================================================================
	//                               All middle dissipators
	// ======================================================================================

	for (int diss_id = 1; diss_id <= integrable_N - 3; diss_id++)
	{
		Eigen::MatrixXcd curr_diss = eye;

		int left_s_id = 1;
		int left_f_id = diss_id - 1;
		for (int spin_id = left_s_id; spin_id <= left_f_id; spin_id++)
		{
			curr_diss = Eigen::kroneckerProduct(curr_diss, eye).eval();
		}

		curr_diss = Eigen::kroneckerProduct(curr_diss, con).eval();

		int right_s_id = diss_id + 2;
		int right_f_id = integrable_N - 1;
		for (int spin_id = right_s_id; spin_id <= right_f_id; spin_id++)
		{
			curr_diss = Eigen::kroneckerProduct(curr_diss, eye).eval();
		}
		for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
		{
			for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
			{
				int index = st_id_1 * md->sys_size + st_id_2;
				md->dissipators[diss_id][index].real = curr_diss(st_id_1, st_id_2).real();
				md->dissipators[diss_id][index].imag = curr_diss(st_id_1, st_id_2).imag();
			}
		}
	}

	// ======================================================================================
	//                               Border dissipator (N, 1)
	// ======================================================================================
	Eigen::MatrixXcd border_diss = Eigen::MatrixXcd::Zero(md->sys_size, md->sys_size);
	for (int state_id = 0; state_id <= md->sys_size - 1; state_id++)
	{
		int state_ind = state_id;
		vector<int> state_bits = convert_int_to_vector_of_bits(state_ind, integrable_N);

		if (state_bits[0] == state_bits[integrable_N - 1])
		{
			if (state_bits[0] == 0)
			{
				border_diss(state_id, state_id) = integrable_tau;
			}
			else
			{
				border_diss(state_id, state_id) = integrable_k;
			}
		}
		else
		{
			vector<int> inv_bits = state_bits;
			inv_bits[0] = state_bits[integrable_N - 1];
			inv_bits[integrable_N - 1] = state_bits[0];
			reverse(inv_bits.begin(), inv_bits.end());
			int row_id = accumulate(inv_bits.rbegin(), inv_bits.rend(), 0, [](int x, int y) { return (x << 1) + y; });
			border_diss(row_id, state_id) = 1;
		}
	}
	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;
			md->dissipators[integrable_N - 1][index].real = border_diss(st_id_1, st_id_2).real();
			md->dissipators[integrable_N - 1][index].imag = border_diss(st_id_1, st_id_2).imag();
		}
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

			hamitlonian_part[index].real = md->hamiltonian[index].real;
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

			md->drv_part[index].real = md->hamiltonian_drv[index].real;
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
			md->hamiltonians_qj[0][index].imag -= (hamitlonian_part[index].real + E_0 * md->hamiltonian_drv[index].real);

			md->hamiltonians_qj[1][index].real += (hamitlonian_part[index].imag + E_1 * 0.0);
			md->hamiltonians_qj[1][index].imag -= (hamitlonian_part[index].real + E_1 * md->hamiltonian_drv[index].real);
		}
	}

	delete[] diss_part;
	delete[] hamitlonian_part;
}

void DimerSyncNewDelBehaviour::init_hamiltonians_qj(AllData* ad) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	double diss_gamma = double(cp->params.find("diss_gamma")->second);
	diss_gamma = diss_gamma / double(md->sys_size - 1);

	MKL_Complex16 diss_gamma_cmplx = { diss_gamma, 0.0 };
	MKL_Complex16 zero_cmplx = { 0.0, 0.0 };

	md->hamiltonians_qj = new MKL_Complex16 * [md->num_ham_qj];
	for (int qj_ham_id = 0; qj_ham_id < md->num_ham_qj; qj_ham_id++)
	{
		md->hamiltonians_qj[qj_ham_id] = new MKL_Complex16[md->sys_size * md->sys_size];
	}

	MKL_Complex16* diss_part = new MKL_Complex16[md->sys_size * md->sys_size];
	MKL_Complex16* hamitlonian_part = new MKL_Complex16[md->sys_size * md->sys_size];
	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			diss_part[index].real = 0.0;
			diss_part[index].imag = 0.0;

			hamitlonian_part[index].real = md->hamiltonian[index].real;
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
	md->drv_part_2 = new MKL_Complex16[md->sys_size * md->sys_size];
	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;
			md->non_drv_part[index].real = hamitlonian_part[index].real;
			md->non_drv_part[index].imag = hamitlonian_part[index].imag;

			md->drv_part[index].real = md->hamiltonian_drv[index].real;
			md->drv_part[index].imag = 0.0;

			md->drv_part_2[index].real = md->hamiltonian_drv_2[index].real;
			md->drv_part_2[index].imag = 0.0;
		}
	}

	double prm_E = double(cp->params.find("dimersync_prm_E")->second);
	double drv_ampl_1 = double(cp->params.find("dimersync_drv_ampl_1")->second);

	double E_0 = prm_E + drv_ampl_1;
	double E_1 = prm_E - drv_ampl_1;

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			md->hamiltonians_qj[0][index].real += (hamitlonian_part[index].imag + E_0 * 0.0);
			md->hamiltonians_qj[0][index].imag -= (hamitlonian_part[index].real + E_0 * md->hamiltonian_drv[index].real);

			md->hamiltonians_qj[1][index].real += (hamitlonian_part[index].imag + E_1 * 0.0);
			md->hamiltonians_qj[1][index].imag -= (hamitlonian_part[index].real + E_1 * md->hamiltonian_drv[index].real);
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

			hamitlonian_part[index].real = md->hamiltonian[index].real;
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
			md->drv_part[index].imag = md->hamiltonian_drv[index].real;
		}
	}

	double E_0 = jcs_drv_ampl;
	double E_1 = 0.0;

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			md->hamiltonians_qj[0][index].real += (hamitlonian_part[index].imag + E_0 * md->hamiltonian_drv[index].real);
			md->hamiltonians_qj[0][index].imag -= (hamitlonian_part[index].real + E_0 * 0.0);

			md->hamiltonians_qj[1][index].real += (hamitlonian_part[index].imag + E_1 * md->hamiltonian_drv[index].real);
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

	double ps_diss_w = double(cp->params.find("ps_diss_w")->second);

	MKL_Complex16 diss_gamma_cav_cmplx = {diss_gamma_cav, 0.0};
	MKL_Complex16 diss_gamma_s = { ps_diss_w, 0.0 };
	MKL_Complex16 zero_cmplx = {0.0, 0.0};

	md->hamiltonians_qj = new MKL_Complex16*[md->num_ham_qj];
	for (int qj_ham_id = 0; qj_ham_id < md->num_ham_qj; qj_ham_id++)
	{
		md->hamiltonians_qj[qj_ham_id] = new MKL_Complex16[md->sys_size * md->sys_size];
	}

	MKL_Complex16 * diss_part_cav = new MKL_Complex16[md->sys_size * md->sys_size];
	MKL_Complex16 * diss_part_s = new MKL_Complex16[md->sys_size * md->sys_size];
	MKL_Complex16 * hamitlonian_part = new MKL_Complex16[md->sys_size * md->sys_size];
	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			diss_part_cav[index].real = 0.0;
			diss_part_cav[index].imag = 0.0;

			diss_part_s[index].real = 0.0;
			diss_part_s[index].imag = 0.0;

			hamitlonian_part[index].real = md->hamiltonian[index].real;
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

	int diss_type = int(cp->params.find("diss_type")->second);
	int num_spins = int(cp->params.find("ps_num_spins")->second);

	if (diss_type == 1)
	{
		for (int s_id = 0; s_id < num_spins; s_id++)
		{
			cblas_zgemm(
				CblasRowMajor,
				CblasConjTrans,
				CblasNoTrans,
				md->sys_size,
				md->sys_size,
				md->sys_size,
				&diss_gamma_s,
				md->dissipators[s_id + 1],
				md->sys_size,
				md->dissipators[s_id + 1],
				md->sys_size,
				&zero_cmplx,
				diss_part_s,
				md->sys_size
			);

			for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
			{
				for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
				{
					int index = st_id_1 * md->sys_size + st_id_2;
					hamitlonian_part[index].real += 0.5 * diss_part_s[index].imag;
					hamitlonian_part[index].imag -= 0.5 * diss_part_s[index].real;
				}
			}
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
			md->drv_part[index].imag = md->hamiltonian_drv[index].real;
		}
	}

	double E_0 = ps_drv_ampl;
	double E_1 = 0.0;

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			md->hamiltonians_qj[0][index].real += (hamitlonian_part[index].imag + E_0 * md->hamiltonian_drv[index].real);
			md->hamiltonians_qj[0][index].imag -= (hamitlonian_part[index].real + E_0 * 0.0);

			md->hamiltonians_qj[1][index].real += (hamitlonian_part[index].imag + E_1 * md->hamiltonian_drv[index].real);
			md->hamiltonians_qj[1][index].imag -= (hamitlonian_part[index].real + E_1 * 0.0);
		}
	}

	delete[] diss_part_cav;
	delete[] diss_part_s;
	delete[] hamitlonian_part;
}

void MBLNewDelBehaviour::init_hamiltonians_qj(AllData * ad) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	double diss_gamma = double(cp->params.find("diss_gamma")->second);

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

			hamitlonian_part[index].real = md->hamiltonian[index].real;
			hamitlonian_part[index].imag = 0.0;

			for (int qj_ham_id = 0; qj_ham_id < md->num_ham_qj; qj_ham_id++)
			{
				md->hamiltonians_qj[qj_ham_id][index].real = 0.0;
				md->hamiltonians_qj[qj_ham_id][index].imag = 0.0;
			}
		}
	}

	for (int diss_id = 0; diss_id < md->num_diss; diss_id++)
	{
		cblas_zgemm(
			CblasRowMajor,
			CblasConjTrans,
			CblasNoTrans,
			md->sys_size,
			md->sys_size,
			md->sys_size,
			&diss_gamma_cmplx,
			md->dissipators[diss_id],
			md->sys_size,
			md->dissipators[diss_id],
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
	}

	md->non_drv_part = new MKL_Complex16[md->sys_size * md->sys_size];
	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;
			md->non_drv_part[index].real = hamitlonian_part[index].real;
			md->non_drv_part[index].imag = hamitlonian_part[index].imag;
		}
	}

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;
			md->hamiltonians_qj[0][index].real += hamitlonian_part[index].imag;
			md->hamiltonians_qj[0][index].imag -= hamitlonian_part[index].real;
		}
	}

	delete[] diss_part;
	delete[] hamitlonian_part;
}

void LndHamNewDelBehaviour::init_hamiltonians_qj(AllData* ad) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	MKL_Complex16 one_cmplx = { 1.0, 0.0 };
	MKL_Complex16 zero_cmplx = { 0.0, 0.0 };

	md->hamiltonians_qj = new MKL_Complex16 * [md->num_ham_qj];
	for (int qj_ham_id = 0; qj_ham_id < md->num_ham_qj; qj_ham_id++)
	{
		md->hamiltonians_qj[qj_ham_id] = new MKL_Complex16[md->sys_size * md->sys_size];
	}

	MKL_Complex16* hamitlonian_part = new MKL_Complex16[md->sys_size * md->sys_size];
	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			hamitlonian_part[index].real = md->hamiltonian[index].real;
			hamitlonian_part[index].imag = md->hamiltonian[index].imag;

			for (int qj_ham_id = 0; qj_ham_id < md->num_ham_qj; qj_ham_id++)
			{
				md->hamiltonians_qj[qj_ham_id][index].real = 0.0;
				md->hamiltonians_qj[qj_ham_id][index].imag = 0.0;
			}
		}
	}

	for (int diss_id = 0; diss_id < md->num_diss; diss_id++)
	{
		ds_mtx diss_part_eigen = md->dissipators_eigen[diss_id].adjoint() * md->dissipators_eigen[diss_id];

		for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
		{
			for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
			{
				int index = st_id_1 * md->sys_size + st_id_2;
				hamitlonian_part[index].real += 0.5 * diss_part_eigen(st_id_1, st_id_2).imag();
				hamitlonian_part[index].imag -= 0.5 * diss_part_eigen(st_id_1, st_id_2).real();
			}
		}
	}

	md->non_drv_part = new MKL_Complex16[md->sys_size * md->sys_size];
	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;
			md->non_drv_part[index].real = hamitlonian_part[index].real;
			md->non_drv_part[index].imag = hamitlonian_part[index].imag;
		}
	}

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;
			md->hamiltonians_qj[0][index].real += hamitlonian_part[index].imag;
			md->hamiltonians_qj[0][index].imag -= hamitlonian_part[index].real;
		}
	}
	delete[] hamitlonian_part;
}

void IntegrableNewDelBehaviour::init_hamiltonians_qj(AllData* ad) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	MKL_Complex16 diss_gamma_cmplx = { 1, 0.0 };
	MKL_Complex16 zero_cmplx = { 0.0, 0.0 };

	md->hamiltonians_qj = new MKL_Complex16 * [md->num_ham_qj];
	for (int qj_ham_id = 0; qj_ham_id < md->num_ham_qj; qj_ham_id++)
	{
		md->hamiltonians_qj[qj_ham_id] = new MKL_Complex16[md->sys_size * md->sys_size];
	}

	MKL_Complex16* diss_part = new MKL_Complex16[md->sys_size * md->sys_size];
	MKL_Complex16* hamitlonian_part = new MKL_Complex16[md->sys_size * md->sys_size];
	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			diss_part[index].real = 0.0;
			diss_part[index].imag = 0.0;

			hamitlonian_part[index].real = md->hamiltonian[index].real;
			hamitlonian_part[index].imag = 0.0;

			for (int qj_ham_id = 0; qj_ham_id < md->num_ham_qj; qj_ham_id++)
			{
				md->hamiltonians_qj[qj_ham_id][index].real = 0.0;
				md->hamiltonians_qj[qj_ham_id][index].imag = 0.0;
			}
		}
	}

	for (int diss_id = 0; diss_id < md->num_diss; diss_id++)
	{
		cblas_zgemm(
			CblasRowMajor,
			CblasConjTrans,
			CblasNoTrans,
			md->sys_size,
			md->sys_size,
			md->sys_size,
			&diss_gamma_cmplx,
			md->dissipators[diss_id],
			md->sys_size,
			md->dissipators[diss_id],
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
	}

	md->non_drv_part = new MKL_Complex16[md->sys_size * md->sys_size];
	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;
			md->non_drv_part[index].real = hamitlonian_part[index].real;
			md->non_drv_part[index].imag = hamitlonian_part[index].imag;
		}
	}

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;
			md->hamiltonians_qj[0][index].real += hamitlonian_part[index].imag;
			md->hamiltonians_qj[0][index].imag -= hamitlonian_part[index].real;
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

	free_random_obs(ad);
}

void DimerSyncNewDelBehaviour::free_hamiltonians(AllData* ad) const
{
	MainData* md = ad->md;

	delete[] md->hamiltonian;
	delete[] md->hamiltonian_drv;
	delete[] md->hamiltonian_drv_2;

	free_random_obs(ad);
}

void JCSNewDelBehaviour::free_hamiltonians(AllData * ad) const
{
	MainData * md = ad->md;

	delete[] md->hamiltonian;
	delete[] md->hamiltonian_drv;
	delete[] md->special;
	delete[] md->special_2;
	delete[] md->special_3;

	free_random_obs(ad);
}

void PSNewDelBehaviour::free_hamiltonians(AllData * ad) const
{
	MainData * md = ad->md;

	delete[] md->hamiltonian;
	delete[] md->hamiltonian_drv;
	delete[] md->special;
	delete[] md->special_2;
	delete[] md->special_3;
	delete[] md->special_4;

	free_random_obs(ad);
}

void MBLNewDelBehaviour::free_hamiltonians(AllData * ad) const
{
	MainData * md = ad->md;

	delete[] md->hamiltonian;
	delete[] md->special;

	free_random_obs(ad);
}

void LndHamNewDelBehaviour::free_hamiltonians(AllData* ad) const
{
	MainData* md = ad->md;

	delete[] md->hamiltonian;
	delete[] md->special;

	free_random_obs(ad);
}

void IntegrableNewDelBehaviour::free_hamiltonians(AllData* ad) const
{
	MainData* md = ad->md;

	delete[] md->hamiltonian;
	delete[] md->special;

	free_random_obs(ad);
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

void DimerSyncNewDelBehaviour::free_dissipators(AllData* ad) const
{
	MainData* md = ad->md;

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

void MBLNewDelBehaviour::free_dissipators(AllData * ad) const
{
	MainData * md = ad->md;

	for (int diss_id = 0; diss_id < md->num_diss; diss_id++)
	{
		delete[] md->dissipators[diss_id];
	}
	delete[] md->dissipators;
}

void LndHamNewDelBehaviour::free_dissipators(AllData* ad) const
{
}

void IntegrableNewDelBehaviour::free_dissipators(AllData* ad) const
{
	MainData* md = ad->md;

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

void DimerSyncNewDelBehaviour::free_hamiltonians_qj(AllData* ad) const
{
	MainData* md = ad->md;

	delete[] md->non_drv_part;
	delete[] md->drv_part;
	delete[] md->drv_part_2;

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

void MBLNewDelBehaviour::free_hamiltonians_qj(AllData * ad) const
{
	MainData * md = ad->md;

	delete[] md->non_drv_part;

	for (int qj_ham_id = 0; qj_ham_id < md->num_ham_qj; qj_ham_id++)
	{
		delete[] md->hamiltonians_qj[qj_ham_id];
	}
	delete[] md->hamiltonians_qj;
}

void LndHamNewDelBehaviour::free_hamiltonians_qj(AllData* ad) const
{
	MainData* md = ad->md;

	delete[] md->non_drv_part;

	for (int qj_ham_id = 0; qj_ham_id < md->num_ham_qj; qj_ham_id++)
	{
		delete[] md->hamiltonians_qj[qj_ham_id];
	}
	delete[] md->hamiltonians_qj;
}

void IntegrableNewDelBehaviour::free_hamiltonians_qj(AllData* ad) const
{
	MainData* md = ad->md;

	delete[] md->non_drv_part;

	for (int qj_ham_id = 0; qj_ham_id < md->num_ham_qj; qj_ham_id++)
	{
		delete[] md->hamiltonians_qj[qj_ham_id];
	}
	delete[] md->hamiltonians_qj;
}




void init_random_obs(AllData* ad)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int num_random_obs = int(cp->params.find("num_random_obs")->second);
	int random_obs_seed = double(cp->params.find("random_obs_seed")->second);
	int random_obs_mns = double(cp->params.find("random_obs_mns")->second);
	int random_obs_type = double(cp->params.find("random_obs_type")->second);

	if (num_random_obs > 0)
	{
		md->random_obs_mtxs.reserve(num_random_obs);
		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			md->random_obs_mtxs.push_back(new MKL_Complex16[md->sys_size * md->sys_size]);
		}

		double* disorder_real = new double[md->sys_size * md->sys_size];
		double* disorder_imag = new double[md->sys_size * md->sys_size];

		MKL_Complex16* x_mtx = new MKL_Complex16[md->sys_size * md->sys_size];

		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			VSLStreamStatePtr stream;
			vslNewStream(&stream, VSL_BRNG_MCG31, 3873247);
			vslLeapfrogStream(stream, random_obs_seed + obs_id, random_obs_mns);

			vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, md->sys_size * md->sys_size, disorder_real, 0.0, 1.0);
			vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, md->sys_size * md->sys_size, disorder_imag, 0.0, 1.0);

			for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
			{
				for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
				{
					int index = st_id_1 * md->sys_size + st_id_2;

					md->random_obs_mtxs[obs_id][index].real = 0.0;
					md->random_obs_mtxs[obs_id][index].imag = 0.0;

					x_mtx[index].real = 0.5 * disorder_real[index];
					x_mtx[index].imag = 0.5 * disorder_imag[index];
				}
			}

			MKL_Complex16 alpha = { 1.0, 0.0 };
			MKL_Complex16 beta = { 0.0, 0.0 };

			cblas_zgemm(
				CblasRowMajor,
				CblasNoTrans,
				CblasConjTrans,
				md->sys_size,
				md->sys_size,
				md->sys_size,
				&alpha,
				x_mtx,
				md->sys_size,
				x_mtx,
				md->sys_size,
				&beta,
				md->random_obs_mtxs[obs_id],
				md->sys_size
			);


			if (random_obs_type > 0)
			{
				double trace = 0.0;

				for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
				{
					int index = st_id_1 * md->sys_size + st_id_1;
					trace += md->random_obs_mtxs[obs_id][index].real;
				}

				if (random_obs_type == 2)
				{
					trace = trace / double(md->sys_size);
				}

				for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
				{
					for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
					{
						int index = st_id_1 * md->sys_size + st_id_2;

						md->random_obs_mtxs[obs_id][index].real = md->random_obs_mtxs[obs_id][index].real / trace;
						md->random_obs_mtxs[obs_id][index].imag = md->random_obs_mtxs[obs_id][index].imag / trace;
					}
				}
			}
		}

		delete[] disorder_real;
		delete[] disorder_imag;

		delete[] x_mtx;
	}
}

void free_random_obs(AllData* ad)
{
	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);
	if (num_random_obs > 0)
	{
		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			delete[] ad->md->random_obs_mtxs[obs_id];
		}
	}
}
