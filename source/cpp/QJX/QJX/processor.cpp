#include "processor.h"

Processor::Processor(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	this->rp = rp;
	this->cp = cp;
	this->md = md;
	this->qjd = qjd;
}

void Processor::set_debug_behaviour(DebugBehavior* db)
{
	this->db = db;
}

void Processor::set_output_behaviour(OuputBehavior* ob)
{
	this->ob = ob;
}

void Processor::set_init_behaviour(InitBehavior* ib)
{
	this->ib = ib;
}

void Processor::set_free_behaviour(FreeBehavior* fb)
{
	this->fb = fb;
}

void Processor::set_qj_init_behaviour(QJInitBehavior* qj_ib)
{
	this->qj_ib = qj_ib;
}

void Processor::set_qj_free_behaviour(QJFreeBehavior* fb)
{
	this->qj_fb = qj_fb;
}

void Processor::set_qj_experiment_behaviour(QJExperimentBehavior* eb)
{
	this->qj_eb = qj_eb;
}

void Processor::process()
{
	ib->init_sizes(rp, cp, md);
	ib->init_hamiltonians(rp, cp, md);
	ib->init_dissipators(rp, cp, md);
	ib->init_hamiltonians_qj(rp, cp, md);

	ob->suffix_param(cp, 4);
	db->save(rp, cp, md);

	qj_ib->init_data(rp, cp, md, qjd);

	qj_eb->trans_process(rp, cp, md, qjd);
	qj_eb->obs_process(rp, cp, md, qjd);

	qj_fb->free_data(rp, cp, md, qjd);

	fb->free_hamiltonians(rp, cp, md);
	fb->free_dissipators(rp, cp, md);
	fb->free_hamiltonians_qj(rp, cp, md);
}